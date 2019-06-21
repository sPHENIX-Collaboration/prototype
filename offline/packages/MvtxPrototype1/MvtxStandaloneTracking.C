#include "MvtxStandaloneTracking.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefUtil.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TDecompSVD.h>
#include <TMath.h>

#include <iostream>
#include <utility>
#include <stdio.h>
#include <map>

MvtxStandaloneTracking::MvtxStandaloneTracking()
  : window_x_(10)
  , window_z_(10)
  , ghostrejection_(false)
  , verbosity_(0)
{

}

MvtxStandaloneTracking::~MvtxStandaloneTracking()
{
}

void
MvtxStandaloneTracking::RunTracking(PHCompositeNode* topNode, MvtxTrackList &trklst, std::vector<int> &lyrs)
{
  // check for appropriate layers
  if ( lyrs.size() < 3 || lyrs.size() > 4 )
  {
    std::cout << PHWHERE << "ERROR: Inappropriate number of input layers - " << lyrs.size() << std::endl;
    return;
  }

  trklst.clear();

  //------
  //--- get cluster container
  //------
  clusters_ = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
  if (!clusters_)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TrkrClusterContainer" << std::endl;
    return;
  }


  //------
  //--- associate clusters
  //------
  AssociateClusters(trklst, lyrs);

  if ( verbosity_ > 0 )
    std::cout << PHWHERE << " Finished associating clusters. N candidates:" << trklst.size() << std::endl;
  //------
  //--- fit tracks
  //------
  for ( unsigned int itrk = 0; itrk < trklst.size(); itrk++)
  {
    TrackFitXY(trklst.at(itrk));
    TrackFitZY(trklst.at(itrk));
  } // itrk

  if ( verbosity_ > 1 )
    PrintTrackCandidates(trklst);


  //------
  // --- choose best tracks
  //------
  if ( ghostrejection_ && trklst.size() > 1 )
    RunGhostRejection(trklst);

  // --- done
  return;
}

void
MvtxStandaloneTracking::AssociateClusters(MvtxTrackList &trklst, std::vector<int> &lyrs)
{

  // --- utility class
  TrkrDefUtil util;


  // --- loop over all clusters in the first desired layer
  TrkrClusterContainer::ConstRange clusrange0 =
    clusters_->GetClusters(TrkrDefs::TRKRID::mvtx_id, lyrs.at(0));
  for ( TrkrClusterContainer::ConstIterator iter0 = clusrange0.first;
        iter0 != clusrange0.second;
        ++iter0)
  {

    // -- loop over all clusters in the second desired layer
    TrkrClusterContainer::ConstRange clusrange1 =
      clusters_->GetClusters(TrkrDefs::TRKRID::mvtx_id, lyrs.at(1));
    for ( TrkrClusterContainer::ConstIterator iter1 = clusrange1.first;
          iter1 != clusrange1.second;
          ++iter1)
    {

      // get clusters
      TrkrCluster* clus0 = iter0->second;
      TrkrCluster* clus1 = iter1->second;

      // calculate slope & interecept in xy plane
      double mxy = CalcSlope(clus0->GetY(), clus0->GetX(), clus1->GetY(), clus1->GetX());
      double bxy = CalcIntecept(clus0->GetY(), clus0->GetX(), mxy);

      // calculate slope & interecept in zy plane
      double mzy = CalcSlope(clus0->GetY(), clus0->GetZ(), clus1->GetY(), clus1->GetZ());
      double bzy = CalcIntecept(clus0->GetY(), clus0->GetZ(), mxy);

      // -- loop over all clusters in the third desired layer
      TrkrClusterContainer::ConstRange clusrange2 =
        clusters_->GetClusters(TrkrDefs::TRKRID::mvtx_id, lyrs.at(2));
      for ( TrkrClusterContainer::ConstIterator iter2 = clusrange2.first;
            iter2 != clusrange2.second;
            ++iter2)
      {

        // check that the projection is within our window
        if ( fabs((iter2->second)->GetX() - CalcProjection((iter2->second)->GetY(), mxy, bxy)) < window_x_ &&
             fabs((iter2->second)->GetZ() - CalcProjection((iter2->second)->GetY(), mzy, bzy)) < window_z_ )
        {

          // make the candidate
          MvtxTrack trk;
          (trk.ClusterList).push_back(iter0->second);
          (trk.ClusterList).push_back(iter1->second);
          (trk.ClusterList).push_back(iter2->second);

          // if there's another layer, require it
          if ( lyrs.size() > 3 )
          {
            TrkrClusterContainer::ConstRange clusrange3 =
              clusters_->GetClusters(TrkrDefs::TRKRID::mvtx_id, lyrs.at(3));
            for ( TrkrClusterContainer::ConstIterator iter3 = clusrange3.first;
                  iter3 != clusrange3.second;
                  ++iter3)
            {

              // check that the projection is within our window
              if ( fabs((iter3->second)->GetX() - CalcProjection((iter3->second)->GetY(), mxy, bxy)) < window_x_ &&
                   fabs((iter3->second)->GetZ() - CalcProjection((iter3->second)->GetY(), mzy, bzy)) < window_z_ )
              {

                (trk.ClusterList).push_back(iter3->second);
                trklst.push_back(trk);
              }
            }
          }
          // else we're done
          else
          {
            trklst.push_back(trk);
          }

        }
      } // clusrange 2
    } // clusrange 1
  } // clusrange 0

}


void
MvtxStandaloneTracking::RunGhostRejection(MvtxTrackList &trklst)
{
  if ( verbosity_ > 0 )
    std::cout << PHWHERE << " Running Ghost Rejection on "
              << trklst.size() << " tracks" << std::endl;

  // --- First, make a map of all cluster keys & track index
  std::multimap<TrkrDefs::cluskey, unsigned int> key_trk_map;

  for (unsigned int itrk = 0; itrk < trklst.size(); itrk++)
  {
    for (unsigned int iclus = 0; iclus < trklst.at(itrk).ClusterList.size(); iclus++)
    {
      TrkrDefs::cluskey ckey = trklst.at(itrk).ClusterList.at(iclus)->GetClusKey();
      key_trk_map.insert(std::make_pair(ckey, itrk));
    } // iclus
  } // itrk

  // --- find clusters associated with more than one track and pick the best track
  std::set<unsigned int> remove_set;
  for ( auto iter = key_trk_map.begin(); iter != key_trk_map.end(); ++iter)
  {
    // get the upper bound for this key
    auto upiter = key_trk_map.upper_bound(iter->first);

    // iterate over common clusters and get the best track
    double chi2_best = 9999.;
    unsigned int idx_best = 0;
    int ntrk = 0;
    for ( auto jter = iter; jter != upiter; ++jter)
    {
      ntrk++;
      double chi2 = trklst.at(jter->second).chi2_xy + trklst.at(jter->second).chi2_zy;
      if ( chi2 < chi2_best && chi2 > 0 )
      {
        chi2_best = chi2;
        idx_best = jter->second;
      }
    }

    // FOR TESTING
    if ( ntrk > 1 && verbosity_ > 1 )
    {
      std::cout << PHWHERE << " Tracks sharing cluster:" << ntrk << std::endl;
      for ( auto jter = iter; jter != upiter; ++jter)
      {
        double chi2 = trklst.at(jter->second).chi2_xy + trklst.at(jter->second).chi2_zy;
        std::cout << "     "
                  << " trk idx: " << jter->second
                  << " chi2:" << chi2
                  << " m_xy:" << trklst.at(jter->second).m_xy;

        if ( jter->second == idx_best)
        {
          std::cout << "  <--- BEST" << std::endl;
        }
        else
        {
          std::cout << std::endl;
        }
      }
    }

    // remove all pairs that aren't the best
    for ( auto jter = iter; jter != upiter; ++jter)
    {
      if ( jter->second != idx_best )
      {
        remove_set.insert(jter->second);
      }
    }
  } // iter


  // --- Now remove tracks from the track list
  // reverse iterate so we don't have to keep track of indeces changed by removal
  if ( verbosity_ > 0 )
  {
    std::cout << PHWHERE << " List of idx to remove:";
    for ( auto it = remove_set.begin(); it != remove_set.end(); ++it)
    {
      std::cout << " " << *it;
    }
    std::cout << std::endl;
  }


  for ( auto rit = remove_set.rbegin(); rit != remove_set.rend(); ++rit)
    trklst.erase(trklst.begin() + *rit);

  // --- done

  if ( verbosity_ > 1 )
  {
    std::cout << PHWHERE << " Remaining tracks:" << std::endl;
    for ( unsigned int i = 0; i < trklst.size(); i++)
    {
      double chi2 = trklst.at(i).chi2_xy + trklst.at(i).chi2_zy;
      std::cout << "    chi2:" << chi2 << " m_xy:" << trklst.at(i).m_xy << std::endl;
    }
  }

  return;
}


TVectorD
MvtxStandaloneTracking::SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L)
{
  // Simple generalized linear least-squares fitter.
  // Solve y(X) = beta'*X + error(L) for beta (vector of regression coefs.)
  // L is the inverse covariance matrix for y.
  // Least-squares solution involves solving X' L X beta = X' L y

  TMatrixD XT(X); XT.T();
  TMatrixD A = XT * L * X;
  TVectorD b = XT * L * y;

  // Now solve A*beta = b using SVD. Decompose A = U S V'.
  TDecompSVD svd(A);
  TMatrixD UT = svd.GetU(); UT.T();

  // Construct Moore-Penrose pseudoinverse of S
  TVectorD s = svd.GetSig();
  TMatrixD Sd(s.GetNrows(), s.GetNrows());
  for (int i = 0; i < s.GetNrows(); i++)
    Sd(i, i) = s(i) > 0 ? 1. / s(i) : 0.;

  TVectorD beta = svd.GetV() * Sd * UT * b;

  return beta;
}

void
MvtxStandaloneTracking::TrackFitXY(MvtxTrack &trk)
{
  // Longitudinal/polar component
  // Fit track using a straight line x(y') = x0 + c*y'.
  // Assigns residuals, z-intercept, and polar angle.

  // m = # measurements; n = # parameters.
  int m = (trk.ClusterList).size(), n = 2;

  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);
  TVectorD x(m);

  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    x(iclus) = clus->GetX();
    y(iclus) = clus->GetY();

    X(iclus, 0) = 1;
    X(iclus, 1) = y(iclus);

    Cinv(iclus, iclus) = 2.0 * sqrt(clus->GetSize(0, 0));
  }

  TVectorD beta = SolveGLS(X, x, Cinv);
  double x0 = beta(0), c = beta(1);

  trk.m_xy = c;
  trk.b_xy = x0;

  double chi2  = 0;
  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    double cx = clus->GetX();
    double cy = clus->GetY();

    double px = cy * c + x0;

    double dx = cx - px;
    trk.dx.push_back(dx);

    chi2 += pow(dx, 2) / clus->GetError(0, 0);
  }
  chi2 /= double(m - 2);
  trk.chi2_xy = chi2;

  return;
}

void
MvtxStandaloneTracking::TrackFitZY(MvtxTrack &trk)
{
  // Longitudinal/polar component
  // Fit track using a straight line z(y') = z0 + c*y'.
  // Assigns residuals, z-intercept, and polar angle.

  // m = # measurements; n = # parameters.
  int m = (trk.ClusterList).size(), n = 2;

  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);
  TVectorD z(m);

  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    z(iclus) = clus->GetZ();
    y(iclus) = clus->GetY();

    X(iclus, 0) = 1;
    X(iclus, 1) = y(iclus);

    Cinv(iclus, iclus) = 2.0 * sqrt(clus->GetSize(2, 2));
  }

  TVectorD beta = SolveGLS(X, z, Cinv);
  double z0 = beta(0), c = beta(1);

  trk.m_zy = c;
  trk.b_zy = z0;

  double chi2 = 0;
  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    double cz = clus->GetZ();
    double cy = clus->GetY();

    double pz = cy * c + z0;

    double dz = cz - pz;
    trk.dz.push_back(dz);

    chi2 += pow(dz, 2) / clus->GetError(2, 2);
  }
  chi2 /= double(m - 2);
  trk.chi2_zy = chi2;

  return;
}

void
MvtxStandaloneTracking::PrintTrackCandidates(MvtxTrackList &trklst)
{
  std::cout << "===================================================" << std::endl;
  std::cout << "== " << PHWHERE << " Found " << trklst.size() << " Track Candidates" << std::endl;
  std::cout << "===================================================" << std::endl;
  for ( unsigned int i = 0; i < trklst.size(); i++)
  {
    std::cout << "== " << i << std::endl;
    for ( unsigned int j = 0; j < trklst.at(i).ClusterList.size(); j++)
    {
      std::cout << "    clus " << j
                << " key:0x" << std::hex << trklst.at(i).ClusterList.at(j)->GetClusKey() << std::dec
                << " (" << trklst.at(i).ClusterList.at(j)->GetX()
                << ", " << trklst.at(i).ClusterList.at(j)->GetY()
                << ", " << trklst.at(i).ClusterList.at(j)->GetZ()
                << ")"
                << " dx:" << trklst.at(i).dx.at(j)
                << " dz:" << trklst.at(i).dz.at(j)
                << std::endl;
    }
    std::cout << "    xy"
              << " m:" << trklst.at(i).m_xy
              << " b:" << trklst.at(i).b_xy
              << " chi2:" << trklst.at(i).chi2_xy
              << std::endl;
    std::cout << "    zy"
              << " m:" << trklst.at(i).m_zy
              << " b:" << trklst.at(i).b_zy
              << " chi2:" << trklst.at(i).chi2_zy
              << std::endl;
  } // i
  std::cout << "===================================================" << std::endl;
}

double
MvtxStandaloneTracking::CalcSlope(double x0, double y0, double x1, double y1)
{
  return (y1 - y0) / (x1 - x0);
}

double
MvtxStandaloneTracking::CalcIntecept(double x0, double y0, double m)
{
  return y0 - x0 * m;
}

double
MvtxStandaloneTracking::CalcProjection(double x, double m, double b)
{
  return m * x + b;
}


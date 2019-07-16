#include "MvtxPrototype2Eval.h"

#include "SegmentationAlpide.h"

#include <mvtx/MvtxDefs.h>
#include <mvtx/MvtxHit.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include <TMatrixF.h>
#include <TVector3.h>
#include <TH2D.h>
#include <TH1D.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

static const float twopi = 2.0 * M_PI;


MvtxPrototype2Eval::MvtxPrototype2Eval(const string &name)
  : SubsysReco(name)
  , _outfile("MvtxQAHisto.root")
  , m_hits(nullptr)
  , m_clusterlist(nullptr)
{
}

int MvtxPrototype2Eval::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  int chipColor[] = {kBlue, kRed, kGreen+2, kMagenta+2};

  h1d_nevents = new TH1D("h1d_nevents", ";Nevents", 1, -.5, .5);
  for (int ilay=0; ilay<NLAYER; ++ilay)
  {
    h2d_hit[ilay] = new TH2D(Form("h2d_hit_lay%d", ilay),"", 9 * NCOL, -.5, (9*NCOL) - .5, NROW, -.5, NROW - .5);
    h2d_clu[ilay] = new TH2D(Form("h2d_clu_lay%d", ilay),"", 9 * NCOL, -.5, (9*NCOL) - .5, NROW, -.5, NROW - .5);

    // difference in mean pixel index per event (from chip 0)
    h1d_clu_diffcol[ilay] = new TH1D(Form("h1d_clu_diffcol_lay%d", ilay), ";col index diff", 1023, -511.5, 511.5);
    h1d_clu_diffcol[ilay]->SetLineWidth(2);
    h1d_clu_diffcol[ilay]->SetLineColor(chipColor[ilay]);

    h1d_clu_diffrow[ilay] = new TH1D(Form("h1d_clu_diffrow_lay%d", ilay), ";row index diff", 1023, -511.5, 511.5);
    h1d_clu_diffrow[ilay]->SetLineWidth(2);
    h1d_clu_diffrow[ilay]->SetLineColor(chipColor[ilay]);
  }

  h1d_clus_size = new TH1D("h1d_clus_size","",50,0.5,50.5);

  return Fun4AllReturnCodes::EVENT_OK;

}

int MvtxPrototype2Eval::InitRun(PHCompositeNode *topNode)
{
  _n_events = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxPrototype2Eval::End(PHCompositeNode *topNode)
{

  cout << "-----MvtxQAHisto::End------(" << _n_events << ")" << endl;

  PHTFileServer::get().cd( _outfile );
  PHTFileServer::get().write( _outfile );

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxPrototype2Eval::process_event(PHCompositeNode *topNode)
{
  _n_events++;
  h1d_nevents->Fill(0);
  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }


  // loop over each MvtxHitSet object (chip)
  TrkrHitSetContainer::ConstRange hitsetrange =
    m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
      hitsetitr != hitsetrange.second;
      ++hitsetitr)
  {

    TrkrHitSet *hitset = hitsetitr->second;

    int layer = TrkrDefs::getLayer(hitset->getHitSetKey());
    int stave = MvtxDefs::getStaveId(hitset->getHitSetKey());
    int chip  = MvtxDefs::getChipId(hitset->getHitSetKey());

    if ( Verbosity() > 10 )
    {
      cout << "layer: " << layer << " stave: " << stave << " chip: " << chip << endl;
    }

    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
        hitr != hitrangei.second;
        ++hitr)
    {

      int col = MvtxDefs::getCol( hitr->first);
      int row = MvtxDefs::getRow( hitr->first);
      int stave_col = chip * SegmentationAlpide::NCols + col;

      //if (!MvtxPrototype2Geom::detectorToStave(chip, row, col, sRow, sCol)){
      //  cout << "ERROR" << endl;
      //  return Fun4AllReturnCodes::ABORTEVENT;
      //}
      h2d_hit[layer]->Fill(stave_col, row);

      //cout << "col: " << col << " row: " << row << endl;

    }
  }

  //Cluster
  TrkrClusterContainer::ConstRange clusrange = m_clusterlist->getClusters();

  int   n_clu_per_layer[NLAYER] = {0};
  float mean_clu_row[NLAYER] = {0};
  float mean_clu_col[NLAYER] = {0};
  for (TrkrClusterContainer::ConstIterator clusitr = clusrange.first;
      clusitr != clusrange.second;
      ++clusitr)
  {
    int layer = TrkrDefs::getLayer(clusitr->first);
    //int stave = MvtxDefs::getStaveId(clusitr->first);
    int chip  = MvtxDefs::getChipId(clusitr->first);

    TrkrCluster *clus = clusitr->second;
    int stave_col = chip * SegmentationAlpide::NCols + (int)clus->getZ();
    n_clu_per_layer[layer]++;
    mean_clu_col[layer] += stave_col;
    mean_clu_row[layer] += clus->getX();
   // if (!MvtxPrototype2Geom::detectorToStave(chip, (int)clus->getX(), (int)clus->getZ(), sRow, sCol)){
   //   cout << "ERROR" << endl;
   //   return Fun4AllReturnCodes::ABORTEVENT;
   // }
    //cout << chip << " " << (int)clus->getX() << " " << (int)clus->getZ() << " " << sRow << " " << sCol << endl;
    h2d_clu[layer]->Fill(stave_col, clus->getX());

    h1d_clus_size->Fill(clus->getAdc());
  }

  for ( int ilay = 0; ilay < NLAYER; ++ilay )
  {
    if ( n_clu_per_layer[ilay] > 0 )
    {
      mean_clu_col[ilay] /= (float)n_clu_per_layer[ilay];
      mean_clu_row[ilay] /= (float)n_clu_per_layer[ilay];
    }
    if ( mean_clu_col[0]  >= 0 && mean_clu_row[0] )
    {
      h1d_clu_diffcol[ilay]->Fill(mean_clu_col[ilay] - mean_clu_col[0]);
      h1d_clu_diffrow[ilay]->Fill(mean_clu_row[ilay] - mean_clu_row[0]);
    }
  }

  // done
  return Fun4AllReturnCodes::EVENT_OK;
}
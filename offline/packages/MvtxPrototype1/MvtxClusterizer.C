#include "MvtxClusterizer.h"

#include "mvtx/MvtxDefUtil.h"
#include "mvtx/MvtxHitSetv1.h"

#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetv1.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>
#include <g4detectors/PHG4CylinderGeom_Siladders.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

#include <TMatrixF.h>
#include <TVector3.h>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
using namespace boost;

#include <iostream>
#include <stdexcept>
#include <cmath>

using namespace std;

static const float twopi = 2.0 * M_PI;

bool MvtxClusterizer::are_adjacent(const pixel lhs,
                                   const pixel rhs )
{

  if ( GetZClustering() )
  {
    // column is first, row is second
    if ( fabs(lhs.first - rhs.first) <= 1 )
    {
      if ( fabs(lhs.second - rhs.second) <= 1 )
      {
        return true;
      }
    }
  }
  else
  {
    if ( fabs(lhs.first - rhs.first) == 0 )
    {
      if ( fabs(lhs.second - rhs.second) <= 1 )
      {
        return true;
      }
    }
  }

  return false;
}


MvtxClusterizer::MvtxClusterizer(const string &name) :
  SubsysReco(name),
  hits_(NULL),
  clusterlist_(NULL),
  make_z_clustering_(true),
  _timer(PHTimeServer::get()->insert_new(name))
{

}

int MvtxClusterizer::InitRun(PHCompositeNode* topNode)
{

  //-----------------
  // Add Cluster Node
  //-----------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVX node if required
  PHCompositeNode* svxNode
    = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "TRKR"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("TRKR");
    dstNode->addNode(svxNode);
  }

  // Create the Cluster node if required
  TrkrClusterContainer *trkrclusters
    = findNode::getClass<TrkrClusterContainer>(dstNode, "TrkrClusterContainer");
  if (!trkrclusters) {
    trkrclusters = new TrkrClusterContainer();
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
      new PHIODataNode<PHObject>(trkrclusters, "TrkrClusterContainer", "PHObject");
    svxNode->addNode(TrkrClusterContainerNode);
  }

  //----------------
  // Report Settings
  //----------------

  if (verbosity > 0)
  {
    cout << "====================== MvtxClusterizer::InitRun() =====================" << endl;
    cout << " Z-dimension Clustering = " << boolalpha << make_z_clustering_ << noboolalpha << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxClusterizer::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();

  // get node containing the digitized hits
  hits_ = findNode::getClass<TrkrHitSetContainer>(topNode, "TrkrHitSetContainer");
  if (!hits_)
  {
    cout << PHWHERE << "ERROR: Can't find node TrkrHitSetContainer" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for clusters
  clusterlist_ = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
  if (!clusterlist_)
  {
    cout << PHWHERE << " ERROR: Can't find TrkrClusterContainer." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
	clusterlist_->Reset();

  // run clustering
	ClusterMvtx(topNode);
  PrintClusters(topNode);

  // done
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxClusterizer::ClusterMvtx(PHCompositeNode *topNode) {

  if (verbosity > 0)
    cout << "Entering MvtxClusterizer::ClusterMvtx " << endl;

  //----------
  // Get Nodes
  //----------

  //-----------
  // Clustering
  //-----------

  // loop over each MvtxHitSet object (chip)
  TrkrHitSetContainer::ConstRange hitsetrange =
    //hits_->GetHitSets(TrkrDefs::TRKRID::mvtx_id);
    hits_->GetHitSets();
  for ( TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
        hitsetitr != hitsetrange.second;
        ++hitsetitr )
  {

    MvtxHitSetv1* hitset = static_cast<MvtxHitSetv1*>(hitsetitr->second);
		if (verbosity>2)
			hitset->identify();

    TrkrDefUtil trkrutil;
    MvtxDefUtil mvtxutil;


    // fill a vector of hits to make things easier
    // D. McGlinchey - this probably isn't necessary. Use iterator directly?
    std::vector<pixel> hitvec;
    MvtxHitSetv1::ConstRange hitrangei = hitset->GetHits();
    for ( MvtxHitSetv1::ConstIterator hitr = hitrangei.first;
          hitr != hitrangei.second;
          ++hitr)
    {
      hitvec.push_back(make_pair(hitr->first, hitr->second));
    }
		if ( verbosity>2 )
			cout << "hitvec.size(): " << hitvec.size() << endl;

    // do the clustering
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    // loop over hits in this chip
    for ( unsigned int i = 0; i < hitvec.size(); i++)
    {
      for ( unsigned int j = 0; j < hitvec.size(); j++)
      {
        if ( are_adjacent(hitvec[i], hitvec[j]) )
          add_edge(i, j, G);
      }
    }


    // Find the connections between the vertices of the graph (vertices are the rawhits,
    // connections are made when they are adjacent to one another)
    vector<int> component(num_vertices(G));

    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]);

    // Loop over the components(hit cells) compiling a list of the
    // unique connected groups (ie. clusters).
    set<int> cluster_ids; // unique components
    multimap<int, pixel> clusters;
    for (unsigned int i = 0; i < component.size(); i++)
    {
      cluster_ids.insert( component[i] );
      clusters.insert( make_pair(component[i], hitvec[i]) );
    }

    // loop over the componenets and make clusters
    for (set<int>::iterator clusiter = cluster_ids.begin();
         clusiter != cluster_ids.end();
         clusiter++ )
    {

      int clusid = *clusiter;
      pair< multimap<int, pixel>::iterator,
            multimap<int, pixel>::iterator> clusrange = clusters.equal_range(clusid);

      multimap<int, pixel>::iterator mapiter = clusrange.first;

      if (verbosity > 2)
        cout << "Filling cluster id " << clusid << endl;

      // make cluster
      TrkrDefs::cluskey ckey = mvtxutil.GenClusKey(hitset->GetHitSetKey(),clusid);
      TrkrClusterv1* clus = static_cast<TrkrClusterv1*>((clusterlist_->FindOrAddCluster(ckey))->second);

      // determine the size of the cluster in phi and z
      set<int> phibins;
      set<int> zbins;

      // determine the cluster position...
      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      unsigned nhits = 0;


      for (mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ )
      {

        // size
        zbins.insert((mapiter->second).first);
        phibins.insert((mapiter->second).second);

        // find the center of the pixel in local coords
        xsum += (mapiter->second).second;
        zsum += (mapiter->second).first + 0.5;

        ++nhits;
      } //mapitr

      float thickness = 50.e-4/28e-4;
      float phisize = phibins.size();
      float zsize = zbins.size();

      double clusx = NAN;
      double clusy = NAN;
      double clusz = NAN;

      clusx = xsum / nhits;
      clusy = ysum / nhits;
      clusz = zsum / nhits;

      clus->SetPosition(0, clusx);
      clus->SetPosition(1, clusy);
      clus->SetPosition(2, clusz);
      clus->SetLocal();
      
      clus->SetAdc(nhits);

      float invsqrt12 = 1.0 / sqrt(12.0);

      TMatrixF DIM(3, 3);
      DIM[0][0] = pow(0.5 * phisize, 2);

      DIM[0][1] = 0.0;
      DIM[0][2] = 0.0;
      DIM[1][0] = 0.0;
      DIM[1][1] = pow(0.5 * thickness, 2);
      DIM[1][2] = 0.0;
      DIM[2][0] = 0.0;
      DIM[2][1] = 0.0;
      DIM[2][2] = pow(0.5 * zsize, 2);

      TMatrixF ERR(3, 3);
      ERR[0][0] = pow(0.5 * phisize * invsqrt12, 2);
      ERR[0][1] = 0.0;
      ERR[0][2] = 0.0;
      ERR[1][0] = 0.0;
      ERR[1][1] = pow(0.5 * thickness * invsqrt12, 2);
      ERR[1][2] = 0.0;
      ERR[2][0] = 0.0;
      ERR[2][1] = 0.0;
      ERR[2][2] = pow(0.5 * zsize * invsqrt12, 2);


      clus->SetSize( 0 , 0 , DIM[0][0] );
      clus->SetSize( 0 , 1 , DIM[0][1] );
      clus->SetSize( 0 , 2 , DIM[0][2] );
      clus->SetSize( 1 , 0 , DIM[1][0] );
      clus->SetSize( 1 , 1 , DIM[1][1] );
      clus->SetSize( 1 , 2 , DIM[1][2] );
      clus->SetSize( 2 , 0 , DIM[2][0] );
      clus->SetSize( 2 , 1 , DIM[2][1] );
      clus->SetSize( 2 , 2 , DIM[2][2] );

      clus->SetError( 0 , 0 , ERR[0][0] );
      clus->SetError( 0 , 1 , ERR[0][1] );
      clus->SetError( 0 , 2 , ERR[0][2] );
      clus->SetError( 1 , 0 , ERR[1][0] );
      clus->SetError( 1 , 1 , ERR[1][1] );
      clus->SetError( 1 , 2 , ERR[1][2] );
      clus->SetError( 2 , 0 , ERR[2][0] );
      clus->SetError( 2 , 1 , ERR[2][1] );
      clus->SetError( 2 , 2 , ERR[2][2] );


      if ( verbosity > 2 )
        clus->identify();

    } // clusitr
  } // hitsetitr

  return;
}


void MvtxClusterizer::PrintClusters(PHCompositeNode * topNode)
{

  if (verbosity >= 1) {

    TrkrClusterContainer *clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
    if (!clusterlist) return;

    cout << "================= MvtxClusterizer::process_event() ====================" << endl;

    cout << " Found and recorded the following " << clusterlist->size() << " clusters: " << endl;

    clusterlist->identify();

    cout << "===========================================================================" << endl;
  }

  return;
}

#include "MvtxPrototype2Clusterizer.h"
#include "MvtxPrototype2Geom.h"

//#include <g4detectors/PHG4CylinderCellGeom.h>
//#include <g4detectors/PHG4CylinderCellGeomContainer.h>
//#include <g4detectors/PHG4CylinderGeom.h>
//#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <mvtx/MvtxDefs.h>
#include <mvtx/MvtxHit.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include <TMatrixF.h>
#include <TVector3.h>

#define BOOST_NO_HASH  // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
using namespace boost;

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

static const float twopi = 2.0 * M_PI;

bool MvtxPrototype2Clusterizer::are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs)
{
  if (GetZClustering())
  {
    // column is first, row is second
    if (fabs( MvtxDefs::getCol(lhs.first) - MvtxDefs::getCol(rhs.first) ) <= 1)
    {
      if (fabs( MvtxDefs::getRow(lhs.first) - MvtxDefs::getRow(rhs.first) ) <= 1)
      {
        return true;
      }
    }
  }
  else
  {
    if (fabs( MvtxDefs::getCol(lhs.first) - MvtxDefs::getCol(rhs.first) ) == 0)
    {
      if (fabs( MvtxDefs::getRow(lhs.first) - MvtxDefs::getRow(rhs.first) ) <= 1)
      {
        return true;
      }
    }
  }

  return false;
}

MvtxPrototype2Clusterizer::MvtxPrototype2Clusterizer(const string &name)
  : SubsysReco(name)
  , m_hits(nullptr)
  , m_clusterlist(nullptr)
  , m_clusterhitassoc(nullptr)
  , m_makeZClustering(true)
{
}

int MvtxPrototype2Clusterizer::InitRun(PHCompositeNode *topNode)
{
  //-----------------
  // Add Cluster Node
  //-----------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVX node if required
  PHCompositeNode *svxNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", "TRKR"));
  if (!svxNode)
  {
    svxNode = new PHCompositeNode("TRKR");
    dstNode->addNode(svxNode);
  }

  // Create the Cluster node if required
  TrkrClusterContainer *trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
      dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
      {
	DetNode = new PHCompositeNode("TRKR");
	dstNode->addNode(DetNode);
      }

    trkrclusters = new TrkrClusterContainer();
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
      new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  TrkrClusterHitAssoc *clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode,"TRKR_CLUSTERHITASSOC");
  if(!clusterhitassoc)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
	{
	  DetNode = new PHCompositeNode("TRKR");
	  dstNode->addNode(DetNode);
	}

      clusterhitassoc = new TrkrClusterHitAssoc();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(clusterhitassoc, "TRKR_CLUSTERHITASSOC", "PHObject");
      DetNode->addNode(newNode);
    }

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== MvtxPrototype2Clusterizer::InitRun() =====================" << endl;
    cout << " Z-dimension Clustering = " << boolalpha << m_makeZClustering << noboolalpha << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxPrototype2Clusterizer::process_event(PHCompositeNode *topNode)
{
  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // run clustering
  ClusterMvtx(topNode);
  PrintClusters(topNode);

  // done
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxPrototype2Clusterizer::ClusterMvtx(PHCompositeNode *topNode)
{

	if (Verbosity() > 0)
		cout << "Entering MvtxPrototype2Clusterizer::ClusterMvtx " << endl;

	//-----------
	// Clustering
	//-----------

	// loop over each MvtxHitSet object (chip)
	TrkrHitSetContainer::ConstRange hitsetrange =
		m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
	for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
			hitsetitr != hitsetrange.second;
			++hitsetitr)
	{

		TrkrHitSet *hitset = hitsetitr->second;

		if(Verbosity() > 10) cout << "MvtxPrototype2Clusterizer found hitsetkey " << hitsetitr->first << endl;

		if (Verbosity() > 10)
			hitset->identify();

		// fill a vector of hits to make things easier
		std::vector <std::pair< TrkrDefs::hitkey, TrkrHit*> > hitvec;

		TrkrHitSet::ConstRange hitrangei = hitset->getHits();
		for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
				hitr != hitrangei.second;
				++hitr)
		{
			hitvec.push_back(make_pair(hitr->first, hitr->second));
		}
		if (Verbosity() > 10)
			cout << "hitvec.size(): " << hitvec.size() << endl;

		// do the clustering
		typedef adjacency_list<vecS, vecS, undirectedS> Graph;
		Graph G;

		// loop over hits in this chip
		for (unsigned int i = 0; i < hitvec.size(); i++)
		{
			for (unsigned int j = 0; j < hitvec.size(); j++)
			{
				if (are_adjacent(hitvec[i], hitvec[j]))
					add_edge(i, j, G);
			}
		}

		// Find the connections between the vertices of the graph (vertices are the rawhits,
		// connections are made when they are adjacent to one another)
		vector<int> component(num_vertices(G));

		// this is the actual clustering, performed by boost
		connected_components(G, &component[0]);

		// Loop over the components(hits) compiling a list of the
		// unique connected groups (ie. clusters).
		set<int> cluster_ids;  // unique components
		//multimap<int, pixel> clusters;
		multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*> >  clusters;
		for (unsigned int i = 0; i < component.size(); i++)
		{
			cluster_ids.insert(component[i]);
			clusters.insert(make_pair(component[i], hitvec[i]));
		}

		// loop over the componenets and make clusters
		for (set<int>::iterator clusiter = cluster_ids.begin(); clusiter != cluster_ids.end(); ++clusiter)
		{
			int clusid = *clusiter;
			pair<multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator,
				multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator>  clusrange = clusters.equal_range(clusid);
			multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator mapiter = clusrange.first;

			if (Verbosity() > 2)
				cout << "Filling cluster id " << clusid << endl;

			// make the cluster directly in the node tree
			TrkrDefs::cluskey ckey = MvtxDefs::genClusKey(hitset->getHitSetKey(), clusid);
			TrkrClusterv1 *clus = static_cast<TrkrClusterv1 *>((m_clusterlist->findOrAddCluster(ckey))->second);

			// we need the geometry object for this layer to get the global positions
			//int layer = TrkrDefs::getLayer(ckey);
			//int stave =  MvtxDefs::getStaveId(ckey);
			//int chip = MvtxDefs::getChipId(ckey);

			// determine the size of the cluster in phi and z
			set<int> phibins;
			set<int> zbins;

			// determine the cluster position...
			// need to define x, y, z coordinate
			double xsum = 0.0;
			double ysum = 0.0;
			double zsum = 0.0;
			unsigned nhits = 0;

			double clusx = NAN;
			double clusy = NAN;
			double clusz = NAN;

			for (mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
			{
				// size
				int col =  MvtxDefs::getCol( (mapiter->second).first);
				int row = MvtxDefs::getRow( (mapiter->second).first);
				zbins.insert(col);
				phibins.insert(row);

				zsum += col;
				xsum += row;
				ysum += 0;

				// add the association between this cluster key and this hitkey to the table
				m_clusterhitassoc->addAssoc(ckey, mapiter->second.first);

				++nhits;
			}//mapiter

			float thickness = 50.e-4/28e-4; //sensor thickness converted to pixel size
			float phisize = phibins.size();
			float zsize = zbins.size();

			// This is the local position
			clusx = xsum/nhits + 0.5;
			clusy = ysum/nhits + 0.5;
			clusz = zsum/nhits + 0.5;
			clus->setAdc(nhits);

			/*
			cout << "new mvtx clusterizer layer: " << layer << " stave: " << stave << " chip: " << chip
				<< " x: " << clusx << " y: " << clusy << " z: " << clusz << endl;
			*/

			clus->setPosition(0, clusx);
			clus->setPosition(1, clusy);
			clus->setPosition(2, clusz);
			clus->setLocal();

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

			clus->setSize( 0 , 0 , DIM[0][0] );
			clus->setSize( 0 , 1 , DIM[0][1] );
			clus->setSize( 0 , 2 , DIM[0][2] );
			clus->setSize( 1 , 0 , DIM[1][0] );
			clus->setSize( 1 , 1 , DIM[1][1] );
			clus->setSize( 1 , 2 , DIM[1][2] );
			clus->setSize( 2 , 0 , DIM[2][0] );
			clus->setSize( 2 , 1 , DIM[2][1] );
			clus->setSize( 2 , 2 , DIM[2][2] );

			clus->setError( 0 , 0 , ERR[0][0] );
			clus->setError( 0 , 1 , ERR[0][1] );
			clus->setError( 0 , 2 , ERR[0][2] );
			clus->setError( 1 , 0 , ERR[1][0] );
			clus->setError( 1 , 1 , ERR[1][1] );
			clus->setError( 1 , 2 , ERR[1][2] );
			clus->setError( 2 , 0 , ERR[2][0] );
			clus->setError( 2 , 1 , ERR[2][1] );
			clus->setError( 2 , 2 , ERR[2][2] );

			if (Verbosity() > 2)
				clus->identify();

		}
	}

	if(Verbosity() > 5)
	{
		// check that the associations were written correctly
		m_clusterhitassoc->identify();
	}

  return;
}

void MvtxPrototype2Clusterizer::PrintClusters(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
  {
    TrkrClusterContainer *clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!clusterlist) return;

    cout << "================= Aftyer MvtxPrototype2Clusterizer::process_event() ====================" << endl;

    cout << " There are " << clusterlist->size() << " clusters recorded: " << endl;

    clusterlist->identify();

    cout << "===========================================================================" << endl;
  }

  return;
}

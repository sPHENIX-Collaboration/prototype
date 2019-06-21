////////////////////////////////////////////////////////////////////////////////
//
// This module is desgined to grab svtx tracks and put truth and cluster
// information into a TTree for GenFit testing
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 1 Apr 2016
//
////////////////////////////////////////////////////////////////////////////////


#include "MvtxQAHisto.h"

#include <phool/phool.h>
#include <phool/getClass.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

#include <trackbase/TrkrDefUtil.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetv1.h>
#include <mvtx/MvtxDefUtil.h>
#include <mvtx/MvtxHitSetv1.h>

#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TVector3.h>

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
MvtxQAHisto::MvtxQAHisto(const string &name):
  SubsysReco( name )
{
  //initialize
  _outfile = "MvtxQAHisto.root";

	beam_x[0] = 826.693;
	beam_x[1] = 822.267;
	beam_x[2] = 818.413;
	beam_x[3] = 830.190;

	beam_y[0] = 158.773;
	beam_y[1] = 167.537;
	beam_y[2] = 181.318;
	beam_y[3] = 190.988;

}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int MvtxQAHisto::Init(PHCompositeNode *topNode)
{

	//recoConsts *rc = recoConsts::instance();
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

	for (int ichip=0; ichip<4; ichip++){
		h2d_hit[ichip] = new TH2F(Form("h2d_hit_chip%d",ichip),"",1024,-0.5,1023.5,512,-0.5,511.5);
		h1d_hit_per_evt[ichip] = new TH1F(Form("h1d_hit_per_evt_chip%d",ichip),"",101,-0.5,100);
		h1d_clus_per_evt[ichip] = new TH1F(Form("h1d_clus_per_evt_chip%d",ichip),"",51,-0.5,50.5);

		h2d_hit_beam[ichip] = new TH2F(Form("h2d_hit_beam_chip%d",ichip),"",1025,-512.5,512.5,513,-256.5,256.5);
		h2d_hit_trk[ichip] = new TH2F(Form("h2d_hit_trk_chip%d",ichip),"",1025,-512.5,512.5,513,-256.5,256.5);
		h2d_clus[ichip] = new TH2F(Form("h2d_clus_chip%d",ichip),"",1024,0,3.0,512,0,1.5);
		h2d_clus_beam[ichip] = new TH2F(Form("h2d_clus_beam_chip%d",ichip),"",1024,-3.0/2,3.0/2,512,-1.5/2,1.5/2);

		h1d_clus_size_x[ichip] = new TH1F(Form("h1d_clus_size_x_chip%d",ichip),"",51,-0.5,50.5);
		h1d_clus_size_z[ichip] = new TH1F(Form("h1d_clus_size_z_chip%d",ichip),"",51,-0.5,50.5);

	}
	h1d_trk_finder_x = new TH1F("h1d_trk_finder_x","",513,-256.5,256.5);
	h1d_trk_finder_z = new TH1F("h1d_trk_finder_z","",1025,-512.5,512.5);
	h2d_trk_finder = new TH2F("h2d_trk_finder","",1025,-512.5,512.5,513,-256.5,256.5);

	h1d_clus_associated = new TH1F("h1d_clus_associated","",11,-0.5,10.5);
	h1d_clus_eff = new TH1F("h1d_clus_eff","",5,-0.5,4.5);

  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int MvtxQAHisto::process_event(PHCompositeNode *topNode)
{
  _event++;
  //if (1)
  if (_event % 1000 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

	TrkrDefUtil trkrdefutil;
	MvtxDefUtil mvtxdefutil;

	int nhit_per_chip[4] = {0};
	int nclus_per_chip[4] = {0};

	double avg_x = 0.0;
	double avg_y = 0.0;
	double norm = 0.0;

	h1d_trk_finder_x->Reset();
	h1d_trk_finder_z->Reset();
	h2d_trk_finder->Reset();

	TrkrHitSetContainer::ConstRange iter_range = hitsetcon->GetHitSets();
	for ( TrkrHitSetContainer::ConstIterator iter = iter_range.first; iter!=iter_range.second; ++iter){

		int ichip = int(mvtxdefutil.GetChipId(iter->first));

		MvtxHitSetv1 *hitset = static_cast<MvtxHitSetv1 *>(iter->second);

		MvtxHitSetv1::ConstRange hit_iter_range = hitset->GetHits();
		for ( MvtxHitSetv1::ConstIterator hit_iter = hit_iter_range.first; hit_iter!=hit_iter_range.second; ++hit_iter){

			int icol = int(hit_iter->first);
			int irow = int(hit_iter->second);

			h2d_hit[ichip]->Fill(icol,irow);

			h1d_trk_finder_z->Fill(icol-beam_x[ichip]);
			h1d_trk_finder_x->Fill(irow-beam_y[ichip]);
			h2d_trk_finder->Fill(icol-beam_x[ichip],irow-beam_y[ichip]);

			nhit_per_chip[ichip]++;

		}//hit_iter
	}//iter

	TrkrClusterContainer::ConstRange clus_iter_range = cluscon->GetClusters();
	for (TrkrClusterContainer::ConstIterator clus_iter = clus_iter_range.first; clus_iter!=clus_iter_range.second; ++clus_iter){

		int ichip  = trkrdefutil.GetLayer(clus_iter->first);

		TrkrCluster *clus = clus_iter->second;

		h2d_clus[ichip]->Fill(clus->GetZ(),clus->GetX());
		h2d_clus_beam[ichip]->Fill(clus->GetZ() - beam_x[ichip]*(28e-4),clus->GetX() - beam_y[ichip]*(28e-4));

		h1d_clus_size_z[ichip]->Fill(clus->GetZSize());
		h1d_clus_size_x[ichip]->Fill(clus->GetPhiSize());

		nclus_per_chip[ichip]++;
	}//

	int nchip_w_clus = 0;

	for (int ichip=0; ichip<4; ichip++){
		h1d_hit_per_evt[ichip]->Fill(nhit_per_chip[ichip]);
		h1d_clus_per_evt[ichip]->Fill(nclus_per_chip[ichip]);

		if ( nclus_per_chip[ichip]>0 ) nchip_w_clus++;
	}

	//if ( nchip_w_clus<3 ) return 0;
	//if ( 1 ) return 0;

	//fast track reconstruction
	/*
	double trk_z = h1d_trk_finder_z->GetBinCenter(h1d_trk_finder_z->GetMaximumBin());
	double trk_x = h1d_trk_finder_x->GetBinCenter(h1d_trk_finder_x->GetMaximumBin());
	double trk_z_w = h1d_trk_finder_z->GetBinContent(h1d_trk_finder_z->GetMaximumBin());
	double trk_x_w = h1d_trk_finder_x->GetBinContent(h1d_trk_finder_x->GetMaximumBin());
	*/

	/*
	double trk_z = 0;
	double trk_x = 0;
	double trk_w = 0.1;
	for (int ix=0; ix<h2d_trk_finder->GetNbinsX(); ix++){
		for (int iy=0; iy<h2d_trk_finder->GetNbinsY(); iy++){
			if ( h2d_trk_finder->GetBinContent(ix+1,iy+1)<trk_w ) continue; 

			trk_z = h1d_trk_finder_z->GetBinCenter(ix+1);
			trk_x = h1d_trk_finder_x->GetBinCenter(iy+1);
			trk_w = h2d_trk_finder->GetBinContent(ix+1,iy+1);

		}
	}

	if ( trk_w<3 ) return 0;

	int nclus_associated = 0;
	int nclus_associated_per_chip[4] = {0};

	//TrkrClusterContainer::ConstRange clus_iter_range = cluscon->GetClusters();
	for (TrkrClusterContainer::ConstIterator clus_iter = clus_iter_range.first; clus_iter!=clus_iter_range.second; ++clus_iter){

		int ichip  = trkrdefutil.GetLayer(clus_iter->first);

		TrkrCluster *clus = clus_iter->second;

		double clusz = clus->GetZ()/(28e-4) - beam_x[ichip];
		double clusx = clus->GetX()/(28e-4) - beam_y[ichip];

		if ( fabs(clusz-trk_z)<5 && fabs(clusx-trk_x)<5 ){
			nclus_associated_per_chip[ichip]++;
			nclus_associated++;
		}

	}//

	h1d_clus_associated->Fill(nclus_associated);

	if ( nclus_associated<3 ) return 0;

	for (int ichip=0; ichip<4; ichip++){
		if ( nclus_associated_per_chip[ichip]>0 ) h1d_clus_eff->Fill(ichip);
	}

	//cout << "N Associated Clus: " << nclus_associated << endl;


	//cout << "Z: " << trk_z << ", X: " << trk_x << ", S: " << trk_w << endl; 

	//TrkrHitSetContainer::ConstRange iter_range = hitsetcon->GetHitSets();
	for ( TrkrHitSetContainer::ConstIterator iter = iter_range.first; iter!=iter_range.second; ++iter){

		int ichip = int(mvtxdefutil.GetChipId(iter->first));

		MvtxHitSetv1 *hitset = static_cast<MvtxHitSetv1 *>(iter->second);

		MvtxHitSetv1::ConstRange hit_iter_range = hitset->GetHits();
		for ( MvtxHitSetv1::ConstIterator hit_iter = hit_iter_range.first; hit_iter!=hit_iter_range.second; ++hit_iter){

			int icol = int(hit_iter->first);
			int irow = int(hit_iter->second);

			h2d_hit_trk[ichip]->Fill(icol-beam_x[ichip]-trk_z,irow-beam_y[ichip]-trk_x);

		}//hit_iter
	}//iter
	*/


  return 0;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int MvtxQAHisto::End(PHCompositeNode *topNode)
{

	cout << "-----MvtxQAHisto::End------" << endl;

  PHTFileServer::get().cd( _outfile );
	PHTFileServer::get().write( _outfile );
  //PHTFileServer::get().close();

  return 0;
}


//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void MvtxQAHisto::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
	//
	//
	cluscon = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
	if (cluscon==0)
	{
		cout << "MvtxQAHisto::Process_Event - TrkrClusterContainer not found" << endl;
		exit(-1);
	}

	hitsetcon = findNode::getClass<TrkrHitSetContainer>(topNode, "TrkrHitSetContainer");
	if (hitsetcon==0)
	{
		cout << "MvtxQAHisto::Process_Event - TrkrHitSetContainer not found" << endl;
		exit(-1);
	}

}




#include "MvtxUnpackPRDF.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>
#include <Event/packet.h>
#include <Event/packet_hbd_fpgashort.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

#include <trackbase/TrkrDefUtil.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetv1.h>
#include <mvtx/MvtxDefUtil.h>
#include <mvtx/MvtxHitSetv1.h>

#include <iostream>
#include <string>
#include <cassert>

using namespace std;

//____________________________________
MvtxUnpackPRDF::MvtxUnpackPRDF() :
    SubsysReco("MvtxUnpackPRDF"),
		/*PHCompositeNode **/ dstNode(NULL),
		/*Event**/_event(NULL),
		/*Packet_hbd_fpgashort**/_packet(NULL),
		/*int*/_nevents(0),
		/*int*/_verbosity(0),
		/*bool*/_first(true)
{
	for (int ichip=0; ichip<NCHIP; ichip++){
		_nevent_per_chip[ichip] = 0;
		_npixel_per_chip[ichip] = 0;
	}
}

//____________________________________
int
MvtxUnpackPRDF::Init(PHCompositeNode *topNode)
{

	cout << "-----MvtxUnpackPRDF::Init-----" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int
MvtxUnpackPRDF::InitRun(PHCompositeNode *topNode)
{

	CreateNodeTree(topNode);

	cout << "-----MvtxUnpackPRDF::Init-----" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int
MvtxUnpackPRDF::process_event(PHCompositeNode *topNode)
{
	_nevents++;
	_event = findNode::getClass<Event>(topNode, "PRDF");
	if (_event==0)
	{
		cout << "MvtxUnpackPRDF::Process_Event - Event not found" << endl;
		return -1;
	}
	if (_event->getEvtType()!=1)
	{
		cout << "MvtxUnpackPRDF::Process_Event - non-data event type " << _event->getEvtType() << endl;
		return -1;
	}

	_hitsetcon = findNode::getClass<TrkrHitSetContainer>(topNode, "TrkrHitSetContainer");
	if (_hitsetcon==0)
	{
		cout << "MvtxUnpackPRDF::Process_Event - TrkrHitSetContainer not found" << endl;
		return -1;
	}
	else{

		MvtxDefUtil *mvtxdef = new MvtxDefUtil();

		for (int ichip=0; ichip<NCHIP; ichip++){
			MvtxHitSetv1 *hitset = new MvtxHitSetv1();
			hitset->SetHitSetKey(mvtxdef->GenHitSetKey(char(ichip),uint8_t(0),uint8_t(ichip)));
			_hitsetcon->AddHitSetSpecifyKey(hitset->GetHitSetKey(),hitset);
			if ( _verbosity>10 )
			{
				cout << "-----MvtxUnpackPRDF::CreateNodeTree, Create hitset " << hitset->GetHitSetKey() << endl;
			}
		}

		if ( _verbosity>10 )
		{
			cout << "-----MvtxUnpackPRDF::process_event Create HitSet" << endl;
			_hitsetcon->identify();
		}
		delete mvtxdef;
	}

	if ( _verbosity>90 )
	{
		cout << "EVENT: " << _nevents << endl;
	}

	MakeHits();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________
void
MvtxUnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
{

  PHNodeIterator nodeItr(topNode);
  //DST node
  dstNode = static_cast<PHCompositeNode *>(nodeItr.findFirst("PHCompositeNode",
                                                              "DST"));
  if (!dstNode)
  {
    cout << "PHComposite node created: DST" << endl;
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

	PHCompositeNode* trkrNode 
		= dynamic_cast<PHCompositeNode*>(nodeItr.findFirst("PHCompositeNode","TRKR"));
	if (!trkrNode) {
		trkrNode = new PHCompositeNode("TRKR");
		dstNode->addNode(trkrNode);
	}

	TrkrHitSetContainer *hitsetcon = findNode::getClass<TrkrHitSetContainer>(dstNode,"TrkrHitSetContainer");
	if (!hitsetcon) {
		hitsetcon = new TrkrHitSetContainer();
		PHIODataNode<PHObject> *hitsetconNode =
			new PHIODataNode<PHObject>(hitsetcon,"TrkrHitSetContainer","PHObject");
		trkrNode->addNode(hitsetconNode);
	}//

  //DATA nodes
	/*
  data_node = static_cast<PHCompositeNode *>(nodeItr.findFirst("PHCompositeNode",
                                                               "RAW_DATA"));
  if (!data_node)
  {
		cout << "PHComposite node created: RAW_DATA" << endl;
    data_node = new PHCompositeNode("RAW_DATA");
    dst_node->addNode(data_node);
  }
	*/

}

//_______________________________________
void
MvtxUnpackPRDF::MakeHits()
{

	int nhits_for_this = 0;

	Packet *p = _event->getPacket(2000);
	if (p) 
	{   


		int chipmax = p->iValue(0, "HIGHEST_CHIP") + 1;
		if ( _verbosity>100 )
		{
			cout << "CHIPMAX: " << chipmax << endl;
		}

		if ( chipmax>NCHIP )
		{
			cout << "WARNING!! Event: " << _nevents << " More chips than expected!" 
				<< " NCHIP:" << NCHIP << " HIGHEST_CHIP:" << chipmax << endl;

			chipmax = NCHIP;
		}

		int excess_data_bytes = p->iValue(0, "EXCESS_DATA_BYTES");
		if ( excess_data_bytes>0 )
		{
			cout << "WARNING!! Event: " << _nevents << " Data found past chip trailer"
				<< " EXCESS_DATA_BYTES: " << excess_data_bytes << endl;
		}

		//float mrow_chip0 = -1;
		//float mcol_chip0 = -1;

		for (int ichip=0; ichip<NCHIP; ichip++)
		{
			int header_found = p->iValue(ichip, "HEADER_FOUND");
			int trailer_found = p->iValue(ichip, "TRAILER_FOUND");
			int bunchcounter = p->iValue(ichip, "BUNCHCOUNTER");
			//int unexpected_bytes = p->iValue(ichip, "UNEXPECTED_BYTES");
			//int readout_flags = p->iValue(ichip, "READOUT_FLAGS");

			if ( header_found==0 || trailer_found==0 )
			{
				if ( _verbosity>100 )
				{
					cout << "WARNING!! Event: " << _nevents << " Missing chip " << ichip << " HEADER_FOUND: " << header_found << " TRAILER_FOUND: " << trailer_found << " BUNCHCOUNTER: " << bunchcounter << endl;
				}
				//continue;
			}
		}//ichip

		MvtxDefUtil *mvtxdef = new MvtxDefUtil();

		for (int ichip=0; ichip<chipmax; ichip++)
		{

			int header_found = p->iValue(ichip, "HEADER_FOUND");
			int trailer_found = p->iValue(ichip, "TRAILER_FOUND");
			if (header_found==0 || trailer_found==0) continue;
			TrkrDefs::hitsetkey mvtx_hitsetkey = mvtxdef->GenHitSetKey(char(ichip),uint8_t(0),uint8_t(ichip));

			for (int irow=0; irow<NROW; irow++)
			{
				int regions = p->iValue(ichip, irow);
				if ( regions==0 ) continue;

				if ( _verbosity>100 )
					cout << "ICHIP: " << ichip << ",IROW: " << irow << ", REGIONS: " << regions << endl;

				for (int iregion=0; iregion<NREGION; iregion++)
				{
					if ( (regions>>iregion)&1 )
					{
						if ( _verbosity>100 )
							cout << "REGIONS: " << regions << ", I: " << iregion << endl;

						int bits = p->iValue(ichip, iregion, irow);
						if ( bits==0 ) continue;

						for (int icol=0; icol<NCOL_PER_REGION; icol++)
						{
							if ( (bits>>icol)&1 )
							{
								_npixel_per_chip[ichip]++;
								int row = irow;
								int col = iregion*NCOL_PER_REGION + icol;

								MvtxHitSetv1 *hitset = static_cast<MvtxHitSetv1 *>(_hitsetcon->FindOrAddHitSet(mvtx_hitsetkey)->second);
								if ( hitset )
								{
									//cout << "ADD HIT, COL: " << col << ", ROW: " << row << endl;
									hitset->AddHit(uint16_t(col),uint16_t(row));
									nhits_for_this++;
								}
								else
								{
									cout << "NO HITSET!" << endl;
								}

							}
						}//icol

					}
				}//iregion
			}//irow
	  
			_nevent_per_chip[ichip]++;

		}//ichip

		if ( _verbosity>10 )
		{
			cout << "-----MvtxUnpackPRDF::process_event Check HitSetCon" << endl;
			_hitsetcon->identify();

			int nhits = 0;
			for (int ichip=0; ichip<NCHIP; ichip++){
				TrkrDefs::hitsetkey mvtx_hitsetkey = mvtxdef->GenHitSetKey(char(ichip),uint8_t(0),uint8_t(ichip));
				MvtxHitSetv1 *hitset = static_cast<MvtxHitSetv1 *>(_hitsetcon->FindHitSet(mvtx_hitsetkey));
				MvtxHitSetv1::ConstRange itr_range = hitset->GetHits();

				for (MvtxHitSetv1::ConstIterator itr = itr_range.first; itr != itr_range.second; ++itr){
					nhits++;
				}
			}
			cout << "N HITS: " << nhits_for_this << ", " << nhits << endl;
		}

		delete mvtxdef;
		delete p;
	}//p

}

//___________________________________
int
MvtxUnpackPRDF::End(PHCompositeNode *topNode)
{
	cout << "-----MvtxUnpackPRDF::End-----" << endl;
	if ( _verbosity>10 )
	{

		MvtxDefUtil *mvtxdef = new MvtxDefUtil();

		cout << "-----MvtxUnpackPRDF::End::PrintSummary-----" << endl;
		for (int ichip=0; ichip<NCHIP; ichip++){

			cout << "CHIP, " << ichip << ", " 
				<< ", N EVENT: " << _nevent_per_chip[ichip] 
				<< ", N HITS: " << _npixel_per_chip[ichip] 
				<< endl; 
		}

		delete mvtxdef;
	}

  return Fun4AllReturnCodes::EVENT_OK;
}


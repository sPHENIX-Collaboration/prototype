#include "MvtxPrototype2UnpackPRDF.h"

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

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
//#include <trackbase/TrkrHitSetv1.h>
#include <mvtx/MvtxDefs.h>
#include <mvtx/MvtxHit.h>
//#include <mvtx/MvtxHitSetv1.h>

#include <iostream>
#include <string>
#include <cassert>

using namespace std;

//____________________________________
MvtxPrototype2UnpackPRDF::MvtxPrototype2UnpackPRDF() :
    SubsysReco("MvtxPrototype2UnpackPRDF"),
		/*PHCompositeNode **/ dstNode(NULL),
		/*Event**/_event(NULL),
		/*Packet_hbd_fpgashort**/_packet(NULL),
		/*int*/_nevents(0),
		/*bool*/_first(true)
{
	for (int istave=0; istave<NSTAVE; istave++){
		for (int ichip=0; ichip<NCHIP; ichip++){
			_nevent_per_chip[istave][ichip] = 0;
			_npixel_per_chip[istave][ichip] = 0;
		}
	}
}

//____________________________________
int
MvtxPrototype2UnpackPRDF::Init(PHCompositeNode *topNode)
{

	cout << "-----MvtxPrototype2UnpackPRDF::Init-----" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int
MvtxPrototype2UnpackPRDF::InitRun(PHCompositeNode *topNode)
{

	CreateNodeTree(topNode);

	cout << "-----MvtxPrototype2UnpackPRDF::Init-----" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int
MvtxPrototype2UnpackPRDF::process_event(PHCompositeNode *topNode)
{
	_nevents++;
	_event = findNode::getClass<Event>(topNode, "PRDF");
	if (_event==0)
	{
		cout << "MvtxPrototype2UnpackPRDF::Process_Event - Event not found" << endl;
		return -1;
	}
	if (_event->getEvtType()!=1)
	{
		cout << "MvtxPrototype2UnpackPRDF::Process_Event - non-data event type " << _event->getEvtType() << endl;
		return -1;
	}

	_hitsetcon = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
	if (_hitsetcon==0)
	{
		cout << "MvtxPrototype2UnpackPRDF::Process_Event - TRKR_HITSET not found" << endl;
		return -1;
	}

	/*
	else{

		MvtxDefUtil *mvtxdef = new MvtxDefUtil();

		for (int ichip=0; ichip<NCHIP; ichip++){
			MvtxHitSetv1 *hitset = new MvtxHitSetv1();
			hitset->SetHitSetKey(mvtxdef->GenHitSetKey(char(ichip),uint8_t(0),uint8_t(ichip)));
			_hitsetcon->AddHitSetSpecifyKey(hitset->GetHitSetKey(),hitset);
			if ( _verbosity>10 )
			{
				cout << "-----MvtxPrototype2UnpackPRDF::CreateNodeTree, Create hitset " << hitset->GetHitSetKey() << endl;
			}
		}

		if ( _verbosity>10 )
		{
			cout << "-----MvtxPrototype2UnpackPRDF::process_event Create HitSet" << endl;
			_hitsetcon->identify();
		}
		delete mvtxdef;
	}
	*/

	if ( Verbosity()>90 )
	{
		cout << "EVENT: " << _nevents << endl;
	}

	MakeHits();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________
void
MvtxPrototype2UnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
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

	TrkrHitSetContainer *hitsetcon = findNode::getClass<TrkrHitSetContainer>(dstNode,"TRKR_HITSET");
	if (!hitsetcon) {
		hitsetcon = new TrkrHitSetContainer();
		PHIODataNode<PHObject> *hitsetconNode =
			new PHIODataNode<PHObject>(hitsetcon,"TRKR_HITSET","PHObject");
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
MvtxPrototype2UnpackPRDF::MakeHits()
{

	//int nhits_for_this = 0;

	Packet *p = _event->getPacket(2000);
	if (p) 
	{   

		for (int iru=1; iru<NMAXRU+1; iru++)
		{

			if ( p->iValue(iru)==-1 )
			{
				cout << "No such RU: " << iru << endl;
				continue;
			}
			if ( p->iValue(iru, "BAD_RUIDS") )
			{
				cout << "Bad RU: " << iru << endl;
				continue;
			}

			for (int ich=1; ich<NMAXRUCHN+1; ich++)
			{

				int chip_id = p->iValue(iru, ich, "CHIP_ID"); 

				//skip masked or dead chip
				if ( chip_id==-1 )
				{
					continue;
				}

				int stave_id_from_map = GetStaveID(iru, ich);
				int chip_id_from_map = GetChipID(iru, ich);

				//skip unused channel 
				if ( stave_id_from_map==-1 || chip_id_from_map==-1 )
				{
					continue;
				}

				//skip inconsistent chip id
				if ( chip_id!=chip_id_from_map )
				{
					cout << "WARNING!! Event: " << _nevents << " Inconsistent chip ID"
						<< " from packet: " << chip_id
						<< " from mapping: " << chip_id_from_map << endl;
				}

				if ( Verbosity()>100 )
					cout << "iru: " << iru << " ich: " << ich << " chip_id from packet: " << chip_id << " chip_id from map: " << GetChipID(iru, ich) << endl;

				int excess_data_bytes = p->iValue(iru, ich, "EXCESS_DATA_BYTES");
				if ( excess_data_bytes>0 )
				{
					cout << "WARNING!! Event: " << _nevents << " Data found past chip trailer"
						<< " EXCESS_DATA_BYTES: " << excess_data_bytes << endl;
				}

				int bad_bytes = p->iValue(iru, ich, "BAD_BYTES");
				if ( bad_bytes>0 )
				{
					cout << "WARNING!! Event: " << _nevents << " Bad data found"
						<< " BAD_BYTES: " << bad_bytes << endl;
				}

				int header_found = p->iValue(iru, ich, "HEADER_FOUND");
				int trailer_found = p->iValue(iru, ich, "TRAILER_FOUND");
				int bunchcounter = p->iValue(iru, ich, "BUNCHCOUNTER");

				if ( header_found==0 || trailer_found==0 )
				{
					if ( Verbosity()>50 )
					{
						cout << "WARNING!! Event: " << _nevents << " Missing RU " << iru << " CH: " << ich << " HEADER_FOUND: " << header_found << " TRAILER_FOUND: " << trailer_found << " BUNCHCOUNTER: " << bunchcounter << endl;
					}
					//continue;
				}

				int nhits_per_ch = p->iValue(iru, ich);

				for (int ihit=0; ihit<nhits_per_ch; ihit++){

					int encoded_hit = p->iValue(iru, ich, ihit);

					int row_hit = DecodeRow(encoded_hit);
					int col_hit = DecodeCol(encoded_hit);

					if ( row_hit>=NROW || col_hit>=NCOL )
					{
						cout << "WARNING!! Event: " << _nevents << " Hit out of window"
							<< " RU: " << iru
							<< " CH: " << ich
							<< " ROW: " << row_hit
							<< " COL: " << col_hit
							<< endl;
						continue;
					}

					_npixel_per_chip[stave_id_from_map][chip_id_from_map]++;

					TrkrDefs::hitsetkey hitsetkey = MvtxDefs::genHitSetKey(0, stave_id_from_map, chip_id_from_map);
					TrkrHitSetContainer::Iterator hitsetit = _hitsetcon->findOrAddHitSet(hitsetkey); 

					// generate the key for this hit
					TrkrDefs::hitkey hitkey = MvtxDefs::genHitKey(uint16_t(col_hit),uint16_t(row_hit));
					// See if this hit already exists
					TrkrHit *hit = nullptr;
					hit = hitsetit->second->getHit(hitkey);
					if ( !hit )
					{
						hit = new MvtxHit();
						hitsetit->second->addHitSpecificKey(hitkey, hit);
					}

				}//ihit

				_nevent_per_chip[stave_id_from_map][chip_id_from_map]++;
			}//ich
		}//iru

		if ( Verbosity()>10 )
		{
			cout << "-----MvtxPrototype2UnpackPRDF::process_event Check HitSetCon" << endl;
			//_hitsetcon->identify();

			for (int istave=0; istave<NSTAVE; istave++){
				for (int ichip=0; ichip<NCHIP; ichip++){
					TrkrDefs::hitsetkey mvtx_hitsetkey = MvtxDefs::genHitSetKey(0,uint8_t(istave),uint8_t(ichip));
					TrkrHitSet *hitset = static_cast<TrkrHitSet*>(_hitsetcon->findHitSet(mvtx_hitsetkey));

					if ( hitset ){
						int nhits = 0;
						TrkrHitSet::ConstRange itr_range = hitset->getHits();
						for (TrkrHitSet::ConstIterator itr = itr_range.first; itr != itr_range.second; ++itr){
							nhits++;
						}
						cout << "Stave: " << istave << " Chip: " << ichip << " Nhits: " << nhits << endl;
					}
					/*
					*/
				}
			}
		}

		/*
		if ( _verbosity>10 )
		{
			cout << "-----MvtxPrototype2UnpackPRDF::process_event Check HitSetCon" << endl;
			_hitsetcon->identify();

			int nhits = 0;
		}
		*/

		//delete mvtxdef;
		delete p;
	}//p

}

//___________________________________
int
MvtxPrototype2UnpackPRDF::GetStaveID(int ruid, int chid)
{

	if ( ruid<=0 || ruid>NMAXRU || chid<=0 || chid>NMAXRUCHN ){
		return -1;
	}

	/*
	const int map_stave[NMAXRU][NMAXRUCHN] = {
		{
			0, 0, 0, 0, 0, 0, 0,
			0, 2, 2, 2, 2, 2, 2,
			2, 2, 2, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 0, -1
		},
		{
			3, 3, 3, 3, 3, 3, 3,
			3, 3, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1
		}
	};
	*/

	//order
	//A105
	//C104
	//C105
	//E103

	const int map_stave[NMAXRU][NMAXRUCHN] = {
		{
			3, 3, 3, 3, 3, 3, 3,
			3, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 3, -1
		},
		{
			2, 2, 2, 2, 2, 2, 2,
			2, 2, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1
		}
	};

	return map_stave[ruid-1][chid-1]; 

}

//___________________________________
int
MvtxPrototype2UnpackPRDF::GetChipID(int ruid, int chid)
{

	if ( ruid<=0 || ruid>NMAXRU || chid<=0 || chid>NMAXRUCHN ){
		return -1;
	}

	const int map_chip[NMAXRU][NMAXRUCHN] = {
		{
			0, 1, 2, 3, 4, 5, 6,
			7, 8, 7, 6, 5, 4, 3,
			2, 1, 0, 8, 7, 6, 5,
			4, 3, 2, 1, 0, 8, -1
		},
		{
			2, 1, 0, 3, 4, 5, 6,
			7, 8, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1
		}
	};

	return map_chip[ruid-1][chid-1];

}

//___________________________________
int
MvtxPrototype2UnpackPRDF::DecodeRow(int val) const
{
	return (val >> 16);
}

//___________________________________
int
MvtxPrototype2UnpackPRDF::DecodeCol(int val) const
{
	return (val & 0xffff);
}

//___________________________________
int
MvtxPrototype2UnpackPRDF::End(PHCompositeNode *topNode)
{
	cout << "-----MvtxPrototype2UnpackPRDF::End-----" << endl;

	if ( Verbosity()>1 )
	{
		cout << "-----MvtxPrototype2UnpackPRDF::End::PrintSummary-----" << endl;
		for (int istave=0; istave<NSTAVE; istave++){
			for (int ichip=0; ichip<NCHIP; ichip++){

				cout << "STAVE: " << istave << ", CHIP: " << ichip << ", " 
					<< ", N EVENT: " << _nevent_per_chip[istave][ichip] 
					<< ", N HITS: " << _npixel_per_chip[istave][ichip] 
					<< endl; 
			}
		}

	}

  return Fun4AllReturnCodes::EVENT_OK;
}


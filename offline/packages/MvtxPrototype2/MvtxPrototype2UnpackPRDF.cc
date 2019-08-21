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

#include <cassert>
#include <iostream>
#include <string>

using namespace std;

map< pair< int, int>, pair< int, int> > MvtxPrototype2UnpackPRDF::s_map_chips =
{
  {{1,1}, {0,0}},
  {{1,2}, {0,1}},
  {{1,3}, {0,2}},
  {{1,4}, {0,3}},
  {{1,5}, {0,4}},
  {{1,6}, {0,5}},
  {{1,7}, {0,6}},
  {{1,8}, {0,7}},
  {{1,9}, {3,8}},
  {{1,10}, {3,7}},
  {{1,11}, {3,6}},
  {{1,12}, {3,5}},
  {{1,13}, {3,4}},
  {{1,14}, {3,3}},
  {{1,15}, {3,2}},
  {{1,16}, {3,1}},
  {{1,17}, {3,0}},
  {{1,18}, {2,8}},
  {{1,19}, {2,7}},
  {{1,20}, {2,6}},
  {{1,21}, {2,5}},
  {{1,22}, {2,4}},
  {{1,23}, {2,3}},
  {{1,24}, {2,2}},
  {{1,25}, {2,1}},
  {{1,26}, {2,0}},
  {{1,27}, {0,8}},
  {{2,1}, {1,2}},
  {{2,2}, {1,1}},
  {{2,3}, {1,0}},
  {{2,4}, {1,3}},
  {{2,5}, {1,4}},
  {{2,6}, {1,5}},
  {{2,7}, {1,6}},
  {{2,8}, {1,7}},
  {{2,9}, {1,8}}
  //order
  //E103
  //C105
  //C104
  //A105
}; //<ruid, ruch> to <stave, chipID>

/**
 * Layers order in geom
 * A105 -> layer 0, stave 0
 * C104 -> layer 1, stave 0
 * C105 -> layer 2, stave 0
 * E103 -> layer 3, stave 0
 * <{stave_from_map_chips, layers_in_geom}
 */
map<int, int> MvtxPrototype2UnpackPRDF::s_map_layers = { {0,0}, {1,1}, {2,2}, {3,3} };

//____________________________________
MvtxPrototype2UnpackPRDF::MvtxPrototype2UnpackPRDF() :
    SubsysReco("MvtxPrototype2UnpackPRDF"),
    /*PHCompositeNode **/ dstNode(NULL),
    /*Event**/_event(NULL),
    /*Packet_hbd_fpgashort**/_packet(NULL),
    /*int*/_nevents(0),
    /*bool*/_first(true)
{
  for (int istave=0; istave<NLAYER; istave++){
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
  _event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event==0)
  {
    cout << "MvtxPrototype2UnpackPRDF::Process_Event - Event not found" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if ( _event->getEvtType() == BEGRUNEVENT )
  {
    cout << __FILE__ << " - run event header found. Aborting "<< endl;
    _nevents = 0;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if ( _event->getEvtType() != DATAEVENT )
  {
    cout << "MvtxPrototype2UnpackPRDF::Process_Event - non-data event type "  \
         << _event->getEvtType() << ". Aborting."<< endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _nevents++;

  _hitsetcon = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if ( _hitsetcon == 0 )
  {
    cout << "MvtxPrototype2UnpackPRDF::Process_Event - TRKR_HITSET not found" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
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
  double significand = _nevents / pow(10, (int) (log10(_nevents)));
  if ( Verbosity() >= VERBOSITY_MORE && fmod(significand,1.0) == 0 && significand <=10  )
  {
    cout << "-----MvtxPrototype2UnpackPRDF::process_event Check HitSetCon" << endl;
    //_hitsetcon->identify();
    for (int istave=0; istave<NLAYER; istave++)
    {
      for (int ichip=0; ichip<NCHIP; ichip++)
      {
        TrkrDefs::hitsetkey mvtx_hitsetkey =
          MvtxDefs::genHitSetKey(uint8_t(istave), 0, uint8_t(ichip));
        TrkrHitSet *hitset =
          static_cast<TrkrHitSet*>(_hitsetcon->findHitSet(mvtx_hitsetkey));

        if ( hitset )
        {
          int nhits = hitset->size();
          cout << "Stave: " << istave
               << " Chip: " << ichip << " Nhits: " << nhits << endl;
        }
      }
    }
  }

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
    bool event_err = false;

    if ( !event_err )
    {
      for (int iru=0; iru<NMAXRU+1; iru++)
      {
/*
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
*/
        if ( p->iValue(iru) != -1 )
        {
          for ( int ich = 0; ich < NMAXRUCHN+1; ich++)
          {
            if ( p->iValue(iru, ich) > 0 )
            {
/*
              //int chip_id = p->iValue(iru, ich, "CHIP_ID");
              //skip masked or dead chip
              if ( chip_id==-1 )
              {
                continue;
              }
*/
              if ( s_map_chips.count({iru,ich}) != 1 )
              {
                cout << PHWHERE << "invalid: (iru " << iru << ", ich " \
                     << ich << ") " << endl;
              }
              else
              {
                std::pair<int,int> chip_pos = s_map_chips[{iru, ich}];
                int stave_id_from_map = chip_pos.first;
                int chip_id_from_map = chip_pos.second;
  
                //skip unused channel
                if ( stave_id_from_map==-1 || chip_id_from_map==-1 )
                {
                  continue;
                }
  /*
                //skip inconsistent chip id
                if ( chip_id!=chip_id_from_map )
                {
                  cout << "WARNING!! Event: " << _nevents << " Inconsistent chip ID"
                  << " from packet: " << chip_id
                  << " from mapping: " << chip_id_from_map << endl;
                }
    
                if ( Verbosity()>100 )
                  cout << "iru: " << iru << " ich: " << ich << " chip_id from packet: " \
                       << chip_id << " chip_id from map: " << chip_id_from_map << endl;
    
                int excess_data_bytes = p->iValue(iru, ich, "EXCESS_DATA_BYTES");
                if ( excess_data_bytes>0 )
                {
                  cout << "WARNING!! Event: " << _nevents << " Data found past chip trailer" \
                       << " EXCESS_DATA_BYTES: " << excess_data_bytes << endl;
                }
  
                int bad_bytes = p->iValue(iru, ich, "BAD_BYTES");
                if ( bad_bytes>0 )
                {
                  cout << "WARNING!! Event: " << _nevents << " Bad data found" \
                       << " BAD_BYTES: " << bad_bytes << endl;
                }
    
                int header_found = p->iValue(iru, ich, "HEADER_FOUND");
                int trailer_found = p->iValue(iru, ich, "TRAILER_FOUND");
                int bunchcounter = p->iValue(iru, ich, "BUNCHCOUNTER");
    
                if ( header_found==0 || trailer_found==0 )
                {
                  if ( Verbosity()>50 )
                  {
                    cout << "WARNING!! Event: " << _nevents << " Missing RU "          \
                         << iru << " CH: " << ich << " HEADER_FOUND: " << header_found \
                         << " TRAILER_FOUND: " << trailer_found << " BUNCHCOUNTER: "   \
                         << bunchcounter << endl;
                  }
                  //continue;
                }
  */
                int nhits_per_ch = p->iValue(iru, ich);
                for (int ihit = 0; ihit < nhits_per_ch; ihit++)
                {
                  int encoded_hit = p->iValue(iru, ich, ihit);
                  int row_hit = DecodeRow(encoded_hit);
                  int col_hit = DecodeCol(encoded_hit);
  
                  if ( row_hit>=NROW || col_hit>=NCOL )
                  {
                    cout << PHWHERE << "WARNING!! Event: " << _nevents \
                         << " Hit out of window"
                         << " RU: " << iru
                         << " CH: " << ich
                         << " ROW: " << row_hit
                         << " COL: " << col_hit
                         << endl;
                    continue;
                  }
  
                  _npixel_per_chip[stave_id_from_map][chip_id_from_map]++;
  
                  auto lay = s_map_layers[stave_id_from_map];
                  TrkrDefs::hitsetkey hitsetkey = \
                    MvtxDefs::genHitSetKey(lay, 0, chip_id_from_map);
                  TrkrHitSetContainer::Iterator hitsetit = \
                    _hitsetcon->findOrAddHitSet(hitsetkey);
  
                  // generate the key for this hit
                  TrkrDefs::hitkey hitkey = \
                    MvtxDefs::genHitKey(uint16_t(col_hit),uint16_t(row_hit));
                  // See if this hit already exists
                  TrkrHit *hit = hitsetit->second->getHit(hitkey);
                  if ( !hit )
                  {
                    hit = new MvtxHit();
                    hitsetit->second->addHitSpecificKey(hitkey, hit);
                  }
                  else
                  {
                    cout << PHWHERE << "WARNING!! Event: " << _nevents
                         << " duplicated Hit "
                         << " STAVE: " << stave_id_from_map
                         << " CHIP: "  << chip_id_from_map
                         << " ROW: "   << row_hit
                         << " COL: "   << col_hit
                         << endl;
                  }

                }//ihit

                _nevent_per_chip[stave_id_from_map][chip_id_from_map]++;
              } //if good ru,ch
            } //if ru,ch have data
          }//ich
        }//if iru has data
      }//iru
    }// !err

    delete p;
  }//p

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

  if ( Verbosity() > VERBOSITY_SOME )
  {
    cout << "-----MvtxPrototype2UnpackPRDF::End::PrintSummary-----" << endl;
    for (int istave=0; istave<NLAYER; istave++){
      for (int ichip=0; ichip<NCHIP; ichip++){

        cout << "LAYER: " << istave << ", CHIP: " << ichip << ", "
          << ", N EVENT: " << _nevent_per_chip[istave][ichip]
          << ", N HITS: " << _npixel_per_chip[istave][ichip]
          << endl;
      }
    }

  }

  return Fun4AllReturnCodes::EVENT_OK;
}


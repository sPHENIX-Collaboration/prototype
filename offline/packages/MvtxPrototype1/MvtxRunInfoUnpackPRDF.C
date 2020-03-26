#include "MvtxRunInfoUnpackPRDF.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packetConstants.h>
#include <Event/packet.h>
#include <Event/packet_hbd_fpgashort.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

#include <ffaobjects/EventHeaderv1.h>

#include <iostream>
#include <string>
#include <cassert>
#include <sstream>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

//____________________________________
MvtxRunInfoUnpackPRDF::MvtxRunInfoUnpackPRDF() :
    SubsysReco("MvtxRunInfoUnpackPRDF"),
		/*PHCompositeNode **/ dstNode(NULL),
		/*Event**/_event(NULL),
		/*Packet_hbd_fpgashort**/_packet(NULL)
{

}

//____________________________________
int
MvtxRunInfoUnpackPRDF::Init(PHCompositeNode *topNode)
{

	cout << "-----MvtxRunInfoUnpackPRDF::Init-----" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int
MvtxRunInfoUnpackPRDF::InitRun(PHCompositeNode *topNode)
{

	CreateNodeTree(topNode);

	cout << "-----MvtxRunInfoUnpackPRDF::Init-----" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int
MvtxRunInfoUnpackPRDF::process_event(PHCompositeNode *topNode)
{
	/*
	_event = findNode::getClass<Event>(topNode, "PRDF");
	if (_event==0)
	{
		cout << "MvtxRunInfoUnpackPRDF::Process_Event - Event not found" << endl;
		return -1;
	}
	*/

	stringstream ss_packet;
	string s_packet;
	char par_name[30], par_sval[30], par_unit[30];
	int par_ival; 

	Event *event = findNode::getClass<Event>(topNode, "PRDF");
	if ( event==NULL )
	{
		if ( Verbosity()>=VERBOSITY_SOME )
			cout << "RunInfoUnpackPRDF::Process_Event - Event not found" << endl;
		return Fun4AllReturnCodes::DISCARDEVENT;
	}

	EventHeaderv1 *eventheader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
	if ( eventheader )
	{
		eventheader->set_RunNumber(event->getRunNumber());
		eventheader->set_EvtSequence(event->getEvtSequence());
		eventheader->set_EvtType(event->getEvtType());
		eventheader->set_TimeStamp(event->getTime());
		if ( verbosity )
		{   
			eventheader->identify();
		}   
	}

	if ( event->getEvtType()!=BEGRUNEVENT )
		return Fun4AllReturnCodes::EVENT_OK;
	else
	{
		if ( verbosity>=VERBOSITY_SOME )
		{
			cout << "RunInfoUnpackPRDF::process_event - with BEGRUNEVENT events ";
			event->identify();
		}

		PHParameters Params("RunInfo");

		Packet *packet = event->getPacket(910);
		if ( !packet )
		{
			cout << "RunInfoUnpackPRDF::process_event - failed to locate packet 910" << endl;
		}else{
			ss_packet.clear();
			packet->dump(ss_packet);

			getline(ss_packet,s_packet);
			sscanf(s_packet.c_str(), "%s = %i %s", par_name, &par_ival, par_unit);
			Params.set_int_param(par_name,par_ival);

			getline(ss_packet,s_packet);
			sscanf(s_packet.c_str(), "%s = %s", par_name, par_sval);
			Params.set_string_param(par_name,par_sval);
		}//910

		packet = event->getPacket(920);
		if ( !packet )
		{
			cout << "RunInfoUnpackPRDF::process_event - failed to locate packet 910" << endl;
			//Params.set_double_param(name, NAN);
		}else{
			ss_packet.clear();
			packet->dump(ss_packet);

			getline(ss_packet,s_packet);
			sscanf(s_packet.c_str(), "%s = %i", par_name, &par_ival);
			Params.set_int_param(par_name,par_ival);

			getline(ss_packet,s_packet);
			sscanf(s_packet.c_str(), "%s = %s", par_name, par_sval);
			Params.set_string_param(par_name,par_sval);

			getline(ss_packet,s_packet);
			sscanf(s_packet.c_str(), "%s = %i %s", par_name, &par_ival, par_unit);
			Params.set_int_param(par_name,par_ival);

			while ( getline(ss_packet,s_packet) )
			{
				sscanf(s_packet.c_str(), "%s = %i", par_name, &par_ival);
				Params.set_int_param(par_name,par_ival);
			}
		}//920

		packet = event->getPacket(930);
		if ( !packet )
		{
			cout << "RunInfoUnpackPRDF::process_event - failed to locate packet 930" << endl;
			//Params.set_double_param(name, NAN);
		}else{
			ss_packet.clear();
			packet->dump(ss_packet);

			while ( getline(ss_packet,s_packet) )
			{
				sscanf(s_packet.c_str(), "%s = %i %s", par_name, &par_ival, par_unit);
				Params.set_int_param(par_name,par_ival);
			}
		}//930

		Params.SaveToNodeTree(topNode, "RUN_INFO");

		if ( verbosity>=VERBOSITY_SOME )
		{
			Params.Print();
		}
			
	}//


  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________
void
MvtxRunInfoUnpackPRDF::CreateNodeTree(PHCompositeNode *topNode)
{

  PHNodeIterator nodeItr(topNode);

  PHCompositeNode *run_node = static_cast<PHCompositeNode *>(nodeItr.findFirst("PHCompositeNode", "RUN"));
  if (!run_node)
  {
    cout << "PHComposite node created: RUN" << endl;
    run_node = new PHCompositeNode("RUN");
    topNode->addNode(run_node);
  }

  PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(run_node,"RUN_INFO");
	if (!nodeparams)
	{
		run_node->addNode(new PHIODataNode<PdbParameterMap>(new PdbParameterMap(),"RUN_INFO"));
	}

  //DST node
  PHCompositeNode *dst_node = static_cast<PHCompositeNode *>(nodeItr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    cout << "PHComposite node created: DST" << endl;
    dst_node = new PHCompositeNode("DST");
    topNode->addNode(dst_node);
  }

  EventHeaderv1 *eventheader = new EventHeaderv1();
  PHObjectNode_t *EventHeaderNode = new PHObjectNode_t(eventheader, "EventHeader", "PHObject");  // contain PHObject
  dst_node->addNode(EventHeaderNode);

}

//___________________________________
int
MvtxRunInfoUnpackPRDF::End(PHCompositeNode *topNode)
{
	cout << "-----MvtxRunInfoUnpackPRDF::End-----" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}


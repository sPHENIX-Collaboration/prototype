#include "EventInfoSummary.h"

#include "RawTower_Prototype3.h"

#include <calobase/RawTower.h>  // for RawTower
#include <calobase/RawTowerContainer.h>

#include <phparameter/PHParameters.h>

#include <pdbcalbase/PdbParameterMap.h>

#include <fun4all/Fun4AllBase.h>  // for Fun4AllB...
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIOData...
#include <phool/PHNodeIterator.h>  // for PHNodeIt...
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#pragma GCC diagnostic pop

#include <cassert>
#include <cmath>  // for NAN
#include <iostream>
#include <string>
#include <utility>  // for pair

// with the exception of std,  please no using namespaces
// it severely muddies the ability to figure out
// where things come from
using namespace std;

//____________________________________
EventInfoSummary::EventInfoSummary()
  : SubsysReco("EventInfoSummary")
  , eventinfo_node_name("EVENT_INFO")
{
}

//_____________________________________
int EventInfoSummary::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int EventInfoSummary::process_event(PHCompositeNode *topNode)
{
  Event *event = findNode::getClass<Event>(topNode, "PRDF");
  if (!event)
  {
    if (Verbosity() >= VERBOSITY_SOME)
      cout << "EventInfoSummary::Process_Event - Event not found" << endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  // search for run info
  if (event->getEvtType() != DATAEVENT)
    return Fun4AllReturnCodes::EVENT_OK;
  else  // DATAEVENT
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      cout << "EventInfoSummary::process_event - with DATAEVENT events ";
      event->identify();
    }

    map<int, Packet *> packet_list;

    PHParameters Params("EventInfo");

    // spill indicator
    {
      RawTowerContainer *TOWER_RAW_SPILL_WARBLER =
          findNode::getClass<RawTowerContainer>(topNode,
                                                "TOWER_RAW_SPILL_WARBLER");
      assert(TOWER_RAW_SPILL_WARBLER);

      RawTower_Prototype3 *raw_tower = dynamic_cast<RawTower_Prototype3 *>(
          TOWER_RAW_SPILL_WARBLER->getTower(0));
      assert(raw_tower);

      boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::variance>> acc;

      for (int i = 0; i < RawTower_Prototype3::NSAMPLES; i++)
      {
        acc(raw_tower->get_signal_samples(i));
      }

      const double warbler_rms = boost::accumulators::variance(acc);
      const bool is_in_spill = warbler_rms > (1000 * 1000);
      Params.set_double_param("beam_SPILL_WARBLER_RMS", warbler_rms);
      Params.set_double_param("beam_Is_In_Spill", is_in_spill);
      Params.set_int_param("beam_Is_In_Spill", is_in_spill);
    }

    // energy sums
    {
      RawTowerContainer *TOWER_CALIB_CEMC =
          findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
      assert(TOWER_CALIB_CEMC);

      RawTowerContainer *TOWER_CALIB_LG_HCALIN =
          findNode::getClass<RawTowerContainer>(topNode,
                                                "TOWER_CALIB_LG_HCALIN");

      RawTowerContainer *TOWER_CALIB_LG_HCALOUT =
          findNode::getClass<RawTowerContainer>(topNode,
                                                "TOWER_CALIB_LG_HCALOUT");

      // process inner HCAL
      if (TOWER_CALIB_CEMC)
      {
        double sum_energy_calib = 0;

        auto range = TOWER_CALIB_CEMC->getTowers();
        for (auto it = range.first; it != range.second; ++it)
        {
          RawTower *tower = it->second;
          assert(tower);

          const int col = tower->get_bineta();
          const int row = tower->get_binphi();

          if (col < 0 or col >= 8)
            continue;
          if (row < 0 or row >= 8)
            continue;

          const double energy_calib = tower->get_energy();
          sum_energy_calib += energy_calib;

        }  //       for (auto it = range.first; it != range.second; ++it)
        Params.set_double_param("CALIB_CEMC_Sum", sum_energy_calib);
      }  // process inner HCAL

      // process inner HCAL
      if (TOWER_CALIB_LG_HCALIN)
      {
        double sum_energy_calib = 0;

        auto range = TOWER_CALIB_LG_HCALIN->getTowers();
        for (auto it = range.first; it != range.second; ++it)
        {
          RawTower *tower = it->second;
          assert(tower);

          const int col = tower->get_bineta();
          const int row = tower->get_binphi();

          if (col < 0 or col >= 4)
            continue;
          if (row < 0 or row >= 4)
            continue;

          const double energy_calib = tower->get_energy();
          sum_energy_calib += energy_calib;

        }  //       for (auto it = range.first; it != range.second; ++it)
        Params.set_double_param("CALIB_LG_HCALIN_Sum", sum_energy_calib);
      }  // process inner HCAL

      // process outer HCAL
      if (TOWER_CALIB_LG_HCALOUT)
      {
        double sum_energy_calib = 0;

        auto range = TOWER_CALIB_LG_HCALOUT->getTowers();
        for (auto it = range.first; it != range.second; ++it)
        {
          RawTower *tower = it->second;
          assert(tower);

          const int col = tower->get_bineta();
          const int row = tower->get_binphi();

          if (col < 0 or col >= 4)
            continue;
          if (row < 0 or row >= 4)
            continue;

          const double energy_calib = tower->get_energy();
          sum_energy_calib += energy_calib;

        }  //       for (auto it = range.first; it != range.second; ++it)

        Params.set_double_param("CALIB_LG_HCALOUT_Sum", sum_energy_calib);
      }  // process outer HCAL
    }

    // generic packets
    for (typ_channel_map::const_iterator it = channel_map.begin();
         it != channel_map.end(); ++it)
    {
      const string &name = it->first;
      const channel_info &info = it->second;

      if (packet_list.find(info.packet_id) == packet_list.end())
      {
        packet_list[info.packet_id] = event->getPacket(info.packet_id);
      }

      Packet *packet = packet_list[info.packet_id];

      if (!packet)
      {
        //          if (Verbosity() >= VERBOSITY_SOME)
        cout << "EventInfoSummary::process_event - failed to locate packet "
             << info.packet_id << " from ";
        event->identify();

        Params.set_double_param(name, NAN);
        continue;
      }

      const int ivalue = packet->iValue(info.offset);

      const double dvalue = ivalue * info.calibration_const;

      if (Verbosity() >= VERBOSITY_SOME)
      {
        cout << "EventInfoSummary::process_event - " << name << " = " << dvalue
             << ", raw = " << ivalue << " @ packet " << info.packet_id
             << ", offset " << info.offset << endl;
      }

      Params.set_double_param(name, dvalue);
    }

    for (map<int, Packet *>::iterator it = packet_list.begin();
         it != packet_list.end(); ++it)
    {
      if (it->second)
        delete it->second;
    }

    Params.SaveToNodeTree(topNode, eventinfo_node_name);

    if (Verbosity() >= VERBOSITY_SOME)
      Params.Print();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________
void EventInfoSummary::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator nodeItr(topNode);

  PHNodeIterator iter(topNode);

  // DST node
  PHCompositeNode *dst_node =
      static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    cout << "PHComposite node created: DST" << endl;
    dst_node = new PHCompositeNode("DST");
    topNode->addNode(dst_node);
  }

  PdbParameterMap *nodeparams =
      findNode::getClass<PdbParameterMap>(dst_node, eventinfo_node_name);
  if (not nodeparams)
  {
    dst_node->addNode(new PHIODataNode<PdbParameterMap>(new PdbParameterMap(),
                                                        eventinfo_node_name));
  }
}

void EventInfoSummary::add_channel(
    const std::string &name,        //! name of the channel
    const int packet_id,            //! packet id
    const unsigned int offset,      //! offset in packet data
    const double calibration_const  //! conversion constant from integer to
                                    //! meaningful value
)
{
  channel_map.insert(
      make_pair(name, channel_info(packet_id, offset, calibration_const)));
}

/*
 * TpcPrototypeUnpacker.cc
 *
 *  Created on: Sep 19, 2018
 *      Author: jinhuang
 */

#include "TpcPrototypeUnpacker.h"

#include "TpcPrototypeDefs.h"

#include <g4tpc/PHG4TpcPadPlane.h>  // for PHG4TpcPadPlane

#include <tpc/TpcDefs.h>
#include <tpc/TpcHit.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phfield/PHFieldConfigv1.h>
#include <phfield/PHFieldConfigv2.h>
#include <phfield/PHFieldUtility.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>
//#include <Event/packetConstants.h>
#include <Event/oncsSubConstants.h>

#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <boost/bimap.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <tuple>

using namespace std;
using namespace TpcPrototypeDefs::FEEv2;

TpcPrototypeUnpacker::TpcPrototypeUnpacker(const std::string& outputfilename)
  : SubsysReco("TpcPrototypeUnpacker")
  , padplane(nullptr)
  , tpcCylinderCellGeom(nullptr)
  , trkrclusters(nullptr)
  , m_outputFileName(outputfilename)
  , m_eventT(nullptr)
  , m_peventHeader(&m_eventHeader)
  , m_nClusters(-1)
  , m_IOClusters(nullptr)
  , m_chanT(nullptr)
  , m_pchanHeader(&m_chanHeader)
  , m_chanData(kSAMPLE_LENGTH, 0)
  , enableClustering(true)
  , m_clusteringZeroSuppression(50)
  , m_nPreSample(5)
  , m_nPostSample(5)
  , m_XRayLocationX(-1)
  , m_XRayLocationY(-1)
  , m_pdfMaker(nullptr)
{
}

TpcPrototypeUnpacker::~TpcPrototypeUnpacker()
{
  if (m_IOClusters)
  {
    m_IOClusters->Clear();
    delete m_IOClusters;
  }
  if (m_pdfMaker)
  {
    delete m_pdfMaker;
  }
  if (padplane)
  {
    delete padplane;
  }
}

int TpcPrototypeUnpacker::ResetEvent(PHCompositeNode* topNode)
{
  m_eventHeader = EventHeader();
  m_padPlaneData.Reset();
  m_clusters.clear();
  m_chanHeader = ChannelHeader();

  m_nClusters = -1;
  assert(m_IOClusters);
  m_IOClusters->Clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPrototypeUnpacker::Init(PHCompositeNode* topNode)
{
  if (padplane)
    padplane->Init(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPrototypeUnpacker::InitRun(PHCompositeNode* topNode)
{
  if (padplane)
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      cout << "TpcPrototypeUnpacker::InitRun - making pad plane"
           << endl;
    }

    //setup the constant field
    const int field_ret = InitField(topNode);
    if (field_ret != Fun4AllReturnCodes::EVENT_OK)
    {
      cout << "TpcPrototypeUnpacker::InitRun- Error - Failed field init with status = " << field_ret << endl;
      return field_ret;
    }

    string seggeonodename = "CYLINDERCELLGEOM_SVTX";  // + detector;
    tpcCylinderCellGeom = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
    if (!tpcCylinderCellGeom)
    {
      tpcCylinderCellGeom = new PHG4CylinderCellGeomContainer();
      PHNodeIterator iter(topNode);
      PHCompositeNode* runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(tpcCylinderCellGeom, seggeonodename.c_str(), "PHObject");
      runNode->addNode(newNode);
    }

    padplane->InitRun(topNode);
    padplane->CreateReadoutGeometry(topNode, tpcCylinderCellGeom);
  }

  // Create the Cluster node if required
  trkrclusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator iter(topNode);

    // Looking for the DST node
    PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      cout << PHWHERE << "DST Node missing, doing nothing." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    PHNodeIterator dstiter(dstNode);
    PHCompositeNode* DetNode =
        dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    trkrclusters = new TrkrClusterContainer();
    PHIODataNode<PHObject>* TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "TpcPrototypeUnpacker::InitRun - Making PHTFileServer " << m_outputFileName
         << endl;

    m_pdfMaker = new SampleFit_PowerLawDoubleExp_PDFMaker();
  }
  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);

  TH1D* h = new TH1D("hNormalization",  //
                     "Normalization;Items;Summed quantity", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event count");
  h->GetXaxis()->SetBinLabel(i++, "Collision count");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit Edep");
  h->GetXaxis()->SetBinLabel(i++, "TPC Pad Hit");
  h->GetXaxis()->SetBinLabel(i++, "TPC Charge e");
  h->GetXaxis()->SetBinLabel(i++, "TPC Charge fC");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  m_eventT = new TTree("eventT", "TPC FEE per-event Tree");
  assert(m_eventT);
  m_eventT->Branch("evthdr", &m_peventHeader);
  m_eventT->Branch("nClusters", &m_nClusters, "nClusters/I");
  m_IOClusters = new TClonesArray("TpcPrototypeUnpacker::ClusterData", 1000);
  m_eventT->Branch("Clusters", &m_IOClusters);

  m_chanT = new TTree("chanT", "TPC FEE per-channel Tree");
  assert(m_chanT);
  m_chanT->Branch("event", &m_eventHeader.event, "event/I");
  m_chanT->Branch("chanhdr", &m_pchanHeader);
  m_chanT->Branch("adc", m_chanData.data(), str(boost::format("adc[%d]/i") % kSAMPLE_LENGTH).c_str());

  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  //  {
  //    const PHG4CylinderCellGeom* layer_geom = seggeo->GetLayerCellGeom(layer);

  //    const string histNameCellHit(boost::str(boost::format{"hCellHit_Layer%1%"} % layer));
  //    const string histNameCellCharge(boost::str(boost::format{"hCellCharge_Layer%1%"} % layer));

  //  }

  //  hm->registerHisto(new TH2D("hLayerCellHit",  //
  //                             "Number of ADC time-bin hit per channel;Layer ID;Hit number",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             300, -.5, 299.5));
  //  hm->registerHisto(new TH2D("hLayerCellCharge",  //
  //                             "Charge integrated over drift window per channel;Layer ID;Charge [fC]",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             1000, 0, 1e7 * eplus / (1e-15 * coulomb)));
  //
  //  hm->registerHisto(new TH2D("hLayerSumCellHit",  //
  //                             "Number of ADC time-bin hit integrated over channels per layer;Layer ID;Hit number",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             10000, -.5, 99999.5));
  //  hm->registerHisto(new TH2D("hLayerSumCellCharge",  //
  //                             "Charge integrated over drift window and channel per layer;Layer ID;Charge [fC]",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             10000, 0, 1000 * 4e6 * eplus / (1e-15 * coulomb)));

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPrototypeUnpacker::End(PHCompositeNode* topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "TpcPrototypeUnpacker::End - write to " << m_outputFileName << endl;
  }
  PHTFileServer::get().cd(m_outputFileName);

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  for (unsigned int i = 0; i < hm->nHistos(); i++)
    hm->getHisto(i)->Write();

  // help index files with TChain
  TTree* T_Index = new TTree("T_Index", "T_Index");
  assert(T_Index);
  T_Index->Write();

  m_eventT->Write();
  m_chanT->Write();

  if (m_pdfMaker)
  {
    delete m_pdfMaker;
    m_pdfMaker = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPrototypeUnpacker::process_event(PHCompositeNode* topNode)
{
  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto("hNormalization"));
  assert(h_norm);

  Event* event = findNode::getClass<Event>(topNode, "PRDF");
  if (event == nullptr)
  {
    if (Verbosity() >= VERBOSITY_SOME)
      cout << "TpcPrototypeUnpacker::Process_Event - Event not found" << endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  if (Verbosity() >= VERBOSITY_SOME)
    event->identify();

  // search for data event
  if (event->getEvtType() == BEGRUNEVENT)
  {
    //    get_motor_loc(event);

    return Fun4AllReturnCodes::EVENT_OK;
  }
  if (event->getEvtType() != DATAEVENT)
    return Fun4AllReturnCodes::DISCARDEVENT;

  m_eventHeader.run = event->getRunNumber();
  m_eventHeader.event = event->getEvtSequence();

  //  m_eventHeader.xray_x = m_XRayLocationX;
  //  m_eventHeader.xray_y = m_XRayLocationY;

  if (m_pdfMaker)
  {
    m_pdfMaker->MakeSectionPage(str(boost::format("ADC signal fits for Run %1% and event %2%") % m_eventHeader.run % m_eventHeader.event));
  }

  static const char* MAX_LENGTH = "MAX_SAMPLES";
  static const char* IS_PRESENT = "IS_PRESENT";
  static const char* BX_COUNTER = "BX";

  Packet* p = event->getPacket(kPACKET_ID);
  if (p == nullptr)
    return Fun4AllReturnCodes::DISCARDEVENT;

  if (Verbosity() >= VERBOSITY_SOME) p->identify();

  if (Verbosity() >= VERBOSITY_EVEN_MORE)
  {
    cout << "TpcPrototypeUnpacker::process_event - package dump" << endl;
    p->dump();
  }

  int max_length[kN_FEES] = {-1};
  for (unsigned int i = 0; i < kN_FEES; i++)
  {
    try
    {
      max_length[i] = p->iValue(i, MAX_LENGTH);
    }
    catch (const std::out_of_range& e)
    {
      max_length[i] = -1;
      break;
    }

    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "TpcPrototypeUnpacker::process_event - max_length[" << i << "]=" << max_length[i] << endl;
    }
  }

  m_eventHeader.bx_counter = 0;
  //  bool first_channel = true;

  for (unsigned int fee = 0; fee < kN_FEES; ++fee)
  {
    try
    {
      if (!p->iValue(fee, IS_PRESENT))
      {
        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << "TpcPrototypeUnpacker::process_event - missing fee[" << fee << "]=" << endl;
        }
        continue;
      }
    }
    catch (const std::out_of_range& e)
    {
      cout << "TpcPrototypeUnpacker::process_event - out_of_range in p->iValue(" << fee << ", IS_PRESENT)" << endl;

      break;
    }

    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "TpcPrototypeUnpacker::process_event - processing fee[" << fee << "]=" << endl;
    }

    for (unsigned int channel = 0; channel < kN_CHANNELS; channel++)
    {
      m_chanHeader = ChannelHeader();
      m_chanHeader.fee_id = fee;
      m_chanHeader.size = p->iValue(fee, channel, "NR_SAMPLES");  // number of words until the next channel (header included). this is the real packet_length

      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << "TpcPrototypeUnpacker::process_event - processing fee["
             << fee << "], chan[" << channel << "] NR_SAMPLES = "
             << p->iValue(fee, channel, "NR_SAMPLES") << endl;
      }

      unsigned int real_t = 0;
      uint32_t old_bx_count = 0;
      uint16_t adc = 0;
      uint32_t bx_count = 0;
      m_chanData.resize(kSAMPLE_LENGTH, 0);

      for (unsigned int t = 0; t < kSAMPLE_LENGTH; t++)
      {
        try
        {
          adc = p->iValue(fee, channel, t);
          bx_count = p->iValue(fee, channel, t, BX_COUNTER);

          if (real_t == 0)
            m_chanHeader.bx_counter = p->iValue(fee, channel, t, BX_COUNTER);
        }
        catch (const std::out_of_range& e)
        {
          adc = 0;
          bx_count = 0;
          continue;
        }

        if (bx_count >= old_bx_count + 128)
        {
          old_bx_count = bx_count;
          real_t = 0;
          m_chanData.resize(kSAMPLE_LENGTH, 0);
        }

        assert(real_t >= 0);
        assert(real_t < kSAMPLE_LENGTH);
        m_chanData[real_t] = adc;
        //          h_raw_wave->Fill(real_t, adc);
        //          h_raw->Fill(real_t, total_chan, adc);
        //          adc_data->push_back(adc);
        real_t++;
      }  //      for (unsigned int t = 3; t < kSAMPLE_LENGTH; t++)

      //      m_chanHeader.packet_type = p->iValue(channel * kPACKET_LENGTH + 2) & 0xffff;  // that's the Elink packet type
      //      m_chanHeader.bx_counter = ((p->iValue(channel * kPACKET_LENGTH + 4) & 0xffff) << 4) | (p->iValue(channel * kPACKET_LENGTH + 5) & 0xffff);
      //      m_chanHeader.sampa_address = (p->iValue(channel * kPACKET_LENGTH + 3) >> 5) & 0xf;
      //      m_chanHeader.sampa_channel = p->iValue(channel * kPACKET_LENGTH + 3) & 0x1f;
      m_chanHeader.fee_channel = channel;

      //            const pair<int, int> pad = SAMPAChan2PadXY(m_chanHeader.fee_channel);

      m_chanHeader.pad_x = TpcR2Map.GetRpos(fee, channel);
      m_chanHeader.pad_y = TpcR2Map.Getphipos(fee, channel);

      //      if (first_channel)
      //      {
      //        first_channel = false;
      //        m_eventHeader.bx_counter = m_chanHeader.bx_counter;
      //      }
      //      else if (m_eventHeader.bx_counter != m_chanHeader.bx_counter)
      //      {
      //        m_eventHeader.bx_counter_consistent = false;
      //
      //        //      printf("TpcPrototypeUnpacker::process_event - ERROR: Malformed packet, event number %i, reason: bx_counter mismatch (expected 0x%x, got 0x%x)\n", m_eventHeader.event, m_eventHeader.bx_counter, m_chanHeader.bx_counter);
      //        //
      //        //      event->identify();
      //        //      p->identify();
      //        //      return Fun4AllReturnCodes::DISCARDEVENT;
      //      }

      //      if (m_chanHeader.fee_channel > 255 || m_chanHeader.sampa_address > 7 || m_chanHeader.sampa_channel > 31)
      //      {
      //        printf("TpcPrototypeUnpacker::process_event - ERROR: Malformed packet, event number %i, reason: bad channel (got %i, sampa_addr: %i, sampa_chan: %i)\n", m_eventHeader.event, m_chanHeader.fee_channel, m_chanHeader.sampa_address, m_chanHeader.sampa_channel);
      //
      //        event->identify();
      //        p->identify();
      //        return Fun4AllReturnCodes::DISCARDEVENT;
      //      }

      //    SampaChannel *chan = fee_data->append(new SampaChannel(fee_channel, bx_counter, packet_type));

      //      assert(m_chanData.size() == kSAMPLE_LENGTH);
      //      fill(m_chanData.begin(), m_chanData.end(), 0);
      //      for (unsigned int sample = 0; sample < kSAMPLE_LENGTH; sample++)
      //      {
      //        //        chan->append(p->iValue(channel * PACKET_LENGTH + 9 + sample) & 0xffff);
      //        uint32_t value = p->iValue(channel * kPACKET_LENGTH + 9 + sample) & 0xffff;
      //        m_chanData[sample] = value;
      //      }
      //
      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << "TpcPrototypeUnpacker::process_event - "
             << "m_chanHeader.m_size = " << int(m_chanHeader.size) << ", "
             << "m_chanHeader.m_packet_type = " << int(m_chanHeader.packet_type) << ", "
             << "m_chanHeader.m_bx_counter = " << int(m_chanHeader.bx_counter) << ", "
             << "m_chanHeader.m_sampa_address = " << int(m_chanHeader.sampa_address) << ", "
             << "m_chanHeader.m_sampa_channel = " << int(m_chanHeader.sampa_channel) << ", "
             << "m_chanHeader.m_fee_channel = " << int(m_chanHeader.fee_channel) << ": "
             << " ";

        for (unsigned int sample = 0; sample < kSAMPLE_LENGTH; sample++)
        {
          cout << "data[" << sample << "] = " << int(m_chanData[sample]) << " ";
        }

        cout << endl;
      }

      //      // fill event data
      //      if (PadPlaneData::IsValidPad(m_chanHeader.pad_x, m_chanHeader.pad_y))
      //      {
      vector<int>& paddata = m_padPlaneData.getPad(m_chanHeader.pad_x, m_chanHeader.pad_y);
      //
      for (unsigned int sample = 0; sample < kSAMPLE_LENGTH; sample++)
      {
        paddata[sample] = int(m_chanData[sample]);
      }
      //
      auto pedestal_max = roughZeroSuppression(paddata);
      m_chanHeader.pedestal = pedestal_max.first;
      m_chanHeader.max = pedestal_max.second;
      //      }
      // output per-channel TTree
      m_chanT->Fill();
    }  //    for (unsigned int channel = 0; channel < kN_CHANNELS; channel++)
  }    //  for (unsigned int fee = 0; fee < kN_FEES; ++fee)

  if (enableClustering)
  {
    int ret = Clustering();

    if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  }
  h_norm->Fill("Event count", 1);
  m_eventT->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPrototypeUnpacker::Clustering()
{
  if (Verbosity())
    cout << __PRETTY_FUNCTION__ << " entry" << endl;

  // find cluster
  m_padPlaneData.Clustering(m_clusteringZeroSuppression, Verbosity() >= VERBOSITY_SOME);
  const multimap<int, PadPlaneData::SampleID>& groups = m_padPlaneData.getGroups();

  // export clusters
  assert(m_clusters.size() == 0);  //already cleared.
  for (const auto& iter : groups)
  {
    const int& i = iter.first;
    const PadPlaneData::SampleID& id = iter.second;
    m_clusters[i].padxs.insert(id.padx);
    m_clusters[i].padys.insert(id.pady);
    m_clusters[i].samples.insert(id.sample);
  }

  // process cluster
  for (auto& iter : m_clusters)
  {
    ClusterData& cluster = iter.second;

    assert(cluster.padxs.size() > 0);
    assert(cluster.padys.size() > 0);
    assert(cluster.samples.size() > 0);

    //expand cluster by +/-1 in y
    if (*(cluster.padys.begin()) - 1 >= 0)
      cluster.padys.insert(*(cluster.padys.begin()) - 1);
    if (*(cluster.padys.rbegin()) + 1 < (int) kMaxPadY)
      cluster.padys.insert(*(cluster.padys.rbegin()) + 1);

    cluster.min_sample = max(0, *cluster.samples.begin() - m_nPreSample);
    cluster.max_sample = min((int) (kSAMPLE_LENGTH) -1, *cluster.samples.rbegin() + m_nPostSample);
    const int n_sample = cluster.max_sample - cluster.min_sample + 1;

    cluster.sum_samples.assign(n_sample, 0);
    for (int pad_x = *cluster.padxs.begin(); pad_x <= *cluster.padxs.rbegin(); ++pad_x)
    {
      cluster.padx_samples[pad_x].assign(n_sample, 0);
    }
    for (int pad_y = *cluster.padys.begin(); pad_y <= *cluster.padys.rbegin(); ++pad_y)
    {
      cluster.pady_samples[pad_y].assign(n_sample, 0);
    }

    for (int pad_x = *cluster.padxs.begin(); pad_x <= *cluster.padxs.rbegin(); ++pad_x)
    {
      for (int pad_y = *cluster.padys.begin(); pad_y <= *cluster.padys.rbegin(); ++pad_y)
      {
        assert(m_padPlaneData.IsValidPad(pad_x, pad_y));

        vector<int>& padsamples = m_padPlaneData.getPad(pad_x, pad_y);

        for (int i = 0; i < n_sample; ++i)
        {
          int adc = padsamples.at(cluster.min_sample + i);
          cluster.sum_samples[i] += adc;
          cluster.padx_samples[pad_x][i] += adc;
          cluster.pady_samples[pad_y][i] += adc;
        }

      }  //    	    for (int pad_y = *cluster.padys.begin(); pad_y<=*cluster.padys.rbegin() ;++pady)

    }  //    for (int pad_x = *cluster.padxs.begin(); pad_x<=*cluster.padxs.rbegin() ;++padx)

    if (m_pdfMaker)
    {
      m_pdfMaker->MakeSectionPage(str(boost::format("Event %1% Cluster %2%: sum all channel fit followed by fit of each pad component") % m_eventHeader.event % iter.first));
    }

    // fit - overal cluster
    map<int, double> parameters_constraints;
    {
      double peak = NAN;
      double peak_sample = NAN;
      double pedstal = NAN;
      map<int, double> parameters_io;
      SampleFit_PowerLawDoubleExp(cluster.sum_samples, peak,
                                  peak_sample, pedstal, parameters_io, Verbosity());

      parameters_constraints[1] = parameters_io[1];
      parameters_constraints[2] = parameters_io[2];
      parameters_constraints[3] = parameters_io[3];
      parameters_constraints[5] = parameters_io[5];
      parameters_constraints[6] = parameters_io[6];

      cluster.peak = peak;
      cluster.peak_sample = peak_sample;
      cluster.pedstal = pedstal;
    }

    // fit - X
    {
      //      double sum_peak = 0;
      //      double sum_peak_padx = 0;
      //      for (int pad_x = *cluster.padxs.begin(); pad_x <= *cluster.padxs.rbegin(); ++pad_x)
      //      {
      //        double peak = NAN;
      //        double peak_sample = NAN;
      //        double pedstal = NAN;
      //        map<int, double> parameters_io(parameters_constraints);
      //
      //        SampleFit_PowerLawDoubleExp(cluster.padx_samples[pad_x], peak,
      //                                    peak_sample, pedstal, parameters_io, Verbosity());
      //
      //        cluster.padx_peaks[pad_x] = peak;
      //        sum_peak += peak;
      //        sum_peak_padx += peak * pad_x;
      //      }
      //      cluster.avg_padx = sum_peak_padx / sum_peak;
      //      cluster.size_pad_x = cluster.padxs.size();

      assert(cluster.padxs.size() == 1);
      cluster.avg_padx = *cluster.padxs.begin();
      cluster.size_pad_x = cluster.padxs.size();
    }

    // fit - Y
    {
      double sum_peak = 0;
      double sum_peak_pady = 0;
      for (int pad_y = *cluster.padys.begin(); pad_y <= *cluster.padys.rbegin(); ++pad_y)
      {
        double peak = NAN;
        double peak_sample = NAN;
        double pedstal = NAN;
        map<int, double> parameters_io(parameters_constraints);

        SampleFit_PowerLawDoubleExp(cluster.pady_samples[pad_y], peak,
                                    peak_sample, pedstal, parameters_io, Verbosity());

        cluster.pady_peaks[pad_y] = peak;
        sum_peak += peak;
        sum_peak_pady += peak * pad_y;
      }
      cluster.avg_pady = sum_peak_pady / sum_peak;
      cluster.size_pad_y = cluster.padys.size();
    }
  }  //   for (auto& iter : m_clusters)

  // sort by energy
  map<double, int> cluster_energy;
  for (auto& iter : m_clusters)
  {
    //reverse energy sorting
    cluster_energy[-iter.second.peak] = iter.first;
  }

  // save clusters
  m_nClusters = 0;
  assert(m_IOClusters);
  for (const auto& iter : cluster_energy)
  {
    ClusterData& cluster = m_clusters[iter.second];

    //output to DST clusters
    int ret = exportDSTCluster(cluster, m_nClusters);
    if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

    // super awkward ways of ROOT filling TClonesArray
    new ((*m_IOClusters)[m_nClusters++]) ClusterData(cluster);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPrototypeUnpacker::exportDSTCluster(ClusterData& cluster, const int iclus)
{
  assert(tpcCylinderCellGeom);
  assert(trkrclusters);

  assert(cluster.avg_padx >= 0);
  assert(cluster.avg_padx < (int) kMaxPadX);
  const uint8_t layer = static_cast<uint8_t>(cluster.avg_padx);
  const PHG4CylinderCellGeom* layergeom = tpcCylinderCellGeom->GetLayerCellGeom(layer);
  const int NPhiBins = layergeom->get_phibins();
  const int NZBins = layergeom->get_zbins();

  // create the cluster entry directly in the node tree
  const TrkrDefs::hitsetkey hit_set_key(TpcDefs::genHitSetKey(layer, 0, 0));
  const TrkrDefs::cluskey ckey = TpcDefs::genClusKey(hit_set_key, iclus);
  TrkrClusterv1* clus = static_cast<TrkrClusterv1*>((trkrclusters->findOrAddCluster(ckey))->second);
  assert(clus);

  //  calculate geometry
  const double radius = layergeom->get_radius();  // returns center of layer
  if (cluster.avg_pady < 0)
  {
    cout << __PRETTY_FUNCTION__ << " WARNING - cluster.avg_pady = " << cluster.avg_pady << " < 0"
         << ", cluster.avg_padx = " << cluster.avg_padx << endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (cluster.avg_pady >= NPhiBins)
  {
    cout << __PRETTY_FUNCTION__ << " WARNING - cluster.avg_pady = " << cluster.avg_pady << " > " << NPhiBins
         << ", cluster.avg_padx = " << cluster.avg_padx << endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }
  const int lowery = floor(cluster.avg_pady);

  if (lowery < 0)
  {
    cout << __PRETTY_FUNCTION__ << " WARNING - cluster.avg_pady = " << cluster.avg_pady << " -> "
         << " lower = " << lowery << " < 0"
         << ", cluster.avg_padx = " << cluster.avg_padx << endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const double clusphi = layergeom->get_phicenter(lowery)                                             //
                         + (layergeom->get_phicenter(lowery + 1) - layergeom->get_phicenter(lowery))  //
                               * (cluster.avg_pady - lowery);

  assert(cluster.min_sample >= 0);
  assert(cluster.min_sample + cluster.peak_sample < NZBins);
  const double clusz = layergeom->get_zcenter(cluster.min_sample)  //
                       + (layergeom->get_zcenter(cluster.min_sample + 1) - layergeom->get_zcenter(cluster.min_sample)) * cluster.peak_sample;

  const double phi_size = cluster.size_pad_y ;                       // * radius * layergeom->get_phistep();
  const double z_size = (cluster.max_sample - cluster.min_sample) ;  // * layergeom->get_zstep();

  static const double phi_err = 170e-4;

  static const double z_err = 1000e-4;

  // Fill in the cluster details
  //================
  clus->setAdc(cluster.peak);
  clus->setPosition(0, radius * cos(clusphi));
  clus->setPosition(1, radius * sin(clusphi));
  clus->setPosition(2, clusz);
  clus->setGlobal();

  //update cluster positions
  cluster.avg_pos_x = clus->getPosition(0);
  cluster.avg_pos_y = clus->getPosition(1);
  cluster.avg_pos_z = clus->getPosition(2);

  TMatrixF DIM(3, 3);
  DIM[0][0] = 0.0;
  DIM[0][1] = 0.0;
  DIM[0][2] = 0.0;
  DIM[1][0] = 0.0;
  DIM[1][1] = phi_size * phi_size;  //cluster_v1 expects polar coordinates covariance
  DIM[1][2] = 0.0;
  DIM[2][0] = 0.0;
  DIM[2][1] = 0.0;
  DIM[2][2] = z_size * z_size;

  TMatrixF ERR(3, 3);
  ERR[0][0] = 0.0;
  ERR[0][1] = 0.0;
  ERR[0][2] = 0.0;
  ERR[1][0] = 0.0;
  ERR[1][1] = phi_err * phi_err;  //cluster_v1 expects rad, arc, z as elementsof covariance
  ERR[1][2] = 0.0;
  ERR[2][0] = 0.0;
  ERR[2][1] = 0.0;
  ERR[2][2] = z_err * z_err;

  TMatrixF ROT(3, 3);
  ROT[0][0] = cos(clusphi);
  ROT[0][1] = -sin(clusphi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(clusphi);
  ROT[1][1] = cos(clusphi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0;
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;

  TMatrixF ROT_T(3, 3);
  ROT_T.Transpose(ROT);

  TMatrixF COVAR_DIM(3, 3);
  COVAR_DIM = ROT * DIM * ROT_T;

  clus->setSize(0, 0, COVAR_DIM[0][0]);
  clus->setSize(0, 1, COVAR_DIM[0][1]);
  clus->setSize(0, 2, COVAR_DIM[0][2]);
  clus->setSize(1, 0, COVAR_DIM[1][0]);
  clus->setSize(1, 1, COVAR_DIM[1][1]);
  clus->setSize(1, 2, COVAR_DIM[1][2]);
  clus->setSize(2, 0, COVAR_DIM[2][0]);
  clus->setSize(2, 1, COVAR_DIM[2][1]);
  clus->setSize(2, 2, COVAR_DIM[2][2]);
  //cout << " covar_dim[2][2] = " <<  COVAR_DIM[2][2] << endl;

  TMatrixF COVAR_ERR(3, 3);
  COVAR_ERR = ROT * ERR * ROT_T;

  clus->setError(0, 0, COVAR_ERR[0][0]);
  clus->setError(0, 1, COVAR_ERR[0][1]);
  clus->setError(0, 2, COVAR_ERR[0][2]);
  clus->setError(1, 0, COVAR_ERR[1][0]);
  clus->setError(1, 1, COVAR_ERR[1][1]);
  clus->setError(1, 2, COVAR_ERR[1][2]);
  clus->setError(2, 0, COVAR_ERR[2][0]);
  clus->setError(2, 1, COVAR_ERR[2][1]);
  clus->setError(2, 2, COVAR_ERR[2][2]);

  if (Verbosity() >= 2)
  {
    cout << __PRETTY_FUNCTION__ << "Dump clusters after TpcPrototypeClusterizer" << endl;
    clus->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

TpcPrototypeUnpacker::PadPlaneData::
    PadPlaneData()
  : m_data(kMaxPadY, vector<vector<int>>(kMaxPadX, vector<int>(kSAMPLE_LENGTH, 0)))
{
}

void TpcPrototypeUnpacker::PadPlaneData::Reset()
{
  for (auto& padrow : m_data)
  {
    for (auto& pad : padrow)
    {
      fill(pad.begin(), pad.end(), 0);
    }
  }

  m_groups.clear();
}

bool TpcPrototypeUnpacker::PadPlaneData::IsValidPad(const int pad_x, const int pad_y)
{
  return (pad_x >= 0) and
         (pad_x < int(kMaxPadX)) and
         (pad_y >= 0) and
         (pad_y < int(kMaxPadY));
}

vector<int>& TpcPrototypeUnpacker::PadPlaneData::getPad(const int pad_x, const int pad_y)
{
  assert(pad_x >= 0);
  assert(pad_x < int(kMaxPadX));
  assert(pad_y >= 0);
  assert(pad_y < int(kMaxPadY));

  return m_data[pad_y][pad_x];
}

std::pair<int, int> TpcPrototypeUnpacker::roughZeroSuppression(std::vector<int>& data)
{
  std::vector<int> sorted_data(data);

  sort(sorted_data.begin(), sorted_data.end());

  const int pedestal = sorted_data[sorted_data.size() / 2];
  const int max = sorted_data.back();

  for (auto& d : data)
    d -= pedestal;

  return make_pair(pedestal, max);
}

bool operator<(const TpcPrototypeUnpacker::PadPlaneData::SampleID& s1, const TpcPrototypeUnpacker::PadPlaneData::SampleID& s2)
{
  if (s1.pady == s2.pady)
  {
    if (s1.padx == s2.padx)
    {
      return s1.sample < s2.sample;
    }
    else
      return s1.padx < s2.padx;
  }
  else
    return s1.pady < s2.pady;
}

//! 3-D Graph clustering based on PHMakeGroups()
void TpcPrototypeUnpacker::PadPlaneData::Clustering(int zero_suppression, bool verbosity)
{
  using namespace boost;
  typedef adjacency_list<vecS, vecS, undirectedS> Graph;
  typedef bimap<Graph::vertex_descriptor, SampleID> VertexList;

  Graph G;
  VertexList vertex_list;

  for (unsigned int pady = 0; pady < kMaxPadY; ++pady)
  {
    for (unsigned int padx = 0; padx < kMaxPadX; ++padx)
    {
      for (unsigned int sample = 0; sample < kSAMPLE_LENGTH; sample++)
      {
        if (m_data[pady][padx][sample] > zero_suppression)
        {
          SampleID id{(int) (pady), (int) (padx), (int) (sample)};
          Graph::vertex_descriptor v = boost::add_vertex(G);
          vertex_list.insert(VertexList::value_type(v, id));

          add_edge(v, v, G);
        }
      }  //      for (unsigned int sample = 0; sample < kSAMPLE_LENGTH; sample++)
    }
  }  //   for (unsigned int pady = 0; pady < kMaxPadY; ++pady)

  // connect 2-D adjacent samples within each pad_x
  vector<SampleID> search_directions;
  search_directions.push_back(SampleID{0, 0, 1});
  //  search_directions.push_back(SampleID{0, 1, 0});
  search_directions.push_back(SampleID{1, 0, 0});

  for (const auto& it : vertex_list.right)
  {
    const SampleID id = it.first;
    const Graph::vertex_descriptor v = it.second;

    for (const SampleID& search_direction : search_directions)
    {
      //      const SampleID next_id = id + search_direction;
      SampleID next_id(id);
      next_id.adjust(search_direction);

      auto next_it = vertex_list.right.find(next_id);
      if (next_it != vertex_list.right.end())
      {
        add_edge(v, next_it->second, G);
      }
    }

  }  //  for (const auto & it : vertex_list)

  // Find the connections between the vertices of the graph (vertices are the rawhits,
  // connections are made when they are adjacent to one another)
  std::vector<int> component(num_vertices(G));
  connected_components(G, &component[0]);

  // Loop over the components(vertices) compiling a list of the unique
  // connections (ie clusters).
  set<int> comps;                // Number of unique components
  assert(m_groups.size() == 0);  // no overwrite

  for (unsigned int i = 0; i < component.size(); i++)
  {
    comps.insert(component[i]);
    m_groups.insert(make_pair(component[i], vertex_list.left.find(vertex(i, G))->second));
  }

  //debug prints
  if (verbosity)
    for (const int& comp : comps)
    {
      cout << "TpcPrototypeUnpacker::PadPlaneData::Clustering - find cluster " << comp << " containing ";
      const auto range = m_groups.equal_range(comp);

      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const SampleID& id = iter->second;
        cout << "adc[" << id.pady << "][" << id.padx << "][" << id.sample << "] = " << m_data[id.pady][id.padx][id.sample] << ", ";
      }
      cout << endl;
    }  //  for (const int& comp : comps)
}

Fun4AllHistoManager*
TpcPrototypeUnpacker::getHistoManager()
{
  static string histname("TpcPrototypeUnpacker_HISTOS");

  Fun4AllServer* se = Fun4AllServer::instance();
  Fun4AllHistoManager* hm = se->getHistoManager(histname);

  if (not hm)
  {
    cout
        << "TpcPrototypeUnpacker::get_HistoManager - Making Fun4AllHistoManager "
        << histname << endl;
    hm = new Fun4AllHistoManager(histname);
    se->registerHistoManager(hm);
  }

  assert(hm);

  return hm;
}

void TpcPrototypeUnpacker::registerPadPlane(PHG4TpcPadPlane* inpadplane)
{
  cout << "Registering padplane " << endl;
  padplane = inpadplane;
  padplane->Detector("TPC");
  padplane->UpdateInternalParameters();
  if (Verbosity())
    cout << __PRETTY_FUNCTION__ << " padplane registered and parameters updated" << endl;

  return;
}

int TpcPrototypeUnpacker::InitField(PHCompositeNode* topNode)
{
  if (Verbosity() >= 1) cout << "TpcPrototypeUnpacker::InitField - create magnetic field setup" << endl;

  unique_ptr<PHFieldConfig> default_field_cfg(nullptr);

  default_field_cfg.reset(new PHFieldConfigv2(0, 0, 0));

  PHField* phfield = PHFieldUtility::GetFieldMapNode(default_field_cfg.get(), topNode, Verbosity() + 1);
  assert(phfield);

  return Fun4AllReturnCodes::EVENT_OK;
}

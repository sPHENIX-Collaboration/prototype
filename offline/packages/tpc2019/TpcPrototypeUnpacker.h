/*
 * TpcPrototypeUnpacker.h
 *
 *  Created on: Sep 19, 2018
 *      Author: jinhuang
 */

#ifndef CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TpcPrototypeUnpacker_H_
#define CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TpcPrototypeUnpacker_H_

#include <fun4all/SubsysReco.h>

#include <TObject.h>

#include <stdint.h>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

class PHCompositeNode;
class Fun4AllHistoManager;
class TTree;
class TClonesArray;
class Event;
class PHG4TpcPadPlane;
class PHG4CylinderCellGeomContainer;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;

namespace TpcPrototypeDefs
{
namespace FEEv2
{
class SampleFit_PowerLawDoubleExp_PDFMaker;
}
}  // namespace TpcPrototypeDefs

class TpcPrototypeUnpacker : public SubsysReco
{
 public:
  TpcPrototypeUnpacker(const std::string &outputfilename =
                           "TpcPrototypeUnpacker.root");
  virtual ~TpcPrototypeUnpacker();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void setClusteringZeroSuppression(int threshold)
  {
    m_clusteringZeroSuppression = threshold;
  }

  void setNPostSample(int nPostSample)
  {
    m_nPostSample = nPostSample;
  }

  void setNPreSample(int nPreSample)
  {
    m_nPreSample = nPreSample;
  }

  void registerPadPlane(PHG4TpcPadPlane *padplane);

  void setEnableClustering(bool b) { enableClustering = b; }
  //! simple event header class for ROOT file IO
  class EventHeader : public TObject
  {
   public:
    int run;
    int event;

    uint32_t bx_counter;
    bool bx_counter_consistent;

    int xray_x;
    int xray_y;

    EventHeader()
      : run(-1)
      , event(-1)
      , bx_counter(0)
      , bx_counter_consistent(true)
      , xray_x(-1)
      , xray_y(-1)
    {
    }

    ClassDef(TpcPrototypeUnpacker::EventHeader, 1)
  };

  //! buffer for full event data
  class PadPlaneData
  {
   public:
    PadPlaneData();
    void Reset();

    struct SampleID
    {
      int pady;
      int padx;
      int sample;

      void adjust(const SampleID &adjustment)
      {
        pady += adjustment.pady;
        padx += adjustment.padx;
        sample += adjustment.sample;
      }
    };

    static bool IsValidPad(const int pad_x, const int pad_y);
    std::vector<int> &getPad(const int pad_x, const int pad_y);
    int getSample(const SampleID &id);

    //! 3-D Graph clustering based on PHMakeGroups()
    void Clustering(int zero_suppression, bool verbosity = false);

#ifndef __CINT__

    const std::vector<std::vector<std::vector<int>>> &getData() const
    {
      return m_data;
    }

    const std::multimap<int, SampleID> &getGroups() const
    {
      return m_groups;
    }

   private:
    //! full event data in index order of m_data[pady][padx][sample]
    std::vector<std::vector<std::vector<int>>> m_data;

    std::multimap<int, SampleID> m_groups;

#endif  // #ifndef __CINT__
  };

  //! buffer for a cluster's data
  class ClusterData : public TObject
  {
   public:
    ClusterData()
      : min_sample(-1)
      , max_sample(-1)
      , peak(NAN)
      , peak_sample(NAN)
      , pedstal(NAN)
      , avg_padx(-1)
      , avg_pady(NAN)
      , size_pad_x(-1)
      , size_pad_y(-1)
      , avg_pos_x(NAN)
      , avg_pos_y(NAN)
      , avg_pos_z(NAN)
    {
    }

    void Clear(Option_t * /*option*/ = "");

    std::set<int> padxs;
    std::set<int> padys;
    std::set<int> samples;

    std::map<int, std::vector<double>> padx_samples;
    std::map<int, std::vector<double>> pady_samples;
    std::vector<double> sum_samples;

    int min_sample;
    int max_sample;

    double peak;
    double peak_sample;
    double pedstal;

    std::map<int, double> padx_peaks;
    std::map<int, double> pady_peaks;

    //! pad coordinate
    int avg_padx;
    double avg_pady;

    //! pad size
    int size_pad_x;
    int size_pad_y;

    //! pad coordinate
    double avg_pos_x;
    double avg_pos_y;
    double avg_pos_z;

    ClassDef(TpcPrototypeUnpacker::ClusterData, 3);
  };

  //! simple channel header class for ROOT file IO
  class ChannelHeader : public TObject
  {
   public:
    int size;
    //! = p->iValue(channel * kPACKET_LENGTH + 2) & 0xffff;  // that's the Elink packet type
    int packet_type;
    //! = ((p->iValue(channel * kPACKET_LENGTH + 4) & 0xffff) << 4) | (p->iValue(channel * kPACKET_LENGTH + 5) & 0xffff);
    int bx_counter;
    //! = (p->iValue(channel * kPACKET_LENGTH + 3) >> 5) & 0xf;
    int sampa_address;
    //! = p->iValue(channel * kPACKET_LENGTH + 3) & 0x1f;
    int sampa_channel;
    //! = (sampa_address << 5) | sampa_channel;
    int fee_channel;

    //! which fee
    int fee_id;

    //! pad coordinate
    int pad_x;
    int pad_y;

    int pedestal;
    int max;

    ChannelHeader()
      : size(-1)
      , packet_type(-1)
      , bx_counter(-1)
      , sampa_address(-1)
      , sampa_channel(-1)
      , fee_channel(-1)
      , fee_id(-1)
      , pad_x(-1)
      , pad_y(-1)
      , pedestal(-1)
      , max(-1)
    {
    }

    ClassDef(TpcPrototypeUnpacker::ChannelHeader, 1)
  };

 private:
  PHG4TpcPadPlane *padplane;
  PHG4CylinderCellGeomContainer *tpcCylinderCellGeom;
  TrkrClusterContainer *trkrclusters;
  int exportDSTCluster(ClusterData &cluster, const int i);
  int InitField(PHCompositeNode *topNode);

#ifndef __CINT__

  // IO stuff

  Fun4AllHistoManager *getHistoManager();

  std::string m_outputFileName;

  TTree *m_eventT;

  EventHeader m_eventHeader;
  EventHeader *m_peventHeader;  //! ->m_eventHeader,  for filling TTree

  int m_nClusters;
  TClonesArray *m_IOClusters;

  TTree *m_chanT;

  ChannelHeader m_chanHeader;
  ChannelHeader *m_pchanHeader;  //! ->m_chanHeader,  for filling TTree
  std::vector<uint32_t> m_chanData;

  // clustering stuff
  PadPlaneData m_padPlaneData;
  std::map<int, ClusterData> m_clusters;

  //! rough zero suppression by subtracting sample medium value
  //! \return pair of pedestal and max
  static std::pair<int, int> roughZeroSuppression(std::vector<int> &data);

  bool enableClustering;
  //! Clustering then prepare IOs
  int Clustering(void);

#endif  // #ifndef __CINT__

  int m_clusteringZeroSuppression;
  int m_nPreSample;
  int m_nPostSample;
  int m_XRayLocationX;
  int m_XRayLocationY;

  TpcPrototypeDefs::FEEv2::SampleFit_PowerLawDoubleExp_PDFMaker *m_pdfMaker;
};

bool operator<(const TpcPrototypeUnpacker::PadPlaneData::SampleID &s1, const TpcPrototypeUnpacker::PadPlaneData::SampleID &s2);

#endif /* CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TpcPrototypeUnpacker_H_ */

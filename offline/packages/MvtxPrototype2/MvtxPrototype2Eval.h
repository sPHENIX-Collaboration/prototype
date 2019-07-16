#ifndef MVTX_MVTXEVAL_H
#define MVTX_MVTXEVAL_H

#include <fun4all/SubsysReco.h>
#include <utility>
#include "MvtxPrototype2Geom.h"

class TrkrHit;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TH2D;
class TH1D;

class MvtxPrototype2Eval : public SubsysReco
{
 public:

  MvtxPrototype2Eval(const std::string &name = "MvtxPrototype2Eval");
  virtual ~MvtxPrototype2Eval() {}

  //! module initialization
  int Init(PHCompositeNode *topNode);

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode);

  void set_filename(const char* file)
  { if(file) _outfile = file; }

 private:

  std::string _outfile;

  TrkrHitSetContainer *m_hits;
  TrkrClusterContainer *m_clusterlist;

  int _n_events;
  TH1D *h1d_nevents;
  TH2D *h2d_hit[NLAYER];
  TH2D *h2d_clu[NLAYER];
  TH1D *h1d_clu_diffrow[NLAYER];
  TH1D *h1d_clu_diffcol[NLAYER];
  TH1D *h1d_clus_size;

};

#endif  // MVTX_MVTXCLUSTERIZER_H

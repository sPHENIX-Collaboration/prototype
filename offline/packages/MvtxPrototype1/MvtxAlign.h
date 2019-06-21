#ifndef __PHG4SVTXCLUSTERIZER_H__
#define __PHG4SVTXCLUSTERIZER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <map>
#include <limits.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefUtil.h>

class TrkrHitSetContainer;

class MvtxAlign : public SubsysReco {

public:

  struct AlignmentPar {
    double dx;
    double dy;
    double dz;
  };


  MvtxAlign(const std::string &name = "MvtxAlign");
  virtual ~MvtxAlign() {}

  //! module initialization
  int Init(PHCompositeNode *topNode) {return 0;}

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode) {return 0;}

  //! put in misalignment by hand
  void AddAlignmentPar(TrkrDefs::hitsetkey key, double dx, double dy, double dz);

  //! print stored misalignments
  void PrintAlignmentPars(std::ostream &os = std::cout) const;

  //! set the directory for alignment parameter files
  void SetAlignmentParFileDir(const std::string fdir) { fdir_ = fdir; }

  //! set boolean to read from file rather than setting by hand
  void SetAlignmentParFromFile(const bool yn) { apff_ = yn; }

private:

  // read the alignment parameters from file based on run number
  int ReadAlignmentParFile();

  // node tree storage pointers
  TrkrClusterContainer* clusters_;

  // storage object for misalignments
  std::map<TrkrDefs::hitsetkey, AlignmentPar> alignmap_;

  // directory for alignment parameter files
  std::string fdir_;
  int runnumber_;
  bool apff_; // alignment par from file

  PHTimeServer::timer _timer;
};

#endif

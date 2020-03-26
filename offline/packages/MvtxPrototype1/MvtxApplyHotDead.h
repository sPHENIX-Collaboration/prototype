#ifndef __PHG4SVTXCLUSTERIZER_H__
#define __PHG4SVTXCLUSTERIZER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <map>
#include <limits.h>

#include <mvtx/MvtxHitSetv1.h>
#include <trackbase/TrkrDefUtil.h>

class TrkrHitSetContainer;

class MvtxApplyHotDead : public SubsysReco {

public:

  typedef std::pair<unsigned int, unsigned int> Pixel;
  typedef std::multimap<TrkrDefs::hitsetkey, Pixel> PixelMap;
  typedef PixelMap::const_iterator ConstIterator;

  MvtxApplyHotDead(const std::string &name = "MvtxApplyHotDead",
                   const std::string file = "");
  virtual ~MvtxApplyHotDead() {}

  //! module initialization
  int Init(PHCompositeNode *topNode) {return 0;}

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode) {return 0;}

  //! print pixels to be masked
  void PrintHotDeadMap(std::ostream& os = std::cout) const;

private:

  //! read in the pixels to mask
  void ReadHotDeadFile();

  // node tree storage pointers
  TrkrHitSetContainer* hits_;

  std::string hdfile_; // hot dead file name
  PixelMap hdmap_;


  PHTimeServer::timer _timer;
};

#endif

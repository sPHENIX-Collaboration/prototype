#ifndef __PHG4SVTXCLUSTERIZER_H__
#define __PHG4SVTXCLUSTERIZER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <map>
#include <limits.h>

#include <mvtx/MvtxHitSetv1.h>

class TrkrHitSetContainer;
class TrkrClusterContainer;

class MvtxClusterizer : public SubsysReco {

public:

  typedef std::pair<unsigned int, unsigned int> pixel;

  MvtxClusterizer(const std::string &name = "MvtxClusterizer");
  virtual ~MvtxClusterizer(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode){return 0;}
  
  //! option to turn off z-dimension clustering
  void SetZClustering(const bool make_z_clustering) {
    make_z_clustering_ = make_z_clustering;
  }
  bool GetZClustering() const {
      return make_z_clustering_;
  }

private:

  bool are_adjacent(const pixel lhs, const pixel rhs );

  void ClusterMvtx(PHCompositeNode *topNode);

  void PrintClusters(PHCompositeNode *topNode);
  
  // node tree storage pointers
  TrkrHitSetContainer* hits_;
  TrkrClusterContainer* clusterlist_;

  // settings
  bool make_z_clustering_;    // z_clustering_option

  PHTimeServer::timer _timer;
};

#endif

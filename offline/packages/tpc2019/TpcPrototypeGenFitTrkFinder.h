/*!
 *  \file		TpcPrototypeGenFitTrkFinder.h
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_TpcPrototypeGenFitTrkFinder_H
#define TRACKRECO_TpcPrototypeGenFitTrkFinder_H

#include <fun4all/SubsysReco.h>

#include <TMatrixFfwd.h>  // for TMatrixF
#include <TVector3.h>     // for TVector3

#include <cstddef>  // for NULL
#include <memory>   // for shared_ptr
#include <string>
#include <vector>

class TClonesArray;

namespace PHGenFit
{
class Track;
} /* namespace PHGenFit */

class SvtxTrack;
namespace PHGenFit
{
class Fitter;
} /* namespace PHGenFit */

class SvtxTrackMap;
class SvtxVertexMap;
class TrkrCluster;
class SvtxVertex;
class PHCompositeNode;
class AssocInfoContainer;
class TrkrClusterContainer;
class TTree;

//! \brief	find tracks
class TpcPrototypeGenFitTrkFinder : public SubsysReco
{
 public:
  //! Default constructor
  TpcPrototypeGenFitTrkFinder(const std::string& name = "TpcPrototypeGenFitTrkFinder", int layers = 16);

  //! dtor
  ~TpcPrototypeGenFitTrkFinder();

  //!Initialization, called for initialization
  int Init(PHCompositeNode*);

  //!Initialization Run, called for initialization of a run
  int InitRun(PHCompositeNode*);

  //!Process Event, called for each event
  int process_event(PHCompositeNode*);

  //!End, write and close files
  int End(PHCompositeNode*);

  //! For evalution
  //! Change eval output filename
  void set_eval_filename(const char* file)
  {
    if (file)
      _eval_outname = file;
  }
  std::string get_eval_filename() const
  {
    return _eval_outname;
  }

  void fill_eval_tree(PHCompositeNode*);
  void init_eval_tree();
  void reset_eval_variables();

  bool is_do_eval() const
  {
    return _do_eval;
  }

  void set_do_eval(bool doEval)
  {
    _do_eval = doEval;
  }

  bool is_do_evt_display() const
  {
    return _do_evt_display;
  }

  void set_do_evt_display(bool doEvtDisplay)
  {
    _do_evt_display = doEvtDisplay;
  }

  const std::string& get_track_fitting_alg_name() const
  {
    return _track_fitting_alg_name;
  }

  void set_track_fitting_alg_name(const std::string& trackFittingAlgName)
  {
    _track_fitting_alg_name = trackFittingAlgName;
  }

  int get_primary_pid_guess() const
  {
    return _primary_pid_guess;
  }

  void set_primary_pid_guess(int primaryPidGuess)
  {
    _primary_pid_guess = primaryPidGuess;
  }

 private:
  //! Event counter
  int _event;

  //! Get all the nodes
  int GetNodes(PHCompositeNode*);

  //!Create New nodes
  int CreateNodes(PHCompositeNode*);

  typedef std::vector<const TrkrCluster*> tracklet_t;

  double getChi2Ndf(const tracklet_t & tracklet);

  PHGenFit::Fitter* _fitter;

  //! KalmanFitterRefTrack, KalmanFitter, DafSimple, DafRef
  std::string _track_fitting_alg_name;

  int nLayer;
  unsigned int minLayer;
  unsigned int maxTracklet;
  int _primary_pid_guess;
  double rphiWindow;
  double ZWindow;

  //! Input Node pointers
  //  PHG4TruthInfoContainer* _truth_container;
  TrkrClusterContainer* _clustermap;
  SvtxTrackMap* _trackmap;
  AssocInfoContainer *_assoc_container;
  SvtxVertexMap* _vertexmap;

  //! Evaluation
  //! switch eval out
  bool _do_eval;

  //! eval output filename
  std::string _eval_outname;

  TTree* _eval_tree;
  int _tca_ntrack;
  TClonesArray* _tca_trackmap;

  TTree* _cluster_eval_tree;
  float _cluster_eval_tree_x;
  float _cluster_eval_tree_y;
  float _cluster_eval_tree_z;
  float _cluster_eval_tree_gx;
  float _cluster_eval_tree_gy;
  float _cluster_eval_tree_gz;

  bool _do_evt_display;
};

#endif

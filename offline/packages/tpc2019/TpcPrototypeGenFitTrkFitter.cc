/*!
 *  \file		TpcPrototypeGenFitTrkFitter.C
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "TpcPrototypeGenFitTrkFitter.h"
#include "TpcPrototypeTrack.h"
#include "TpcPrototypeCluster.h"

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack_v1.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackState.h>  // for SvtxTrackState
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>  // for SvtxVertexMap
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex_v1.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/Measurement.h>  // for Measurement
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>
#include <phgenfit/Track.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <phfield/PHFieldUtility.h>
#include <phgeom/PHGeomUtility.h>

#include <GenFit/AbsMeasurement.h>  // for AbsMeasurement
#include <GenFit/EventDisplay.h>    // for EventDisplay
#include <GenFit/Exception.h>       // for Exception
#include <GenFit/GFRaveConverters.h>
#include <GenFit/GFRaveTrackParameters.h>  // for GFRaveTrackParame...
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>
#include <GenFit/TrackPoint.h>  // for TrackPoint

//Rave
#include <rave/ConstantMagneticField.h>
#include <rave/VacuumPropagator.h>  // for VacuumPropagator
#include <rave/VertexFactory.h>

#include <TClonesArray.h>
#include <TMatrixDSymfwd.h>  // for TMatrixDSym
#include <TMatrixFfwd.h>     // for TMatrixF
#include <TMatrixT.h>        // for TMatrixT, operator*
#include <TMatrixTSym.h>     // for TMatrixTSym
#include <TMatrixTUtils.h>   // for TMatrixTRow
#include <TRotation.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVectorDfwd.h>  // for TVectorD
#include <TVectorT.h>     // for TVectorT

#include <cassert>
#include <cmath>  // for sqrt, NAN
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

class PHField;
class TGeoManager;
namespace genfit
{
class AbsTrackRep;
}

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

#define WILD_FLOAT -9999.

#define _DEBUG_MODE_ 0

//#define _DEBUG_

using namespace std;

class PHRaveVertexFactory
{
 public:
  //! ctor
  PHRaveVertexFactory(const int Verbosity())
  {
    rave::ConstantMagneticField mfield(0., 0., 0.);  // RAVE use Tesla
    _factory = new rave::VertexFactory(mfield, rave::VacuumPropagator(),
                                       "default", Verbosity());

    IdGFTrackStateMap_.clear();
  }

  //! dotr
  ~PHRaveVertexFactory()
  {
    clearMap();

    delete _factory;
  }

  void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
                    const std::vector<genfit::Track*>& tracks, const bool use_beamspot = false)
  {
    clearMap();

    try
    {
      genfit::RaveToGFVertices(vertices,
                               _factory->create(
                                   genfit::GFTracksToTracks(tracks, NULL,
                                                            IdGFTrackStateMap_, 0),
                                   use_beamspot),
                               IdGFTrackStateMap_);
    }
    catch (genfit::Exception& e)
    {
      std::cerr << e.what();
    }
  }

  void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
                    const std::vector<genfit::Track*>& tracks,
                    std::vector<genfit::MeasuredStateOnPlane*>& GFStates,
                    const bool use_beamspot = false)
  {
    clearMap();

    try
    {
      genfit::RaveToGFVertices(vertices,
                               _factory->create(
                                   genfit::GFTracksToTracks(tracks, &GFStates,
                                                            IdGFTrackStateMap_, 0),
                                   use_beamspot),
                               IdGFTrackStateMap_);
    }
    catch (genfit::Exception& e)
    {
      std::cerr << e.what();
    }
  }

 private:
  void clearMap()
  {
    for (unsigned int i = 0; i < IdGFTrackStateMap_.size(); ++i)
      delete IdGFTrackStateMap_[i].state_;

    IdGFTrackStateMap_.clear();
  }

  std::map<int, genfit::trackAndState> IdGFTrackStateMap_;

  rave::VertexFactory* _factory;
};

/*
 * Constructor
 */
TpcPrototypeGenFitTrkFitter::TpcPrototypeGenFitTrkFitter(const string& name)
  : SubsysReco(name)
  , _flags(NONE)
  , _output_mode(TpcPrototypeGenFitTrkFitter::MakeNewNode)
  , _over_write_svtxtrackmap(true)
  , _over_write_svtxvertexmap(true)
  , _fit_primary_tracks(false)
  , _use_truth_vertex(false)
  , _fitter(NULL)
  , _track_fitting_alg_name("DafRef")
  , _primary_pid_guess(2212)
  , _fit_min_pT(0.1)
  , _vertex_min_ndf(20)
  , _vertex_finder(NULL)
  , _vertexing_method("avf-smoothing:1")
  //  , _truth_container(NULL)
  , _clustermap(NULL)
  , _trackmap(NULL)
  , _vertexmap(NULL)
  , _trackmap_refit(NULL)
  , _primary_trackmap(NULL)
  , _vertexmap_refit(NULL)
  , _do_eval(false)
  , _eval_outname("TpcPrototypeGenFitTrkFitter_eval.root")
  , _eval_tree(NULL)
  , _tca_ntrack(-1)
  , _tca_particlemap(NULL)
  , _tca_vtxmap(NULL)
  , _tca_trackmap(NULL)
  , _tca_tpctrackmap(nullptr)
  , _tca_vertexmap(NULL)
  , _tca_trackmap_refit(NULL)
  , _tca_primtrackmap(NULL)
  , _tca_vertexmap_refit(NULL)
  , _do_evt_display(false)
{
  Verbosity(0);

  _event = 0;

  _cluster_eval_tree = NULL;
  _cluster_eval_tree_x = WILD_FLOAT;
  _cluster_eval_tree_y = WILD_FLOAT;
  _cluster_eval_tree_z = WILD_FLOAT;
  _cluster_eval_tree_gx = WILD_FLOAT;
  _cluster_eval_tree_gy = WILD_FLOAT;
  _cluster_eval_tree_gz = WILD_FLOAT;
}

/*
 * Init
 */
int TpcPrototypeGenFitTrkFitter::Init(PHCompositeNode* topNode)
{
  //	CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Init run
 */
int TpcPrototypeGenFitTrkFitter::InitRun(PHCompositeNode* topNode)
{
  CreateNodes(topNode);

  TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
  PHField* field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

  //_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
  _fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
                                          field, _track_fitting_alg_name,
                                          "RKTrackRep", _do_evt_display);

  if (!_fitter)
  {
    cerr << PHWHERE << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _fitter->set_verbosity(Verbosity());

  //LogDebug(genfit::FieldManager::getInstance()->getFieldVal(TVector3(0, 0, 0)).Z());

  _vertex_finder = new genfit::GFRaveVertexFactory(Verbosity());
  //_vertex_finder->setMethod("kalman-smoothing:1"); //! kalman-smoothing:1 is the defaul method
  _vertex_finder->setMethod(_vertexing_method.data());
  //_vertex_finder->setBeamspot();

  //_vertex_finder = new PHRaveVertexFactory(Verbosity());

  if (!_vertex_finder)
  {
    cerr << PHWHERE << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (_do_eval)
  {
    if (Verbosity() >= 1)
      cout << PHWHERE << " Openning file: " << _eval_outname << endl;
    PHTFileServer::get().open(_eval_outname, "RECREATE");
    init_eval_tree();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
/*
 * process_event():
 *  Call user instructions for every event.
 *  This function contains the analysis structure.
 *
 */
int TpcPrototypeGenFitTrkFitter::process_event(PHCompositeNode* topNode)
{
  _event++;

  if (Verbosity() > 1)
    std::cout << PHWHERE << "Events processed: " << _event << std::endl;
  //	if (_event % 1000 == 0)
  //		cout << PHWHERE << "Events processed: " << _event << endl;

  //cout << "Start TpcPrototypeGenFitTrkFitter::process_event" << endl;

  GetNodes(topNode);

  //! stands for Refit_GenFit_Tracks
  vector<genfit::Track*> rf_gf_tracks;
  rf_gf_tracks.clear();

  vector<std::shared_ptr<PHGenFit::Track> > rf_phgf_tracks;
  //  rf_phgf_tracks.clear();

  map<unsigned int, unsigned int> svtxtrack_genfittrack_map;

  if (_trackmap_refit)
    _trackmap_refit->empty();

  // _trackmap is SvtxTrackMap from the node tree
  for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* svtx_track = iter->second;
    if (Verbosity() > 50)
    {
      cout << "   process SVTXTrack " << iter->first << endl;
      svtx_track->identify();
    }
    if (!svtx_track)
      continue;
    if (!(svtx_track->get_pt() > _fit_min_pT))
      continue;

    //! stands for Refit_PHGenFit_Track
    std::shared_ptr<PHGenFit::Track> rf_phgf_track = ReFitTrack(topNode, svtx_track);

    if (rf_phgf_track)
    {
      svtxtrack_genfittrack_map[svtx_track->get_id()] =
          rf_phgf_tracks.size();
      rf_phgf_tracks.push_back(rf_phgf_track);
      if (rf_phgf_track->get_ndf() > _vertex_min_ndf)
        rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());
    }
  }

  /*
	 * add tracks to event display
	 * needs to make copied for smart ptrs will be destroied even
	 * there are still references in TGeo::EventView
	 */
  if (_do_evt_display)
  {
    //search for unused clusters

    assert(_clustermap);
    //unused clusters
    set<TrkrDefs::cluskey> cluster_ids;
    //all clusters
    auto cluster_range = _clustermap->getClusters(TrkrDefs::tpcId);
    for (auto iter = cluster_range.first; iter != cluster_range.second; ++iter)
    {
      cluster_ids.insert(iter->first);
    }
    // minus used clusters
    for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
         ++iter)
    {
      SvtxTrack* svtx_track = iter->second;
      for (auto iter = svtx_track->begin_cluster_keys();
           iter != svtx_track->end_cluster_keys(); ++iter)
      {
        cluster_ids.erase(*iter);
      }
    }

    //    // add tracks
    vector<genfit::Track*> copy;
    for (genfit::Track* t : rf_gf_tracks)
    {
      copy.push_back(new genfit::Track(*t));
    }
    //

    //add clusters
    for (const auto cluster_id : cluster_ids)
    {
      const TrkrCluster* cluster =
          _clustermap->findCluster(cluster_id);
      assert(cluster);

      std::shared_ptr<PHGenFit::Track> cluster_holder = DisplayCluster(cluster);
      if (cluster_holder.get())
        copy.push_back(new genfit::Track(*cluster_holder->getGenFitTrack()));
    }

    _fitter->getEventDisplay()->addEvent(copy);
  }

  for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();)
  {
    std::shared_ptr<PHGenFit::Track> rf_phgf_track = NULL;

    if (svtxtrack_genfittrack_map.find(iter->second->get_id()) != svtxtrack_genfittrack_map.end())
    {
      unsigned int itrack =
          svtxtrack_genfittrack_map[iter->second->get_id()];
      rf_phgf_track = rf_phgf_tracks[itrack];
    }

    if (rf_phgf_track)
    {
      //FIXME figure out which vertex to use.
      SvtxVertex* vertex = NULL;
      if (_over_write_svtxvertexmap)
      {
        if (_vertexmap->size() > 0)
          vertex = _vertexmap->get(0);
        //cout << PHWHERE << "        will use vertex " << vertex->get_x() << "  " << vertex->get_y() << "  " << vertex->get_z() << endl;
      }
      else
      {
        if (_vertexmap_refit->size() > 0)
          vertex = _vertexmap_refit->get(0);
      }

      std::shared_ptr<SvtxTrack> rf_track = MakeSvtxTrack(iter->second, rf_phgf_track,
                                                          vertex);
      if (!rf_track)
      {
        //if (_output_mode == OverwriteOriginalNode)
        if (_over_write_svtxtrackmap)
        {
          auto key = iter->first;
          ++iter;
          _trackmap->erase(key);
          continue;
        }
      }

      //			delete vertex;//DEBUG

      //			rf_phgf_tracks.push_back(rf_phgf_track);
      //			rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());

      if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
        if (_trackmap_refit)
        {
          _trackmap_refit->insert(rf_track.get());
          //					delete rf_track;
        }

      if (_over_write_svtxtrackmap || _output_mode == DebugMode)
      {
        *(dynamic_cast<SvtxTrack_v1*>(iter->second)) =
            *(dynamic_cast<SvtxTrack_v1*>(rf_track.get()));
        //				delete rf_track;
      }
    }
    else
    {
      if (_over_write_svtxtrackmap)
      {
        auto key = iter->first;
        ++iter;
        _trackmap->erase(key);
        continue;
      }
    }

    ++iter;
  }
  // Need to keep tracks if _do_evt_display
  if (!_do_evt_display)
  {
    rf_phgf_tracks.clear();
  }

  if (_do_eval)
  {
    fill_eval_tree(topNode);
  }
#ifdef _DEBUG_
  cout << __LINE__ << endl;
#endif
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * End
 */
int TpcPrototypeGenFitTrkFitter::End(PHCompositeNode* topNode)
{
  if (_do_eval)
  {
    if (Verbosity() >= 1)
      cout << PHWHERE << " Writing to file: " << _eval_outname << endl;
    PHTFileServer::get().cd(_eval_outname);
    _eval_tree->Write();
    _cluster_eval_tree->Write();
  }

  if (_do_evt_display)
  {
    //    if(opts[i] == 'A') drawAutoScale_ = true;
    //    if(opts[i] == 'B') drawBackward_ = true;
    //    if(opts[i] == 'D') drawDetectors_ = true;
    //    if(opts[i] == 'E') drawErrors_ = true;
    //    if(opts[i] == 'F') drawForward_ = true;
    //    if(opts[i] == 'H') drawHits_ = true;
    //    if(opts[i] == 'M') drawTrackMarkers_ = true;
    //    if(opts[i] == 'P') drawPlanes_ = true;
    //    if(opts[i] == 'S') drawScaleMan_ = true;
    //    if(opts[i] == 'T') drawTrack_ = true;
    //    if(opts[i] == 'X') drawSilent_ = true;
    //    if(opts[i] == 'G') drawGeometry_ = true;
    _fitter->getEventDisplay()->setOptions("ADEHPTG");

    _fitter->displayEvent();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * dtor
 */
TpcPrototypeGenFitTrkFitter::~TpcPrototypeGenFitTrkFitter()
{
  delete _fitter;
  delete _vertex_finder;
}

/*
 * fill_eval_tree():
 */
void TpcPrototypeGenFitTrkFitter::fill_eval_tree(PHCompositeNode* topNode)
{
  //! Make sure to reset all the TTree variables before trying to set them.
  reset_eval_variables();

  int i = 0;
  for (SvtxTrackMap::ConstIter itr = _trackmap->begin();
       itr != _trackmap->end(); ++itr)
  {
    //    new ((*_tca_trackmap)[i])(SvtxTrack_v1)(
    //        *dynamic_cast<SvtxTrack_v1*>(itr->second));

    new ((*_tca_tpctrackmap)[i])(TpcPrototypeTrack)(*(MakeTpcPrototypeTrack(itr->second)));

    i++;
  }
  _tca_ntrack = i;

  i = 0;
  if (_vertexmap)
    for (SvtxVertexMap::ConstIter itr = _vertexmap->begin();
         itr != _vertexmap->end(); ++itr)
      new ((*_tca_vertexmap)[i++])(SvtxVertex_v1)(
          *dynamic_cast<SvtxVertex_v1*>(itr->second));

  if (_trackmap_refit)
  {
    i = 0;
    for (SvtxTrackMap::ConstIter itr = _trackmap_refit->begin();
         itr != _trackmap_refit->end(); ++itr)
      new ((*_tca_trackmap_refit)[i++])(SvtxTrack_v1)(
          *dynamic_cast<SvtxTrack_v1*>(itr->second));
  }

  if (_fit_primary_tracks)
  {
    i = 0;
    for (SvtxTrackMap::ConstIter itr = _primary_trackmap->begin();
         itr != _primary_trackmap->end(); ++itr)
      new ((*_tca_primtrackmap)[i++])(SvtxTrack_v1)(
          *dynamic_cast<SvtxTrack_v1*>(itr->second));
  }

  if (_vertexmap_refit)
  {
    i = 0;
    for (SvtxVertexMap::ConstIter itr = _vertexmap_refit->begin();
         itr != _vertexmap_refit->end(); ++itr)
      new ((*_tca_vertexmap_refit)[i++])(SvtxVertex_v1)(
          *dynamic_cast<SvtxVertex_v1*>(itr->second));
  }

  _eval_tree->Fill();

  return;
}

/*
 * init_eval_tree
 */
void TpcPrototypeGenFitTrkFitter::init_eval_tree()
{
  if (!_tca_particlemap)
    _tca_particlemap = new TClonesArray("PHG4Particlev2");
  if (!_tca_vtxmap)
    _tca_vtxmap = new TClonesArray("PHG4VtxPointv1");

  if (!_tca_trackmap)
    _tca_trackmap = new TClonesArray("SvtxTrack_v1");
  if (!_tca_tpctrackmap)
    _tca_tpctrackmap = new TClonesArray("TpcPrototypeTrack");
  if (!_tca_vertexmap)
    _tca_vertexmap = new TClonesArray("SvtxVertex_v1");
  if (!_tca_trackmap_refit)
    _tca_trackmap_refit = new TClonesArray("SvtxTrack_v1");
  if (_fit_primary_tracks)
    if (!_tca_primtrackmap)
      _tca_primtrackmap = new TClonesArray("SvtxTrack_v1");
  if (!_tca_vertexmap_refit)
    _tca_vertexmap_refit = new TClonesArray("SvtxVertex_v1");

  //! create TTree
  _eval_tree = new TTree("T", "TpcPrototypeGenFitTrkFitter Evaluation");
  _eval_tree->Branch("nTrack", &_tca_ntrack, "nTrack/I");

  //  _eval_tree->Branch("PrimaryParticle", _tca_particlemap);
  //  _eval_tree->Branch("TruthVtx", _tca_vtxmap);

  _eval_tree->Branch("SvtxTrack", _tca_trackmap);
  _eval_tree->Branch("TPCTrack", _tca_tpctrackmap);
  //  _eval_tree->Branch("SvtxVertex", _tca_vertexmap);
  //  _eval_tree->Branch("SvtxTrackRefit", _tca_trackmap_refit);
  if (_fit_primary_tracks)
    _eval_tree->Branch("PrimSvtxTrack", _tca_primtrackmap);
  //  _eval_tree->Branch("SvtxVertexRefit", _tca_vertexmap_refit);

  _cluster_eval_tree = new TTree("cluster_eval", "cluster eval tree");
  _cluster_eval_tree->Branch("x", &_cluster_eval_tree_x, "x/F");
  _cluster_eval_tree->Branch("y", &_cluster_eval_tree_y, "y/F");
  _cluster_eval_tree->Branch("z", &_cluster_eval_tree_z, "z/F");
  _cluster_eval_tree->Branch("gx", &_cluster_eval_tree_gx, "gx/F");
  _cluster_eval_tree->Branch("gy", &_cluster_eval_tree_gy, "gy/F");
  _cluster_eval_tree->Branch("gz", &_cluster_eval_tree_gz, "gz/F");
}

/*
 * reset_eval_variables():
 *  Reset all the tree variables to their default values.
 *  Needs to be called at the start of every event
 */
void TpcPrototypeGenFitTrkFitter::reset_eval_variables()
{
  _tca_particlemap->Clear();
  _tca_vtxmap->Clear();

  _tca_ntrack = -1;
  _tca_trackmap->Clear();
  _tca_tpctrackmap->Clear();
  _tca_vertexmap->Clear();
  _tca_trackmap_refit->Clear();
  if (_fit_primary_tracks)
    _tca_primtrackmap->Clear();
  _tca_vertexmap_refit->Clear();

  _cluster_eval_tree_x = WILD_FLOAT;
  _cluster_eval_tree_y = WILD_FLOAT;
  _cluster_eval_tree_z = WILD_FLOAT;
  _cluster_eval_tree_gx = WILD_FLOAT;
  _cluster_eval_tree_gy = WILD_FLOAT;
  _cluster_eval_tree_gz = WILD_FLOAT;
}

int TpcPrototypeGenFitTrkFitter::CreateNodes(PHCompositeNode* topNode)
{
  // create nodes...
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVTX node
  PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst(
      "PHCompositeNode", "SVTX"));
  if (!tb_node)
  {
    tb_node = new PHCompositeNode("SVTX");
    dstNode->addNode(tb_node);
    if (Verbosity() > 0)
      cout << "SVTX node added" << endl;
  }

  if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
  {
    _trackmap_refit = new SvtxTrackMap_v1;
    PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
        _trackmap_refit, "SvtxTrackMapRefit", "PHObject");
    tb_node->addNode(tracks_node);
    if (Verbosity() > 0)
      cout << "Svtx/SvtxTrackMapRefit node added" << endl;
  }

  if (_fit_primary_tracks)
  {
    _primary_trackmap = new SvtxTrackMap_v1;
    PHIODataNode<PHObject>* primary_tracks_node =
        new PHIODataNode<PHObject>(_primary_trackmap, "PrimaryTrackMap",
                                   "PHObject");
    tb_node->addNode(primary_tracks_node);
    if (Verbosity() > 0)
      cout << "Svtx/PrimaryTrackMap node added" << endl;
  }

  if (!(_over_write_svtxvertexmap))
  {
    _vertexmap_refit = new SvtxVertexMap_v1;
    PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
        _vertexmap_refit, "SvtxVertexMapRefit", "PHObject");
    tb_node->addNode(vertexes_node);
    if (Verbosity() > 0)
      cout << "Svtx/SvtxVertexMapRefit node added" << endl;
  }
  else if (!findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap"))
  {
    _vertexmap = new SvtxVertexMap_v1;
    PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
        _vertexmap, "SvtxVertexMap", "PHObject");
    tb_node->addNode(vertexes_node);
    if (Verbosity() > 0)
      cout << "Svtx/SvtxVertexMap node added" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int TpcPrototypeGenFitTrkFitter::GetNodes(PHCompositeNode* topNode)
{
  //DST objects
  //Truth container
  //  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
  //                                                                "G4TruthInfo");
  //  if (!_truth_container && _event < 2)
  //  {
  //    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
  //         << endl;
  //    return Fun4AllReturnCodes::ABORTEVENT;
  //  }

  // Input Svtx Clusters
  //_clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_clustermap && _event < 2)
  {
    cout << PHWHERE << " TRKR_CLUSTER node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Svtx Tracks
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_trackmap && _event < 2)
  {
    cout << PHWHERE << " SvtxTrackMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Svtx Vertices
  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertexmap && _event < 2)
  {
    cout << PHWHERE << " SvtxVertexrMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Output Svtx Tracks
  if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
  {
    _trackmap_refit = findNode::getClass<SvtxTrackMap>(topNode,
                                                       "SvtxTrackMapRefit");
    if (!_trackmap_refit && _event < 2)
    {
      cout << PHWHERE << " SvtxTrackMapRefit node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // Output Primary Svtx Tracks
  if (_fit_primary_tracks)
  {
    _primary_trackmap = findNode::getClass<SvtxTrackMap>(topNode,
                                                         "PrimaryTrackMap");
    if (!_primary_trackmap && _event < 2)
    {
      cout << PHWHERE << " PrimaryTrackMap node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // Output Svtx Vertices
  if (!(_over_write_svtxvertexmap))
  {
    _vertexmap_refit = findNode::getClass<SvtxVertexMap>(topNode,
                                                         "SvtxVertexMapRefit");
    if (!_vertexmap_refit && _event < 2)
    {
      cout << PHWHERE << " SvtxVertexMapRefit node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::shared_ptr<PHGenFit::Track> TpcPrototypeGenFitTrkFitter::DisplayCluster(const TrkrCluster* cluster)
{
  assert(cluster);

  if (Verbosity() >= 1)
  {
    cout << __PRETTY_FUNCTION__ << ": process cluster: ";
    cluster->identify();
  }

  // prepare seed
  TVector3 seed_mom(100, 0, 0);
  TVector3 seed_pos(0, 0, 0);
  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = 100.;
    }
  }

  // Create measurements
  std::vector<PHGenFit::Measurement*> measurements;

  TrkrDefs::cluskey cluster_key = cluster->getClusKey();

  TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));

  seed_mom.SetPhi(pos.Phi());
  seed_mom.SetTheta(pos.Theta());

  seed_pos = pos;

  //TODO use u, v explicitly?
  TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);

  for (int radial_shift = -1; radial_shift <= 1; ++radial_shift)
  {
    const double new_radial = pos.Perp() + radial_shift * 0.1;
    TVector3 new_pos(pos);
    new_pos.SetPerp(new_radial);

    //------------------------------
    // new

    // Replace n for the silicon subsystems

    // get the trkrid
    int layer = TrkrDefs::getLayer(cluster_key);

    // end new
    //-----------------

    PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(new_pos, n,
                                                                  cluster->getRPhiError(), cluster->getZError());

    if (Verbosity() > 50)
    {
      cout << "Add meas layer " << layer << " cluskey " << cluster_key
           << endl
           << " pos.X " << pos.X() << " pos.Y " << pos.Y() << " pos.Z " << pos.Z()
           << "  n.X " << n.X() << " n.Y " << n.Y()
           << " RPhiErr " << cluster->getRPhiError()
           << " ZErr " << cluster->getZError()
           << endl;
    }
    measurements.push_back(meas);
  }

  /*!
   * mu+: -13
   * mu-: 13
   * pi+: 211
   * pi-: -211
   * e-:  11
   * e+:  -11
   */
  //TODO Add multiple TrackRep choices.
  //int pid = 211;
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(-13);
  std::shared_ptr<PHGenFit::Track> track(new PHGenFit::Track(rep, seed_pos, seed_mom,
                                                             seed_cov));

  //TODO unsorted measurements, should use sorted ones?
  track->addMeasurements(measurements);

  //  if (measurements.size()==1) return track;

  /*!
   *  Fit the track
   *  ret code 0 means 0 error or good status
   */
  if (_fitter->processTrack(track.get(), false) != 0)
  {
    if (Verbosity() >= 1)
    {
      LogWarning("Track fitting failed");
      cout << __PRETTY_FUNCTION__ << ": failed cluster build: track->getChisq() " << track->get_chi2() << " get_ndf " << track->get_ndf()
           << " mom.X " << track->get_mom().X()
           << " mom.Y " << track->get_mom().Y()
           << " mom.Z " << track->get_mom().Z()
           << endl;
    }
    //delete track;
    return nullptr;
  }

  if (Verbosity() > 50)
    cout << " track->getChisq() " << track->get_chi2() << " get_ndf " << track->get_ndf()
         << " mom.X " << track->get_mom().X()
         << " mom.Y " << track->get_mom().Y()
         << " mom.Z " << track->get_mom().Z()
         << endl;

  return track;
}

//struct CompMeasurementByR {
//  bool operator() (PHGenFit::Measurement *m1,PHGenFit::Measurement *m2) {
//	  float x1 = m1->getMeasurement()
//
//	  return (i<j);}
//} myobject;

/*
 * fit track with SvtxTrack as input seed.
 * \param intrack Input SvtxTrack
 * \param invertex Input Vertex, if fit track as a primary vertex
 */
//PHGenFit::Track* TpcPrototypeGenFitTrkFitter::ReFitTrack(PHCompositeNode *topNode, const SvtxTrack* intrack,
std::shared_ptr<PHGenFit::Track> TpcPrototypeGenFitTrkFitter::ReFitTrack(PHCompositeNode* topNode, const SvtxTrack* intrack,
                                                                         const SvtxVertex* invertex)
{
  //std::shared_ptr<PHGenFit::Track> empty_track(NULL);
  if (!intrack)
  {
    cerr << PHWHERE << " Input SvtxTrack is NULL!" << endl;
    return NULL;
  }

  // prepare seed
  TVector3 seed_mom(100, 0, 0);
  TVector3 seed_pos(0, 0, 0);
  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = 100.;
    }
  }

  // Create measurements
  std::vector<PHGenFit::Measurement*> measurements;

  //! 1000 is a arbitrary number for now
  const double vertex_chi2_over_dnf_cut = 1000;
  const double vertex_cov_element_cut = 10000;  //arbitrary cut cm*cm

  if (invertex and invertex->size_tracks() > 1 and invertex->get_chisq() / invertex->get_ndof() < vertex_chi2_over_dnf_cut)
  {
    TVector3 pos(invertex->get_x(), invertex->get_y(), invertex->get_z());
    TMatrixDSym cov(3);
    cov.Zero();
    bool is_vertex_cov_sane = true;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
      {
        cov(i, j) = invertex->get_error(i, j);

        if (i == j)
        {
          if (!(invertex->get_error(i, j) > 0 and invertex->get_error(i, j) < vertex_cov_element_cut))
            is_vertex_cov_sane = false;
        }
      }

    if (is_vertex_cov_sane)
    {
      PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(
          pos, cov);
      measurements.push_back(meas);
      if (Verbosity() >= 2)
      {
        meas->getMeasurement()->Print();
      }
    }
  }

  // sort clusters with radius before fitting
  if (Verbosity() > 20) intrack->identify();
  std::map<float, TrkrDefs::cluskey> m_r_cluster_id;
  for (auto iter = intrack->begin_cluster_keys();
       iter != intrack->end_cluster_keys(); ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    TrkrCluster* cluster = _clustermap->findCluster(cluster_key);
    float x = cluster->getPosition(0);
    float y = cluster->getPosition(1);
    float r = sqrt(x * x + y * y);
    m_r_cluster_id.insert(std::pair<float, TrkrDefs::cluskey>(r, cluster_key));
    int layer_out = TrkrDefs::getLayer(cluster_key);
    if (Verbosity() > 20) cout << "    Layer " << layer_out << " cluster " << cluster_key << " radius " << r << endl;
  }

  for (auto iter = m_r_cluster_id.begin();
       iter != m_r_cluster_id.end();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = iter->second;
    TrkrCluster* cluster = _clustermap->findCluster(cluster_key);
    if (!cluster)
    {
      LogError("No cluster Found!");
      continue;
    }

    TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));

    seed_mom.SetPhi(pos.Phi());
    seed_mom.SetTheta(pos.Theta());

    //TODO use u, v explicitly?
    TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);

    //------------------------------
    // new

    // Replace n for the silicon subsystems

    // get the trkrid
    unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
    int layer = TrkrDefs::getLayer(cluster_key);

    if (trkrid == TrkrDefs::mvtxId)
    {
      assert(0);
    }
    else if (trkrid == TrkrDefs::inttId)
    {
      assert(0);
    }
    // end new
    //-----------------

    PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
                                                                  cluster->getRPhiError(), cluster->getZError());

    if (Verbosity() > 50)
    {
      cout << "Add meas layer " << layer << " cluskey " << cluster_key
           << endl
           << " pos.X " << pos.X() << " pos.Y " << pos.Y() << " pos.Z " << pos.Z()
           << "  n.X " << n.X() << " n.Y " << n.Y()
           << " RPhiErr " << cluster->getRPhiError()
           << " ZErr " << cluster->getZError()
           << endl;
    }
    measurements.push_back(meas);
  }

  /*!
	 * mu+:	-13
	 * mu-:	13
	 * pi+:	211
	 * pi-:	-211
	 * e-:	11
	 * e+:	-11
	 */
  //TODO Add multiple TrackRep choices.
  //int pid = 211;
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
  std::shared_ptr<PHGenFit::Track> track(new PHGenFit::Track(rep, seed_pos, seed_mom,
                                                             seed_cov));
  track->addMeasurements(measurements);

  //  if (measurements.size()==1) return track;

  /*!
	 *  Fit the track
	 *  ret code 0 means 0 error or good status
	 */
  if (_fitter->processTrack(track.get(), false) != 0)
  {
    if (Verbosity() >= 1)
    {
      LogWarning("Track fitting failed");
      cout << __PRETTY_FUNCTION__ << " track->getChisq() " << track->get_chi2() << " get_ndf " << track->get_ndf()
           << " mom.X " << track->get_mom().X()
           << " mom.Y " << track->get_mom().Y()
           << " mom.Z " << track->get_mom().Z()
           << endl;
    }
    //delete track;
    return nullptr;
  }

  if (Verbosity() > 50)
    cout << " track->getChisq() " << track->get_chi2() << " get_ndf " << track->get_ndf()
         << " mom.X " << track->get_mom().X()
         << " mom.Y " << track->get_mom().Y()
         << " mom.Z " << track->get_mom().Z()
         << endl;

  return track;
}

shared_ptr<TpcPrototypeTrack> TpcPrototypeGenFitTrkFitter::MakeTpcPrototypeTrack(const SvtxTrack* svtxtrack)
{
  assert(svtxtrack);
  if (Verbosity() >= 1)
  {
    cout << __PRETTY_FUNCTION__ << " refit ";
    svtxtrack->identify();
  }

  shared_ptr<TpcPrototypeTrack> track(new TpcPrototypeTrack);
  track->event = _event;
  track->trackID = svtxtrack->get_id();
  track->chisq = svtxtrack->get_chisq();
  track->ndf = svtxtrack->get_ndf();
  track->px = svtxtrack->get_px();
  track->py = svtxtrack->get_py();
  track->pz = svtxtrack->get_pz();
  track->x = svtxtrack->get_x();
  track->y = svtxtrack->get_y();
  track->z = svtxtrack->get_z();

  track->nCluster = 0;

  // sort intput cluters to one per layer
  assert(_clustermap);
  vector<TrkrCluster*> clusterLayer(track->nLayer, nullptr);
  for (auto iter = svtxtrack->begin_cluster_keys();
       iter != svtxtrack->end_cluster_keys(); ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    TrkrCluster* cluster = _clustermap->findCluster(cluster_key);
    int layer = TrkrDefs::getLayer(cluster_key);

    if (Verbosity())
    {
      cout << __PRETTY_FUNCTION__ << " - layer sorting cluster ";
      cluster->identify();
    }

    assert(layer >= 0);
    assert(layer < track->nLayer);
    assert(cluster);

    if (clusterLayer[layer])
    {
      cout << __PRETTY_FUNCTION__ << " -WARNING- more than one cluster at layer " << layer << ": " << endl;
      clusterLayer[layer]->identify();
      cluster->identify();
    }
    else
    {
      clusterLayer[layer] = cluster;
      track->nCluster += 1;
    }
  }

  if (track->nCluster < 4)
    return track;

  // refit by "excluding" each cluster
  // prepare seed
  TVector3 seed_mom(svtxtrack->get_px(), svtxtrack->get_py(), svtxtrack->get_pz());
  seed_mom.SetMag(120);

  TVector3 seed_pos(svtxtrack->get_x(), svtxtrack->get_y(), svtxtrack->get_z());
  TMatrixDSym seed_cov(6);
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = 100.;
    }
  }
  for (int layerStudy = 0; layerStudy < track->nLayer; ++layerStudy)
  {
    //scale uncertainty for the cluster in study with this factor
    const static double errorScaleFactor = 100;

    if (clusterLayer[layerStudy] == nullptr) continue;

    // Create measurements
    std::vector<PHGenFit::Measurement*> measurements;

    int indexStudy = -1;
    int currentIndex = 0;
    for (const auto cluster : clusterLayer)
    {
      if (!cluster) continue;

      assert(cluster);
      TrkrDefs::cluskey cluster_key = cluster->getClusKey();
      int layer = TrkrDefs::getLayer(cluster_key);

      const double scale_this_layer = (layer == layerStudy) ? errorScaleFactor : 1;

      TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));

      //      seed_mom.SetPhi(pos.Phi());
      //      seed_mom.SetTheta(pos.Theta());

      //TODO use u, v explicitly?
      TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);

      PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
                                                                    cluster->getRPhiError() * scale_this_layer,
                                                                    cluster->getZError() * scale_this_layer);

      measurements.push_back(meas);

      if (layer == layerStudy) indexStudy = currentIndex;
      ++currentIndex;
    }  //    for (const auto cluster : clusterLayer)
    assert(indexStudy >= 0);
    assert(indexStudy < track->nLayer);

    // do the fit
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
    std::shared_ptr<PHGenFit::Track> phgf_track(new PHGenFit::Track(rep, seed_pos, seed_mom,
                                                                    seed_cov));
    phgf_track->addMeasurements(measurements);
    if (_fitter->processTrack(phgf_track.get(), false) != 0)
    {
      if (Verbosity() >= 1)
      {
        LogWarning("Track fitting failed");
        cout << __PRETTY_FUNCTION__ << " track->getChisq() " << phgf_track->get_chi2() << " get_ndf " << phgf_track->get_ndf()
             << " mom.X " << phgf_track->get_mom().X()
             << " mom.Y " << phgf_track->get_mom().Y()
             << " mom.Z " << phgf_track->get_mom().Z()
             << endl;
      }
      //delete track;
      continue;
    }

    //propagate to state
    {
      std::shared_ptr<const genfit::MeasuredStateOnPlane> gf_state = NULL;
      try
      {
        const genfit::Track* gftrack = phgf_track->getGenFitTrack();
        assert(gftrack);
        const genfit::AbsTrackRep* rep = gftrack->getCardinalRep();
        assert(rep);
        genfit::TrackPoint* trpoint = gftrack->getPointWithMeasurementAndFitterInfo(indexStudy, gftrack->getCardinalRep());
        assert(trpoint);
        genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(trpoint->getFitterInfo(rep));
        assert(kfi);
        //gf_state = std::shared_ptr <genfit::MeasuredStateOnPlane> (const_cast<genfit::MeasuredStateOnPlane*> (&(kfi->getFittedState(true))));
        const genfit::MeasuredStateOnPlane* temp_state = &(kfi->getFittedState(true));
        gf_state = std::shared_ptr<genfit::MeasuredStateOnPlane>(new genfit::MeasuredStateOnPlane(*temp_state));
      }
      catch (...)
      {
        if (Verbosity() > 1)
          LogWarning("Exrapolation failed!");
        continue;
      }
      if (!gf_state)
      {
        if (Verbosity() > 1)
          LogWarning("Exrapolation failed!");
        continue;
      }
      TVector3 extra_pos(gf_state->getPos());

      const TrkrCluster* cluster(clusterLayer[layerStudy]);
      assert(cluster);
      TVector3 cluster_pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));

      TVector3 n_dir(cluster->getPosition(0), cluster->getPosition(1), 0);
      n_dir.SetMag(1);
      TVector3 z_dir(0, 0, 1);
      TVector3 azimuth_dir(z_dir.Cross(n_dir));
      TVector3 pos_diff = cluster_pos - extra_pos;
      const double n_residual = pos_diff.Dot(n_dir);
      const double z_residual = pos_diff.Dot(z_dir);
      const double azimuth_residual = pos_diff.Dot(azimuth_dir);

      if (Verbosity() >= 1)
      {
        cout << __PRETTY_FUNCTION__;
        cout << " - layer " << layerStudy << " at index " << indexStudy;
        cout << " cluster @ " << cluster_pos.x() << ", " << cluster_pos.y() << ", " << cluster_pos.z();
        cout << " extrapolate to " << extra_pos.x() << ", " << extra_pos.y() << ", " << extra_pos.z();
        cout << ", n_residual = " << n_residual;
        cout << ", z_residual = " << z_residual;
        cout << ", azimuth_residual = " << azimuth_residual;
        cout << endl;

        cout << "Cluster data which provide much more detailed information on the raw signals: "<<endl;

        // in case TpcPrototypeCluster specific functionality is needed, first to a conversion and check
        const TpcPrototypeCluster * prototype_cluster = dynamic_cast<const TpcPrototypeCluster *>(cluster);
        assert(prototype_cluster);

        prototype_cluster->identify();
      }
      assert(abs(n_residual) < 1e-4);  //same layer check

      (track->clusterKey)[layerStudy] = cluster->getClusKey();
      (track->clusterlayer)[layerStudy] = layerStudy;
      (track->clusterid)[layerStudy] = TrkrDefs::getClusIndex(cluster->getClusKey());

      (track->clusterX)[layerStudy] = cluster->getPosition(0);
      (track->clusterY)[layerStudy] = cluster->getPosition(1);
      (track->clusterZ)[layerStudy] = cluster->getPosition(2);
      (track->clusterE)[layerStudy] = cluster->getAdc();
      (track->clusterSizePhi)[layerStudy] = cluster->getPhiSize();
      (track->clusterResidualPhi)[layerStudy] = azimuth_residual;
      (track->clusterProjectionPhi)[layerStudy] = extra_pos.Phi();
      (track->clusterResidualZ)[layerStudy] = z_residual;
    }  //propagate to state

  }  // refit by "excluding" each cluster  for (int layer = 0; layer < track->nLayer; ++layer)

  return track;
}

/*
 * Make SvtxTrack from PHGenFit::Track and SvtxTrack
 */
//SvtxTrack* TpcPrototypeGenFitTrkFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
std::shared_ptr<SvtxTrack> TpcPrototypeGenFitTrkFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
                                                                      const std::shared_ptr<PHGenFit::Track>& phgf_track, const SvtxVertex* vertex)
{
  double chi2 = phgf_track->get_chi2();
  double ndf = phgf_track->get_ndf();

  TVector3 vertex_position(0, 0, 0);
  TMatrixF vertex_cov(3, 3);
  double dvr2 = 0;
  double dvz2 = 0;

  //  if (_use_truth_vertex)
  //  {
  //    PHG4VtxPoint* first_point = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
  //    vertex_position.SetXYZ(first_point->get_x(), first_point->get_y(), first_point->get_z());
  //    if (Verbosity() > 1)
  //    {
  //      cout << PHWHERE << "Using: truth vertex: {" << vertex_position.X() << ", " << vertex_position.Y() << ", " << vertex_position.Z() << "} " << endl;
  //    }
  //  }
  //  else if (vertex)
  //  {
  //    vertex_position.SetXYZ(vertex->get_x(), vertex->get_y(),
  //                           vertex->get_z());
  //    dvr2 = vertex->get_error(0, 0) + vertex->get_error(1, 1);
  //    dvz2 = vertex->get_error(2, 2);
  //
  //    for (int i = 0; i < 3; i++)
  //      for (int j = 0; j < 3; j++)
  //        vertex_cov[i][j] = vertex->get_error(i, j);
  //  }

  //genfit::MeasuredStateOnPlane* gf_state_beam_line_ca = NULL;
  std::shared_ptr<genfit::MeasuredStateOnPlane> gf_state_beam_line_ca = NULL;
  try
  {
    gf_state_beam_line_ca = std::shared_ptr<genfit::MeasuredStateOnPlane>(phgf_track->extrapolateToLine(vertex_position,
                                                                                                        TVector3(0., 0., 1.)));
  }
  catch (...)
  {
    if (Verbosity() >= 2)
      LogWarning("extrapolateToLine failed!");
  }
  if (!gf_state_beam_line_ca) return NULL;

  /*!
	 *  1/p, u'/z', v'/z', u, v
	 *  u is defined as momentum X beam line at POCA of the beam line
	 *  v is alone the beam line
	 *  so u is the dca2d direction
	 */

  double u = gf_state_beam_line_ca->getState()[3];
  double v = gf_state_beam_line_ca->getState()[4];

  double du2 = gf_state_beam_line_ca->getCov()[3][3];
  double dv2 = gf_state_beam_line_ca->getCov()[4][4];
  //cout << PHWHERE << "        u " << u << " v " << v << " du2 " << du2 << " dv2 " << dv2 << " dvr2 " << dvr2 << endl;
  //delete gf_state_beam_line_ca;

  //const SvtxTrack_v1* temp_track = static_cast<const SvtxTrack_v1*> (svtx_track);
  //	SvtxTrack_v1* out_track = new SvtxTrack_v1(
  //			*static_cast<const SvtxTrack_v1*>(svtx_track));
  std::shared_ptr<SvtxTrack_v1> out_track = std::shared_ptr<SvtxTrack_v1>(new SvtxTrack_v1(*static_cast<const SvtxTrack_v1*>(svtx_track)));

  out_track->set_dca2d(u);
  out_track->set_dca2d_error(sqrt(du2 + dvr2));

  std::shared_ptr<genfit::MeasuredStateOnPlane> gf_state_vertex_ca = NULL;
  try
  {
    gf_state_vertex_ca = std::shared_ptr<genfit::MeasuredStateOnPlane>(phgf_track->extrapolateToPoint(vertex_position));
  }
  catch (...)
  {
    if (Verbosity() >= 2)
      LogWarning("extrapolateToPoint failed!");
  }
  if (!gf_state_vertex_ca)
  {
    //delete out_track;
    return NULL;
  }

  TVector3 mom = gf_state_vertex_ca->getMom();
  TVector3 pos = gf_state_vertex_ca->getPos();
  TMatrixDSym cov = gf_state_vertex_ca->get6DCov();

  //	genfit::MeasuredStateOnPlane* gf_state_vertex_ca =
  //			phgf_track->extrapolateToLine(vertex_position,
  //					TVector3(0., 0., 1.));

  u = gf_state_vertex_ca->getState()[3];
  v = gf_state_vertex_ca->getState()[4];

  du2 = gf_state_vertex_ca->getCov()[3][3];
  dv2 = gf_state_vertex_ca->getCov()[4][4];

  double dca3d = sqrt(u * u + v * v);
  double dca3d_error = sqrt(du2 + dv2 + dvr2 + dvz2);

  out_track->set_dca(dca3d);
  out_track->set_dca_error(dca3d_error);

  //
  // in: X, Y, Z; out; r: n X Z, Z X r, Z

  float dca3d_xy = NAN;
  float dca3d_z = NAN;
  float dca3d_xy_error = NAN;
  float dca3d_z_error = NAN;

  try
  {
    TMatrixF pos_in(3, 1);
    TMatrixF cov_in(3, 3);
    TMatrixF pos_out(3, 1);
    TMatrixF cov_out(3, 3);

    TVectorD state6(6);      // pos(3), mom(3)
    TMatrixDSym cov6(6, 6);  //

    gf_state_vertex_ca->get6DStateCov(state6, cov6);

    TVector3 vn(state6[3], state6[4], state6[5]);

    // mean of two multivariate gaussians Pos - Vertex
    pos_in[0][0] = state6[0] - vertex_position.X();
    pos_in[1][0] = state6[1] - vertex_position.Y();
    pos_in[2][0] = state6[2] - vertex_position.Z();

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        cov_in[i][j] = cov6[i][j] + vertex_cov[i][j];
      }
    }

    // vn is momentum vector, pos_in is position vector (of what?)
    pos_cov_XYZ_to_RZ(vn, pos_in, cov_in, pos_out, cov_out);

    if (Verbosity() > 50)
    {
      cout << " vn.X " << vn.X() << " vn.Y " << vn.Y() << " vn.Z " << vn.Z() << endl;
      cout << " pos_in.X " << pos_in[0][0] << " pos_in.Y " << pos_in[1][0] << " pos_in.Z " << pos_in[2][0] << endl;
      cout << " pos_out.X " << pos_out[0][0] << " pos_out.Y " << pos_out[1][0] << " pos_out.Z " << pos_out[2][0] << endl;
    }

    dca3d_xy = pos_out[0][0];
    dca3d_z = pos_out[2][0];
    dca3d_xy_error = sqrt(cov_out[0][0]);
    dca3d_z_error = sqrt(cov_out[2][2]);
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("DCA calculationfailed!");
  }

  out_track->set_dca3d_xy(dca3d_xy);
  out_track->set_dca3d_z(dca3d_z);
  out_track->set_dca3d_xy_error(dca3d_xy_error);
  out_track->set_dca3d_z_error(dca3d_z_error);

  //if(gf_state_vertex_ca) delete gf_state_vertex_ca;

  out_track->set_chisq(chi2);
  out_track->set_ndf(ndf);
  out_track->set_charge(phgf_track->get_charge());

  out_track->set_px(mom.Px());
  out_track->set_py(mom.Py());
  out_track->set_pz(mom.Pz());

  out_track->set_x(pos.X());
  out_track->set_y(pos.Y());
  out_track->set_z(pos.Z());

  for (int i = 0; i < 6; i++)
  {
    for (int j = i; j < 6; j++)
    {
      out_track->set_error(i, j, cov[i][j]);
    }
  }

  const genfit::Track* gftrack = phgf_track->getGenFitTrack();
  const genfit::AbsTrackRep* rep = gftrack->getCardinalRep();
  for (unsigned int id = 0; id < gftrack->getNumPointsWithMeasurement(); ++id)
  {
    genfit::TrackPoint* trpoint = gftrack->getPointWithMeasurementAndFitterInfo(id, gftrack->getCardinalRep());

    if (!trpoint)
    {
      if (Verbosity() > 1)
        LogWarning("!trpoint");
      continue;
    }

    genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(trpoint->getFitterInfo(rep));
    if (!kfi)
    {
      if (Verbosity() > 1)
        LogWarning("!kfi");
      continue;
    }

    std::shared_ptr<const genfit::MeasuredStateOnPlane> gf_state = NULL;
    try
    {
      //gf_state = std::shared_ptr <genfit::MeasuredStateOnPlane> (const_cast<genfit::MeasuredStateOnPlane*> (&(kfi->getFittedState(true))));
      const genfit::MeasuredStateOnPlane* temp_state = &(kfi->getFittedState(true));
      gf_state = std::shared_ptr<genfit::MeasuredStateOnPlane>(new genfit::MeasuredStateOnPlane(*temp_state));
    }
    catch (...)
    {
      if (Verbosity() > 1)
        LogWarning("Exrapolation failed!");
    }
    if (!gf_state)
    {
      if (Verbosity() > 1)
        LogWarning("Exrapolation failed!");
      continue;
    }
    genfit::MeasuredStateOnPlane temp;
    float pathlength = -phgf_track->extrapolateToPoint(temp, vertex_position, id);

    std::shared_ptr<SvtxTrackState> state = std::shared_ptr<SvtxTrackState>(new SvtxTrackState_v1(pathlength));
    state->set_x(gf_state->getPos().x());
    state->set_y(gf_state->getPos().y());
    state->set_z(gf_state->getPos().z());

    state->set_px(gf_state->getMom().x());
    state->set_py(gf_state->getMom().y());
    state->set_pz(gf_state->getMom().z());

    //gf_state->getCov().Print();

    for (int i = 0; i < 6; i++)
    {
      for (int j = i; j < 6; j++)
      {
        state->set_error(i, j, gf_state->get6DCov()[i][j]);
      }
    }

    out_track->insert_state(state.get());
  }

  return out_track;
}

/*
 * Fill SvtxVertexMap from GFRaveVertexes and Tracks
 */
bool TpcPrototypeGenFitTrkFitter::FillSvtxVertexMap(
    const std::vector<genfit::GFRaveVertex*>& rave_vertices,
    const std::vector<genfit::Track*>& gf_tracks)
{
  if (_over_write_svtxvertexmap)
  {
    _vertexmap->clear();
  }

  for (unsigned int ivtx = 0; ivtx < rave_vertices.size(); ++ivtx)
  {
    genfit::GFRaveVertex* rave_vtx = rave_vertices[ivtx];

    if (!rave_vtx)
    {
      cerr << PHWHERE << endl;
      return false;
    }

    std::shared_ptr<SvtxVertex> svtx_vtx(new SvtxVertex_v1());

    svtx_vtx->set_chisq(rave_vtx->getChi2());
    svtx_vtx->set_ndof(rave_vtx->getNdf());
    svtx_vtx->set_position(0, rave_vtx->getPos().X());
    svtx_vtx->set_position(1, rave_vtx->getPos().Y());
    svtx_vtx->set_position(2, rave_vtx->getPos().Z());

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        svtx_vtx->set_error(i, j, rave_vtx->getCov()[i][j]);

    for (unsigned int i = 0; i < rave_vtx->getNTracks(); i++)
    {
      //TODO Assume id's are sync'ed between _trackmap_refit and gf_tracks, need to change?
      const genfit::Track* rave_track =
          rave_vtx->getParameters(i)->getTrack();
      for (unsigned int j = 0; j < gf_tracks.size(); j++)
      {
        if (rave_track == gf_tracks[j])
          svtx_vtx->insert_track(j);
      }
    }

    if (_over_write_svtxvertexmap)
    {
      if (_vertexmap)
      {
        _vertexmap->insert_clone(svtx_vtx.get());
      }
      else
      {
        LogError("!_vertexmap");
      }
    }
    else
    {
      if (_vertexmap_refit)
      {
        _vertexmap_refit->insert_clone(svtx_vtx.get());
      }
      else
      {
        LogError("!_vertexmap_refit");
      }
    }

    //		if (Verbosity() >= 2) {
    //			cout << PHWHERE << endl;
    //			svtx_vtx->Print();
    //			_vertexmap_refit->Print();
    //		}

    //delete svtx_vtx;
  }  //loop over RAVE vertices

  return true;
}

bool TpcPrototypeGenFitTrkFitter::pos_cov_uvn_to_rz(const TVector3& u, const TVector3& v,
                                                    const TVector3& n, const TMatrixF& pos_in, const TMatrixF& cov_in,
                                                    TMatrixF& pos_out, TMatrixF& cov_out) const
{
  if (pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
    return false;
  }

  if (cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
    return false;
  }

  TVector3 Z_uvn(u.Z(), v.Z(), n.Z());
  TVector3 up_uvn = TVector3(0., 0., 1.).Cross(Z_uvn);  // n_uvn X Z_uvn

  if (up_uvn.Mag() < 0.00001)
  {
    if (Verbosity() > 0) LogWarning("n is parallel to z");
    return false;
  }

  // R: rotation from u,v,n to n X Z, nX(nXZ), n
  TMatrixF R(3, 3);
  TMatrixF R_T(3, 3);

  try
  {
    // rotate u along z to up
    float phi = -atan2(up_uvn.Y(), up_uvn.X());
    R[0][0] = cos(phi);
    R[0][1] = -sin(phi);
    R[0][2] = 0;
    R[1][0] = sin(phi);
    R[1][1] = cos(phi);
    R[1][2] = 0;
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;

    R_T.Transpose(R);
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("Can't get rotation matrix");

    return false;
  }

  pos_out.ResizeTo(3, 1);
  cov_out.ResizeTo(3, 3);

  pos_out = R * pos_in;
  cov_out = R * cov_in * R_T;

  return true;
}

bool TpcPrototypeGenFitTrkFitter::get_vertex_error_uvn(const TVector3& u,
                                                       const TVector3& v, const TVector3& n, const TMatrixF& cov_in,
                                                       TMatrixF& cov_out) const
{
  /*!
	 * Get matrix that rotates frame (u,v,n) to (x,y,z)
	 * or the matrix that rotates vector defined in (x,y,z) to defined (u,v,n)
	 */

  TMatrixF R = get_rotation_matrix(u, v, n);
  //
  //	LogDebug("TpcPrototypeGenFitTrkFitter::get_vertex_error_uvn::R = ");
  //	R.Print();
  //	cout<<"R.Determinant() = "<<R.Determinant()<<"\n";

  if (!(abs(R.Determinant() - 1) < 0.01))
  {
    if (Verbosity() > 0)
      LogWarning("!(abs(R.Determinant()-1)<0.0001)");
    return false;
  }

  if (R.GetNcols() != 3 || R.GetNrows() != 3)
  {
    if (Verbosity() > 0)
      LogWarning("R.GetNcols() != 3 || R.GetNrows() != 3");
    return false;
  }

  if (cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3)
  {
    if (Verbosity() > 0)
      LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
    return false;
  }

  TMatrixF R_T(3, 3);

  R_T.Transpose(R);

  cov_out.ResizeTo(3, 3);

  cov_out = R * cov_in * R_T;

  return true;
}

bool TpcPrototypeGenFitTrkFitter::pos_cov_XYZ_to_RZ(
    const TVector3& n, const TMatrixF& pos_in, const TMatrixF& cov_in,
    TMatrixF& pos_out, TMatrixF& cov_out) const
{
  if (pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
    return false;
  }

  if (cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3)
  {
    if (Verbosity() > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
    return false;
  }

  // produces a vector perpendicular to both the momentum vector and beam line - i.e. in the direction of the dca_xy
  // only the angle of r will be used, not the magnitude
  TVector3 r = n.Cross(TVector3(0., 0., 1.));
  if (r.Mag() < 0.00001)
  {
    if (Verbosity() > 0) LogWarning("n is parallel to z");
    return false;
  }

  // R: rotation from u,v,n to n X Z, nX(nXZ), n
  TMatrixF R(3, 3);
  TMatrixF R_T(3, 3);

  try
  {
    // rotate u along z to up
    float phi = -atan2(r.Y(), r.X());
    R[0][0] = cos(phi);
    R[0][1] = -sin(phi);
    R[0][2] = 0;
    R[1][0] = sin(phi);
    R[1][1] = cos(phi);
    R[1][2] = 0;
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;

    R_T.Transpose(R);
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("Can't get rotation matrix");

    return false;
  }

  pos_out.ResizeTo(3, 1);
  cov_out.ResizeTo(3, 3);

  pos_out = R * pos_in;
  cov_out = R * cov_in * R_T;

  return true;
}

/*!
 * Get 3D Rotation Matrix that rotates frame (x,y,z) to (x',y',z')
 * Default rotate local to global, or rotate vector in global to local representation
 */
TMatrixF TpcPrototypeGenFitTrkFitter::get_rotation_matrix(const TVector3 x,
                                                          const TVector3 y, const TVector3 z, const TVector3 xp, const TVector3 yp,
                                                          const TVector3 zp) const
{
  TMatrixF R(3, 3);

  TVector3 xu = x.Unit();
  TVector3 yu = y.Unit();
  TVector3 zu = z.Unit();

  const float max_diff = 0.01;

  if (!(
          abs(xu * yu) < max_diff and
          abs(xu * zu) < max_diff and
          abs(yu * zu) < max_diff))
  {
    if (Verbosity() > 0)
      LogWarning("input frame error!");
    return R;
  }

  TVector3 xpu = xp.Unit();
  TVector3 ypu = yp.Unit();
  TVector3 zpu = zp.Unit();

  if (!(
          abs(xpu * ypu) < max_diff and
          abs(xpu * zpu) < max_diff and
          abs(ypu * zpu) < max_diff))
  {
    if (Verbosity() > 0)
      LogWarning("output frame error!");
    return R;
  }

  /*!
	 * Decompose x',y',z' in x,y,z and call them u,v,n
	 * Then the question will be rotate the standard X,Y,Z to u,v,n
	 */

  TVector3 u(xpu.Dot(xu), xpu.Dot(yu), xpu.Dot(zu));
  TVector3 v(ypu.Dot(xu), ypu.Dot(yu), ypu.Dot(zu));
  TVector3 n(zpu.Dot(xu), zpu.Dot(yu), zpu.Dot(zu));

  try
  {
    std::shared_ptr<TRotation> rotation(new TRotation());
    //TRotation *rotation = new TRotation();

    //! Rotation that rotate standard (X, Y, Z) to (u, v, n)
    rotation->RotateAxes(u, v, n);

    R[0][0] = rotation->XX();
    R[0][1] = rotation->XY();
    R[0][2] = rotation->XZ();
    R[1][0] = rotation->YX();
    R[1][1] = rotation->YY();
    R[1][2] = rotation->YZ();
    R[2][0] = rotation->ZX();
    R[2][1] = rotation->ZY();
    R[2][2] = rotation->ZZ();
    //
    //		LogDebug("TpcPrototypeGenFitTrkFitter::get_rotation_matrix: TRotation:");
    //		R.Print();
    //		cout<<"R.Determinant() = "<<R.Determinant()<<"\n";

    //delete rotation;

    //		TMatrixF ROT1(3, 3);
    //		TMatrixF ROT2(3, 3);
    //		TMatrixF ROT3(3, 3);
    //
    //		// rotate n along z to xz plane
    //		float phi = -atan2(n.Y(), n.X());
    //		ROT1[0][0] = cos(phi);
    //		ROT1[0][1] = -sin(phi);
    //		ROT1[0][2] = 0;
    //		ROT1[1][0] = sin(phi);
    //		ROT1[1][1] = cos(phi);
    //		ROT1[1][2] = 0;
    //		ROT1[2][0] = 0;
    //		ROT1[2][1] = 0;
    //		ROT1[2][2] = 1;
    //
    //		// rotate n along y to z
    //		TVector3 n1(n);
    //		n1.RotateZ(phi);
    //		float theta = -atan2(n1.X(), n1.Z());
    //		ROT2[0][0] = cos(theta);
    //		ROT2[0][1] = 0;
    //		ROT2[0][2] = sin(theta);
    //		ROT2[1][0] = 0;
    //		ROT2[1][1] = 1;
    //		ROT2[1][2] = 0;
    //		ROT2[2][0] = -sin(theta);
    //		ROT2[2][1] = 0;
    //		ROT2[2][2] = cos(theta);
    //
    //		// rotate u along z to x
    //		TVector3 u2(u);
    //		u2.RotateZ(phi);
    //		u2.RotateY(theta);
    //		float phip = -atan2(u2.Y(), u2.X());
    //		ROT3[0][0] = cos(phip);
    //		ROT3[0][1] = -sin(phip);
    //		ROT3[0][2] = 0;
    //		ROT3[1][0] = sin(phip);
    //		ROT3[1][1] = cos(phip);
    //		ROT3[1][2] = 0;
    //		ROT3[2][0] = 0;
    //		ROT3[2][1] = 0;
    //		ROT3[2][2] = 1;
    //
    //		// R: rotation from u,v,n to (z X n), v', z
    //		R = ROT3 * ROT2 * ROT1;
    //
    //		R.Invert();
    //		LogDebug("TpcPrototypeGenFitTrkFitter::get_rotation_matrix: Home Brew:");
    //		R.Print();
    //		cout<<"R.Determinant() = "<<R.Determinant()<<"\n";
  }
  catch (...)
  {
    if (Verbosity() > 0)
      LogWarning("Can't get rotation matrix");

    return R;
  }

  return R;
}

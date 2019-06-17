/*!
 *  \file		TpcPrototypeGenFitTrkFinder.C
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "TpcPrototypeGenFitTrkFinder.h"
#include "TpcPrototypeTrack.h"

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack_v1.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackState.h>  // for SvtxTrackState
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>  // for SvtxVertexMap
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex_v1.h>

#include <mvtx/MvtxDefs.h>

#include <intt/InttDefs.h>

#include <g4detectors/PHG4CylinderGeom.h>  // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

//
#include <intt/CylinderGeomIntt.h>

#include <mvtx/CylinderGeom_Mvtx.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>  // for PHG4VtxPoint
#include <g4main/PHG4VtxPointv1.h>

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
#include <TMath.h>           // for ATan2
#include <TMatrixDSymfwd.h>  // for TMatrixDSym
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

/*
 * Constructor
 */
TpcPrototypeGenFitTrkFinder::TpcPrototypeGenFitTrkFinder(const string& name, int layers)
  : SubsysReco(name)
  , _fitter(NULL)
  , _track_fitting_alg_name("KalmanFitter")
  , nLayer(layers)
  , _primary_pid_guess(2212)
  , rphiWindow(2)
  , ZWindow(2)
  , _clustermap(NULL)
  , _trackmap(NULL)
  , _vertexmap(NULL)
  , _do_eval(false)
  , _eval_outname("TpcPrototypeGenFitTrkFinder_eval.root")
  , _eval_tree(NULL)
  , _tca_ntrack(-1)
  , _tca_trackmap(NULL)
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
int TpcPrototypeGenFitTrkFinder::Init(PHCompositeNode* topNode)
{
  //	CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Init run
 */
int TpcPrototypeGenFitTrkFinder::InitRun(PHCompositeNode* topNode)
{
  CreateNodes(topNode);

  TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
  PHField* field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

  //_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
  _fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
                                          field, _track_fitting_alg_name,
                                          "RKTrackRep", _do_evt_display);
  _fitter->set_verbosity(Verbosity());

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
int TpcPrototypeGenFitTrkFinder::process_event(PHCompositeNode* topNode)
{
  _event++;

  if (Verbosity() > 1)
    std::cout << PHWHERE << "Events processed: " << _event << std::endl;
  //	if (_event % 1000 == 0)
  //		cout << PHWHERE << "Events processed: " << _event << endl;

  //cout << "Start TpcPrototypeGenFitTrkFinder::process_event" << endl;

  GetNodes(topNode);

  vector<std::shared_ptr<PHGenFit::Track> > rf_phgf_tracks;
  //  rf_phgf_tracks.clear();

  // Need to keep tracks if _do_evt_display
  if (!_do_evt_display)
  {
    rf_phgf_tracks.clear();
  }

  typedef vector<const TrkrCluster*> tracklet_t;
  set<tracklet_t> tracklets;

  //progressive search track finding
  for (int layer = 0; layer < nLayer; ++layer)
  {
    set<tracklet_t> new_tracklets(tracklets);
    auto cluster_range = _clustermap->getClusters(TrkrDefs::tpcId, layer);

    for (auto iter = cluster_range.first; iter != cluster_range.second; ++iter)
    {
      const TrkrCluster* cluster = iter->second;
      assert(cluster);

      TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));
      TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);

      TVector3 n_dir(cluster->getPosition(0), cluster->getPosition(1), 0);
      n_dir.SetMag(1);
      TVector3 z_dir(0, 0, 1);
      TVector3 azimuth_dir(z_dir.Cross(n_dir));

      bool matched = false;
      for (const auto tracklet : tracklets)
      {
        assert(tracklet.size() >= 1);
        const TrkrCluster* last_cluster = tracklet.back();
        assert(last_cluster);

        TVector3 last_pos(last_cluster->getPosition(0), last_cluster->getPosition(1), last_cluster->getPosition(2));

        TVector3 pos_diff = pos - last_pos;
        const double n_residual = pos_diff.Dot(n_dir);
        const double z_residual = pos_diff.Dot(z_dir);
        const double azimuth_residual = pos_diff.Dot(azimuth_dir);

        assert(n_residual > 0);
        const double z_diff = z_residual / n_residual;
        const double azimuth_diff = azimuth_residual / n_residual;
        if (abs(azimuth_diff) < rphiWindow and abs(z_diff) < ZWindow)
        {  //have match

          if (Verbosity())
          {
            cout << __PRETTY_FUNCTION__ << " adding cluster at layer " << layer << " to tracklet length " << tracklet.size() << endl;
          }

          auto iter_old_tracklet =
              new_tracklets.find(tracklet);
          if (iter_old_tracklet != new_tracklets.end())
            new_tracklets.erase(iter_old_tracklet);  //remove the shorter tracklet

          tracklet_t new_tracklet(tracklet);
          new_tracklet.push_back(cluster);
          new_tracklets.insert(new_tracklet);  // insert the longer tracklet
          matched = true;
        }

      }  //     for (auto& tracklet : tracklets)

      if (not matched)
      {  //new track seed from unused cluster
        new_tracklets.insert(tracklet_t(1, cluster));
        if (Verbosity())
        {
          cout << __PRETTY_FUNCTION__ << " init tracket with cluster at layer " << layer << endl;
        }
      }
    }  //      for (auto iter = cluster_range.first; iter != cluster_range.second; ++iter)

    tracklets = new_tracklets;
  }  //  for (int layer = 0; layer < nLayer; ++layer)

  if (Verbosity())
  {
    cout << __PRETTY_FUNCTION__ << "print trackets: ";
    for (const auto& tracklet : tracklets)
    {
      cout << "\t"
           << "size = " << tracklet.size() << endl;
    }
  }

  if (_do_eval)
  {
    fill_eval_tree(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * End
 */
int TpcPrototypeGenFitTrkFinder::End(PHCompositeNode* topNode)
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
TpcPrototypeGenFitTrkFinder::~TpcPrototypeGenFitTrkFinder()
{
  delete _fitter;
}

/*
 * fill_eval_tree():
 */
void TpcPrototypeGenFitTrkFinder::fill_eval_tree(PHCompositeNode* topNode)
{
  //! Make sure to reset all the TTree variables before trying to set them.
  reset_eval_variables();

  int i = 0;
  for (SvtxTrackMap::ConstIter itr = _trackmap->begin();
       itr != _trackmap->end(); ++itr)
  {
    new ((*_tca_trackmap)[i])(SvtxTrack_v1)(
        *dynamic_cast<SvtxTrack_v1*>(itr->second));
  }
  _tca_ntrack = i;

  _eval_tree->Fill();

  return;
}

/*
 * init_eval_tree
 */
void TpcPrototypeGenFitTrkFinder::init_eval_tree()
{
  if (!_tca_trackmap)
    _tca_trackmap = new TClonesArray("SvtxTrack_v1");

  //! create TTree
  _eval_tree = new TTree("T", "TpcPrototypeGenFitTrkFinder Evaluation");
  _eval_tree->Branch("nTrack", &_tca_ntrack, "nTrack/I");

  //  _eval_tree->Branch("PrimaryParticle", _tca_particlemap);
  //  _eval_tree->Branch("TruthVtx", _tca_vtxmap);

  _eval_tree->Branch("SvtxTrack", _tca_trackmap);

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
void TpcPrototypeGenFitTrkFinder::reset_eval_variables()
{
  _tca_ntrack = -1;
  _tca_trackmap->Clear();

  _cluster_eval_tree_x = WILD_FLOAT;
  _cluster_eval_tree_y = WILD_FLOAT;
  _cluster_eval_tree_z = WILD_FLOAT;
  _cluster_eval_tree_gx = WILD_FLOAT;
  _cluster_eval_tree_gy = WILD_FLOAT;
  _cluster_eval_tree_gz = WILD_FLOAT;
}

int TpcPrototypeGenFitTrkFinder::CreateNodes(PHCompositeNode* topNode)
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

  _trackmap = new SvtxTrackMap_v1;
  PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
      _trackmap, "SvtxTrackMap", "PHObject");
  tb_node->addNode(tracks_node);
  if (Verbosity() > 0)
    cout << "Svtx/SvtxTrackMap node added" << endl;

  _vertexmap = new SvtxVertexMap_v1;
  PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
      _vertexmap, "SvtxVertexMap", "PHObject");
  tb_node->addNode(vertexes_node);
  if (Verbosity() > 0)
    cout << "Svtx/SvtxVertexMap node added" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int TpcPrototypeGenFitTrkFinder::GetNodes(PHCompositeNode* topNode)
{
  // Input Svtx Clusters
  //_clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_clustermap && _event < 2)
  {
    return Fun4AllReturnCodes::ABORTRUN;
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

  return Fun4AllReturnCodes::EVENT_OK;
}

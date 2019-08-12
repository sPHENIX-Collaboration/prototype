/*!
 *  \file		TpcPrototypeGenFitTrkFinder.C
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "TpcPrototypeGenFitTrkFinder.h"

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackreco/AssocInfoContainer.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
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

#include <GenFit/EventDisplay.h>    // for EventDisplay
#include <GenFit/RKTrackRep.h>

#include <TClonesArray.h>
#include <TMatrixDSymfwd.h>  // for TMatrixDSym
#include <TMatrixTSym.h>     // for TMatrixTSym
#include <TMatrixTUtils.h>   // for TMatrixTRow
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>                              // for find
#include <set>                                    // for set

#include <cassert>
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
namespace PHGenFit { class Measurement; }

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
  , _fitter(nullptr)
  , _track_fitting_alg_name("KalmanFitter")
  , nLayer(layers)
  , minLayer(8)
  , maxTracklet(50)
  , _primary_pid_guess(2212)
  , rphiWindow(2)
  , ZWindow(2)
  , _clustermap(nullptr)
  , _trackmap(nullptr)
  , _assoc_container(nullptr)
  , _vertexmap(nullptr)
  , _do_eval(false)
  , _eval_outname("TpcPrototypeGenFitTrkFinder_eval.root")
  , _eval_tree(nullptr)
  , _tca_ntrack(-1)
  , _tca_trackmap(nullptr)
  , _do_evt_display(false)
{
  Verbosity(0);

  _event = 0;

  _cluster_eval_tree = nullptr;
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

  assert(_clustermap);
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

          if (Verbosity() >= 2)
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

          if (new_tracklets.size() >= maxTracklet)
          {
            if (Verbosity())
              cout << __PRETTY_FUNCTION__ << " skipping rest tracklet at layer " << layer
                   << " due to tracklet count " << new_tracklets.size() << " > " << maxTracklet << endl;
            break;
          }
        }  //have match

      }  //     for (auto& tracklet : tracklets)

      if (new_tracklets.size() >= maxTracklet)
      {
        if (Verbosity())
          cout << __PRETTY_FUNCTION__ << " skipping rest clusters at layer " << layer
               << " due to tracklet count " << new_tracklets.size() << " > " << maxTracklet << endl;
        break;
      }

      if (not matched)
      {  //new track seed from unused cluster
        new_tracklets.insert(tracklet_t(1, cluster));
        if (Verbosity() >= 2)
        {
          cout << __PRETTY_FUNCTION__ << " init tracket with cluster at layer " << layer << endl;
        }
      }  //      if (not matched)

    }  //      for (auto iter = cluster_range.first; iter != cluster_range.second; ++iter)

    tracklets = new_tracklets;
  }  //  for (int layer = 0; layer < nLayer; ++layer)

  if (Verbosity())
  {
    cout << __PRETTY_FUNCTION__ << "print initial trackets: ";
    for (const auto& tracklet : tracklets)
    {
      cout << ","
           << "size = " << tracklet.size();
    }
    cout << endl;
  }

  // track quality map
  multimap<double, tracklet_t> quality_tracklets_map;
  for (const auto& tracklet : tracklets)
  {
    double chi2Ndf = getChi2Ndf(tracklet);
    if (chi2Ndf > 0)
      quality_tracklets_map.insert(std::pair<double, tracklet_t>(chi2Ndf, tracklet));
  }
  if (Verbosity())
  {
    cout << __PRETTY_FUNCTION__ << "print fitted trackets: ";
    for (const auto& tracklet : quality_tracklets_map)
    {
      cout << ","
           << "size = " << tracklet.second.size() << " chi2/ndf = " << tracklet.first;
    }
    cout << endl;
  }  //  if (Verbosity())

  // deghost
  for (auto iter_current = quality_tracklets_map.begin(); iter_current != quality_tracklets_map.end(); ++iter_current)
  {
    const tracklet_t& track_current = iter_current->second;

    auto iter_check = iter_current;
    ++iter_check;
    for (; iter_check != quality_tracklets_map.end();)
    {
      const tracklet_t& track_check = iter_check->second;

      unsigned int identical_cluster = 0;

      for (const auto cluster : track_current)
      {
        assert(cluster);

        if (find(track_check.begin(), track_check.end(), cluster) != track_check.end())
          ++identical_cluster;
      }

      if (identical_cluster >= track_current.size() / 3 or identical_cluster >= track_check.size() / 3)
      {
        if (Verbosity())
        {
          cout << __PRETTY_FUNCTION__ << " found " << identical_cluster << "-shared ghost track size" << track_check.size() << " chi2ndf" << iter_check->first  //
               << " from base track size" << track_current.size() << " chi2ndf" << iter_current->first << endl;
        }

        auto iter_tmp = iter_check;
        ++iter_check;

        quality_tracklets_map.erase(iter_tmp);
      }
      else
      {
        if (Verbosity())
        {
          cout << __PRETTY_FUNCTION__ << "low " << identical_cluster << "-shared track size" << track_check.size() << " chi2ndf" << iter_check->first  //
               << " from base track size" << track_current.size() << " chi2ndf" << iter_current->first << endl;
        }
        ++iter_check;
      }
    }
  }
  if (Verbosity())
  {
    cout << __PRETTY_FUNCTION__ << "print deghosted trackets: ";
    for (const auto& tracklet : quality_tracklets_map)
    {
      cout << ","
           << "size = " << tracklet.second.size() << " chi2/ndf = " << tracklet.first;
    }
    cout << endl;
  }  //  if (Verbosity())

  //output
  assert(_trackmap);
  assert(_assoc_container);
  for (const auto& quality_tracklet : quality_tracklets_map)
  {
    const tracklet_t& tracklet = quality_tracklet.second;
    assert(tracklet.front());
    assert(tracklet.back());
    TVector3 pos_front(tracklet.front()->getPosition(0), tracklet.front()->getPosition(1), tracklet.front()->getPosition(2));
    TVector3 pos_back(tracklet.back()->getPosition(0), tracklet.back()->getPosition(1), tracklet.back()->getPosition(2));

    TVector3 seed_mom(pos_back - pos_front);
    seed_mom.SetMag(120);

    TVector3 seed_pos(pos_front);

    std::unique_ptr<SvtxTrack_v1> svtx_track(new SvtxTrack_v1());

    svtx_track->set_id(_trackmap->size());

    // dummy values, set px to make it through the minimum pT cut
    svtx_track->set_px(seed_mom.x());
    svtx_track->set_py(seed_mom.y());
    svtx_track->set_pz(seed_mom.z());
    svtx_track->set_x(seed_pos.x());
    svtx_track->set_y(seed_pos.y());
    svtx_track->set_z(seed_pos.z());
    for (const TrkrCluster* cluster : tracklet)
    {
      svtx_track->insert_cluster_key(cluster->getClusKey());
      _assoc_container->SetClusterTrackAssoc(cluster->getClusKey(), svtx_track->get_id());
    }

    _trackmap->insert(svtx_track.get());
  }

  if (_do_eval)
  {
    fill_eval_tree(topNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

double TpcPrototypeGenFitTrkFinder::getChi2Ndf(const tracklet_t& tracklet)
{
  assert(_clustermap);

  if (tracklet.size() < minLayer)
  {
    if (Verbosity())
    {
      cout << __PRETTY_FUNCTION__ << " drop short tracklet size = " << tracklet.size() << " < " << minLayer << endl;
    }
    return -1;
  }

  // prepare seed
  assert(tracklet.front());
  assert(tracklet.back());
  TVector3 pos_front(tracklet.front()->getPosition(0), tracklet.front()->getPosition(1), tracklet.front()->getPosition(2));
  TVector3 pos_back(tracklet.back()->getPosition(0), tracklet.back()->getPosition(1), tracklet.back()->getPosition(2));

  TVector3 seed_mom(pos_back - pos_front);
  seed_mom.SetMag(120);

  TVector3 seed_pos(pos_front);
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
  for (const auto cluster : tracklet)
  {
    assert(cluster);
    if (Verbosity() >= 2)
    {
      TrkrDefs::cluskey cluster_key = cluster->getClusKey();
      int layer = TrkrDefs::getLayer(cluster_key);

      cout << __PRETTY_FUNCTION__ << "add cluster on layer " << layer << ": ";
      cluster->identify();
    }
    TVector3 pos(cluster->getPosition(0), cluster->getPosition(1), cluster->getPosition(2));

    //TODO use u, v explicitly?
    TVector3 n(cluster->getPosition(0), cluster->getPosition(1), 0);

    PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
                                                                  cluster->getRPhiError(),
                                                                  cluster->getZError());

    measurements.push_back(meas);

  }  //    for (const auto cluster : clusterLayer)

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
    return -1;
  }

  if (Verbosity() >= 1)
  {
    cout << __PRETTY_FUNCTION__ << " track->getChisq() " << phgf_track->get_chi2() << " get_ndf " << phgf_track->get_ndf()
         << " mom.X " << phgf_track->get_mom().X()
         << " mom.Y " << phgf_track->get_mom().Y()
         << " mom.Z " << phgf_track->get_mom().Z()
         << endl;
  }
  return phgf_track->get_chi2() / phgf_track->get_ndf();
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

  _assoc_container = new AssocInfoContainer;
  PHIODataNode<PHObject>* assoc_node = new PHIODataNode<PHObject>(
      _assoc_container, "AssocInfoContainer", "PHObject");
  tb_node->addNode(assoc_node);
  if (Verbosity() > 0)
    cout << "Svtx/AssocInfoContainer node added" << endl;
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

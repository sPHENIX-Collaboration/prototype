#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <fun4all/Fun4AllServer.h>

#include <g4eval/SvtxEvaluator.h>
#include <g4eval/TrkrEvaluator.h>

#include <g4detectors/PHG4BlockSubsystem.h>
#include <g4main/PHG4Reco.h>

#include <g4tpc/PHG4TpcDigitizer.h>
#include <g4tpc/PHG4TpcElectronDrift.h>
#include <g4tpc/PHG4TpcPadPlane.h>
#include <g4tpc/PHG4TpcPadPlaneReadout.h>
#include <g4tpc/PHG4TpcSubsystem.h>
#include <trackreco/PHGenFitTrkFitter.h>
#include <trackreco/PHGenFitTrkProp.h>
#include <trackreco/PHHoughSeeding.h>
#include <trackreco/PHInitVertexing.h>
#include <trackreco/PHTrackSeeding.h>
#include <trackreco/PHTruthTrackSeeding.h>
#include <trackreco/PHTruthVertexing.h>

#include <g4eval/PHG4DSTReader.h>

#include <g4histos/G4HitNtuple.h>

#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>

#include <phool/recoConsts.h>

#include <tpc2019/TpcPrototypeClusterizer.h>
#include <tpc2019/TpcPrototypeGenFitTrkFitter.h>

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4tpc.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libg4tpc.so)
R__LOAD_LIBRARY(libg4intt.so)
R__LOAD_LIBRARY(libg4mvtx.so)
R__LOAD_LIBRARY(libg4hough.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libtpc2019.so)
R__LOAD_LIBRARY(libtrack_reco.so)
#endif

int n_tpc_layer_inner = 16;
int tpc_layer_rphi_count_inner = 1152;
int n_tpc_layer_mid = 16;
int n_tpc_layer_outer = 16;
int n_gas_layer = n_tpc_layer_inner + n_tpc_layer_mid + n_tpc_layer_outer;

int Fun4All_G4_TPC(int nEvents = 1, bool eventDisp = false, int verbosity = 1)
{
  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors");
  gSystem->Load("libg4testbench");
  gSystem->Load("libg4histos");
  gSystem->Load("libg4eval.so");
  gSystem->Load("libqa_modules");
  gSystem->Load("libg4tpc");
  gSystem->Load("libtrack_io.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libtpc2019.so");
  gSystem->Load("libg4eval.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4hough.so");
  gSystem->Load("libtrack_reco.so");

  bool bh_on = false;  // the surrounding boxes need some further thinking
  bool dstreader = true;
  bool dstoutput = false;

  const double TPCDriftLength = 40;

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();
  // only set this if you want a fixed random seed to make
  // results reproducible for testing
  // rc->set_IntFlag("RANDOMSEED",12345678);

  // simulated setup sits at eta=1, theta=40.395 degrees
  double theta = 90 - 5;
  double phi = 180 + 360 / 12 / 2;
  // shift in x with respect to midrapidity setup
  double add_place_z = -TPCDriftLength * .5;
  // Test beam generator
  PHG4SimpleEventGenerator *gen = new PHG4SimpleEventGenerator();
  gen->add_particles("proton", 1);  // mu-,e-,anti_proton,pi-
  gen->set_vertex_distribution_mean(0.0, 0.0, add_place_z);
  gen->set_vertex_distribution_width(0.0, .7, .7);  // Rough beam profile size @ 16 GeV measured by Abhisek
  gen->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus,
                                        PHG4SimpleEventGenerator::Gaus,
                                        PHG4SimpleEventGenerator::Gaus);  // Gauss beam profile
  double angle = theta * TMath::Pi() / 180.;
  double eta = -1. * TMath::Log(TMath::Tan(angle / 2.));
  gen->set_eta_range(eta - 0.001, eta + 0.001);                                          // 1mrad angular divergence
  gen->set_phi_range(TMath::Pi() * phi / 180 - 0.001, TMath::Pi() * phi / 180 + 0.001);  // 1mrad angular divergence
  const double momentum = 120;
  gen->set_p_range(momentum, momentum, momentum * 2e-2);  // 2% momentum smearing
  se->registerSubsystem(gen);

  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field(0);
  //  g4Reco->SetPhysicsList("QGSP_BERT_HP"); // uncomment this line to enable the high-precision neutron simulation physics list, QGSP_BERT_HP

  //----------------------------------------
  // TPC G4
  //----------------------------------------

  PHG4TpcSubsystem *tpc = new PHG4TpcSubsystem("TPC");
  tpc->SetActive();
  tpc->SuperDetector("TPC");
  tpc->set_double_param("steplimits", 1);
  tpc->set_double_param("tpc_length", TPCDriftLength *2);  // 2x 40 cm drift
  // By default uses "sPHENIX_TPC_Gas", defined in PHG4Reco. That is 90:10 Ne:C4
  tpc->SetAbsorberActive();
  g4Reco->registerSubsystem(tpc);

  if (bh_on)
  {
    // BLACKHOLE, box surrounding the prototype to check for leakage
    PHG4BlockSubsystem *bh[5];
    // surrounding outer hcal
    // top
    bh[0] = new PHG4BlockSubsystem("bh1", 1);
    bh[0]->set_double_param("size_x", 270.);
    bh[0]->set_double_param("size_y", 0.01);
    bh[0]->set_double_param("size_z", 165.);
    bh[0]->set_double_param("place_x", 270. / 2.);
    bh[0]->set_double_param("place_y", 125. / 2.);
    // bottom
    bh[1] = new PHG4BlockSubsystem("bh2", 2);
    bh[1]->set_double_param("size_x", 270.);
    bh[1]->set_double_param("size_y", 0.01);
    bh[1]->set_double_param("size_z", 165.);
    bh[1]->set_double_param("place_x", 270. / 2.);
    bh[1]->set_double_param("place_y", -125. / 2.);
    // right side
    bh[2] = new PHG4BlockSubsystem("bh3", 3);
    bh[2]->set_double_param("size_x", 200.);
    bh[2]->set_double_param("size_y", 125.);
    bh[2]->set_double_param("size_z", 0.01);
    bh[2]->set_double_param("place_x", 200. / 2.);
    bh[2]->set_double_param("place_z", 165. / 2.);
    // left side
    bh[3] = new PHG4BlockSubsystem("bh4", 4);
    bh[3]->set_double_param("size_x", 270.);
    bh[3]->set_double_param("size_y", 125.);
    bh[3]->set_double_param("size_z", 0.01);
    bh[3]->set_double_param("place_x", 270. / 2.);
    bh[3]->set_double_param("place_z", -165. / 2.);
    // back
    bh[4] = new PHG4BlockSubsystem("bh5", 5);
    bh[4]->set_double_param("size_x", 0.01);
    bh[4]->set_double_param("size_y", 125.);
    bh[4]->set_double_param("size_z", 165.);
    bh[4]->set_double_param("place_x", 270.);
    for (int i = 0; i < 5; i++)
    {
      bh[i]->BlackHole();
      bh[i]->SetActive();
      bh[i]->SuperDetector("BlackHole");
      bh[i]->OverlapCheck(true);
      g4Reco->registerSubsystem(bh[i]);
    }
  }
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  se->registerSubsystem(g4Reco);

  //=========================
  // setup Tpc readout for filling cells
  // g4tpc/PHG4TpcElectronDrift uses
  // g4tpc/PHG4TpcPadPlaneReadout
  //=========================

  PHG4TpcPadPlane *padplane = new PHG4TpcPadPlaneReadout();

  // The pad plane readout default is set in PHG4TpcPadPlaneReadout
  // We may want to change the number of inner layers, and can do that here
  padplane->set_int_param("tpc_minlayer_inner", 0);  // sPHENIX layer number of first Tpc readout layer
  padplane->set_int_param("ntpc_layers_inner", n_tpc_layer_inner);
  padplane->set_int_param("ntpc_phibins_inner", tpc_layer_rphi_count_inner);

  //only build mid-layer and to size
  padplane->set_int_param("ntpc_layers_inner", 0);
  padplane->set_int_param("ntpc_layers_outer", 0);
  padplane->set_double_param("maxdriftlength", TPCDriftLength);
  padplane->set_int_param("ntpc_phibins_mid", 16*8*12);

  PHG4TpcElectronDrift *edrift = new PHG4TpcElectronDrift();
  edrift->Detector("TPC");
  // fudge factors to get drphi 150 microns (in mid and outer Tpc) and dz 500 microns cluster resolution
  // They represent effects not due to ideal gas properties and ideal readout plane behavior
  // defaults are 0.12 and 0.15, they can be changed here to get a different resolution
  edrift->set_double_param("added_smear_trans", 0.12);
  edrift->set_double_param("added_smear_long", 0.15);
  edrift->registerPadPlane(padplane);
  edrift->Verbosity(verbosity);;
  se->registerSubsystem(edrift);

  // Tpc
  //====
  PHG4TpcDigitizer *digitpc = new PHG4TpcDigitizer();
  digitpc->SetTpcMinLayer(0);
  double ENC = 670.0;  // standard
  digitpc->SetENC(ENC);
  double ADC_threshold = 4.0 * ENC;
  digitpc->SetADCThreshold(ADC_threshold);  // 4 * ENC seems OK

  cout << " Tpc digitizer: Setting ENC to " << ENC << " ADC threshold to " << ADC_threshold << endl;

  se->registerSubsystem(digitpc);

  // For the Tpc
  //==========
  TpcPrototypeClusterizer *tpcclusterizer = new TpcPrototypeClusterizer();
  tpcclusterizer->Verbosity(verbosity);;
  se->registerSubsystem(tpcclusterizer);

  //-------------
  // Tracking
  //------------

  // This should be true for everything except testing wirh truth track seeding!
  const bool use_track_prop = false;
  if (use_track_prop)
  {
    //--------------------------------------------------
    // Normal track seeding and propagation
    //--------------------------------------------------

    // for now, we cheat to get the initial vertex for the full track reconstruction case
    PHInitVertexing *init_vtx = new PHTruthVertexing("PHTruthVertexing");
    init_vtx->Verbosity(verbosity);;
    se->registerSubsystem(init_vtx);

    // find seed tracks using a subset of TPC layers
    PHTrackSeeding *track_seed = new PHHoughSeeding("PHHoughSeeding", 0, 0, n_gas_layer);
    track_seed->Verbosity(verbosity);;
    se->registerSubsystem(track_seed);

    // Find all clusters associated with each seed track
    PHGenFitTrkProp *track_prop = new PHGenFitTrkProp("PHGenFitTrkProp", 0, 0, n_gas_layer);
    track_prop->Verbosity(verbosity);;
    se->registerSubsystem(track_prop);
  }
  else
  {
    //--------------------------------------------------
    // Track finding using truth information
    //--------------------------------------------------

    PHInitVertexing *init_vtx = new PHTruthVertexing("PHTruthVertexing");
    init_vtx->Verbosity(verbosity);;
    se->registerSubsystem(init_vtx);

    // For each truth particle, create a track and associate clusters with it using truth information
    PHTruthTrackSeeding *pat_rec = new PHTruthTrackSeeding("PHTruthTrackSeeding");
    pat_rec->Verbosity(verbosity);;
    pat_rec->set_min_clusters_per_track(10);
    se->registerSubsystem(pat_rec);
  }
  //
  //------------------------------------------------
  // Fitting of tracks using Kalman Filter
  //------------------------------------------------

  TpcPrototypeGenFitTrkFitter *kalman = new TpcPrototypeGenFitTrkFitter();
  kalman->Verbosity(verbosity);;
  kalman->set_do_evt_display(eventDisp);
  kalman->set_do_eval(true);
  se->registerSubsystem(kalman);

  //----------------
  // Tracking evaluation
  //----------------

  SvtxEvaluator *eval;
  eval = new SvtxEvaluator("SVTXEVALUATOR", "G4TPC_eval.root", "SvtxTrackMap", 0, 0, n_gas_layer);
  eval->do_cluster_eval(true);
  eval->do_g4hit_eval(true);
  eval->do_hit_eval(true);  // enable to see the hits that includes the chamber physics...
  eval->do_gpoint_eval(true);
  eval->do_eval_light(false);
  eval->scan_for_embedded(false);  // take all tracks if false - take only embedded tracks if true
  eval->Verbosity(verbosity);;
  se->registerSubsystem(eval);

  //----------------------
  // save a comprehensive  evaluation file
  //----------------------
  if (dstreader)
  {
    PHG4DSTReader *ana = new PHG4DSTReader(string("G4TPC_DSTReader.root"));
    ana->set_save_particle(true);
    ana->set_load_all_particle(true);
    ana->set_load_active_particle(true);
    ana->set_save_vertex(true);

    ana->AddNode("ABSORBER_TPC");
    ana->AddNode("TPC");
    if (bh_on)
      ana->AddNode("BlackHole");  // add a G4Hit node

    se->registerSubsystem(ana);
  }

  if (dstoutput)
  {
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", "G4TPC.root");
    se->registerOutputManager(out);
  }

  Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
  se->registerInputManager(in);
  if (nEvents <= 0)
  {
    return 0;
  }

  gSystem->ListLibraries();
  se->run(nEvents);

  se->End();

  //   std::cout << "All done" << std::endl;
  delete se;
  //   return 0;
  gSystem->Exit(0);
  return 0;
}

// for using QuickTest to check if macro loads
void RunLoadTest() {}

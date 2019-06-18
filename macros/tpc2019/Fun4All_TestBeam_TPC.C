#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllPrdfInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <fun4all/Fun4AllServer.h>

#include <g4eval/SvtxEvaluator.h>
#include <g4eval/TrkrEvaluator.h>

#include <g4detectors/PHG4BlockSubsystem.h>
#include <g4main/PHG4Reco.h>
#include <phgeom/PHGeomFileImport.h>

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
#include <tpc2019/TpcPrototypeGenFitTrkFinder.h>
#include <tpc2019/TpcPrototypeGenFitTrkFitter.h>
#include <tpc2019/TpcPrototypeUnpacker.h>

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

int n_tpc_layer_inner = 0;
int tpc_layer_rphi_count_inner = 1152;
int n_tpc_layer_mid = 16;
int n_tpc_layer_outer = 0;
int n_gas_layer = n_tpc_layer_inner + n_tpc_layer_mid + n_tpc_layer_outer;

int Fun4All_TestBeam_TPC(int nEvents = 2, int nSkip = 260,
                         const string &input_file = "data/tpc_cosmics_00000146-0000.evt",
                         bool eventDisp = false, int verbosity = 1)
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
  bool dstoutput = false;

  const double TPCDriftLength = 40;

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(2);
  recoConsts *rc = recoConsts::instance();
  // only set this if you want a fixed random seed to make
  // results reproducible for testing
  rc->set_IntFlag("RANDOMSEED", 12345678);

  //swap out with the test beam geometry for the analysis stage
  PHGeomFileImport *import = new PHGeomFileImport("TpcPrototypeGeometry.gdml");
  se->registerSubsystem(import);

  PHG4TpcPadPlane *padplane = new PHG4TpcPadPlaneReadout();

  // The pad plane readout default is set in PHG4TpcPadPlaneReadout
  //only build mid-layer and to size
  padplane->set_int_param("tpc_minlayer_inner", 0);
  padplane->set_int_param("ntpc_layers_inner", 0);
  padplane->set_int_param("ntpc_layers_outer", 0);
  padplane->set_double_param("maxdriftlength", TPCDriftLength);
  padplane->set_int_param("ntpc_phibins_mid", 16 * 8 * 12);

  TpcPrototypeUnpacker *tpcfee = new TpcPrototypeUnpacker((input_file) + string("_TpcPrototypeUnpacker.root"));
  //  tpcfee->Verbosity(TPCFEETestRecov1::VERBOSITY_SOME);
  tpcfee->Verbosity(TpcPrototypeUnpacker::VERBOSITY_MORE);

  tpcfee->setNPreSample(5);
  tpcfee->setNPostSample(7);

  se->registerSubsystem(tpcfee);

  //
  //  // For the Tpc
  //  //==========
  //  TpcPrototypeClusterizer *tpcclusterizer = new TpcPrototypeClusterizer();
  //  tpcclusterizer->Verbosity(verbosity);
  //  ;
  //  se->registerSubsystem(tpcclusterizer);
  //
  //  //-------------
  //  // Tracking
  //  //------------
  //
  ////  // Find all clusters associated with each seed track
  //  TpcPrototypeGenFitTrkFinder *finder = new TpcPrototypeGenFitTrkFinder();
  //  finder->Verbosity(verbosity);
  //  finder->set_do_evt_display(eventDisp);
  //  finder->set_do_eval(true);
  //  finder->set_eval_filename(input_file + "_TpcPrototypeGenFitTrkFinder.root");
  //  se->registerSubsystem(finder);
  ////
  ////  //
  ////  //------------------------------------------------
  ////  // Fitting of tracks using Kalman Filter
  ////  //------------------------------------------------
  ////
  //  TpcPrototypeGenFitTrkFitter *kalman = new TpcPrototypeGenFitTrkFitter();
  //  kalman->Verbosity(verbosity);
  //  kalman->set_do_evt_display(eventDisp);
  //  kalman->set_eval_filename(input_file + "_TpcPrototypeGenFitTrkFinder.root");
  //  kalman->set_do_eval(true);
  //  se->registerSubsystem(kalman);

  if (dstoutput)
  {
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", input_file + "_DST.root");
    se->registerOutputManager(out);
  }

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
  in->fileopen(input_file);
  se->registerInputManager(in);

  gSystem->ListLibraries();
  se->skip(nSkip);
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

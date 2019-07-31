#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllPrdfInputManager.h>

#include <phool/recoConsts.h>

#include <mvtxprototype2/MvtxPrototype2Geom.h>
#include <mvtxprototype2/MvtxPrototype2UnpackPRDF.h>
#include <mvtxprototype2/MvtxPrototype2Clusterizer.h>
#include <mvtxprototype2/MvtxPrototype2Align.h>
#include <anamvtxprototype2/AnaMvtxBeamTest2019.h>

#include <TSystem.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libMvtxPrototype2.so)
R__LOAD_LIBRARY(libAnaMvtxBeamTest2019.so)
#endif

#include <string>
#include <iostream>
using namespace std;

void Fun4All_TestBeam(
		int nEvents = 10,
		const string &input_file = "calib-00000648-0000.prdf")
{
  gSystem->Load("libfun4all");
  gSystem->Load("libMvtxPrototype2.so");
  gSystem->Load("libAnaMvtxBeamTest2019.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(Fun4AllServer::VERBOSITY_SOME);

	char ifile[input_file.length()+1]; // + 1 for the \0 which marks the end of string
	strcpy(ifile, input_file.c_str());
	strtok(ifile,"-");
	int runnumber = atoi(strtok(0,"-"));
	int segnumber = atoi(strtok(strtok(0,"-"),"."));
	cout << "runnumber: " << runnumber << " segment " << segnumber << endl;

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",runnumber);

	/*
	MvtxRunInfoUnpackPRDF *unpack_run = new MvtxRunInfoUnpackPRDF();
	unpack_run->Verbosity(0);
	se->registerSubsystem(unpack_run);
	*/

	//Unpack
	MvtxPrototype2UnpackPRDF *unpack = new MvtxPrototype2UnpackPRDF();
	unpack->Verbosity(1);
	se->registerSubsystem(unpack);

	/*
	MvtxApplyHotDead *hotdead = new MvtxApplyHotDead("MvtxApplyHotDead",Form("hotmap/hotmap_testbeam_%08d.txt",runnumber));
	hotdead->Verbosity(0);
	hotdead->PrintHotDeadMap();
	se->registerSubsystem(hotdead);
	*/

	MvtxPrototype2Clusterizer *clus = new MvtxPrototype2Clusterizer();
	clus->Verbosity(0);
	se->registerSubsystem(clus);

  bool do_align = true;
  if ( do_align )
  {
    MvtxPrototype2Align *align = new MvtxPrototype2Align("MvtxAlign");
    //align->SetAlignmentParFileDir("beamcenter/");
    align->SetAlignmentParFileDir("./");
    align->SetAlignmentParFileName("beamcenter_mean_0848_0849_0851.txt");
    //align->PrintAlignmentPars();
    align->Verbosity(0);
    se->registerSubsystem(align);
  }

  bool do_eval  = true;
  bool do_standalone_tracking = do_eval && true;
  bool do_ghost_rej = do_standalone_tracking && true;
  if ( do_eval )
  {

    AnaMvtxBeamTest2019 *eval = new AnaMvtxBeamTest2019();
	  eval->set_out_filename(Form("MvtxBT2019Eval-%08d-%04d.root",
                           runnumber,segnumber+200));
    eval->set_do_tracking(do_standalone_tracking);
    eval->getStandaloneTracking()->SetGhostRejection(do_ghost_rej);
    eval->getStandaloneTracking()->Verbosity(0);
    eval->Verbosity(0);
	  se->registerSubsystem(eval);

  }

/*
  if ( do_tracking )
  {
    //  // Find all clusters associated with each seed track
    TpcPrototypeGenFitTrkFinder *finder = new TpcPrototypeGenFitTrkFinder();
    finder->Verbosity(verbosity);
    finder->set_do_evt_display(eventDisp);
    finder->set_do_eval(false);
    finder->set_eval_filename(input_file + "_TpcPrototypeGenFitTrkFinder.root");
    se->registerSubsystem(finder);
    ////
    ////  //
    //  //------------------------------------------------
    //  // Fitting of tracks using Kalman Filter
    //  //------------------------------------------------
    //
    TpcPrototypeGenFitTrkFitter *kalman = new TpcPrototypeGenFitTrkFitter();
    kalman->Verbosity(verbosity);
    kalman->set_do_evt_display(eventDisp);
    kalman->set_eval_filename(input_file + "_TpcPrototypeGenFitTrkFitter.root");
    kalman->set_do_eval(true);
    se->registerSubsystem(kalman);
  }
*/

  bool do_dst = false;
  if ( do_dst )
  {
	  Fun4AllDstOutputManager *out_Manager = new Fun4AllDstOutputManager("DSTOUT",  input_file + "_DST.root" );
	  se->registerOutputManager(out_Manager);
  }

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
  in->fileopen(input_file);
  se->registerInputManager(in);

  se->run(nEvents);

  se->End();
}

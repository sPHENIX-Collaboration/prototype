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
#include <mvtxprototype2/MvtxPrototype2Eval.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libMvtxPrototype2.so)
#endif

#include <string>
#include <iostream>
using namespace std;

void Fun4All_TestBeam(
		int nEvents = 10,
		const char *input_file = "calib-00000648-0000.prdf",
		const char *output_file = "DST-calib-00000648-0000.root",
    bool do_align = true
		)
{
  gSystem->Load("libfun4all");
  gSystem->Load("libMvtxPrototype2.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(Fun4AllServer::VERBOSITY_SOME);

	char ifile[strlen(input_file)+1]; // + 1 for the \0 which marks the end of string
	strcpy(ifile,input_file);
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

  if ( do_align )
  {
    MvtxPrototype2Align *align = new MvtxPrototype2Align("MvtxAlign");
    //align->SetAlignmentParFileDir("beamcenter/");
    align->SetAlignmentParFileDir("./");
    align->SetAlignmentParFileName("beamcenter_mean_0848_0849_0851.txt");
    align->PrintAlignmentPars();
    align->Verbosity(0);
    se->registerSubsystem(align);
  }

	MvtxPrototype2Eval *qa = new MvtxPrototype2Eval();
	qa->set_filename(Form("MvtxPrototype2Eval-%08d-%04d.root",runnumber,segnumber+200));
	se->registerSubsystem(qa);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
  in->fileopen(input_file);
  se->registerInputManager(in);

	//Fun4AllDstOutputManager *out_Manager = new Fun4AllDstOutputManager("DSTOUT",output_file);
	//se->registerOutputManager(out_Manager);

  se->run(nEvents);

  se->End();
}

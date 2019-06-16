// $Id: $

/*!
 * \file ConstructGeometry.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <fun4all/Fun4AllServer.h>

#include <g4eval/SvtxEvaluator.h>
#include <g4eval/TrkrEvaluator.h>

#include <g4detectors/PHG4CylinderSubsystem.h>
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

#include <phgeom/PHGeomUtility.h>
#include <phool/recoConsts.h>

#include <TMath.h>
#include <TSystem.h>

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

void ConstructGeometry()
{
  const double placementR = 0.5 * (40 + 60);
  const double rotaitonZ = TMath::Pi() + TMath::Pi() * 2 / 12. / 2.;
  const double driftLength = 40;
  const double cageRadius = 20;

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
  double theta = 90;
  double phi = 180 + 360 / 12 / 2;
  // shift in x with respect to midrapidity setup
  double add_place_z = -driftLength * .5;
  // Test beam generator
  PHG4SimpleEventGenerator *gen = new PHG4SimpleEventGenerator();
  gen->add_particles("proton", 1);  // mu-,e-,anti_proton,pi-
  gen->set_vertex_distribution_mean(0.0, 0.0, add_place_z);
  gen->set_vertex_distribution_width(0.0, .0, .0);  // Rough beam profile size @ 16 GeV measured by Abhisek
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

  {
    PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem("TPC_GasVol", 0);
    cyl->set_double_param("length", driftLength);
    cyl->set_double_param("place_x", placementR * TMath::Cos(rotaitonZ));
    cyl->set_double_param("place_y", placementR * TMath::Sin(rotaitonZ));
    cyl->set_double_param("place_z", -driftLength / 2);
    cyl->set_double_param("radius", 0.0);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_string_param("material", "sPHENIX_TPC_Gas");
    cyl->set_double_param("thickness", cageRadius);
    cyl->SuperDetector("TPC");
    cyl->SetActive();
    cyl->OverlapCheck(1);
    g4Reco->registerSubsystem(cyl);
  }
  {
    PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem("TPC_FieldCage", 1);
    cyl->set_double_param("length", driftLength);
    cyl->set_double_param("place_x", placementR * TMath::Cos(rotaitonZ));
    cyl->set_double_param("place_y", placementR * TMath::Sin(rotaitonZ));
    cyl->set_double_param("place_z", -driftLength / 2);
    cyl->set_double_param("radius", cageRadius + 1e-4);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_string_param("material", "G4_Cu");
    cyl->set_double_param("thickness", 0.00347);
    cyl->SuperDetector("TPC_Support");
    cyl->SetActive();
    cyl->OverlapCheck(1);
    g4Reco->registerSubsystem(cyl);
  }

  for (int sign = -1; sign <= 1; sign += 2)
  {
    const double endcap_thickness = 0.5;

    PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem("TPC_EndCap", 4 + sign);
    cyl->set_double_param("length", endcap_thickness);
    cyl->set_double_param("place_x", placementR * TMath::Cos(rotaitonZ));
    cyl->set_double_param("place_y", placementR * TMath::Sin(rotaitonZ));
    cyl->set_double_param("place_z", -driftLength / 2 + sign * ((driftLength / 2) + endcap_thickness / 2 + 1e-4));
    cyl->set_double_param("radius", 0.0);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_string_param("material", "G4_Al");
    cyl->set_double_param("thickness", cageRadius);
    cyl->SuperDetector("TPC_Support");
    cyl->SetActive();
    cyl->OverlapCheck(1);
    g4Reco->registerSubsystem(cyl);
  }

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  se->registerSubsystem(g4Reco);

  Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
  se->registerInputManager(in);

  char cmd[100];
  g4Reco->InitRun(se->topNode());
  g4Reco->ApplyDisplayAction();
  sprintf(cmd, "/control/execute vis.mac");
  g4Reco->ApplyCommand(cmd);

  se->run(1);
  g4Reco->ApplyCommand("/vis/scene/add/axes 0 0 0 50 cm");
  g4Reco->ApplyCommand("/vis/viewer/zoom 1");

  PHGeomUtility::ExportGeomtry(se->topNode(), "TpcPrototypeGeometry.root");
  PHGeomUtility::ExportGeomtry(se->topNode(), "TpcPrototypeGeometry.gdml");
}

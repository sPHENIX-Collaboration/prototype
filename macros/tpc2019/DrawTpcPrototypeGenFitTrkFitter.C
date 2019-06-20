#include "SaveCanvas.C"
#include "sPhenixStyle.C"

#include <TChain.h>
#include <TCut.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TH3.h>
#include <TPolyLine.h>
#include <TTree.h>

#include <TFile.h>

#include <TColor.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include <TDirectory.h>
#include <TMath.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>
#include <TVirtualFitter.h>

#include <cmath>
#include <iostream>

using namespace std;

// ROOT6 disabled assert. Well....
#ifdef assert
#undef assert
#endif
#define assert(exp)                                                                             \
  {                                                                                             \
    if (!(exp))                                                                                 \
    {                                                                                           \
      cout << "Assert (" << #exp << ") failed at " << __FILE__ << " line " << __LINE__ << endl; \
      exit(1);                                                                                  \
    }                                                                                           \
  }

//! utility function to
void useLogBins(TAxis *axis)
{
  assert(axis);
  assert(axis->GetXmin() > 0);
  assert(axis->GetXmax() > 0);

  const int bins = axis->GetNbins();

  Axis_t from = log10(axis->GetXmin());
  Axis_t to = log10(axis->GetXmax());
  Axis_t width = (to - from) / bins;
  vector<Axis_t> new_bins(bins + 1);

  for (int i = 0; i <= bins; i++)
  {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins.data());
}

TGraphErrors *
FitProfile(const TH2 *h2)
{
  TProfile *p2 = h2->ProfileX();

  int n = 0;
  double x[1000];
  double ex[1000];
  double y[1000];
  double ey[1000];

  for (int i = 1; i <= h2->GetNbinsX(); i++)
  {
    TH1D *h1 = h2->ProjectionY(Form("htmp_%d", rand()), i, i);

    if (h1->GetSum() < 30)
    {
      cout << "FitProfile - ignore bin " << i << endl;
      continue;
    }
    else
    {
      //      cout << "FitProfile - fit bin " << i << endl;
    }

    TF1 fgaus("fgaus", "gaus", -p2->GetBinError(i) * 4,
              p2->GetBinError(i) * 4);

    //    TF1 f2(Form("dgaus"), "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]",
    //           -p2->GetBinError(i) * 4, p2->GetBinError(i) * 4);

    fgaus.SetParameter(1, p2->GetBinContent(i));
    fgaus.SetParameter(2, p2->GetBinError(i));

    h1->Fit(&fgaus, "MQ");

    TF1 f2(Form("dgaus"), "gaus  + [3]",
           -fgaus.GetParameter(2) * 1.5, fgaus.GetParameter(2) * 1.5);

    f2.SetParameters(fgaus.GetParameter(0), fgaus.GetParameter(1),
                     fgaus.GetParameter(2), 0);

    h1->Fit(&f2, "MQR");

    //      new TCanvas;
    //      h1->Draw();
    //      fgaus.Draw("same");
    //      break;

    x[n] = p2->GetBinCenter(i);
    ex[n] = (p2->GetBinCenter(2) - p2->GetBinCenter(1)) / 2;
    y[n] = fgaus.GetParameter(1);
    ey[n] = fgaus.GetParError(1);

    //      p2->SetBinContent(i, fgaus.GetParameter(1));
    //      p2->SetBinError(i, fgaus.GetParameter(2));

    n++;
    delete h1;
  }

  return new TGraphErrors(n, x, y, ex, ey);
}

TFile *_file0 = NULL;
TString description;
TTree *T(nullptr);

void Resolution(const TCut &cut = "Iteration$ >= 5 && Iteration$ <= 10 && TPCTrack.nCluster>=12 && TPCTrack.clusterSizePhi > 3.5")
{
  TH2 *hresidual_phi = new TH2F("hresidual_phi", "hresidual_phi", 50, -2.9, -2.885, 60, -.2, .2);
  T->Draw("TPCTrack.clusterResidualPhi:TPCTrack.clusterProjectionPhi>>hresidual_phi",
          cut, "goff");
  hresidual_phi->SetTitle(";Global Phi [rad];Phi Residual [cm]");

  TH2 *hresidual_phi_highres = new TH2F("hresidual_phi_highres", "hresidual_phi", 500, -2.9, -2.885, 1000, -.2, .2);
  T->Draw("TPCTrack.clusterResidualPhi:TPCTrack.clusterProjectionPhi>>hresidual_phi_highres",
          cut, "goff");
  hresidual_phi_highres->SetTitle(";Global Phi [rad];Phi Residual [cm]");

  TCanvas *c1 = new TCanvas("Resolution", "Resolution", 1800, 860);
  c1->Divide(3, 1);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TGraphErrors *gcenter = FitProfile(hresidual_phi);

  hresidual_phi->Draw("colz");
  gcenter->Draw("p");
  //  hresidual_phi

  //  TF1 *fPeriod = new TF1("fPeriod",
  //                         "[0] * sin(2*pi*x/(2*pi/12/8/16)) + [1] * sin(2*pi*x/(2*pi/12/8/16)) + [2] + [3]*x + [4]*x*x + [5] * sin(4*pi*x/(2*pi/12/8/16)) + [6] * sin(4*pi*x/(2*pi/12/8/16))+ [7] * sin(3*pi*x/(2*pi/12/8/16)) + [8] * sin(3*pi*x/(2*pi/12/8/16))+ [9] * sin(pi*x/(5*pi/12/8/16)) + [10] * sin(pi*x/(5*pi/12/8/16)) + [11] * x*x*x", -3, 0);
  TF1 *fPeriod = new TF1("fPeriod",
                         "([0]+[1]*x) * sin(2*pi*x/(2*pi/12/8/16)) + ([1]+[2]*x) * cos(2*pi*x/(2*pi/12/8/16)) + [3] + [4]*x + [5]*x*x + ([6]+[7]*x) * sin(4*pi*x/(2*pi/12/8/16)) + ([8]+[9]*x) * cos(4*pi*x/(2*pi/12/8/16))+ ([10] + [11]*x) * sin(3*pi*x/(2*pi/12/8/16)) + ([12]+[13]*x) * cos(3*pi*x/(2*pi/12/8/16))+ ([14]+[15]*x) * sin(pi*x/(5*pi/12/8/16)) + ([16]+[17]*x) * cos(pi*x/(5*pi/12/8/16)) + ([18]+[19]*x) * x*x*x", -3, 0);
  gcenter->Fit(fPeriod, "M");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH2 *hresidual_phi_cor = (TH2 *) hresidual_phi->Clone("hresidual_phi_cor");
  hresidual_phi_cor->Reset();
  hresidual_phi_cor->SetTitle(";Global Phi [rad];Modulation-Corrected Phi Residual [cm]");

  double sum = 0;
  for (int binx = 1; binx <= hresidual_phi_highres->GetNbinsX(); ++binx)
    for (int biny = 1; biny <= hresidual_phi_highres->GetNbinsY(); ++biny)
    {
      const double x = hresidual_phi_highres->GetXaxis()->GetBinCenter(binx);
      const double y = hresidual_phi_highres->GetYaxis()->GetBinCenter(biny);
      const double value = hresidual_phi_highres->GetBinContent(binx, biny);

      sum += value;
      //      const double y_cor = y;
      const double y_cor = y - fPeriod->Eval(x);

      hresidual_phi_cor->Fill(x, y_cor, value);
    }
  cout << __PRETTY_FUNCTION__ << " sum = " << sum << endl;

  hresidual_phi_cor->Draw("colz");

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  TH1 *hresidual_cor = hresidual_phi_cor->ProjectionY("hresidual_cor", 1, hresidual_phi_cor->GetNbinsX());

  hresidual_cor->Draw();

  TH1 *hresidual = hresidual_phi->ProjectionY("hresidual", 1, hresidual_phi->GetNbinsX());
  //  hresidual->SetMarkerColor(kRed - 3);
  hresidual->SetLineColor(kRed - 3);
  hresidual->Draw("same");

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString("_DrawJet_") + TString(c1->GetName()), kFALSE);
}

void TrackQA()
{
  TCanvas *c1 = new TCanvas("TrackQA", "TrackQA", 1800, 860);
  c1->Divide(3, 2);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *hnTrack = new TH1F("hnTrack", ";# track per event", 17, -.5, 16.5);
  T->Draw("nTrack>>hnTrack");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *hCluster = new TH1F("hCluster", ";# cluster per track", 16, .5, 16.5);
  T->Draw("TPCTrack.nCluster>>hCluster");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *hclusterSizePhi = new TH1F("hclusterSizePhi", ";Cluster phi size", 16, .5, 16.5);
  T->Draw("TPCTrack.clusterSizePhi>>hclusterSizePhi");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *hclusterProjectionPhi = new TH1F("hclusterProjectionPhi", ";Cluster phi [rad]", 100, -3.5, -2.5);
  T->Draw("TPCTrack.clusterProjectionPhi>>hclusterProjectionPhi");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *hAngle = new TH1F("hAngle", ";Horizontal angle [degree]", 100, -30, 30);
  T->Draw("SvtxTrack.get_pz()/SvtxTrack.get_px()/pi*180>>hAngle");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *hresidualRough = new TH1F("hresidualRough", ";Rought phi residual [cm]", 1000, -1, 1);
  T->Draw("TPCTrack.clusterResidualPhi>>hresidualRough", "Iteration$ >= 5 && Iteration$ <= 10 && TPCTrack.nCluster>=10");

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString("_DrawJet_") + TString(c1->GetName()), kFALSE);
}

void TrackDistortion(const TCut &cut = "TPCTrack.nCluster>=10")
{
  TCanvas *c1 = new TCanvas("TrackDistortion", "TrackDistortion", 1800, 860);
  c1->Divide(2, 1);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH2 *hPhiDistortion = new TH2F("hPhiDistortion", ";Pad Layers;Phi Residual [cm]", 16, -.5, 15.5, 200, -3, 3);
  T->Draw("TPCTrack.clusterResidualPhi:Iteration$>>hPhiDistortion", cut, "goff");
  TGraph *gPhiDistortion = FitProfile(hPhiDistortion);
  hPhiDistortion->Draw("colz");
  gPhiDistortion->Draw("p");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH2 *hZDistortion = new TH2F("hZDistortion", ";Pad Layers;Z Residual [cm]", 16, -.5, 15.5, 200, -3, 3);
  T->Draw("TPCTrack.clusterResidualZ:Iteration$>>hZDistortion", cut, "goff");
  TGraph *gZDistortion = FitProfile(hZDistortion);
  hZDistortion->Draw("colz");
  gZDistortion->Draw("p");

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString("_DrawJet_") + TString(c1->GetName()), kFALSE);
}

void DrawTpcPrototypeGenFitTrkFitter(
    const TString infile = "data/tpc_beam/tpc_beam_00000292-0000.evt_TpcPrototypeGenFitTrkFitter.root", const TString desc = "Run 292"  //
)
{
  gSystem->Load("libtpc2019.so");
  //
  SetsPhenixStyle();
  TVirtualFitter::SetDefaultFitter("Minuit2");
  gStyle->SetLegendTextSize(0);
  //
  if (!_file0)
  {
    TString chian_str = infile;
    chian_str.ReplaceAll("ALL", "*");

    TChain *t = new TChain("T");
    const int n = t->Add(chian_str);

    cout << "Loaded " << n << " root files with eventT in " << chian_str << endl;
    assert(n > 0);

    T = t;

    _file0 = new TFile;
    _file0->SetName(infile);
  }
  //
  description = desc;

  TrackQA();
  TrackDistortion();
  //  Resolution();
  Resolution("Iteration$ >= 7 && Iteration$ <= 8 && TPCTrack.nCluster>=14 && TPCTrack.clusterSizePhi > 3.5");
}

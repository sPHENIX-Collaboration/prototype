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

void Resolution(const TCut &cut = "Iteration$ >= 5 && Iteration$ <= 10 && TPCTrack.nCluster>=12 && TPCTrack.clusterSizePhi > 3.5",
                const double phi_start = -2.905, const double phi_end = -2.885

)
{
  TH2 *hresidual_phi = new TH2F("hresidual_phi", "hresidual_phi", 60, phi_start, phi_end, 60, -.2, .2);
  T->Draw("TPCTrack.clusterResidualPhi:TPCTrack.clusterProjectionPhi>>hresidual_phi",
          cut, "goff");
  hresidual_phi->SetTitle(";Global Phi [rad];Phi Residual [cm]");

  TH2 *hresidual_phi_highres = new TH2F("hresidual_phi_highres", "hresidual_phi", 500, phi_start, phi_end, 1000, -.2, .2);
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

  TF1 *hresidual_fgaus = new TF1("hresidual_fgaus", "gaus", -.2, .2);
  //
  TF1 *hresidual_dgaus = new TF1(Form("hresidual_dgaus"), "gaus + [3]*exp(-0.5*((x-[1])/[4])**2) + [5]",
                                 -.2, .2);
  //
  hresidual_fgaus->SetParameter(1, hresidual->GetMean());
  hresidual_fgaus->SetParameter(2, hresidual->GetRMS());
  //
  hresidual_cor->Fit(hresidual_fgaus, "MQ");
  //
  //  TF1 f2(Form("dgaus"), "gaus  + [3]",
  //         -fgaus.GetParameter(2) * 1.5, fgaus.GetParameter(2) * 1.5);
  //
  hresidual_dgaus->SetParameters(hresidual_fgaus->GetParameter(0),
                                 hresidual_fgaus->GetParameter(1),
                                 hresidual_fgaus->GetParameter(2) / 2,
                                 hresidual_fgaus->GetParameter(0) / 20,
                                 hresidual_fgaus->GetParameter(2) * 10, 0);
  //
  hresidual_cor->Fit(hresidual_dgaus, "M");
  hresidual_dgaus->SetLineColor(kBlue + 2);
  const double resolution = TMath::Min(hresidual_dgaus->GetParameter(2), hresidual_dgaus->GetParameter(4));

  cout << __PRETTY_FUNCTION__ << " hresidual_fgaus->GetParameter(2) = " << hresidual_dgaus->GetParameter(2)
       << "  hresidual_fgaus->GetParameter(4) = " << hresidual_dgaus->GetParameter(4) << " resolution = " << resolution << endl;
  //
  //  //      new TCanvas;
  //  //      h1->Draw();
  //  //      fgaus.Draw("same");
  //  //      break;
  //
  //  x[n] = p2->GetBinCenter(i);
  //  ex[n] = (p2->GetBinCenter(2) - p2->GetBinCenter(1)) / 2;
  //  y[n] = fgaus.GetParameter(1);
  //  ey[n] = fgaus.GetParError(1);

  gPad->SetTopMargin(.3);
  TLegend *leg = new TLegend(.03, .7, .98, .99, description + ": 1-removed residual");
  leg->AddEntry(hresidual, "Raw residual", "l");
  leg->AddEntry(hresidual_cor, "Position dependence corected residual", "p");
  leg->AddEntry(hresidual_dgaus, Form("resolution = %.0f #mum", resolution * 1e4), "l");
  leg->Draw();

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString(c1->GetName()), kFALSE);
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
  T->Draw("atan(TPCTrack.pz/ TPCTrack.px)/pi*180>>hAngle");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *hAngleV = new TH1F("hAngleV", ";Vertical angle [degree]", 100, -30, 30);
  T->Draw("(atan(TPCTrack.py/ TPCTrack.px) - 2*pi/12/2 )/pi*180>>hAngleV");
  //  TH1 *hresidualRough = new TH1F("hresidualRough", ";Rought phi residual [cm]", 1000, -1, 1);
  //  T->Draw("TPCTrack.clusterResidualPhi>>hresidualRough", "Iteration$ >= 5 && Iteration$ <= 10 && TPCTrack.nCluster>=10");

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString(c1->GetName()), kFALSE);
}

void Track3D()
{
  TCanvas *c1 = new TCanvas("Track3D", "Track3D", 1800, 900);
  c1->Divide(2, 1);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogx();
  //  p->DrawFrame(0, 0, 10, 2, ";Transverse Momentum, p_{T} [GeV/c];Nuclear Modification Factor, R_{AA}");
  //
  TH3F *h3ClusterOverlay = new TH3F("h3ClusterOverlay", "h3ClusterOverlay", 128, -40, 10, 64, 35, 65, 128, -15, 15);
  //
  T->SetAlias("PhiCenter", "pi/12 + pi");  // center line azimuthal angle for TPC sector 0
  T->Draw("TPCTrack.clusterX*cos(PhiCenter + pi/2) + TPCTrack.clusterY*sin(PhiCenter + pi/2):TPCTrack.clusterX*cos(PhiCenter) + TPCTrack.clusterY*sin(PhiCenter):TPCTrack.clusterZ>>h3ClusterOverlay", "", "BOX2");
  //               "TPCTrack.clusterE", "BOX2");
  h3ClusterOverlay->SetTitle(";Drift Direction [cm];Pad Row Direction [cm];Azimuth Direction [cm]");
  h3ClusterOverlay->SetLineWidth(0);
  //
  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogx();
  //  p->DrawFrame(0, 0, 10, 2, ";Transverse Momentum, p_{T} [GeV/c];Nuclear Modification Factor, R_{AA}");

  TH1 *h2ClusterOverlay = h3ClusterOverlay->Project3D("zx");
  h2ClusterOverlay->Draw("colz");
  //    eventT->Draw("Clusters.avg_pady:Clusters.min_sample+Clusters.peak_sample>>h2ClusterOverlay",
  //                 "Clusters.peak", "colz");
  //    h2ClusterOverlay->SetTitle(";Time [0-127*50ns];Pads [0-127]");
  //  h2ClusterOverlay->SetLineWidth(0);

  p->SetTopMargin(.9);
  TLegend *leg = new TLegend(.15, .9, .95, .99, description + ": accumulated clusters on tracks");
  leg->Draw();

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString(c1->GetName()), false);
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
             TString(_file0->GetName()) + TString(c1->GetName()), kFALSE);
}

Double_t langaufun(Double_t *x, Double_t *par)
{
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;  // (2 pi)^(-1/2)
  Double_t mpshift = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;  // number of convolution steps
  Double_t sc = 5.0;    // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp - xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for (i = 1.0; i <= np / 2; i++)
  {
    xx = xlow + (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);

    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1 *his,
               const vector<double> &fitrange,
               const vector<double> &startvalues,
               const vector<double> &parlimitslo, const vector<double> &parlimitshi
               //    ,
               //    const vector<double> &fitparams, const vector<double> &fiterrors,
               //    const vector<double> &ChiSqr, Int_t *NDF
)
{
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[4]  reasonable start values for the fit
  //   parlimitslo[4]  lower parameter limits
  //   parlimitshi[4]  upper parameter limits
  //   fitparams[4]    returns the final fit parameters
  //   fiterrors[4]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf

  Int_t i;
  Char_t FunName[100];

  sprintf(FunName, "Fitfcn_%s", his->GetName());

  TF1 *ffitold = (TF1 *) gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName, langaufun, fitrange[0], fitrange[1], 4);
  ffit->SetParameters(startvalues.data());
  ffit->SetParNames("Width", "MP", "Area", "GSigma");

  for (i = 0; i < 4; i++)
  {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  his->Fit(FunName, "RB0");  // fit within specified range, use ParLimits, do not plot

  //   ffit->GetParameters(fitparams);    // obtain fit parameters
  //   for (i=0; i<4; i++) {
  //      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  //   }
  //   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  //   NDF[0] = ffit->GetNDF();           // obtain ndf

  return (ffit);  // return fit function
}

void TrackClusterEnergy(const TCut &cut = "TPCTrack.nCluster>=10")
{
  TCanvas *c1 = new TCanvas("TrackClusterEnergy", "TrackClusterEnergy", 1200, 860);
  //  c1->Divide(2, 1);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  p->SetLogy();

  TH1 *hClusterEnergy = new TH1F("hClusterEnergy", ";Cluster Energy on good track [ADU];Count / bin", 2500, 0, 5000);
  T->Draw("TPCTrack.clusterE>>hClusterEnergy", cut, "goff");
  //  TGraph *gPhiDistortion = FitProfile(hPhiDistortion);

  TF1 *fit = langaufit(hClusterEnergy,
                       {0, 5000},
                       {200, 500, (double) T->GetEntries(cut) * 10, 200},
                       {0, 0, 0, 0},
                       {5000, 5000, 1e10, 5000});
  fit->SetLineColor(kBlue + 2);

  hClusterEnergy->Draw();
  fit->Draw("same");
  //    gPhiDistortion->Draw("p");

  TLegend *leg = new TLegend(.3, .7, .95, .9, description + ": Cluster Energy on good track");
  leg->AddEntry(hClusterEnergy, TString("Data: ") + cut.GetTitle(), "l");
  leg->AddEntry(fit,
                Form("Langdau * Gauss Fit, #mu= %.0f ADU", fit->GetParameter(1)), "l");
  leg->Draw();

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString(c1->GetName()), kFALSE);
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
  Track3D();
  TrackClusterEnergy();
  //  Resolution();
  Resolution("Iteration$ >= 5 && Iteration$ <= 10 && TPCTrack.nCluster>=12 && TPCTrack.clusterSizePhi > 7");
}

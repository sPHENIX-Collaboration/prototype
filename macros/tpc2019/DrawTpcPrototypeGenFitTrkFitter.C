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

void Resolution()
{
  TH2 *hresidual_phi = new TH2F("hresidual_phi", "hresidual_phi", 50, -2.905, -2.884, 100, -.2, .2);
  T->Draw("TPCTrack.clusterResidualPhi:TPCTrack.clusterProjectionPhi>>hresidual_phi",
          "Iteration$ >= 5 && Iteration$ <= 10 && TPCTrack.nCluster>=12 && abs(TPCTrack.clusterSizePhi - 4)<1", "goff");
  hresidual_phi->SetTitle(";Global Phi [rad];Phi Residual [cm]");

  TH2 *hresidual_phi_highres = new TH2F("hresidual_phi", "hresidual_phi", 500, -2.905, -2.884, 1000, -.2, .2);
  T->Draw("TPCTrack.clusterResidualPhi:TPCTrack.clusterProjectionPhi>>hresidual_phi",
          "Iteration$ >= 5 && Iteration$ <= 10 && TPCTrack.nCluster>=12 && abs(TPCTrack.clusterSizePhi - 4)<1", "goff");
  hresidual_phi->SetTitle(";Global Phi [rad];Phi Residual [cm]");

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

  TF1 *fPeriod = new TF1("fPeriod",
                         "[0] * sin(2*pi*x/(2*pi/12/8/16)) + [1] * sin(2*pi*x/(2*pi/12/8/16)) + [2] + [3]*x + [4]*x*x + [5] * sin(4*pi*x/(2*pi/12/8/16)) + [6] * sin(4*pi*x/(2*pi/12/8/16))+ [7] * sin(3*pi*x/(2*pi/12/8/16)) + [8] * sin(3*pi*x/(2*pi/12/8/16))+ [9] * sin(pi*x/(5*pi/12/8/16)) + [10] * sin(pi*x/(5*pi/12/8/16)) + [11] * x*x*x", -3, 0);
  gcenter->Fit(fPeriod, "M");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH2 *hresidual_phi_cor = (TH2 *) hresidual_phi->Clone("hresidual_phi_cor");
  hresidual_phi_cor->Reset();
  hresidual_phi_cor->SetTitle(";Global Phi [rad];Corrected Phi Residual [cm]");

  for (int binx = 1; binx <= hresidual_phi_highres->GetNbinsX(); ++binx)
    for (int biny = 1; biny <= hresidual_phi_highres->GetNbinsY(); ++biny)
    {
      const double x = hresidual_phi_highres->GetXaxis()->GetBinCenter(binx);
      const double y = hresidual_phi_highres->GetYaxis()->GetBinCenter(biny);

      //      const double y_cor = y;
      const double y_cor = y - fPeriod->Eval(x);

      hresidual_phi_cor->Fill(x, y_cor,
                              hresidual_phi_highres->GetBinContent(binx, biny));
    }

  hresidual_phi_cor->Draw("colz");

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  TH1 *hresidual_cor = hresidual_phi_cor->ProjectionY("hresidual_cor");

  hresidual_cor->Draw();

  TH1 *hresidual = hresidual_phi->ProjectionY("hresidual");

  hresidual_cor->Draw("same");

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

  Resolution();
}

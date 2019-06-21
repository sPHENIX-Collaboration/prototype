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

TFile *_file0 = NULL;
TString description;
TTree *eventT(nullptr);
TTree *chanT(nullptr);

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
    y[n] = fgaus.GetParameter(2);
    ey[n] = fgaus.GetParError(2);

    //      p2->SetBinContent(i, fgaus.GetParameter(1));
    //      p2->SetBinError(i, fgaus.GetParameter(2));

    n++;
    delete h1;
  }

  return new TGraphErrors(n, x, y, ex, ey);
}

void Clusters3D()
{
  TCanvas *c1 = new TCanvas("Clusters3D", "Clusters3D", 1200, 900);
  c1->Divide(1, 1);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogx();
  //  p->DrawFrame(0, 0, 10, 2, ";Transverse Momentum, p_{T} [GeV/c];Nuclear Modification Factor, R_{AA}");

  TH3F *h3ClusterOverlay = new TH3F("h3ClusterOverlay", "h3ClusterOverlay", 128, -.5, 127.5, 16, -.5, 15.5, 128, -.5, 127.5);

  eventT->Draw("Clusters.avg_pady:Clusters.avg_padx:Clusters.min_sample+Clusters.peak_sample>>h3ClusterOverlay",
               "Clusters.peak", "BOX2");
  h3ClusterOverlay->SetTitle("Run274 Event 1-100 overlay of Clusters;Time [0-127*50ns];Rows [0-15];Pads [0-127]");
  h3ClusterOverlay->SetLineWidth(0);

  TLegend *leg = new TLegend(.05, .9, .95, .99, description + ": accumulated clusters");
  //  leg->AddEntry("", "Au+Au #sqrt{s_{NN}}=200GeV, 24B 0-10%C", "");
  //  leg->AddEntry("", "p+p #sqrt{s}=200GeV, 200B M.B.", "");
  leg->Draw();

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString("_") + TString("Clusters3D") + TString("_") + TString(c1->GetName()), false);
}

void ClusterQA()
{
  TCanvas *c1 = new TCanvas("ClusterQA", "ClusterQA", 1800, 860);
  c1->Divide(3, 1);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogx();
  //  p->DrawFrame(0, 0, 10, 2, ";Transverse Momentum, p_{T} [GeV/c];Nuclear Modification Factor, R_{AA}");

  TH1 *CLusterEnergyAvg = new TH1F("CLusterEnergyAvg", ";<cluster energy> [ADU];Count", 100, -0, 2000);
  eventT->Draw("Sum$(Clusters[].peak)/nClusters>>CLusterEnergyAvg");
  CLusterEnergyAvg->Fit("gaus");
  TF1 *fgauss = (TF1 *) CLusterEnergyAvg->GetListOfFunctions()->At(0);
  assert(fgauss);
  fgauss->SetLineColor(kBlue + 1);

  TLegend *leg = new TLegend(.2, .8, .99, .99, description);
  //    leg->AddEntry("", "Au+Au #sqrt{s_{NN}}=200GeV, 24B 0-10%C", "");
  leg->AddEntry(fgauss, Form("Energy = %.1f +/- %.1f", fgauss->GetParameter(1), fgauss->GetParError(1)), "l");
  leg->Draw();

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogx();
  //  p->DrawFrame(0, 0, 10, 2, ";Transverse Momentum, p_{T} [GeV/c];Nuclear Modification Factor, R_{AA}");

  TH1 *CLusterEnergyAll = new TH1F("CLusterEnergyAll", ";cluster energy [ADU];Count", 100, -0, 2000);
  eventT->Draw("(Clusters[].peak)>>CLusterEnergyAll");
  //  CLusterEnergyAvg->Fit("gaus");

  p = (TPad *) c1->cd(idx++);
  c1->Update();
  //  p->SetLogx();
  //  p->DrawFrame(0, 0, 10, 2, ";Transverse Momentum, p_{T} [GeV/c];Nuclear Modification Factor, R_{AA}");

  TH1 *nClusters = new TH1F("nClusters", ";Cluster count;Count", 151, -.5, 150.5);
  eventT->Draw("nClusters>>nClusters");
  //  CLusterEnergyAvg->Fit("gaus");

  SaveCanvas(c1,
             TString(_file0->GetName()) + TString("_") + TString("Clusters3D") + TString("_") + TString(c1->GetName()), false);
}

void DrawTpcPrototypeUnpacker(
    const TString infile = "data/tpc_beam/tpc_beam_00000292-0000.evt_TpcPrototypeUnpacker.root", const TString desc = "Run 292"  //
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

    TChain *t = new TChain("eventT");
    const int n = t->Add(chian_str);

    cout << "Loaded " << n << " root files with eventT in " << chian_str << endl;
    assert(n > 0);

    eventT = t;

    t = new TChain("chanT");
    const int n1 = t->Add(chian_str);

    cout << "Loaded " << n1 << " root files with chanT in " << chian_str << endl;
    assert(n1 > 0);

    chanT = t;

    _file0 = new TFile;
    _file0->SetName(infile);
  }
  //
  description = desc;

  Clusters3D();
  ClusterQA();
}

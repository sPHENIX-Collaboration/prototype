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

string basePath;
const static double x0 = 6 - (2 + 7. / 16);  // inch

pair<TGraph *, TGraph *> getResolutions(vector<pair<double, string>> datafiles, Color_t color = kRed + 1)
{
  vector<double> x_cm;
  vector<double> x;
  vector<double> res;
  vector<double> res_err;
  vector<double> res2;
  vector<double> res2_err;

  for (const auto &data : datafiles)
  {
    cout << __PRETTY_FUNCTION__ << " "
         << " processing x = " << data.first << " from " << data.second << endl;

    TFile *f = TFile::Open((TString(basePath + data.second.c_str()) + "_TpcPrototypeGenFitTrkFitter.root_DrawJet_Resolution.root"));
    assert(f);
    assert(f->IsOpen());
    TH1 *hresidual_cor = (TH1 *) f->GetObjectChecked("hresidual_cor", "TH1");
    assert(hresidual_cor);
    TF1 *hresidual_dgaus = (TF1 *) hresidual_cor->GetListOfFunctions()->At(0);
    assert(hresidual_dgaus);

    if (hresidual_dgaus->GetParameter(2) < hresidual_dgaus->GetParameter(4))
    {
      res.push_back(hresidual_dgaus->GetParameter(2) * 1e4);
      res_err.push_back(hresidual_dgaus->GetParError(2) * 1e4);
    }
    else
    {
      res.push_back(hresidual_dgaus->GetParameter(4) * 1e4);
      res_err.push_back(hresidual_dgaus->GetParError(4) * 1e4);
    }
    x.push_back(data.first);
    x_cm.push_back((data.first - x0) * 2.54);
    res2.push_back(res[res.size() - 1] * res[res.size() - 1]);
    res2_err.push_back(2 * res[res.size() - 1] * res_err[res.size() - 1]);

    cout << "x = " << x[res.size() - 1] << " x_cm = " << x_cm[res.size() - 1] << " res = " << res[res.size() - 1] << " res_err = " << res_err[res.size() - 1] << endl;
  }

  TGraphErrors *g = new TGraphErrors(x.size(), x_cm.data(), res.data(), 0, res_err.data());
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->Print();

  TGraphErrors *g2 = new TGraphErrors(x.size(), x.data(), res2.data(), 0, res2_err.data());
  g2->SetMarkerColor(color);
  g2->SetLineColor(color);
  g2->Print();

  return pair<TGraph *, TGraph *>(g, g2);
}

void DrawTpcPrototypeGenFitTrkFitter_Summary(  //
    string basepath = "/phenix/u/jinhuang/links/sPHENIX_work/TPC/fnal_June2019/SimPadPlaneIter4/")
{
  //  gSystem->Load("libtpc2019.so");
  //
  basePath = basepath;
  const double max_res = 300;

  SetsPhenixStyle();
  TVirtualFitter::SetDefaultFitter("Minuit2");
  gStyle->SetLegendTextSize(0);

  TCanvas *c1 = new TCanvas("DrawTpcPrototypeGenFitTrkFitter_Summary", "DrawTpcPrototypeGenFitTrkFitter_Summary", 1800, 700);
  c1->Divide(2, 1);
  int idx = 1;
  TPad *p = nullptr;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  TH1 *frame = p->DrawFrame(0, 0, 20, max_res * max_res);
  frame->SetTitle(";Horizontal Position [in];Resolution^{2} [#mum^{2}]");
  frame->GetYaxis()->SetTitleOffset(1.75);
  TLine *l = new TLine(x0, 0, x0, max_res * max_res);
  l->Draw();

  TLegend *leg1 = new TLegend(.35, .7, .9, .9, "2019 TPC Beam Test Preview");
  TLegend *leg2 = new TLegend(.3, .2, .9, .5, "2019 TPC Beam Test Preview");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  p->DrawFrame(0, 0, 50, max_res)->SetTitle(";Distance to readout, L [cm];Resolution [#mum]");

  {
    Color_t color = kRed + 2;
    const TString name = "Scan2: Gap 300V, GEM 370V";

    pair<TGraph *, TGraph *> scan2 = getResolutions(
        {{6, "tpc_beam_00000288-0000.evt"},
         {6, "tpc_beam_00000292-0000.evt"},
         {10, "tpc_beam_00000293-0000.evt"},
         {14, "tpc_beam_00000294-0000.evt"},
         {18, "tpc_beam_00000295-0000.evt"}},
        color);

    c1->cd(1);
    scan2.second->Draw("p");
    TF1 *fline_scan2 = new TF1("fline_scan2", "pol1", x0, 16);
    scan2.second->Fit(fline_scan2, "MR0");
    fline_scan2->SetLineColor(color);
    fline_scan2->Draw("same");

    leg1->AddEntry(scan2.second, name, "p");

    c1->cd(2);
    scan2.first->Draw("p");
    TF1 *fdiff_scan2 = new TF1("fline_scan2", "sqrt([0]*[0] + [1]*[1]*x)", 0, 30);
    scan2.first->Fit(fdiff_scan2, "MR0");
    fdiff_scan2->SetLineColor(color);
    fdiff_scan2->Draw("same");

    leg2->AddEntry(scan2.first, name, "p");
    leg2->AddEntry(fdiff_scan2, Form("#sqrt{(%.1f #mum)^{2} + (%.1f #mum/#sqrt{cm})^{2} L}", abs(fdiff_scan2->GetParameter(0)), fdiff_scan2->GetParameter(1)), "l");

    c1->Update();
  }

  {
    Color_t color = kBlue + 2;
    const TString name = "Scan3: Gap 150V, GEM 380V";

    pair<TGraph *, TGraph *> scan3 = getResolutions(
        {{6, "tpc_beam_00000297-0000.evt"},
         {6, "tpc_beam_00000298-0000.evt"},
         {10, "tpc_beam_00000299-0000.evt"},
         {14, "tpc_beam_00000300-0000.evt"},
         {18, "tpc_beam_00000301-0000.evt"}},
        color);

    c1->cd(1);
    scan3.second->Draw("p");
    TF1 *fline_scan3 = new TF1("fline_scan3", "pol1", x0, 16);
    scan3.second->Fit(fline_scan3, "MR0");
    fline_scan3->SetLineColor(color);
    fline_scan3->Draw("same");
    leg1->AddEntry(scan3.second, name, "p");

    c1->cd(2);
    scan3.first->Draw("p");
    TF1 *fdiff_scan3 = new TF1("fline_scan3", "sqrt([0]*[0] + [1]*[1]*x)", 0, 30);
    scan3.first->Fit(fdiff_scan3, "MR0");
    fdiff_scan3->SetLineColor(color);
    fdiff_scan3->Draw("same");

    leg2->AddEntry(scan3.first, name, "p");
    leg2->AddEntry(fdiff_scan3, Form("#sqrt{(%.1f #mum)^{2} + (%.1f #mum/#sqrt{cm})^{2} L}", abs(fdiff_scan3->GetParameter(0)), fdiff_scan3->GetParameter(1)), "l");
    c1->Update();
  }

  if (0)
  {
    Color_t color = kGreen + 2;

    pair<TGraph *, TGraph *> scan4 = getResolutions(
        {{6, "tpc_beam_00000357-0000.evt"},
         {6, "tpc_beam_00000358-0000.evt"},
         {10, "tpc_beam_00000359-0000.evt"},
         {10, "tpc_beam_00000360-0000.evt"},
         {13, "tpc_beam_00000356-0000.evt"},
         {14, "tpc_beam_00000361-0000.evt"},
         {14, "tpc_beam_00000362-0000.evt"},
         {18, "tpc_beam_00000363-0000.evt"},
         {18, "tpc_beam_00000366-0000.evt"}},
        color);

    c1->cd(1);
    scan4.second->Draw("p");
    TF1 *fline_scan4 = new TF1("fline_scan4", "pol1", x0, 16);
    scan4.second->Fit(fline_scan4, "MR0");
    fline_scan4->SetLineColor(color);
    fline_scan4->Draw("same");

    c1->cd(2);
    scan4.first->Draw("p");
    TF1 *fdiff_scan4 = new TF1("fline_scan4", "sqrt([0]*[0] + [1]*[1]*x)", 0, 30);
    scan4.first->Fit(fdiff_scan4, "MR0");
    fdiff_scan4->SetLineColor(color);
    fdiff_scan4->Draw("same");

    c1->Update();
  }

  c1->cd(1);
  leg1->Draw();
  c1->cd(2);
  leg2->Draw();
  SaveCanvas(c1,
             TString(basePath) + TString(c1->GetName()), true);
}

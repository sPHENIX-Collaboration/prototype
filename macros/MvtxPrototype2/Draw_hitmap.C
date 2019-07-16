#include "Style.h"

void Draw_hitmap(int runnum=732, const bool bSAVE=false){

  gStyle->SetOptStat(0);

  const int nstave = 4;
  const int nchip = 9;

  TFile *infile = new TFile(Form("MvtxPrototype2Eval-%08d-0200.root",runnum),"read");
  assert(infile);

  TH2D *h2d_hit[nstave];

  for (int is=0; is<nstave; is++){
    h2d_hit[is] = (TH2D*)infile->Get(Form("h2d_hit_stave%d",is));
    assert(h2d_hit[is])
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 300*nchip, 150*nstave);
  c1->Divide(1,nstave);

  for (int istave=0; istave<nstave; istave++){
    c1->cd(istave + 1);
    SetPadStyle();
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.05);

    htmp = (TH1D*)h2d_hit[istave];
    SetHistoStyle("COL","ROW","", 0.07, 0.06);
    h2d_hit[istave]->Draw("box");

    TLatex *tex = new TLatex(100,400,Form("Stave%d", istave));
    tex->SetTextSize(0.07);
    tex->Draw();
  }

  if (bSAVE){
    c1->cd();
    c1->SaveAs(Form("./20190503/hitmap_%08d.pdf", runnum));
  }
}

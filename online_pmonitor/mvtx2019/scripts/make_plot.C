make_plot(TString filename){
gSystem->Load("libmvtx.so");

c1 = new TCanvas("c1","c1",1200,900);

char filein[500];
char fileout[500];

sprintf(fileout,"%s.pdf[",filename.Data());
c1->Print(fileout);

sprintf(filein,"reference_datafiles/%s",filename.Data());

pfileopen(filein);

prun(5000);
c1->Divide(2,2);
c1->cd(1);
hchip->Scale(1.0/hnevent->Integral());
hchip->Draw();
c1->cd(2);
gStyle->SetOptStat(110011);
hnhit_chip_0->Draw();
hnhit_chip_0->GetXaxis()->SetRange(0,50);
hnhit_chip_1->Draw("same");
gPad->SetLogy(1);
c1->cd(3);
h2d_chip_0->Draw("colz");
c1->cd(4);
h2d_chip_1->Draw("colz");

sprintf(fileout,"%s.pdf",filename.Data());
c1->Print(fileout);

//c1->Divide(2,1);
c1->Clear();
c1->Divide(1,2);
c1->cd(1);
hnevent->Draw();

c1->cd(2);
//hnevent->Draw();
hhittime_chip_0->Draw();
hhittime_chip_1->Draw("same");
//c1->cd();

sprintf(fileout,"%s.pdf",filename.Data());
c1->Print(fileout);

sprintf(fileout,"%s.pdf]",filename.Data());
c1->Print(fileout);
}


void Draw_hitmap(int runnum=732, const bool bSAVE=false){

	gStyle->SetOptStat(0);

	const int nstave = 4;
	const int nchip = 9;

	TFile *infile = new TFile(Form("MvtxPrototype2Eval-%08d-0200.root",runnum),"read");

	TH2D *h2d_hit[nstave][nchip];

	for (int is=0; is<nstave; is++){
		for (int ic=0; ic<nchip; ic++){

			h2d_hit[is][ic] = (TH2D*)infile->Get(Form("h2d_hit_stave%d_chip%d",is,ic));

		}
	}

	TCanvas *c1 = new TCanvas("c1","c1",300*nchip,150*nstave);
	c1->Divide(nchip, nstave);

	for (int istave=0; istave<nstave; istave++){
		for (int ichip=0; ichip<nchip; ichip++){

			c1->cd(nchip*istave + ichip + 1);
			SetPadStyle();
			gPad->SetRightMargin(0.01);
			gPad->SetTopMargin(0.05);

			htmp = (TH1D*)h2d_hit[istave][ichip];
			SetHistoStyle("COL","ROW","",0.07,0.06);
			h2d_hit[istave][ichip]->Draw("box");

			TLatex *tex = new TLatex(100,400,Form("Stave%d, Chip%d",istave,ichip));
			tex->SetTextSize(0.07);
			tex->Draw();

		}
	}

	if ( bSAVE ){
		c1->cd();
		c1->SaveAs(Form("./20190503/hitmap_%08d.pdf",runnum));
	}


}

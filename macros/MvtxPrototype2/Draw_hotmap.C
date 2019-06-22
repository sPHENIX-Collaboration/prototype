#include "/phenix/u/shlim/Style.h"

void Draw_hotmap(const int runnumber=28, const int segnumber=100){

	gStyle->SetOptStat(0);

	const int nColor[4] = {1, 2, 4, 8};

	bool bWRITE = true;
	double cut_val = 2.0e-4;

	if ( runnumber==41 || runnumber==42 || runnumber==44 || runnumber==49 || runnumber==50 || runnumber==51 || runnumber==52  
			|| runnumber==78 || runnumber==79 || runnumber==82 ) cut_val = 9.0e-4;
	else if ( runnumber==43 || runnumber==77 || (runnumber>=126 && runnumber<=127) ) cut_val = 2.0e-3;
	else if ( runnumber==46 || runnumber==47 || runnumber==48 || runnumber==98 ) cut_val = 1.0e-3;
	else if ( runnumber==53 || runnumber==54 || (runnumber>=56 && runnumber<=65) || runnumber==97 || runnumber==116 ) cut_val = 4.0e-4;
	else if ( runnumber==121 ) cut_val = 1e10;

	TFile *infile = new TFile(Form("MvtxQAHisto-%08d-%04d.root",runnumber,segnumber),"read");

	TH2F *h2d_hit[4];
	TH2F *h2d_hit_trk[4];
	TH1F *h1d_hit_per_evt[4];
	TH1F *h1d_hit_x[4];
	TH1F *h1d_hit_y[4];
	TH1F *h1d_hit_trk_x[4];
	TH1F *h1d_hit_trk_y[4];

	TH1F *h1d_rate[4];
	TH2F *h2d_hit_cut[4];

	int nbad[4] = {0};

	for (int ichip=0; ichip<4; ichip++){
		h2d_hit[ichip] = (TH2F*)infile->Get(Form("h2d_hit_chip%d",ichip));
		h2d_hit_trk[ichip] = (TH2F*)infile->Get(Form("h2d_hit_trk_chip%d",ichip));
		h1d_hit_per_evt[ichip] = (TH1F*)infile->Get(Form("h1d_hit_per_evt_chip%d",ichip));

		h1d_hit_x[ichip] = (TH1F*)h2d_hit[ichip]->ProjectionX(Form("h1d_hit_x_chip%d",ichip));
		h1d_hit_y[ichip] = (TH1F*)h2d_hit[ichip]->ProjectionY(Form("h1d_hit_y_chip%d",ichip));
		h1d_hit_x[ichip]->SetLineColor(1);
		h1d_hit_y[ichip]->SetLineColor(1);

		h1d_hit_trk_x[ichip] = (TH1F*)h2d_hit_trk[ichip]->ProjectionX(Form("h1d_hit_trk_x_chip%d",ichip));
		h1d_hit_trk_y[ichip] = (TH1F*)h2d_hit_trk[ichip]->ProjectionY(Form("h1d_hit_trk_y_chip%d",ichip));
		h1d_hit_trk_x[ichip]->SetLineColor(1);
		h1d_hit_trk_y[ichip]->SetLineColor(1);

		h1d_rate[ichip] = new TH1F(Form("h1d_rate_chip%d",ichip),"",5000,0.000005,0.1);
	}

	h1d_hit_norm = (TH1F*)infile->Get("h1d_hit_norm");

	ofstream fpar;

	if ( bWRITE ){
		char fname[300];
		sprintf(fname,"hotmap/hotmap_testbeam_%08d.txt",runnumber);

		fpar.open(fname);
	}

	TCanvas *c1 = new TCanvas("c1","c1",150*1.1*4,150*4);
	c1->Divide(4,4);

	for (int ichip=0; ichip<4; ichip++){

		//double norm = h1d_hit_per_evt[ichip]->Integral(2,h1d_hit_per_evt[ichip]->GetNbinsX());
		double norm = h1d_hit_norm->Integral();
		cout << "CHIP: " << ichip << ", EVT: " << norm << endl;

		c1->cd(4*ichip+1);
		SetPadStyle();
		gPad->SetLogz();

		htmp = (TH1F*)h2d_hit[ichip];
		SetHistoStyle();
		htmp->GetXaxis()->SetLabelSize(0.08);

		htmp->Draw("col2");

		h2d_hit_cut[ichip] = (TH2F*)h2d_hit[ichip]->Clone(Form("h2d_hit_cut_%d",ichip));

		for (int ix=0; ix<h1d_hit_x[ichip]->GetNbinsX(); ix++){
			for (int iy=0; iy<h1d_hit_y[ichip]->GetNbinsX(); iy++){

				if ( h2d_hit[ichip]->GetBinContent(ix+1,iy+1)<1 ) continue;

				double val = h2d_hit[ichip]->GetBinContent(ix+1,iy+1)/norm;
				h1d_rate[ichip]->Fill(val);

				if ( val>cut_val ){
					h2d_hit_cut[ichip]->SetBinContent(ix+1,iy+1,0);
					//cout << ichip << " 0 " << ichip << " " << ix << " " << iy << endl;
				}else{
					//cout << val << endl;
				}

				if ( bWRITE ){
					if ( val>cut_val ){
						fpar << ichip << " 0 " << ichip << " " << ix << " " << iy << endl;
						nbad[ichip]++;
					}
				}//bWRITE

			}
		}


		c1->cd(4*ichip+2);
		SetPadStyle();
		gPad->SetLogy();
		gPad->SetLogx();

		htmp = (TH1F*)h1d_rate[ichip];
		SetHistoStyle();
		htmp->SetMinimum(0.5);
		htmp->GetXaxis()->SetLabelSize(0.08);
		htmp->Draw();

		TLine *line = new TLine(cut_val, 0.5, cut_val, htmp->GetMaximum());
		line->SetLineWidth(2);
		line->Draw();

		c1->cd(4*ichip+3);
		SetPadStyle();

		TH1F *h1 = (TH1F*)h2d_hit_cut[ichip]->ProjectionX(Form("h1d_hit_cut_x_%d",ichip));
		htmp = (TH1F*)h1;
		SetHistoStyle();
		htmp->GetXaxis()->SetLabelSize(0.08);
		h1->Draw("");

		c1->cd(4*ichip+4);
		SetPadStyle();

		h1 = (TH1F*)h2d_hit_cut[ichip]->ProjectionY(Form("h1d_hit_cut_y_%d",ichip));
		htmp = (TH1F*)h1;
		SetHistoStyle();
		htmp->GetXaxis()->SetLabelSize(0.08);
		h1->Draw("");

		cout << "N BAD CH: " << nbad[ichip] << endl;

	}//

	if ( bWRITE )
		fpar.close();

	return;

	TCanvas *c3 = new TCanvas("c3","c3",200*1.1*2,200*2);
	c3->Divide(2,2);

	for (int ichip=0; ichip<4; ichip++){
		c3->cd(ichip+1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1F*)gPad->DrawFrame(0,5,10,1.5*h1d_clus_size_x[ichip]->GetMaximum());
		SetHistoStyle();

		h1d_clus_size_x[ichip]->Draw("same");
		h1d_clus_size_z[ichip]->Draw("same");


	}


	/*
	TCanvas *c3 = new TCanvas("c3","c3",150*1.1*4,150*2);
	c3->Divide(4,2);

	for (int ichip=0; ichip<4; ichip++){

		c3->cd(ichip+1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1F*)h1d_hit_trk_x[ichip];
		SetHistoStyle();

		TF1 *fx = new TF1("fx","gaus",-5,5);
		htmp->Fit(fx,"R0Q");

		double sum = 0;
		for (int ii=0; ii<21; ii++){
			sum += fx->Eval(-10+ii);
		}

		cout << "Noise R x: " << (htmp->Integral()-sum)/htmp->Integral() << endl;

		htmp->SetAxisRange(-50,50);
		htmp->Draw("");
		fx->Draw("same");

		c3->cd(ichip+5);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1F*)h1d_hit_trk_y[ichip];
		SetHistoStyle();

		TF1 *fy = new TF1("fy","gaus",-5,5);
		htmp->Fit(fy,"R0Q");

		sum = 0;
		for (int ii=0; ii<21; ii++){
			sum += fy->Eval(-10+ii);
		}

		cout << "Noise R y: " << (htmp->Integral()-sum)/htmp->Integral() << endl;

		htmp->SetAxisRange(-50,50);
		htmp->Draw("");
		fy->Draw("same");

		//cout << "Noise R y: " << (htmp->GetEntries()-fy->Integral(-10,10))/htmp->GetEntries() << endl;
	}//

	int nevents = h1d_clus_associated->Integral(h1d_clus_associated->FindBin(3),h1d_clus_associated->GetNbinsX());

	h1d_clus_eff->Sumw2();
	h1d_clus_eff->Scale(1./nevents);

	TCanvas *c5 = new TCanvas("c5","c5",1.1*400,400);
	SetPadStyle();

	htmp = (TH1F*)gPad->DrawFrame(-0.5,0.95,3.5,1.02);
	SetHistoStyle();

	h1d_clus_eff->SetLineWidth(2);
	h1d_clus_eff->SetLineColor(1);
	h1d_clus_eff->Draw("same");

	*/


}

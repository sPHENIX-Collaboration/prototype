#include "/phenix/u/shlim/Style.h"

void Draw(const int runnumber=28, const int segnumber=200){

	int bPOST_ALIGN = false;
	int bPRE_ALIGN = false;

	float z_min = 200;
	float z_max = 800;
	float x_min = 200;
	float x_max = 500;

	gStyle->SetOptStat(0);

	const int nColor[4] = {1, 2, 4, 8};

	TFile *infile = new TFile(Form("MvtxQAHisto-%08d-%04d.root",runnumber,segnumber),"read");

	TH2F *h2d_hit[4];
	TH2F *h2d_hit_trk[4];
	TH1F *h1d_hit_per_evt[4];
	TH1F *h1d_hit_x[4];
	TH1F *h1d_hit_y[4];
	TH1F *h1d_hit_trk_x[4];
	TH1F *h1d_hit_trk_y[4];

	TH1F *h1d_clus_size_x[4];
	TH1F *h1d_clus_size_z[4];
	TH2F *h2d_clus[4];
	TH1F *h1d_clus_per_evt[4];

	TH1F *h1d_clus_res_x[4];
	TH1F *h1d_clus_res_z[4];
	TProfile *h1p_clus_res_x[4];
	TProfile *h1p_clus_res_z[4];

	TF1 *f1p_clus_res_x[4];
	TF1 *f1p_clus_res_z[4];


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

		h2d_clus[ichip] = (TH2F*)infile->Get(Form("h2d_clus_chip%d",ichip));

		h1d_clus_size_x[ichip] = (TH1F*)infile->Get(Form("h1d_clus_size_x_chip%d",ichip));
		h1d_clus_size_z[ichip] = (TH1F*)infile->Get(Form("h1d_clus_size_z_chip%d",ichip));

		h1d_clus_per_evt[ichip] = (TH1F*)infile->Get(Form("h1d_clus_per_evt_chip%d",ichip));

		h1d_clus_size_x[ichip]->SetLineColor(1);
		h1d_clus_size_z[ichip]->SetLineColor(2);

		h1d_clus_res_x[ichip] = (TH1F*)infile->Get(Form("h1d_clus_res_x_chip%d",ichip));
		h1d_clus_res_z[ichip] = (TH1F*)infile->Get(Form("h1d_clus_res_z_chip%d",ichip));

		h1d_clus_res_x[ichip]->SetLineColor(nColor[ichip]);
		h1d_clus_res_z[ichip]->SetLineColor(nColor[ichip]);

		h1p_clus_res_x[ichip] = (TProfile*)infile->Get(Form("h1p_clus_res_x_chip%d",ichip));
		h1p_clus_res_z[ichip] = (TProfile*)infile->Get(Form("h1p_clus_res_z_chip%d",ichip));

		h1p_clus_res_x[ichip]->SetLineColor(nColor[ichip]);
		h1p_clus_res_z[ichip]->SetLineColor(nColor[ichip]);
	}

	TH1F *h1d_clus_eff = (TH1F*)infile->Get("h1d_clus_eff");
	TH1F *h1d_clus_associated = (TH1F*)infile->Get("h1d_clus_associated");

	TCanvas *c0 = new TCanvas("c0","c0",1.1*2*400,400);
	c0->Divide(2,1);

	c0->cd(1);
	SetPadStyle();
	gPad->SetLogy();

	htmp = (TH1F*)gPad->DrawFrame(0,1,50,1.5*h1d_hit_per_evt[0]->GetMaximum());
	SetHistoStyle();

	for (int ichip=0; ichip<4; ichip++){
		h1d_hit_per_evt[ichip]->SetLineColor(nColor[ichip]);
		h1d_hit_per_evt[ichip]->SetLineWidth(2);
		h1d_hit_per_evt[ichip]->Draw("same");
	}

	c0->cd(2);
	SetPadStyle();
	gPad->SetLogy();

	htmp = (TH1F*)gPad->DrawFrame(0,1,10,1.5*h1d_clus_per_evt[0]->GetMaximum());
	SetHistoStyle();

	for (int ichip=0; ichip<4; ichip++){
		h1d_clus_per_evt[ichip]->SetLineColor(nColor[ichip]);
		h1d_clus_per_evt[ichip]->SetLineWidth(2);
		h1d_clus_per_evt[ichip]->Draw("same");
	}

	//return;
	//
	//
	ofstream fpre;
	if ( bPRE_ALIGN ){
		char fname[300];
		sprintf(fname,"beamcenter/beamcenter_%08d.txt",runnumber);
		fpre.open(fname);
	}

	TCanvas *c1 = new TCanvas("c1","c1",150*1.1*3,150*4);
	c1->Divide(3,4);

	for (int ichip=0; ichip<4; ichip++){

		c1->cd(3*ichip+1);
		SetPadStyle();
		gPad->SetLogz();

		htmp = (TH1F*)h2d_hit[ichip];
		SetHistoStyle();

		htmp->Draw("col2");

		c1->cd(3*ichip+2);
		SetPadStyle();

		htmp = (TH1F*)h1d_hit_x[ichip];
		SetHistoStyle();
		htmp->Draw();

		TF1 *fx = new TF1("fx","gaus",z_min,z_max);
		htmp->Fit(fx,"R0Q");
		//fx->Draw("same");

		fx->SetRange(fx->GetParameter(1)-fx->GetParameter(2), fx->GetParameter(1)+fx->GetParameter(2));
		htmp->Fit(fx,"R0Q");
		fx->Draw("same");

		/*
		TF1 *fxx = new TF1("fxx","gaus(0)+gaus(3)",0,1024);
		fxx->SetParameter(0,fx->GetParameter(0));
		fxx->SetParameter(1,fx->GetParameter(1));
		fxx->SetParameter(2,fx->GetParameter(2));
		fxx->SetParameter(3,0.1*fx->GetParameter(0));
		fxx->FixParameter(4,fx->GetParameter(1));
		fxx->SetParameter(5,10.0*fx->GetParameter(2));
		htmp->Fit(fxx,"R0");
		fxx->Draw("same");
		*/

		c1->cd(3*ichip+3);
		SetPadStyle();

		htmp = (TH1F*)h1d_hit_y[ichip];
		SetHistoStyle();
		htmp->Draw();

		TF1 *fy = new TF1("fy","gaus",x_min,x_max);
		htmp->Fit(fy,"R0Q");
		//fy->Draw("same");
		//
		fy->SetRange(fy->GetParameter(1)-fy->GetParameter(2), fy->GetParameter(1)+fy->GetParameter(2));
		htmp->Fit(fy,"R0Q");
		fy->Draw("same");

		/*
		TF1 *fyy = new TF1("fyy","gaus(0)+gaus(3)",0,1024);
		fyy->SetParameter(0,fy->GetParameter(0));
		fyy->SetParameter(1,fy->GetParameter(1));
		fyy->SetParameter(2,fy->GetParameter(2));
		fyy->SetParameter(3,0.1*fy->GetParameter(0));
		fyy->SetParameter(4,fy->GetParameter(1));
		fyy->SetParameter(5,10.0*fy->GetParameter(2));
		htmp->Fit(fyy,"R0Q");
		fyy->Draw("same");
		*/

		double offset_xx = fx->GetParameter(1);
		double offset_yy = fy->GetParameter(1);
		/*
		if ( fabs(fxx->GetParameter(0))<fabs(fxx->GetParameter(3)) ){
			offset_xx = fxx->GetParameter(4);
		}

		if ( fabs(fyy->GetParameter(0))<fabs(fyy->GetParameter(3)) ){
			offset_yy = fyy->GetParameter(4);
		}
		*/

		if ( bPRE_ALIGN ){
			fpre << ichip << " 0 " << ichip << " " << offset_yy << " " << -2.2*ichip/(28e-4) << " " << offset_xx << endl;
			cout << ichip << " 0 " << ichip << " " << offset_yy << " " << -2.2*ichip/(28e-4) << " " << offset_xx << endl;
		}


		/*
		cout 
			<< "CHIP: " << ichip 
			<< ", x: " << fx->GetParameter(1) 
			<< ", sx: " << fx->GetParameter(2) 
			<< ", y: " << fy->GetParameter(1) 
			<< ", sy: " << fy->GetParameter(2) 
			<< endl;
		*/

	}//

	if ( bPRE_ALIGN ){
		fpre.close();
	}

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

	TCanvas *c4 = new TCanvas("c4","c4",300*1.1*2,300*2);
	c4->Divide(2,2);

	c4->cd(1);
	SetPadStyle();
	gPad->SetLogy();

	htmp = (TH1F*)gPad->DrawFrame(-6,1,6,1.5*h1d_clus_res_z[0]->GetMaximum());
	SetHistoStyle("dz (pixel, 28 #mum)");

	for (int ichip=0; ichip<4; ichip++){
		h1d_clus_res_z[ichip]->SetLineWidth(2);
		h1d_clus_res_z[ichip]->Draw("same");
	}

	c4->cd(2);
	SetPadStyle();
	//gPad->SetLogy();

	htmp = (TH1F*)gPad->DrawFrame(0,-2.5,1024,+2.5);
	SetHistoStyle();
	SetHistoStyle("Column index","dz (pixel, 28 #mum)");

	for (int ichip=0; ichip<4; ichip++){
		h1p_clus_res_z[ichip]->Draw("same");

		f1p_clus_res_z[ichip] = new TF1(Form("f1p_clus_res_z_%d",ichip),"pol1",z_min,z_max);
		h1p_clus_res_z[ichip]->Fit(f1p_clus_res_z[ichip],"R0Q");
		f1p_clus_res_z[ichip]->SetLineColor(nColor[ichip]);
		if ( bPOST_ALIGN )
			f1p_clus_res_z[ichip]->Draw("same");
	}

	c4->cd(3);
	SetPadStyle();
	gPad->SetLogy();

	htmp = (TH1F*)gPad->DrawFrame(-6,1,6,1.5*h1d_clus_res_x[0]->GetMaximum());
	SetHistoStyle("dx (pixel, 28 #mum)");

	for (int ichip=0; ichip<4; ichip++){
		h1d_clus_res_x[ichip]->SetLineWidth(2);
		h1d_clus_res_x[ichip]->Draw("same");
	}

	c4->cd(4);
	SetPadStyle();
	//gPad->SetLogy();

	htmp = (TH1F*)gPad->DrawFrame(0,-2.5,512,+2.5);
	SetHistoStyle("Row index","dx (pixel, 28 #mum)");

	for (int ichip=0; ichip<4; ichip++){
		h1p_clus_res_x[ichip]->Draw("same");

		f1p_clus_res_x[ichip] = new TF1(Form("f1p_clus_res_x_%d",ichip),"pol1",x_min,x_max);
		h1p_clus_res_x[ichip]->Fit(f1p_clus_res_x[ichip],"R0Q");
		f1p_clus_res_x[ichip]->SetLineColor(nColor[ichip]);
		if ( bPOST_ALIGN )
			f1p_clus_res_x[ichip]->Draw("same");
	}

	if ( bPOST_ALIGN ){
		char fname[300];
		sprintf(fname,"beamcenter/beamcenter_%08d.txt",runnumber);
		ifstream fpar;
		fpar.open(fname);
		int tmp_layer, tmp_stave, tmp_chip;
		double tmp_xx, tmp_yy, tmp_zz;

		double par_xx[4], par_yy[4], par_zz[4];

		while ( fpar >> tmp_layer >> tmp_stave >> tmp_chip >> tmp_xx >> tmp_yy >> tmp_zz ){
			par_xx[tmp_chip] = tmp_xx;
			par_yy[tmp_chip] = tmp_yy;
			par_zz[tmp_chip] = tmp_zz;
		}
		fpar.close();

		ofstream fpar_out;
		sprintf(fname,"beamcenter/beamcenter_%08d.txt",runnumber);
		fpar_out.open(fname);

		for (int ichip=0; ichip<4; ichip++){
			cout 
				<< ichip << " " 
				<< "0 " 
				<< ichip << " " 
				<< par_xx[ichip] + (f1p_clus_res_x[ichip]->Eval(par_xx[0]) - f1p_clus_res_x[0]->Eval(par_xx[0])) << " "
				<< par_yy[ichip] << " "
				<< par_zz[ichip] + (f1p_clus_res_z[ichip]->Eval(par_zz[0]) - f1p_clus_res_z[0]->Eval(par_zz[0]))
				<< endl;

			fpar_out
				<< ichip << " " 
				<< "0 " 
				<< ichip << " " 
				<< par_xx[ichip] + (f1p_clus_res_x[ichip]->Eval(par_xx[0]) - f1p_clus_res_x[0]->Eval(par_xx[0])) << " "
				<< par_yy[ichip] << " "
				<< par_zz[ichip] + (f1p_clus_res_z[ichip]->Eval(par_zz[0]) - f1p_clus_res_z[0]->Eval(par_zz[0]))
				<< endl;
		}
		fpar_out.close();
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

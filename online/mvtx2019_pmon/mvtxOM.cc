// Updated by Xiaochun He on May 31, 2019 following Martin Purschke's
// suggestion correction
//
#include "mvtxOM.h"

#include <mvtx/MvtxDefs.h>
#include <pmonitor/pmonitor.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>

#include <fstream>
#include <iostream>
#include <map>

#define IDMVTXV1_MAXRUID       4
#define IDMVTXV1_MAXRUCHN      28

int init_done = 0;

using namespace std;

const int NSTAVE = 4;
const bool chip_expected[4] = {true, true, true, true};
string stave_name[4] = {"E103", "C105", "C104", "A105"};

vector<TLine*> chip_edges, dead_chip_forward, dead_chip_backward;

string outHitLocations = "/home/maps/meeg/felix/daq/felix_rcdaq/online_monitoring/hitLocations/locations.txt";
ofstream write_outHitLocations(outHitLocations.c_str());

map<pair<int,int>,pair<int,int>> chipmap = {
  {{1,1}, {0,0}},
  {{1,2}, {0,1}},
  {{1,3}, {0,2}},
  {{1,4}, {0,3}},
  {{1,5}, {0,4}},
  {{1,6}, {0,5}},
  {{1,7}, {0,6}},
  {{1,8}, {0,7}},
  {{1,9}, {3,8}},
  {{1,10}, {3,7}},
  {{1,11}, {3,6}},
  {{1,12}, {3,5}},
  {{1,13}, {3,4}},
  {{1,14}, {3,3}},
  {{1,15}, {3,2}},
  {{1,16}, {3,1}},
  {{1,17}, {3,0}},
  {{1,18}, {2,8}},
  {{1,19}, {2,7}},
  {{1,20}, {2,6}},
  {{1,21}, {2,5}},
  {{1,22}, {2,4}},
  {{1,23}, {2,3}},
  {{1,24}, {2,2}},
  {{1,25}, {2,1}},
  {{1,26}, {2,0}},
  {{1,27}, {0,8}},
  {{2,1}, {1,2}},
  {{2,2}, {1,1}},
  {{2,3}, {1,0}},
  {{2,4}, {1,3}},
  {{2,5}, {1,4}},
  {{2,6}, {1,5}},
  {{2,7}, {1,6}},
  {{2,8}, {1,7}},
  {{2,9}, {1,8}}
}; //<ruid, ruchn> to <stave, chipID>

int mvtx_evnts;
int mvtx_verbose = 0;
int mvtx_refresh = -1;
int max_npixels = 512*1024; // cut to suppress plotting events with lots of hits
int mvtx_run = 0;
const bool flip_yaxis = false;

//-- histograms filled in event loop
TH1F* hnevnt; // number of events
TH1F* hchip;  // hits vs chip
TH1F* hwarn;
TH1F* herr;
TH1F* hnhit_chip[NSTAVE]; // nhit distribution for each chip
TH2F* h2d_chip[NSTAVE]; // 2D pixel hits for each chip
TH1F* hdiffrow_chip[NSTAVE];
TH1F* hdiffcol_chip[NSTAVE];

//-- variables for online monitoring plotting
TCanvas* com; // canvas
TPad* phitrate;
TPad* pnhit;
TPad* p2d[NSTAVE];
TPad* pmean;
TPad* pdiffcol;
TPad* pdiffrow;
TH1F* haxis_chip; // axis for hits/event vs chip
TH1F* haxis_nhit; // axis for hits/event distribution for each chip
TH1F* haxis_2d;   // axis for 2D pixel hits
TH1F* haxis_diff; // axis for diff pixel index (row or col)
TLatex* lnevents;
TLatex* ldiffcol[NSTAVE];
TLatex* ldiffrow[NSTAVE];
TLatex* lnhitmean[NSTAVE];
TGaxis* reversedaxis;

//-- analysis histograms
TH1F* hhitrate_chip;
TH1F* hwarnrate_chip;
TH1F* herrrate_chip;
TH1F* hnhit_chip_norm[NSTAVE];
TH2F* h2d_chip_norm[NSTAVE];
TGraphAsymmErrors* gmpos_chip[NSTAVE];
TH1F* hdiffrow_chip_norm[NSTAVE];
TH1F* hdiffcol_chip_norm[NSTAVE];
TH1F* hhittime_chip[NSTAVE];

// The following line caused some problem of running ROOT6 in macro
//TF1* fg = new TF1("fg", "gaus", 0, 1024);

TF1* fg;

//-- some constants for different chips
int chipColor[] = {kBlue, kRed, kGreen+2, kMagenta+2};
int chipMarker[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross};

// Show fit to beam center
TCanvas* cBeamCenter = nullptr;
const bool show_beam_fit = true;

unsigned short decode_row(int hit)
{
    return hit >> 16;
}

unsigned short decode_col(int hit)
{
    return hit & 0xffff;
}


//============================================================//

void set_verbose(int v)
{
    mvtx_verbose = v;
}

//============================================================//

void set_refresh(int r)
{
    mvtx_refresh = r;
}

//============================================================//

int pinit()
{

  // Added by Xiaochun He following Martin's recommendation
  fg = 0;

  if (init_done) return 1;
  init_done = 1;

  // reset event counter
  mvtx_evnts = 0;

  hnevnt = new TH1F("hnevent", "", 1, -0.5, 0.5);

  // hits vs chip index
  hchip = new TH1F("hchip", ";chip index;N pixels", 4, -0.5, 3.5);
  hchip->SetLineWidth(2);
  hchip->SetLineColor(kBlue);

  // warnings vs chip index
  hwarn = new TH1F("hwarn", ";chip index; warnings", 4, -0.5, 3.5);
  hwarn->SetLineWidth(2);
  hwarn->SetLineColor(kYellow+2);

  // errors vs chip index
  herr = new TH1F("herr", ";chip index; errors", 4, -0.5, 3.5);
  herr->SetLineWidth(2);
  herr->SetLineColor(kRed);

  char name[500];
  for (int i = 0; i < NSTAVE; i++)
    {
      // total hits/event distribution in each chip
      sprintf(name, "hnhit_chip_%i", i);
      hnhit_chip[i] = new TH1F(name, ";N pixels / event; events", 201, -0.5, 200.5);
      hnhit_chip[i]->SetLineWidth(2);
      hnhit_chip[i]->SetLineColor(chipColor[i]);

      // time distribution
      sprintf(name, "hhittime_chip_%i", i);
      hhittime_chip[i] = new TH1F(name, ";N pixels / 100 events; event number", 100, 0, 10000);
      hhittime_chip[i]->SetLineWidth(2);
      hhittime_chip[i]->SetLineColor(chipColor[i]);

      // 2D pixel hit dsitribution for each chip
      sprintf(name, "h2d_chip_%i", i);
      h2d_chip[i] = new TH2F(name, ";col;row;hits", 9216, -0.5, 9215.5, 512, -0.5, 511.5);

      // mean pixel index for each chip (aggregated over events)
      sprintf(name, "gmpos_chip_%i", i);
      gmpos_chip[i] = new TGraphAsymmErrors();
      gmpos_chip[i]->SetName(name);
      gmpos_chip[i]->SetMarkerStyle(chipMarker[i]);
      gmpos_chip[i]->SetMarkerColor(chipColor[i]);
      gmpos_chip[i]->SetLineColor(chipColor[i]);

      // difference in mean pixel index per event (from chip 0)
      sprintf(name, "hdiffrow_chip_%i", i);
      hdiffrow_chip[i] = new TH1F(name, ";row index diff", 1023, -511.5, 511.5);
      hdiffrow_chip[i]->SetLineWidth(2);
      hdiffrow_chip[i]->SetLineColor(chipColor[i]);

      sprintf(name, "hdiffcol_chip_%i", i);
      hdiffcol_chip[i] = new TH1F(name, ";col index diff", 1023, -511.5, 511.5);
      hdiffcol_chip[i]->SetLineWidth(2);
      hdiffcol_chip[i]->SetLineColor(chipColor[i]);
    }


  // --------------------------------------------
  // -- Setup canvas for Online Monitoring plots
  // --------------------------------------------

  // -- setup axis
  haxis_chip = new TH1F("haxis_chip",
			";chip index; N pixels / event",
			NSTAVE, -0.5, NSTAVE-0.5);
  haxis_chip->SetMinimum(0);
  haxis_chip->SetMaximum(500);

  haxis_nhit = new TH1F("haxis_nhit",
			";N pixels / event; per event",
			51, -0.5, 50.5);
  haxis_nhit->SetMinimum(0);
  haxis_nhit->SetMaximum(1.0);

  haxis_2d = new TH1F("haxis_2d",
		      ";row;col",
		      512, -0.5, 511.5);
  haxis_2d->SetMinimum(-0.5);
  haxis_2d->SetMaximum(511.5);
  if (flip_yaxis)
    {
      haxis_2d->GetYaxis()->SetLabelOffset(999);
      haxis_2d->GetYaxis()->SetTickLength(0);
    }

  haxis_diff = new TH1F("haxis_diff",
			";mean idx diff",
			1023, -511.5, 511.5);
  haxis_diff->SetMinimum(0);
  haxis_diff->SetMaximum(1);

  return 0;

}

//============================================================//

int process_event (Event * e)
{

  // Added by Xiaochu He following Martin's recommedation
  if (!fg) fg = new TF1("fg", "gaus", 0, 1024);

  if ( e->getEvtType() == BEGRUNEVENT )
    {
      mvtx_run = e->getRunNumber();
      mvtx_evnts = 0;
      reset_histos();
      OM();
      return 0;
    }
  if ( e->getEvtType() != DATAEVENT )
    return 0;


  Packet *p = e->getPacket(2000);
  if (p)
    {

      bool evnt_err = false;


      /*
	int bad_ruchns = p->iValue(0, "BAD_RUCHNS");
	if ( bad_ruchns > 0 )
	{
	if ( mvtx_verbose > 0 )
	cout << "WARNING!! Event: " << mvtx_evnts << " Invalid RU channel IDs (really bad data)!"
	<< " BAD_RUCHNS:" << bad_ruchns << endl;
	}

	int bad_chipids = p->iValue(0, "BAD_CHIPIDS");
	if ( bad_chipids > 0 )
	{
	if ( mvtx_verbose > 0 )
	cout << "WARNING!! Event: " << mvtx_evnts << " Invalid chip IDs (bad data)!"
	<< " BAD_CHIPIDS:" << bad_chipids << endl;
	}

	int chipmax = p->iValue(0, "HIGHEST_CHIP") + 1;
	if ( chipmax > NSTAVE )
	{
	if ( mvtx_verbose > 1 )
	cout << "WARNING!! Event: " << mvtx_evnts << " More chips than expected!"
	<< " NSTAVE:" << NSTAVE << " HIGHEST_CHIP:" << chipmax << endl;

	chipmax = NSTAVE;
	}
	if ( mvtx_verbose > 2 )
	{
	cout << "Event:" << mvtx_evnts << " chipmax:" << chipmax << endl;
	}
	float mrow_chip0 = -1;
	float mcol_chip0 = -1;

	int excess_data_bytes = p->iValue(0, "EXCESS_DATA_BYTES");
	if ( excess_data_bytes>0 )
	{
	if ( mvtx_verbose > 1 )
	cout << "WARNING!! Event: " << mvtx_evnts << " Data found past chip trailer"
	<< " EXCESS_DATA_BYTES: " << excess_data_bytes << endl;
	}
	bool evnt_err = false;
	for ( int ichip = 0; ichip < NSTAVE; ichip++)
	{
	int header_found = p->iValue(ichip, "HEADER_FOUND");
	int trailer_found = p->iValue(ichip, "TRAILER_FOUND");
	int bunchcounter = p->iValue(ichip, "BUNCHCOUNTER");
	int unexpected_bytes = p->iValue(ichip, "UNEXPECTED_BYTES");
	int readout_flags = p->iValue(ichip, "READOUT_FLAGS");
        //cout << "HEADER_FOUND: " << header_found << " TRAILER_FOUND: " << trailer_found << " BUNCHCOUNTER: " << bunchcounter << endl;

        bool has_warning = false;
        bool has_error = false;
        if (chip_expected[ichip])
        {
        if (header_found==0 || trailer_found==0)
        {
        if ( mvtx_verbose > 1 )
        cout << "WARNING!! Event: " << mvtx_evnts << " Missing chip " << ichip << " HEADER_FOUND: " << header_found << " TRAILER_FOUND: " << trailer_found << " BUNCHCOUNTER: " << bunchcounter << endl;
        has_warning = true;
        }
        if ( (header_found + trailer_found) == 1 )
        {
        if ( mvtx_verbose > 1 )
        cout << "ERROR!! Event: " << mvtx_evnts << " header and trailer have different states for chip " << ichip << " HEADER_FOUND: " << header_found << " TRAILER_FOUND: " << trailer_found << " BUNCHCOUNTER: " << bunchcounter << endl;
        has_error = true;
        evnt_err = true;
        }
        }
        else
        {
        if (header_found!=0 || trailer_found!=0)
        {
            if ( mvtx_verbose > 1 )
                cout << "WARNING!! Event: " << mvtx_evnts << " Unexpected chip " << ichip << " HEADER_FOUND: " << header_found << " TRAILER_FOUND: " << trailer_found << " BUNCHCOUNTER: " << bunchcounter << endl;
            has_warning = true;
        }
    }
    if (unexpected_bytes!=0)
    {
        if ( mvtx_verbose > 0 )
            cout << "WARNING!! Event: " << mvtx_evnts << " chip " << ichip << " UNEXPECTED_BYTES: " << unexpected_bytes << endl;
        has_warning = true;
    }
    if (readout_flags > 0)
    {
        if ( mvtx_verbose > 1 )
            cout << "WARNING!! Event: " << mvtx_evnts << " chip " << ichip << " READOUT_FLAGS: " << hex << readout_flags << dec << endl;
        has_warning = true;
    }

    if ( has_warning )
        hwarn->Fill(ichip);
    if ( has_error )
        herr->Fill(ichip);
    }
    */

    int npixels[NSTAVE] = {0};
    double mrow[NSTAVE] = {0};
    double mcol[NSTAVE] = {0};
    double mrow_refstave = -1;
    double mcol_refstave = -1;

    if ( !evnt_err ) {
        for (int ruid=0; ruid<IDMVTXV1_MAXRUID+1; ruid++)
        {
            if (p->iValue(ruid)!=-1)
            {
                for ( int ruchn = 0; ruchn < IDMVTXV1_MAXRUCHN+1; ruchn++)
                {
                    if (p->iValue(ruid,ruchn)>0)
                    {
                        for (int i=0;i<p->iValue(ruid,ruchn);i++)
                        {
                            int hit = p->iValue(ruid,ruchn,i);
                            int irow = decode_row(hit);
                            int icol = decode_col(hit);
                            //cout << "(ruid " << ruid << ", ruchn " << ruchn << ") ";
                            //cout << "(row " << irow << ", col " << icol << ") ";
                            if (chipmap.count({ruid,ruchn}) != 1) {
                                cout << "invalid: (ruid " << ruid << ", ruchn " << ruchn << ") " << endl;
                            } else {
                                pair<int, int> chiplocation = chipmap[{ruid,ruchn}];
                                int istave = chiplocation.first;
                                int ichip = chiplocation.second;
                                //cout << "(stave " << istave << ", chip " << ichip << ") ";
                                npixels[istave]++;
                                mrow[istave]+=irow;
                                mcol[istave]+=icol+1024*ichip;
                                if (flip_yaxis)
                                    h2d_chip[istave]->Fill(icol+1024*ichip,irow);
                                else
                                {
                                    h2d_chip[istave]->Fill(icol+1024*ichip,511-irow);
                                    if(ichip == 4){
                                        write_outHitLocations<<e->getEvtSequence()<<", "<< istave<<", "<< icol<<", "<<irow<<endl;
                                        //printf("There is a hit on stave %i at x = %i, y = %i for event %i\n", istave, icol, irow, e->getEvtSequence());
                                        }
                                 }
                            }
                        }
                        //cout << endl;
                    }
                }
            }
        }
        for (int istave=0;istave<NSTAVE;istave++) {
            if (npixels[istave] > 0  && npixels[istave] < max_npixels)
            {
                mrow[istave] /= (float)npixels[istave];
                mcol[istave] /= (float)npixels[istave];

                hchip->Fill(istave, npixels[istave]);
                hhittime_chip[istave]->Fill(mvtx_evnts,npixels[istave]);

                if ( mvtx_verbose > 2 )
                {
                    cout << "    chip:" << istave << " npixels:" << npixels[istave]
                        << " mean:(" << mcol[istave] << ", " << mrow[istave] << ")" << endl;
                }
            }
        }
        if ( npixels[NSTAVE-1] > 0 && npixels[NSTAVE-1] < max_npixels )
        {
            mrow_refstave = mrow[NSTAVE-1];
            mcol_refstave = mcol[NSTAVE-1];
            for (int istave=0;istave<NSTAVE;istave++) {
                if (npixels[istave] > 0  && npixels[istave] < max_npixels)
                {
                    if ( mrow_refstave >= 0 && mcol_refstave >= 0 )
                    {
                        hdiffrow_chip[istave]->Fill(mrow[istave] - mrow_refstave);
                        hdiffcol_chip[istave]->Fill(mcol[istave] - mcol_refstave);
                    }

                }
            }
        }
    }

    /*

    // fill regardless of the number of hits
    hnhit_chip[ichip]->Fill(npixels);

    // only fill for nonzero hits
    if (npixels > 0  && npixels < max_npixels)
    {
    mrow /= (float)npixels;
    mcol /= (float)npixels;

    if ( ichip == NSTAVE - 1 )
    {
    mrow_chip0 = mrow;
    mcol_chip0 = mcol;
    }

    if ( mrow_chip0 >= 0 && mcol_chip0 >= 0 )
    {
    hdiffrow_chip[ichip]->Fill(mrow - mrow_chip0);
    hdiffcol_chip[ichip]->Fill(mcol - mcol_chip0);
    }

    hchip->Fill(ichip, npixels);
    hhittime_chip[ichip]->Fill(mvtx_evnts,npixels);

    if ( mvtx_verbose > 2 )
    {
    cout << "    chip:" << ichip << " npixels:" << npixels
    << " mean:(" << mcol << ", " << mrow << ")" << endl;
    }
    }
    } // ichip
    */
    hnevnt->Fill(0);
    delete p;
}

if ( mvtx_refresh > 0 &&  mvtx_evnts%mvtx_refresh == 0 ) OM();

if ( mvtx_evnts%100000 == 0 )
    cout << " processing event " << mvtx_evnts << endl;

    mvtx_evnts++;
    return 0;
    }

//============================================================//

int process_histos(float thresh)
{
    int sum = 0;
    int row,col;
    for (int i = 0; i < NSTAVE; i++)
    {
        cout << endl;
        cout << "-- stave " << i << endl;
        int tot = 0;
        for (int ix = 1; ix <= h2d_chip[i]->GetNbinsX(); ix++)
        {
            for (int iy = 1; iy <= h2d_chip[i]->GetNbinsY(); iy++)
            {
                if ( h2d_chip[i]->GetBinContent(ix, iy) > thresh )
                {
                    tot++;
                    if (flip_yaxis)
                    {
                        row = 1023-(h2d_chip[i]->GetYaxis()->GetBinCenter(iy));
                    } else
                    {
                        row = h2d_chip[i]->GetYaxis()->GetBinCenter(iy);
                    }
                    col = h2d_chip[i]->GetXaxis()->GetBinCenter(ix);
                    cout << " row:" << row
                        << " col:" << col
                        << endl;
                }
            }
        }
        cout << "Total: " << tot << endl;
        sum += tot;
    }

    return sum;
}

//============================================================//

int mask_pixels(float thresh)
{
    bool flip_yaxis = true; //Only uncomment if we need to flip the yAxis

    string file_name_tb1 = "/home/maps/git/RUv1_Test_sync2018-08/software/py/masklist_testbench1.txt";
    string file_name_tb2 = "/home/maps/git/RUv1_Test_sync2018-08/software/py/masklist_testbench2.txt";
    ofstream write_mask_file_tb1(file_name_tb1.c_str());
    ofstream write_mask_file_tb2(file_name_tb2.c_str());

    write_mask_file_tb1<<"#Connector, ChipID, Col, Row"<<endl;
    write_mask_file_tb2<<"#Connector, ChipID, Col, Row"<<endl;

    int sum = 0;
    int chipid, row, col, tot, chip_tot, prev_tot;
    Int_t stave_map[4] = {0, 4, 1, 2};
    for (int i = 0; i < NSTAVE; i++)
    {
        tot = 0;
        prev_tot = 0;
        for (int ix = 1; ix <= h2d_chip[i]->GetNbinsX(); ix++)
        {
            for (int iy = 1; iy <= h2d_chip[i]->GetNbinsY(); iy++)
            {
                col = h2d_chip[i]->GetXaxis()->GetBinCenter(ix);
                chipid = col/1024;
                col = col%1024;

                if ( h2d_chip[i]->GetBinContent(ix, iy) > thresh )
                {
                    tot++;
                    if (flip_yaxis)
                    {
                        row = 511-(h2d_chip[i]->GetYaxis()->GetBinCenter(iy));
                    }
                    else
                    {
                        row = h2d_chip[i]->GetYaxis()->GetBinCenter(iy);
                    }
                    if (i != 1) { write_mask_file_tb1 << stave_map[i] << "," << chipid  << "," << col << "," << row << endl; }
                    else { write_mask_file_tb2 << stave_map[i] << "," << chipid  << "," << col << "," << row << endl; }
                }
            }
            chip_tot = tot - prev_tot;
            if (col == 1023) {prev_tot = tot; printf("Total pixels masked on stave %i, chip %i: %i\n", i, chipid, chip_tot);}
        }
        sum += tot;
    }
    printf("Total pixels masked: %i\n", sum);
    write_mask_file_tb1.close();
    write_mask_file_tb2.close();
    printf("New pixel mask has been created\n");

    return sum;
}

//============================================================//

int analysis()
{
    //TH1D* hhitrate_chip;
    //TH1D* hnhit_chip_norm[NSTAVE];
    //TH2D* h2d_chip_norm[NSTAVE];
    //TH1F* hdiffrow_chip_norm[NSTAVE];
    //TH1F* hdiffcol_chip_norm[NSTAVE];
    char name[500];

    //-- Check if we've already initialized the histograms
    //   if not, initialize
    if ( !hhitrate_chip )
    {
        hhitrate_chip = (TH1F*) hchip->Clone("hhitrate_chip");
        hwarnrate_chip = (TH1F*) hwarn->Clone("hwarnrate_chip");
        herrrate_chip = (TH1F*) herr->Clone("herrrate_chip");

        for (int ichip = 0; ichip < NSTAVE; ichip++)
        {
            sprintf(name, "hnhit_chip_norm_%i", ichip);
            hnhit_chip_norm[ichip] = (TH1F*) hnhit_chip[ichip]->Clone(name);

            sprintf(name, "h2d_chip_norm_%i", ichip);
            h2d_chip_norm[ichip] = (TH2F*) h2d_chip[ichip]->Clone(name);

            sprintf(name, "hdiffrow_chip_norm_%i", ichip);
            hdiffrow_chip_norm[ichip] = (TH1F*) hdiffrow_chip[ichip]->Clone(name);

            sprintf(name, "hdiffcol_chip_norm_%i", ichip);
            hdiffcol_chip_norm[ichip] = (TH1F*) hdiffcol_chip[ichip]->Clone(name);
        } // ichip
    }

    //-- current number of events
    double nevents = hnevnt->Integral();

    if ( nevents <= 0 )
        return 1;


    //-- hitrate vs chip index
    hhitrate_chip->Reset();
    hwarnrate_chip->Reset();
    herrrate_chip->Reset();
    for (int ib = 1; ib <= hchip->GetNbinsX(); ib++)
    {
        // hits
        float bc = hchip->GetBinContent(ib);
        hhitrate_chip->SetBinContent(ib, bc / nevents);
        hhitrate_chip->SetBinError(ib, sqrt(bc) / nevents);

        // warnings
        bc = hwarn->GetBinContent(ib);
        hwarnrate_chip->SetBinContent(ib, bc / nevents);
        hwarnrate_chip->SetBinError(ib, sqrt(bc) / nevents);

        // errors
        bc = herr->GetBinContent(ib);
        herrrate_chip->SetBinContent(ib, bc / nevents);
        herrrate_chip->SetBinError(ib, sqrt(bc) / nevents);
    } // ib

    if (show_beam_fit) { // create canvas for beamcenter plots
      cBeamCenter = new TCanvas("cBeamCenter", "cBeamCenter", 1350, 900);
      cBeamCenter->Divide(3,4);
    }

    //-- per chip histograms
    for (int ichip = 0; ichip < NSTAVE; ichip++)
    {
        //-- nhit distribution for each chip
        hnhit_chip_norm[ichip]->Reset();
        for (int ib = 1; ib <= hnhit_chip[ichip]->GetNbinsX(); ib++)
        {
            float bc = hnhit_chip[ichip]->GetBinContent(ib);
            hnhit_chip_norm[ichip]->SetBinContent(ib, bc / nevents);
            hnhit_chip_norm[ichip]->SetBinError(ib, sqrt(bc) / nevents);
        } // ib

        //-- 2D distribution for each chip
        h2d_chip_norm[ichip]->Reset();
        for (int ix = 1; ix <= h2d_chip[ichip]->GetNbinsX(); ix++)
            for (int iy = 1; iy <= h2d_chip[ichip]->GetNbinsY(); iy++)
            {
                float bc = h2d_chip[ichip]->GetBinContent(ix, iy);
                h2d_chip_norm[ichip]->SetBinContent(ix, iy, bc / nevents);
            }

        //-- mean pixel index for each chip
        double m, s;
        TH1D* h;

        if ( cBeamCenter )
        {
          cBeamCenter->cd(3*ichip+1);
          gPad->SetLogz();
          TH2D* h2d_beam = (TH2D*)h2d_chip_norm[ichip]->Clone("h2d_beam");
          h2d_beam->GetXaxis()->SetTitle("Stave_Cols");
          h2d_beam->GetYaxis()->SetTitle("Stave_Rows");
          h2d_beam->Draw("colz");
          TLatex lt;
          lt.SetNDC();
          lt.SetTextAlign(22);
          lt.SetTextSize(1.5*lt.GetTextSize());
          lt.SetTextColor(1);
          lt.DrawText(0.5,0.96, stave_name[ichip].c_str());

          switch ( ichip )
          {
            case 0 :
              dead_chip_forward[0]->Draw();
              dead_chip_backward[0]->Draw();
              break;

            case 2 :
              dead_chip_forward[0]->Draw(); dead_chip_backward[0]->Draw();
              dead_chip_forward[3]->Draw(); dead_chip_backward[3]->Draw();
              dead_chip_forward[6]->Draw(); dead_chip_backward[6]->Draw();
              dead_chip_forward[7]->Draw(); dead_chip_backward[7]->Draw();
              dead_chip_forward[8]->Draw(); dead_chip_backward[8]->Draw();
              break;

            case 3 :
              dead_chip_forward[6]->Draw(); dead_chip_backward[6]->Draw();
              break;

            default :
             break;
          }
        }

        // cols projection
        h = (TH1D*) h2d_chip_norm[ichip]->ProjectionX();
        double m_col = h->GetMean();
        double r_col = h->GetRMS();
        fg->SetParameters(h->GetBinContent(h->GetMaximumBin()), m_col, r_col);
        //YCM: initial range to full stave range
        const int lastBinX = 9 * 1024 - 1;
        const int lastBinY = 511;
        fg->SetRange(0, lastBinX);
        h->Fit(fg, "RQ0N");
        m = fg->GetParameter(1);
        s = fg->GetParameter(2);
        //cout << " x1 mean:" << m << " sig:" << s << endl;
        m = m - 2.5*s > lastBinX ? lastBinX : m;
        fg->SetRange(m - 2.5*s, m + 2.5*s);
        fg->SetParameters(fg->GetParameter(0), m, s);
        h->Fit(fg, "RQ0N");
        m = fg->GetParameter(1);
        s = fg->GetParameter(2);
        //cout << " x2 mean:" << m << " sig:" << s << endl;
        if ( m > lastBinX )
        {
            s = s - (m - lastBinX);
            m = lastBinX;
        }
        if ( m < 0 )
        {
            s = s - m;
            m = 0;
        }
        //cout << " (m:" << m << " s:" << s << ")"<< endl;

        // save
        m_col = m;
        r_col = s;
        if ( cBeamCenter )
        {
          cBeamCenter->cd( 3 * ichip + 2 );
          h->GetXaxis()->SetTitle("Stave_Cols");
          h->GetYaxis()->SetTitle("nClusters");
          h->Draw();
          TF1* f_col = (TF1*)fg->Clone("fx");
          f_col->Draw("same");
          TLatex lt;
          lt.SetNDC();
          lt.SetTextAlign(22);
          lt.SetTextSize(1.5*lt.GetTextSize());
          lt.SetTextColor(2);
          lt.DrawLatex(.3, .7, Form("#mu_{COL} = %.2f #pm %.2f ", m_col, r_col));
          lt.DrawLatex(.3, .6, Form("offset = %.2f", m_col));
        }
        else
          delete h;

        // row projection
        h = (TH1D*) h2d_chip_norm[ichip]->ProjectionY();
        double m_row = h->GetMean();
        double r_row = h->GetRMS();
        fg->SetParameters(h->GetBinContent(h->GetMaximumBin()), m_row, r_row);
        fg->SetRange(0, 1023);
        h->Fit(fg, "RQ0N");
        m = fg->GetParameter(1);
        s = fg->GetParameter(2);
        //cout << " y2 mean:" << m << " sig:" << s << endl;
        fg->SetRange(m - 2.5*s, m + 2.5*s);
        fg->SetParameters(fg->GetParameter(0), m, s);
        h->Fit(fg, "RQ0N");
        m = fg->GetParameter(1);
        s = fg->GetParameter(2);
        //cout << " y2 mean:" << m << " sig:" << s << endl;
        if ( m > lastBinY )
        {
            s = s - (m - lastBinY);
            m = lastBinY;
        }
        if ( m < 0 )
        {
            s = s - m;
            m = 0;
        }
        //cout << " (m:" << m << " s:" << s << ")"<< endl;

        // save
        m_row = m;
        r_row = s;
        if ( cBeamCenter )
        {
          cBeamCenter->cd( 3 * ichip + 3 );
          h->GetXaxis()->SetTitle("Stave_Cols");
          h->GetYaxis()->SetTitle("nClusters");
          h->Draw();
          TF1* f_row = (TF1*)fg->Clone("fy");
          f_row->Draw("same");
          TLatex lt;
          lt.SetNDC();
          lt.SetTextAlign(22);
          lt.SetTextSize(1.5*lt.GetTextSize());
          lt.SetTextColor(2);
          lt.DrawLatex(.3, .7, Form("#mu_{ROW}  = %.2f #pm %.2f ", m_row, r_row));
          lt.DrawLatex(.3, .6, Form("offset = %.2f", lastBinY - m_row));
        }
        else
          delete h;

        gmpos_chip[ichip]->SetPoint(0, m_col, m_row);
        gmpos_chip[ichip]->SetPointError(0, r_col, r_col, r_row, r_row);

        //-- difference in mean pixel index per event for each chip
        hdiffrow_chip_norm[ichip]->Reset();
        hdiffcol_chip_norm[ichip]->Reset();
        for (int ib = 1; ib <= hdiffrow_chip[ichip]->GetNbinsX(); ib++)
        {
            hdiffrow_chip_norm[ichip]->SetBinContent(ib, hdiffrow_chip[ichip]->GetBinContent(ib));
            hdiffcol_chip_norm[ichip]->SetBinContent(ib, hdiffcol_chip[ichip]->GetBinContent(ib));
        }

        float introw = hdiffrow_chip_norm[ichip]->Integral();
        if ( introw > 0 )
            hdiffrow_chip_norm[ichip]->Scale(1./introw);

        float intcol = hdiffcol_chip_norm[ichip]->Integral();
        if ( intcol > 0 )
            hdiffcol_chip_norm[ichip]->Scale(1./intcol);
    } // ichip

    return 0;
}

//============================================================//

int OM()
{
    //-- run analysis
    analysis();



    gStyle->SetOptStat(0);
    char name[500];
    //-- setup canvas
    static bool initialized = false;

    if ( !initialized )
    {

        //com3 = new TCanvas("com","MVTX online monitoring", 1600, 800);
        com = new TCanvas("com","MVTX online monitoring", 2000, 800);
        com->SetMargin(0, 0, 0, 0);

        // hitrate vs chip index
        //phitrate = new TPad("phitrate", "hit rate", 0.0, 0.45, 0.29, 0.9);
        phitrate = new TPad("phitrate", "hit rate", 0.0, 0.60, 0.29, 0.89);
        phitrate->SetMargin(0.15, 0.05, 0.15, 0.05);
        phitrate->SetTicks(1, 1);
        phitrate->Draw();
/*
        // nhit distributions
        com->cd();
        pnhit = new TPad("pnhit", "n hit", 0.29, 0.60, 0.50, 0.90);
        pnhit->SetMargin(0.12, 0.02, 0.15, 0.05);
        pnhit->SetTicks(1, 1);
        pnhit->Draw();

        // mean hit positions
        com->cd();
        pmean = new TPad("pmean", "mean", 0.00, 0.00, 0.29, 0.45);
        pmean->SetMargin(0.15, 0.05, 0.15, 0.10);
        pmean->SetTicks(1, 1);
        pmean->Draw();
*/
        // diff col index
        com->cd();
        //pdiffcol = new TPad("pdiffcol", "pdiffcol", 0.29, 0.30, 0.50, 0.60);
        pdiffcol = new TPad("pdiffcol", "pdiffcol", 0.00, 0.30, 0.29, 0.59);
        //pdiffcol->SetMargin(0.12, 0.02, 0.12, 0.10);
        pdiffcol->SetMargin(0.15, 0.05, 0.15, 0.05);
        pdiffcol->SetTicks();
        pdiffcol->Draw();

        // diff row index
        com->cd();
        //pdiffrow = new TPad("pdiffrow", "pdiffrow", 0.29, 0.00, 0.50, 0.30);
        pdiffrow = new TPad("pdiffrow", "pdiffrow", 0.00, 0.00, 0.29, 0.29);
        //pdiffrow->SetMargin(0.12, 0.02, 0.12, 0.10);
        pdiffrow->SetMargin(0.15, 0.05, 0.15, 0.05);
        pdiffrow->SetTicks();
        pdiffrow->Draw();

        // 2D hit distributions
        for (int i = 0; i < NSTAVE; i++)
        {
            com->cd();
            sprintf(name, "p2d_%i", i);
            //p2d[i] = new TPad(name, name,
            //        0.50, 0.9-0.22*i,
            //        1.00, 0.9-0.22*(i+1) );
            p2d[i] = new TPad(name, name,
                    0.30, 0.9-0.22*i,
                    1.00, 0.9-0.22*(i+1) );
            //p2d[i]->SetMargin(0.15, 0.15, 0.15, 0.15);
            p2d[i]->SetMargin(0.01, 0.02, 0.25, 0.25);
            p2d[i]->SetTicks(1, 1);
            p2d[i]->Draw();
        } // i

        //-- info
        com->cd();
        sprintf(name, "Number of Events: %.0f", hnevnt->Integral());
        lnevents = new TLatex(0.5, 0.95, name);
        lnevents->SetNDC();
        lnevents->SetTextAlign(22);
        lnevents->Draw();

        TLatex lt;
        lt.SetNDC();
        lt.SetTextAlign(22);

        //-- hitrate
        phitrate->cd();
        //gPad->SetLogy();
        haxis_chip->GetYaxis()->SetRangeUser(0,5);
        haxis_chip->Draw();
        hhitrate_chip->Draw("hist e same");
        hwarnrate_chip->Draw("hist e same");
        herrrate_chip->Draw("hist e same");

        //-- nhit distributions
        //pnhit->cd();
        //gPad->SetLogy();
        haxis_nhit->GetYaxis()->SetRangeUser(1e-4, 1);
        //haxis_nhit->DrawCopy();
        for (int i = 0; i < NSTAVE; i++)
        {
            hnhit_chip_norm[i]->Draw("hist e same");

            sprintf(name, "chip %i = %.2f", i, hnhit_chip_norm[i]->GetMean());
            lnhitmean[i] = new TLatex(0.8, 0.85 - 0.05*i, name);
            lnhitmean[i]->SetNDC();
            lnhitmean[i]->SetTextColor(chipColor[i]);
            lnhitmean[i]->Draw("same");
        }

        //-- mean positions
        //pmean->cd();
        //haxis_2d->DrawCopy();
        //gPad->Update();
/*        if (flip_yaxis)
        {
            reversedaxis = new TGaxis(gPad->GetUxmin(),
                    gPad->GetUymax(),
                    gPad->GetUxmin()-0.001,
                    gPad->GetUymin(),
                    -0.5,
                    1023.5,
                    510, "+");
            reversedaxis->SetLabelOffset(-0.03);
            reversedaxis->Draw();
        }
        for (int i = 0; i < NSTAVE; i++)
            gmpos_chip[i]->Draw("P");
        lt.DrawLatex(0.5, 0.96, "Beam position & profile");
*/
        //-- diff col
        pdiffcol->cd();
        gPad->SetLogy();
        haxis_diff->GetYaxis()->SetRangeUser(1e-4, 1);
        //haxis_diff->GetXaxis()->SetRangeUser(-25, 25);
        haxis_diff->Draw();
        for (int i = 0; i < NSTAVE - 1; i++)
        {
            hdiffcol_chip_norm[i]->Draw("hist same");

            sprintf(name, "chip %i = %+.1f", i, hdiffcol_chip_norm[i]->GetMean());
            ldiffcol[i] = new TLatex(0.8, 0.90 - 0.05*i, name);
            ldiffcol[i]->SetNDC();
            ldiffcol[i]->SetTextColor(chipColor[i]);
            ldiffcol[i]->Draw("same");

        }
        lt.DrawLatex(0.5, 0.96, "Col");

        //-- diff row
        pdiffrow->cd();
        gPad->SetLogy();
        haxis_diff->GetYaxis()->SetRangeUser(1e-4, 1);
        //haxis_diff->GetXaxis()->SetRangeUser(-25, 25);
        haxis_diff->Draw();
        for (int i = 0; i < NSTAVE - 1; i++)
        {
            hdiffrow_chip_norm[i]->Draw("hist same");

            sprintf(name, "chip %i = %+.1f", i, hdiffrow_chip_norm[i]->GetMean());
            ldiffrow[i] = new TLatex(0.8, 0.90 - 0.05*i, name);
            ldiffrow[i]->SetNDC();
            ldiffrow[i]->SetTextColor(chipColor[i]);
            ldiffrow[i]->Draw("same");

        }
        lt.DrawLatex(0.5, 0.96, "Row");

        for (int i = 1; i < 10; ++i)
        {
            chip_edges.push_back(new TLine((1024*i)-1, 0, (1024*i)-1, 511));
            dead_chip_forward.push_back(new TLine((1024*(i-1))-1, 0, (1024*i)-1, 511));
            dead_chip_backward.push_back(new TLine((1024*(i-1))-1, 511, (1024*i)-1, 0));
        }
        //-- 2D distributions
        for (int i = 0; i < NSTAVE; i++)
        {
            p2d[i]->cd();
            gPad->SetLogz();
            //haxis_2d->DrawCopy();
            gPad->Update();
            if (flip_yaxis)
            {
                reversedaxis = new TGaxis(gPad->GetUxmin(),
                        gPad->GetUymax(),
                        gPad->GetUxmin()-0.001,
                        gPad->GetUymin(),
                        -0.5,
                        511.5,
                        510, "+");
                reversedaxis->SetLabelOffset(-0.03);
                //reversedaxis->Draw();
            }
            h2d_chip_norm[i]->Draw("colz same");
            h2d_chip_norm[i]->GetZaxis()->SetRangeUser(1e-6,1);
            for (int j = 0; j < 8; ++j){chip_edges[j]->Draw();}
            sprintf(name, "Stave %s", stave_name[i].c_str());
            lt.SetTextColor(chipColor[i]);
            lt.DrawLatex(0.5, 0.96, name);
        } // i

        p2d[0]->cd();
        dead_chip_forward[0]->Draw(); dead_chip_backward[0]->Draw();
        p2d[2]->cd();
        dead_chip_forward[0]->Draw(); dead_chip_backward[0]->Draw();
        dead_chip_forward[3]->Draw(); dead_chip_backward[3]->Draw();
        dead_chip_forward[6]->Draw(); dead_chip_backward[6]->Draw();
        dead_chip_forward[7]->Draw(); dead_chip_backward[7]->Draw();
        dead_chip_forward[8]->Draw(); dead_chip_backward[8]->Draw();
        p2d[3]->cd();
        dead_chip_forward[6]->Draw(); dead_chip_backward[6]->Draw();

        initialized = true;
    }

    //-- Update
    com->cd();
    sprintf(name, "Run %i, Number of Events: %.0f", mvtx_run, hnevnt->Integral());
    lnevents->SetText(0.5, 0.95, name);
    for (int i = 0; i < NSTAVE; i++)
    {
        sprintf(name, "stave %s = %.2f", stave_name[i].c_str(), hnhit_chip_norm[i]->GetMean());
        lnhitmean[i]->SetText(0.75, 0.88 - 0.05*i, name);
    }
    for (int i = 0; i < NSTAVE - 1; i++)
    {
        if ( hdiffcol_chip_norm[i]->Integral() > 0 )
        {
            double dcol = hdiffcol_chip_norm[i]->GetBinCenter(hdiffcol_chip_norm[i]->GetMaximumBin());
            sprintf(name, "stave %s = %+.0f", stave_name[i].c_str(), dcol);
        }
        else
            sprintf(name, "");
        ldiffcol[i]->SetText(0.75, 0.83 - 0.05*i, name);

        if ( hdiffrow_chip_norm[i]->Integral() > 0 )
        {
            double drow = hdiffrow_chip_norm[i]->GetBinCenter(hdiffrow_chip_norm[i]->GetMaximumBin());
            sprintf(name, "stave %s = %+.0f", stave_name[i].c_str(), drow);
        }
        else
            sprintf(name, "");
        ldiffrow[i]->SetText(0.75, 0.83 - 0.05*i, name);
    }
    com->Modified();
    phitrate->Modified();
    //pnhit->Modified();
    //pmean->Modified();
    pdiffcol->Modified();
    pdiffrow->Modified();
    for (int i = 0; i < NSTAVE; i++)
        p2d[i]->Modified();
    com->Update();

    return 0;
}

//============================================================//
// Print plots on canvas to a png file

int print_canvas()
{
  char outputFileName[128];
  sprintf(outputFileName, "Run%0.4i.png", mvtx_run);
  com->Print(outputFileName);
  if (show_beam_fit) {
    cBeamCenter->Print(Form("Beam%0.4d.png", mvtx_run));
  }
  return 0;
}
//============================================================//

int print_status()
{
    analysis();
    cout << "==================================================" << endl;
    cout << "==================================================" << endl;
    cout << "== Number of events: " << hnevnt->Integral() << endl;
    for (int i = 0; i < NSTAVE; i++)
    {
        // only print chips with hits
        if ( !(hchip->GetBinContent(i + 1) > 0) )
            continue;

        cout << "==================================================" << endl;
        cout << "==== Chip " << i << endl;
        cout << "== events: " << hnhit_chip[i]->Integral() << endl;
        cout << "== Total hits: " << hchip->GetBinContent(i+1)
            << " (" << hhitrate_chip->GetBinContent(i+1) << ")"
            << endl;
        cout << "== Total Warnings: " << hwarn->GetBinContent(i+1)
            << " (" << hwarnrate_chip->GetBinContent(i+1) << ")"
            << endl;
        cout << "== Total Errors: " << herr->GetBinContent(i+1)
            << " (" << herrrate_chip->GetBinContent(i+1) << ")"
            << endl;
        cout << "== <hits/event>: " << hnhit_chip_norm[i]->GetMean() << endl;
        cout << "== RMS(hits/event): " << hnhit_chip_norm[i]->GetRMS() << endl;

        double mcol, mrow;
        gmpos_chip[i]->GetPoint(0, mcol, mrow);
        cout << "== <col>: " << mcol << endl;
        cout << "== <row>: " << mrow << endl;

        if ( i > 0 )
        {
            double dcol = hdiffcol_chip_norm[i]->GetMaximumBin();
            dcol = hdiffcol_chip_norm[i]->GetBinCenter(dcol);
            cout << "== <dcol>: " << dcol
                << " (" << hdiffcol_chip_norm[i]->GetMean() << ")"
                << endl;

            double drow = hdiffrow_chip_norm[i]->GetMaximumBin();
            drow = hdiffrow_chip_norm[i]->GetBinCenter(drow);
            cout << "== <drow>: " << drow
                << " (" << hdiffrow_chip_norm[i]->GetMean() << ")"
                << endl;
        }
    }
    cout << "==================================================" << endl;
    cout << "==================================================" << endl;
}

//============================================================//

void reset_histos()
{
    hnevnt->Reset();
    hchip->Reset();
    hwarn->Reset();
    herr->Reset();
    for (int i = 0; i < NSTAVE; i++)
    {
        hnhit_chip[i]->Reset();
        h2d_chip[i]->Reset();
        hdiffrow_chip[i]->Reset();
        hdiffcol_chip[i]->Reset();
    } // i
}

//============================================================//

void get_alignment()
{
  string fname(Form("beamcenter/beamcenter_%08d.txt", mvtx_run));
  ofstream fout(fname);

  for (int istave=0; istave<NSTAVE; ++istave)
  {
    double m_col, m_row;
    gmpos_chip[istave]->GetPoint(0, m_col, m_row);
    m_row = 511 - m_row; //chip row flipped wrt histo axis

    fout << istave << " 0 " << fixed << setprecision(2) << m_row << " " <<  0 << " " << m_col << endl;
    cout << istave << " 0 " << fixed << setprecision(2) << m_row << " " <<  0 << " " << m_col << endl;
  }

  fout.close();
  return;
}

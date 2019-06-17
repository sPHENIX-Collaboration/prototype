#include "TpcPrototypeDefs.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using namespace std;

namespace TpcPrototypeDefs
{
//! TPC v1 FEE test stand decoder
namespace FEEv2
{
SampleFit_PowerLawDoubleExp_PDFMaker::SampleFit_PowerLawDoubleExp_PDFMaker()
{
  gStyle->SetOptFit(1111);

  m_canvas = new TCanvas("SampleFit_PowerLawDoubleExp_PDFMaker", "SampleFit_PowerLawDoubleExp_PDFMaker");
  m_pavedtext = new TPaveText(.05, .1, .95, .8);

  m_pavedtext->AddText("SampleFit_PowerLawDoubleExp Fit output");
  m_pavedtext->AddText("A double-component power-law exponential fit of time-dependent ADC pulses.");
  m_pavedtext->AddText("Magenta curve is the sum of the two component, the red and blue curves.");
  m_pavedtext->AddText("Red dot denote the max points");
  m_pavedtext->Draw();

  m_canvas->Print("SampleFit_PowerLawDoubleExp.pdf(");  //open multiplage PDF
}
SampleFit_PowerLawDoubleExp_PDFMaker::~SampleFit_PowerLawDoubleExp_PDFMaker()
{
  if (m_pavedtext) delete m_pavedtext;
  if (m_canvas) delete m_canvas;

  m_canvas = new TCanvas("SampleFit_PowerLawDoubleExp_PDFMaker", "SampleFit_PowerLawDoubleExp_PDFMaker");
  m_pavedtext = new TPaveText(.05, .1, .95, .8);

  m_pavedtext->AddText("SampleFit_PowerLawDoubleExp Fit output");
  m_pavedtext->AddText("End of pages");
  m_pavedtext->Draw();

  m_canvas->Print("SampleFit_PowerLawDoubleExp.pdf)");  //close multiplage PDF
}

void SampleFit_PowerLawDoubleExp_PDFMaker::MakeSectionPage(const string &title)
{
  if (m_pavedtext) delete m_pavedtext;
  if (m_canvas) delete m_canvas;

  m_canvas = new TCanvas("SampleFit_PowerLawDoubleExp_PDFMaker", "SampleFit_PowerLawDoubleExp_PDFMaker");

  m_pavedtext = new TPaveText(.05, .1, .95, .8);

  m_pavedtext->AddText(title.c_str());
  m_pavedtext->Draw();

  m_canvas->Print("SampleFit_PowerLawDoubleExp.pdf");
}

bool SampleFit_PowerLawDoubleExp(        //
    const std::vector<double> &samples,  //
    double &peak,                        //
    double &peak_sample,                 //
    double &pedestal,                    //
    std::map<int, double> &parameters_io,
    const int verbosity)
{
  static const int n_parameter = 7;

  // inital guesses
  int peakPos = 0.;

  //  assert(samples.size() == n_samples);
  const int n_samples = samples.size();

  TGraph gpulse(n_samples);
  for (int i = 0; i < n_samples; i++)
  {
    (gpulse.GetX())[i] = i;

    (gpulse.GetY())[i] = samples[i];
  }

  //Saturation correction - Abhisek
  //  for (int ipoint = 0; ipoint < gpulse.GetN(); ipoint++)
  //    if ((gpulse.GetY())[ipoint] >= ((1 << 10) - 10)  // drop point if touching max or low limit on ADCs
  //        or (not isnormal((gpulse.GetY())[ipoint])))
  //    {
  //      gpulse.RemovePoint(ipoint);
  //      ipoint--;
  //    }

  pedestal = gpulse.GetY()[0];  //(double) PEDESTAL;
  double peakval = pedestal;
  const double risetime = 1.5;

  for (int iSample = 0; iSample < n_samples - risetime * 3; iSample++)
  {
    if (abs(gpulse.GetY()[iSample] - pedestal) > abs(peakval - pedestal))
    {
      peakval = gpulse.GetY()[iSample];
      peakPos = iSample;
    }
  }
  peakval -= pedestal;

  if (verbosity)
  {
    cout << "SampleFit_PowerLawDoubleExp - "
         << "pedestal = " << pedestal << ", "
         << "peakval = " << peakval << ", "
         << "peakPos = " << peakPos << endl;
  }

  // build default value
  struct default_values_t
  {
    default_values_t(double default_value, double min_value, double max_value)
      : def(default_value)
      , min(min_value)
      , max(max_value)
    {
    }
    double def;
    double min;
    double max;
  };

  vector<default_values_t> default_values(n_parameter, default_values_t(numeric_limits<double>::signaling_NaN(), numeric_limits<double>::signaling_NaN(), numeric_limits<double>::signaling_NaN()));

  default_values[0] = default_values_t(peakval * .7, peakval * -1.5, peakval * 1.5);
  default_values[1] = default_values_t(peakPos - risetime, peakPos - 3 * risetime, peakPos + risetime);
  default_values[2] = default_values_t(5., 1, 10.);
  default_values[3] = default_values_t(risetime, risetime * .2, risetime * 10);
  default_values[4] = default_values_t(pedestal, pedestal - abs(peakval), pedestal + abs(peakval));
  //  default_values[5] = default_values_t(0.3, 0, 1);
  //  default_values[6] = default_values_t(5, risetime * .2, risetime * 10);
  default_values[5] = default_values_t(0, 0, 0);  // disable 2nd component
  default_values[6] = default_values_t(risetime, risetime, risetime);

  // fit function
  TF1 fits("f_SignalShape_PowerLawDoubleExp", SignalShape_PowerLawDoubleExp, 0., n_samples, n_parameter);
  fits.SetParNames("Amplitude", "Sample Start", "Power", "Peak Time 1", "Pedestal", "Amplitude ratio", "Peak Time 2");

  for (int i = 0; i < n_parameter; ++i)
  {
    if (parameters_io.find(i) == parameters_io.end())
    {
      fits.SetParameter(i, default_values[i].def);

      if (default_values[i].min < default_values[i].max)
      {
        fits.SetParLimits(i, default_values[i].min, default_values[i].max);
      }
      else
      {
        fits.FixParameter(i, default_values[i].def);
      }

      if (verbosity)
      {
        cout << "SampleFit_PowerLawDoubleExp - parameter [" << i << "]: "
             << "default value = " << default_values[i].def
             << ", min value = " << default_values[i].min
             << ", max value = " << default_values[i].max << endl;
      }
    }
    else
    {
//      fits.SetParLimits(i, parameters_io[i], parameters_io[i]);
      fits.SetParameter(i, parameters_io[i]);
      fits.FixParameter(i, parameters_io[i]);

      if (verbosity)
      {
        cout << "SampleFit_PowerLawDoubleExp - parameter [" << i << "]: fixed to " << parameters_io[i] << endl;
      }
    }
  }

  if (verbosity <= 1)
    gpulse.Fit(&fits, "QRN0W", "goff", 0., (double) n_samples);
  else
    gpulse.Fit(&fits, "RN0VW+", "goff", 0., (double) n_samples);

  // store results
  pedestal = fits.GetParameter(4);

  const double peakpos1 = fits.GetParameter(3);
  const double peakpos2 = fits.GetParameter(6);
  double max_peakpos = fits.GetParameter(1) + (peakpos1 > peakpos2 ? peakpos1 : peakpos2);
  if (max_peakpos > n_samples - 1) max_peakpos = n_samples - 1;

  if (fits.GetParameter(0) > 0)
    peak_sample = fits.GetMaximumX(fits.GetParameter(1), max_peakpos);
  else
    peak_sample = fits.GetMinimumX(fits.GetParameter(1), max_peakpos);

  peak = fits.Eval(peak_sample) - pedestal;

  if (verbosity)
  {
    static int id = 0;
    ++id;

    string c_name(string("SampleFit_PowerLawDoubleExp_") + to_string(id));

    TCanvas *canvas = new TCanvas(
        c_name.c_str(), c_name.c_str());
    canvas->Update();

    TGraph *g_plot = static_cast<TGraph *>(gpulse.DrawClone("ap*l"));
    g_plot->SetTitle((string("ADC data and fit #") + to_string(id) + string(";Sample number;ADC value")).c_str());

    fits.SetLineColor(kMagenta);
    fits.DrawClone("same");
    fits.Print();

    TF1 f1("f_SignalShape_PowerLawExp1", SignalShape_PowerLawExp, 0., n_samples, 5);
    f1.SetParameters(
        fits.GetParameter(0) * (1 - fits.GetParameter(5)) / pow(fits.GetParameter(3), fits.GetParameter(2)) * exp(fits.GetParameter(2)),
        fits.GetParameter(1),
        fits.GetParameter(2),
        fits.GetParameter(2) / fits.GetParameter(3),
        fits.GetParameter(4));
    f1.SetLineColor(kBlue);
    f1.DrawClone("same");

    TF1 f2("f_SignalShape_PowerLawExp2", SignalShape_PowerLawExp, 0., n_samples, 5);
    f2.SetParameters(
        fits.GetParameter(0) * fits.GetParameter(5) / pow(fits.GetParameter(6), fits.GetParameter(2)) * exp(fits.GetParameter(2)),
        fits.GetParameter(1),
        fits.GetParameter(2),
        fits.GetParameter(2) / fits.GetParameter(6),
        fits.GetParameter(4));
    f2.SetLineColor(kRed);
    f2.DrawClone("same");

    TGraph g_max(1);

    g_max.GetX()[0] = peak_sample;
    g_max.GetY()[0] = peak + pedestal;

    g_max.SetMarkerStyle(kFullCircle);
    g_max.SetMarkerSize(2);
    g_max.SetMarkerColor(kRed);

    static_cast<TGraph *>(g_max.DrawClone("p"));

    canvas->Update();

    //    if (id == 1)
    //    {
    //      canvas->Print("SampleFit_PowerLawDoubleExp.pdf(");
    //    }
    canvas->Print("SampleFit_PowerLawDoubleExp.pdf");
  }

  for (int i = 0; i < n_parameter; ++i)
  {
    parameters_io[i] = fits.GetParameter(i);
  }

  if (verbosity)
  {
    cout << "SampleFit_PowerLawDoubleExp - "
         << "peak_sample = " << peak_sample << ", "
         << "max_peakpos = " << max_peakpos << ", "
         << "fits.GetParameter(1) = " << fits.GetParameter(1) << ", "
         << "peak = " << peak << ", "
         << "pedestal = " << pedestal << endl;
  }

  return true;
}

double
SignalShape_PowerLawExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;
  //double  signal = (-1)*par[0]*pow((x[0]-par[1]),par[2])*exp(-(x[0]-par[1])*par[3]);
  double signal = par[0] * pow((x[0] - par[1]), par[2]) * exp(-(x[0] - par[1]) * par[3]);
  return pedestal + signal;
}

double
SignalShape_PowerLawDoubleExp(double *x, double *par)
{
  double pedestal = par[4];
  //                        + ((x[0] - 1.5 * par[1]) > 0) * par[5];  // quick fix on exting tails on the signal function
  if (x[0] < par[1])
    return pedestal;
  //double  signal = (-1)*par[0]*pow((x[0]-par[1]),par[2])*exp(-(x[0]-par[1])*par[3]);
  //  peak / pow(fits.GetParameter(2) / fits.GetParameter(3), fits.GetParameter(2)) * exp(fits.GetParameter(2)) = fits.GetParameter(0);  // exact peak height is (p0*Power(p2/p3,p2))/Power(E,p2)
  //  fits.GetParameter(2) / peak_shift =  fits.GetParameter(3);  // signal peak time

  double signal =                                                                                         //
      par[0]                                                                                              //
      * pow((x[0] - par[1]), par[2])                                                                      //
      * (((1. - par[5]) / pow(par[3], par[2]) * exp(par[2])) * exp(-(x[0] - par[1]) * (par[2] / par[3]))  //
         + (par[5] / pow(par[6], par[2]) * exp(par[2])) * exp(-(x[0] - par[1]) * (par[2] / par[6]))       //
        );
  return pedestal + signal;
}


}  // namespace FEEv2

}  // namespace TpcPrototypeDefs

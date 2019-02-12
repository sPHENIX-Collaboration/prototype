#include "RawTower_Prototype2.h"

#include "PROTOTYPE2_FEM.h"

#include <calobase/RawTowerDefs.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include <cassert>


using namespace std;

RawTower_Prototype2::RawTower_Prototype2() :
    towerid(~0), // initialize all bits on
    energy(0), time(NAN), HBD_channel(-1)
{
  fill_n(signal_samples,NSAMPLES,-9999);
}

RawTower_Prototype2::RawTower_Prototype2(const RawTower & tower)
{
  towerid = (tower.get_id());
  energy = (tower.get_energy());
  time = (tower.get_time());
  HBD_channel = -1;
  fill_n(signal_samples,NSAMPLES,-9999);
}

RawTower_Prototype2::RawTower_Prototype2(RawTowerDefs::keytype id) :
    towerid(id), energy(0), time(NAN), HBD_channel(-1)
{
  fill_n(signal_samples,NSAMPLES,-9999);
}

RawTower_Prototype2::RawTower_Prototype2(const unsigned int icol, const unsigned int irow) :
    towerid(0), energy(0), time(NAN), HBD_channel(-1)
{
  towerid = RawTowerDefs::encode_towerid(RawTowerDefs::NONE, icol, irow);
  fill_n(signal_samples,NSAMPLES,-9999);
}

RawTower_Prototype2::RawTower_Prototype2(const RawTowerDefs::CalorimeterId caloid,
    const unsigned int ieta, const unsigned int iphi) :
    towerid(0), energy(0), time(NAN), HBD_channel(-1)
{
  towerid = RawTowerDefs::encode_towerid(caloid, ieta, iphi);
  fill_n(signal_samples,NSAMPLES,-9999);
}

void
RawTower_Prototype2::Reset()
{
  energy = 0;
  time = NAN;
  fill_n(signal_samples,NSAMPLES,-9999);
}

int
RawTower_Prototype2::isValid() const
{
  return get_energy() != 0;
}

void
RawTower_Prototype2::identify(std::ostream& os) const
{
  os << "RawTower_Prototype2: etabin: " << get_bineta() << ", phibin: " << get_binphi()
      << " energy=" << get_energy() << std::endl;
}

void
RawTower_Prototype2::set_signal_samples(int i, RawTower_Prototype2::signal_type sig)
{
  assert(i>=0);
  assert(i<NSAMPLES);
  signal_samples[i] = sig;
}

RawTower_Prototype2::signal_type
RawTower_Prototype2::get_signal_samples(int i) const
{
  assert(i>=0);
  assert(i<NSAMPLES);
  return signal_samples[i];
}

double
RawTower_Prototype2::get_energy_power_law_exp(int verbosity)
{
  double peak = NAN;
  double peak_sample = NAN;
  double pedstal = NAN;

  vector<double> vec_signal_samples;
  for (int i = 0; i < NSAMPLES; i++)
    {
      vec_signal_samples.push_back(signal_samples[i]);
    }

  PROTOTYPE2_FEM::
  SampleFit_PowerLawExp(vec_signal_samples, peak, peak_sample, pedstal, verbosity);

  return peak;
}

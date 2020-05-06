/**
 * @file trackbase/TpcPrototypeCluster.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TpcPrototypeCluster
 */
#include "TpcPrototypeCluster.h"

#include <cmath>
#include <utility>  // for swap

namespace
{
// square convenience function
template <class T>
inline constexpr T square(const T& x)
{
  return x * x;
}

// get unique index in cov. matrix array from i and j
inline unsigned int covarIndex(unsigned int i, unsigned int j)
{
  if (i > j) std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

// rotate size or covariant matrix to polar coordinates and return the phi component
template <float (TpcPrototypeCluster::*accessor)(unsigned int, unsigned int) const>
float rotate(const TpcPrototypeCluster* cluster)
{
  const auto phi = -std::atan2(cluster->getY(), cluster->getX());
  const auto cosphi = std::cos(phi);
  const auto sinphi = std::sin(phi);

  return square(sinphi) * (cluster->*accessor)(0, 0) +
         square(cosphi) * (cluster->*accessor)(1, 1) +
         2. * cosphi * sinphi * (cluster->*accessor)(0, 1);
}

}  // namespace

TpcPrototypeCluster::TpcPrototypeCluster()
  : m_cluskey(TrkrDefs::CLUSKEYMAX)
  , m_isGlobal(true)
  , m_adc(0xFFFFFFFF)
  , min_sample(-1)
  , max_sample(-1)
  , min_pad_azimuth(-1)
  , max_pad_azimuth(-1)
  , peak(NAN)
  , peak_sample(NAN)
  , pedstal(NAN)
  , avg_pad_radial(-1)
  , avg_pad_azimuth(NAN)
  , size_pad_radial(-1)
  , size_pad_azimuth(-1)
  , delta_azimuth_bin(NAN)
  , delta_z(NAN)
{
  for (int i = 0; i < 3; ++i) m_pos[i] = NAN;

  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      TpcPrototypeCluster::setSize(i, j, NAN);
      TpcPrototypeCluster::setError(i, j, NAN);
    }
  }
}

void TpcPrototypeCluster::identify(std::ostream& os) const
{
  os << "---TpcPrototypeCluster--------------------" << std::endl;
  os << "clusid: " << getClusKey() << std::dec << std::endl;

  os << " (x,y,z) =  (" << getPosition(0);
  os << ", " << getPosition(1) << ", ";
  os << getPosition(2) << ") cm";
  if (m_isGlobal)
    os << " - global coordinates" << std::endl;
  else
    os << " - local coordinates" << std::endl;

  os << " adc = " << getAdc() << std::endl;

  os << " size phi = " << getPhiSize();
  os << " cm, size z = " << getZSize() << " cm" << std::endl;

  os << "         ( ";
  os << getSize(0, 0) << " , ";
  os << getSize(0, 1) << " , ";
  os << getSize(0, 2) << " )" << std::endl;
  os << "  size = ( ";
  os << getSize(1, 0) << " , ";
  os << getSize(1, 1) << " , ";
  os << getSize(1, 2) << " )" << std::endl;
  os << "         ( ";
  os << getSize(2, 0) << " , ";
  os << getSize(2, 1) << " , ";
  os << getSize(2, 2) << " )" << std::endl;

  os << "         ( ";
  os << getError(0, 0) << " , ";
  os << getError(0, 1) << " , ";
  os << getError(0, 2) << " )" << std::endl;
  os << "  err  = ( ";
  os << getError(1, 0) << " , ";
  os << getError(1, 1) << " , ";
  os << getError(1, 2) << " )" << std::endl;
  os << "         ( ";
  os << getError(2, 0) << " , ";
  os << getError(2, 1) << " , ";
  os << getError(2, 2) << " )" << std::endl;

  //  std::set<int> pad_radials;
  os << "  pad_radials  = ( ";
  for (auto& element : pad_radials) os << element << " ";
  os << ")" << std::endl;

  //  std::set<int> pad_azimuths;
  os << "  pad_azimuths  = ( ";
  for (auto& element : pad_azimuths) os << element << " ";
  os << ")" << std::endl;
  //  std::set<int> samples;
  os << "  samples  = ( ";
  for (auto& element : samples) os << element << " ";
  os << ")" << std::endl;
  //
  //  std::map<int, std::vector<double>> pad_radial_samples;
  os << "  pad_radial_samples  = [ " << std::endl;
  for (auto& pair : pad_radial_samples)
  {
    os << pair.first << "->( ";
    for (auto& element : pair.second) os << element << " ";
    os << ")" << std::endl;
  }
  os << "]" << std::endl;
  //  std::map<int, std::vector<double>> pad_azimuth_samples;
  os << "  pad_azimuth_samples  = [ " << std::endl;
  for (auto& pair : pad_azimuth_samples)
  {
    os << pair.first << "->( ";
    for (auto& element : pair.second) os << element << " ";
    os << ")" << std::endl;
  }
  os << "]" << std::endl;
  //  std::vector<double> sum_samples;
  os << "  sum_samples  = ( ";
  for (auto& element : sum_samples) os << element << " ";
  os << ")" << std::endl;
  //
  //  int min_sample;
  os << "  min_sample  =  " << min_sample << std::endl;
  //  int max_sample;
  os << "  max_sample  =  " << max_sample << std::endl;
  //  int min_pad_azimuth;
  os << "  min_pad_azimuth  =  " << min_pad_azimuth << std::endl;
  //  int max_pad_azimuth;
  os << "  max_pad_azimuth  =  " << max_pad_azimuth << std::endl;
  //
  //  double peak;
  os << "  peak  =  " << peak << std::endl;
  //  double peak_sample;
  os << "  peak_sample  =  " << peak_sample << std::endl;
  //  double pedstal;
  os << "  pedstal  =  " << pedstal << std::endl;
  //
  //  //    std::map<int, double> pad_radial_peaks; // radial always have size = 1 in this analysis
  //  std::map<int, double> pad_azimuth_peaks;
  os << "  pad_azimuth_peaks  = ( ";
  for (auto& element : pad_azimuth_peaks) os << element.first << "->" << element.second << " ";
  os << ")" << std::endl;
  //
  //  //! pad coordinate
  //  int avg_pad_radial;
  os << "  avg_pad_radial  =  " << avg_pad_radial << std::endl;
  //  double avg_pad_azimuth;
  os << "  avg_pad_azimuth  =  " << avg_pad_azimuth << std::endl;
  //
  //  //! cluster size in units of pad bins
  //  int size_pad_radial;
  os << "  size_pad_radial  =  " << size_pad_radial << std::endl;
  //  int size_pad_azimuth;
  os << "  size_pad_azimuth  =  " << size_pad_azimuth << std::endl;
  //
  //  //! pad bin size
  //  //! phi size per pad in rad
  //  double delta_azimuth_bin;
  os << "  delta_azimuth_bin  =  " << delta_azimuth_bin << std::endl;
  //  //! z size per ADC sample bin
  //  double delta_z;
  os << "  delta_z  =  " << delta_z << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int TpcPrototypeCluster::isValid() const
{
  if (m_cluskey == TrkrDefs::CLUSKEYMAX) return 0;
  for (int i = 0; i < 3; ++i)
  {
    if (isnan(getPosition(i))) return 0;
  }
  if (m_adc == 0xFFFFFFFF) return 0;
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (isnan(getSize(i, j))) return 0;
      if (isnan(getError(i, j))) return 0;
    }
  }

  return 1;
}

void TpcPrototypeCluster::setSize(unsigned int i, unsigned int j, float value)
{
  m_size[covarIndex(i, j)] = value;
  return;
}

float TpcPrototypeCluster::getSize(unsigned int i, unsigned int j) const
{
  return m_size[covarIndex(i, j)];
}

void TpcPrototypeCluster::setError(unsigned int i, unsigned int j, float value)
{
  m_err[covarIndex(i, j)] = value;
  return;
}

float TpcPrototypeCluster::getError(unsigned int i, unsigned int j) const
{
  return m_err[covarIndex(i, j)];
}

float TpcPrototypeCluster::getPhiSize() const
{
  return 2 * std::sqrt(rotate<&TpcPrototypeCluster::getSize>(this));
}

float TpcPrototypeCluster::getZSize() const
{
  return 2. * sqrt(getSize(2, 2));
}

float TpcPrototypeCluster::getPhiError() const
{
  const float rad = std::sqrt(square(m_pos[0]) + square(m_pos[1]));
  if (rad > 0) return getRPhiError() / rad;
  return 0;
}

float TpcPrototypeCluster::getRPhiError() const
{
  return std::sqrt(rotate<&TpcPrototypeCluster::getError>(this));
}

float TpcPrototypeCluster::getZError() const
{
  return std::sqrt(getError(2, 2));
}

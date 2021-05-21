/**
 * @file trackbase/TpcPrototypeCluster.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Version 1 of TrkrCluster
 */
#ifndef TRACKBASE_TRKRCLUSTERV1_H
#define TRACKBASE_TRKRCLUSTERV1_H

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>

class PHObject;

/**
 * @brief Version 1 of TrkrCluster
 *
 * Note - D. McGlinchey June 2018:
 *   CINT does not like "override", so ignore where CINT
 *   complains. Should be checked with ROOT 6 once
 *   migration occurs.
 */
class TpcPrototypeCluster : public TrkrCluster
{
 public:
  //! ctor
  TpcPrototypeCluster();

  //!dtor
  virtual ~TpcPrototypeCluster() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset() {}
  virtual int isValid() const;
  virtual PHObject* CloneMe() const { return new TpcPrototypeCluster(*this); }
  virtual void setClusKey(TrkrDefs::cluskey id) { m_cluskey = id; }
  virtual TrkrDefs::cluskey getClusKey() const { return m_cluskey; }
  //
  // cluster position
  //
  virtual float getX() const { return m_pos[0]; }
  virtual void setX(float x) { m_pos[0] = x; }
  virtual float getY() const { return m_pos[1]; }
  virtual void setY(float y) { m_pos[1] = y; }
  virtual float getZ() const { return m_pos[2]; }
  virtual void setZ(float z) { m_pos[2] = z; }
  virtual float getPosition(int coor) const { return m_pos[coor]; }
  virtual void setPosition(int coor, float xi) { m_pos[coor] = xi; }
  virtual void setGlobal() { m_isGlobal = true; }
  virtual void setLocal() { m_isGlobal = false; }
  virtual bool isGlobal() { return m_isGlobal; }
  //
  // cluster info
  //
  virtual unsigned int getAdc() const { return m_adc; }
  virtual void setAdc(unsigned int adc) { m_adc = adc; }
  virtual float getSize(unsigned int i, unsigned int j) const;        //< get cluster dimension covar
  virtual void setSize(unsigned int i, unsigned int j, float value);  //< set cluster dimension covar

  virtual float getError(unsigned int i, unsigned int j) const;        //< get cluster error covar
  virtual void setError(unsigned int i, unsigned int j, float value);  //< set cluster error covar

  //
  // convenience interface
  //
  virtual float getPhiSize() const;
  virtual float getZSize() const;

  virtual float getRPhiError() const;
  virtual float getPhiError() const;
  virtual float getZError() const;

  // prototype specific interface.

  double getAvgPadAzimuth() const
  {
    return avg_pad_azimuth;
  }

  void setAvgPadAzimuth(double avgPadAzimuth)
  {
    avg_pad_azimuth = avgPadAzimuth;
  }

  int getAvgPadRadial() const
  {
    return avg_pad_radial;
  }

  void setAvgPadRadial(int avgPadRadial)
  {
    avg_pad_radial = avgPadRadial;
  }

  double getDeltaAzimuthBin() const
  {
    return delta_azimuth_bin;
  }

  void setDeltaAzimuthBin(double deltaAzimuthBin)
  {
    delta_azimuth_bin = deltaAzimuthBin;
  }

  double getDeltaZ() const
  {
    return delta_z;
  }

  void setDeltaZ(double deltaZ)
  {
    delta_z = deltaZ;
  }

  int getMaxPadAzimuth() const
  {
    return max_pad_azimuth;
  }

  void setMaxPadAzimuth(int maxPadAzimuth)
  {
    max_pad_azimuth = maxPadAzimuth;
  }

  int getMaxSample() const
  {
    return max_sample;
  }

  void setMaxSample(int maxSample)
  {
    max_sample = maxSample;
  }

  int getMinPadAzimuth() const
  {
    return min_pad_azimuth;
  }

  void setMinPadAzimuth(int minPadAzimuth)
  {
    min_pad_azimuth = minPadAzimuth;
  }

  int getMinSample() const
  {
    return min_sample;
  }

  void setMinSample(int minSample)
  {
    min_sample = minSample;
  }

  const std::map<int, double>& getPadAzimuthPeaks() const
  {
    return pad_azimuth_peaks;
  }

  void setPadAzimuthPeaks(const std::map<int, double>& padAzimuthPeaks)
  {
    pad_azimuth_peaks = padAzimuthPeaks;
  }

  const std::map<int, std::vector<double>>& getPadAzimuthSamples() const
  {
    return pad_azimuth_samples;
  }

  void setPadAzimuthSamples(const std::map<int, std::vector<double>>& padAzimuthSamples)
  {
    pad_azimuth_samples = padAzimuthSamples;
  }

  const std::set<int>& getPadAzimuths() const
  {
    return pad_azimuths;
  }

  void setPadAzimuths(const std::set<int>& padAzimuths)
  {
    pad_azimuths = padAzimuths;
  }

  const std::map<int, std::vector<double>>& getPadRadialSamples() const
  {
    return pad_radial_samples;
  }

  void setPadRadialSamples(const std::map<int, std::vector<double>>& padRadialSamples)
  {
    pad_radial_samples = padRadialSamples;
  }

  const std::set<int>& getPadRadials() const
  {
    return pad_radials;
  }

  void setPadRadials(const std::set<int>& padRadials)
  {
    pad_radials = padRadials;
  }

  double getPeak() const
  {
    return peak;
  }

  void setPeak(double peak)
  {
    this->peak = peak;
  }

  double getPeakSample() const
  {
    return peak_sample;
  }

  void setPeakSample(double peakSample)
  {
    peak_sample = peakSample;
  }

  double getPedstal() const
  {
    return pedstal;
  }

  void setPedstal(double pedstal)
  {
    this->pedstal = pedstal;
  }

  const std::set<int>& getSamples() const
  {
    return samples;
  }

  void setSamples(const std::set<int>& samples)
  {
    this->samples = samples;
  }

  int getSizePadAzimuth() const
  {
    return size_pad_azimuth;
  }

  void setSizePadAzimuth(int sizePadAzimuth)
  {
    size_pad_azimuth = sizePadAzimuth;
  }

  int getSizePadRadial() const
  {
    return size_pad_radial;
  }

  void setSizePadRadial(int sizePadRadial)
  {
    size_pad_radial = sizePadRadial;
  }

  const std::vector<double>& getSumSamples() const
  {
    return sum_samples;
  }

  void setSumSamples(const std::vector<double>& sumSamples)
  {
    sum_samples = sumSamples;
  }

 protected:
  TrkrDefs::cluskey m_cluskey;  //< unique identifier within container
  float m_pos[3];               //< mean position x,y,z
  bool m_isGlobal;              //< flag for coord sys (true = global)
  unsigned int m_adc;           //< cluster sum adc (D. McGlinchey - Do we need this?)
  float m_size[6];              //< size covariance matrix (packed storage) (+/- cm^2)
  float m_err[6];               //< covariance matrix: rad, arc and z

  // prototype specifics

  std::set<int> pad_radials;
  std::set<int> pad_azimuths;
  std::set<int> samples;

  std::map<int, std::vector<double>> pad_radial_samples;
  std::map<int, std::vector<double>> pad_azimuth_samples;
  std::vector<double> sum_samples;

  int min_sample;
  int max_sample;
  int min_pad_azimuth;
  int max_pad_azimuth;

  double peak;
  double peak_sample;
  double pedstal;

  //    std::map<int, double> pad_radial_peaks; // radial always have size = 1 in this analysis
  std::map<int, double> pad_azimuth_peaks;

  //! pad coordinate
  int avg_pad_radial;
  double avg_pad_azimuth;

  //! cluster size in units of pad bins
  int size_pad_radial;
  int size_pad_azimuth;

  //! pad bin size
  //! phi size per pad in rad
  double delta_azimuth_bin;
  //! z size per ADC sample bin
  double delta_z;

  ClassDefOverride(TpcPrototypeCluster, 1)
};

#endif  //TRACKBASE_TRKRCLUSTERV1_H

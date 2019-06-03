// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PROTOTYPE2_RAWTOWERTEMPERATURE_H
#define PROTOTYPE2_RAWTOWERTEMPERATURE_H

#include <calobase/RawTower.h>

#include <calobase/RawTowerDefs.h>

#include <ctime>
#include <iostream>
#include <vector>

class RawTower_Temperature : public RawTower
{
 public:
  RawTower_Temperature();
  RawTower_Temperature(const unsigned int icol, const unsigned int irow);
  RawTower_Temperature(RawTowerDefs::keytype id);
  virtual ~RawTower_Temperature() {}

  void Reset();
  int isValid() const { return get_nr_entries(); }
  void identify(std::ostream &os = std::cout) const;
  void print(std::ostream &os = std::cout) const;

  int get_column() const { return RawTowerDefs::decode_index1(towerid); }
  int get_row() const { return RawTowerDefs::decode_index2(towerid); }

  void set_id(RawTowerDefs::keytype id) { towerid = id; }
  RawTowerDefs::keytype get_id() const { return towerid; }

  int get_nr_entries() const { return temperatures.size(); }

  int add_entry(const int eventnr, const time_t t, const float temp)
  {
    eventnumbers.push_back(eventnr);
    times.push_back(t);
    temperatures.push_back(temp);
    return get_nr_entries();
  }

  float get_temperature_from_entry(const unsigned int entry) const
  {
    if (entry >= temperatures.size())
      return -1;
    return temperatures[entry];
  }

  time_t get_time_from_entry(const unsigned int entry) const
  {
    if (entry >= times.size())
      return 0;  // 1970...
    return times[entry];
  }

  int get_eventnumber_from_entry(const unsigned int entry) const
  {
    if (entry >= eventnumbers.size())
      return -1;
    return eventnumbers[entry];
  }

  float get_temperature_from_time(const time_t t) const;

  //---Raw data
  // access------------------------------------------------------------

 protected:
  RawTowerDefs::keytype towerid;

  //! Temperature readings
  //! since we do not have more than 100 entries per run typically,
  //! we trade efficiency for some simplicity and just use some vectors.
  //!
  std::vector<int> eventnumbers;
  std::vector<time_t> times;
  std::vector<float> temperatures;

  ClassDef(RawTower_Temperature, 1)
};

#endif

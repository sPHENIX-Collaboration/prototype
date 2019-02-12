#ifndef PROTOTYPE2_GENERICUNPACKPRDFF_H
#define PROTOTYPE2_GENERICUNPACKPRDFF_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class Event;
class RawTowerContainer;

class GenericUnpackPRDF : public SubsysReco {
public:
  GenericUnpackPRDF(const std::string &detector);

  int InitRun(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  void CreateNodeTree(PHCompositeNode *topNode);

  //! add stuff to be unpacked
  void add_channel(const int packet_id, //! packet id
                   const int channel,   //! channel in packet
                   const int tower_id   //! output tower id
  );

private:
  std::string _detector;

  //! packet_id, channel number to define a hbd_channel
  typedef std::pair<int, int> hbd_channel_typ;

  //! list of hbd_channel -> channel id which is also tower id
  typedef std::map<hbd_channel_typ, int> hbd_channel_map;

  hbd_channel_map _hbd_channel_map;

  // output -> Towers
  RawTowerContainer *_towers;
};

#endif

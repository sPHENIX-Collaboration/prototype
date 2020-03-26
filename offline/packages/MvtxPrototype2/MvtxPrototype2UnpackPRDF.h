// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTX_P2_UnpackPRDF_H
#define MVTX_P2_UnpackPRDF_H

#include "MvtxPrototype2Geom.h"

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

#include <map>

class Event;
class Packet;
class Packet_hbd_fpgashort;
class TrkrHitSetContainer;

class MvtxPrototype2UnpackPRDF : public SubsysReco
{
public:
  MvtxPrototype2UnpackPRDF();

  int
  Init(PHCompositeNode *topNode);

  int
  InitRun(PHCompositeNode *topNode);

  int
  process_event(PHCompositeNode *topNode);

  int
  End(PHCompositeNode *topNode);

  void
  CreateNodeTree(PHCompositeNode *topNode);

  void MakeHits();

  int DecodeRow(int val) const;
  int DecodeCol(int val) const;

private:

  PHCompositeNode* dstNode;

  Event* _event;
  Packet_hbd_fpgashort* _packet;
  TrkrHitSetContainer *_hitsetcon;

  static std::map <std::pair<int,int>,std::pair<int,int>> s_map_chips; //<ruid, ruchn> to <stave, chipID>
  static std::map <int, int> s_map_layers; ///< <stave, layer> stave to layer index

  int _nevents;
  int _verbosity;
  bool _first;

  int _nevent_per_chip[NLAYER][NCHIP];
  int _npixel_per_chip[NLAYER][NCHIP];

};

#endif //**MvtxPrototype2UnpackPRDFF**//

#ifndef __MvtxPrototype2UnpackPRDFF__
#define __MvtxPrototype2UnpackPRDFF__

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include "MvtxPrototype2Geom.h"

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

	int GetStaveID(int ruid, int chid);

	int GetChipID(int ruid, int chid);

	int DecodeRow(int val) const;
	int DecodeCol(int val) const;

private:

	PHCompositeNode* dstNode;

  Event* _event;
  Packet_hbd_fpgashort* _packet;
	TrkrHitSetContainer *_hitsetcon;

  map<pair<int,int>,pair<int,int>> mMapChips;
  int _nevents;
	int _verbosity;
	bool _first;

	int _nevent_per_chip[NSTAVE][NCHIP];
	int _npixel_per_chip[NSTAVE][NCHIP];

};

#endif //**MvtxPrototype2UnpackPRDFF**//

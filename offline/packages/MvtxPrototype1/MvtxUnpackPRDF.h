#ifndef __MvtxUnpackPRDFF__
#define __MvtxUnpackPRDFF__

//* Unpacks raw HCAL PRDF files *//
//Abhisek Sen

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

#define NCHIP 4
#define NROW 512
#define NREGION 32
#define NCOL_PER_REGION 32

class Event;
class Packet;
class Packet_hbd_fpgashort;
class TrkrHitSetContainer;

class MvtxUnpackPRDF : public SubsysReco
{
public:
  MvtxUnpackPRDF();

  int
  Init(PHCompositeNode *topNode);

  int
  InitRun(PHCompositeNode *topNode);

  int
  process_event(PHCompositeNode *topNode);

  int
  End(PHCompositeNode *topNode);

	void
	Verbosity(int v) {_verbosity = v;}

  void
  CreateNodeTree(PHCompositeNode *topNode);

	void MakeHits();

private:

	PHCompositeNode* dstNode;

  Event* _event;
  Packet_hbd_fpgashort* _packet;
	TrkrHitSetContainer *_hitsetcon;

  int _nevents;
	int _verbosity;
	bool _first;

	int _nevent_per_chip[NCHIP];
	int _npixel_per_chip[NCHIP];

};

#endif //**MvtxUnpackPRDFF**//

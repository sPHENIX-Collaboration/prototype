#ifndef __MvtxRunInfoUnpackPRDFF__
#define __MvtxRunInfoUnpackPRDFF__

//* RunInfoUnpacks raw HCAL PRDF files *//

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

class Event;
class Packet;
class Packet_hbd_fpgashort;

class MvtxRunInfoUnpackPRDF : public SubsysReco
{
public:
  MvtxRunInfoUnpackPRDF();

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

private:

	PHCompositeNode* dstNode;

  Event* _event;
  Packet_hbd_fpgashort* _packet;

};

#endif //**MvtxRunInfoUnpackPRDFF**//

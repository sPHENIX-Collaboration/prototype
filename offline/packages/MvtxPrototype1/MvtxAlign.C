#include "MvtxAlign.h"

#include <mvtx/MvtxDefUtil.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

MvtxAlign::MvtxAlign(const string &name ) :
  SubsysReco(name),
  clusters_(NULL),
  fdir_(""),
  runnumber_(0),
  apff_(true),
  _timer(PHTimeServer::get()->insert_new(name))
{

}

int MvtxAlign::InitRun(PHCompositeNode* topNode)
{
  // get the run number
  if ( apff_ )
  {
    recoConsts *rc = recoConsts::instance();
    runnumber_ = rc->get_IntFlag("RUNNUMBER");
    cout << PHWHERE << " Runnumber: " << runnumber_ << endl;
    ReadAlignmentParFile();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxAlign::process_event(PHCompositeNode *topNode)
{

  _timer.get()->restart();

  //------
  //--- get cluster container
  //------
  clusters_ = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
  if (!clusters_)
  {
    cout << PHWHERE << "ERROR: Can't find node TrkrClusterContainer" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //------
  //-- shift clusters
  //------
  TrkrDefUtil util;

  TrkrClusterContainer::ConstRange clusrange = clusters_->GetClusters();
  for ( TrkrClusterContainer::ConstIterator iter = clusrange.first;
        iter != clusrange.second;
        ++iter)
  {

    // get the key and check if we need to realign this cluster
    TrkrDefs::cluskey ckey = (iter->second)->GetClusKey();
    TrkrDefs::hitsetkey hkey = util.GetHitSetKeyFromClusKey(ckey);
    auto aligniter = alignmap_.find(hkey);

    if ( aligniter != alignmap_.end() )
    {
      // get a non-const pointer to the cluster
      TrkrClusterv1* clus = static_cast<TrkrClusterv1*>(clusters_->FindCluster(ckey));

      if ( verbosity > 1 )
      {
        cout << " applying alignment to " << endl;
        clus->identify();
      }

      // apply the alignment
      clus->SetX(clus->GetX() + (aligniter->second).dx);
      clus->SetY(clus->GetY() + (aligniter->second).dy);
      clus->SetZ(clus->GetZ() + (aligniter->second).dz);

      if ( verbosity > 1 )
        clus->identify();
    }

  }


  //------
  //--- done
  //------
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}


void MvtxAlign::AddAlignmentPar(TrkrDefs::hitsetkey key, double dx, double dy, double dz)
{

  // check that the misalignment doesn't already exist
  auto iter = alignmap_.find(key);
  if ( iter != alignmap_.end() )
  {
    cout << PHWHERE << " WARNING: overwriting existing misalignment for"
         << " key:0x" << hex << key << dec << endl;

    (iter->second).dx = dx;
    (iter->second).dy = dy;
    (iter->second).dz = dz;
  }
  else
  {
    AlignmentPar mis;
    mis.dx = dx;
    mis.dy = dy;
    mis.dz = dz;

    alignmap_.insert(make_pair(key, mis));

    if ( verbosity > 0 )
    {
      cout << PHWHERE << " Added alignment pars:"
           << " key:0x" << hex << key << dec
           << " dx:" << dx
           << " dy:" << dy
           << " dz:" << dz
           << endl;
    }
  }

}



void MvtxAlign::PrintAlignmentPars(std::ostream& os) const
{
  os << "=============================================================" << endl;
  os << "== " << PHWHERE << "==" << endl;
  os << "=============================================================" << endl;

  for ( auto iter = alignmap_.begin();
        iter != alignmap_.end();
        ++iter)
  {
    os << " key:0x" << hex << iter->first << dec
       << " dx:" << (iter->second).dx
       << " dy:" << (iter->second).dy
       << " dz:" << (iter->second).dz
       << endl;
  }

  os << "=============================================================" << endl;
}


int MvtxAlign::ReadAlignmentParFile()
{

  char fname[500];
  sprintf(fname, "%sbeamcenter_%08d.txt", fdir_.c_str(), runnumber_);
  ifstream fin;
  fin.open(fname);
  if ( fin.is_open() )
  {
    cout << PHWHERE << " Reading alignment parameters from " << fname << endl;
    MvtxDefUtil util;

    string line;
    static bool is_first = true;
    float bc0_x = 0;
    float bc0_y = 0;
    float bc0_z = 0;
    while ( getline(fin, line) )
    {
      int lyr, stave, chip;
      float bcx, bcy, bcz;

      sscanf(line.c_str(), "%d %d %d %f %f %f", &lyr, &stave, &chip, &bcx, &bcy, &bcz);

      if ( is_first )
      {
        bc0_x = bcx;
        bc0_y = bcy;
        bc0_z = bcz;
        is_first = false;
      }

      double dx = bc0_x - bcx;
      double dy = bc0_y - bcy;
      double dz = bc0_z - bcz;

      AddAlignmentPar(
        util.GenHitSetKey(lyr, stave, chip),
        double(dx), double(dy), double(dz));
    }
  }
  else
  {
    cout << PHWHERE << "ERROR: Unable to open file " << fname << endl;
    return 1;
  }

  if ( verbosity > 0 )
    PrintAlignmentPars();

  return 0;
}


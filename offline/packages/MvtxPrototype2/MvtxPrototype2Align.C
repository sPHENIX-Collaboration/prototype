#include "MvtxPrototype2Align.h"
#include "SegmentationAlpide.h"

#include <mvtx/MvtxDefs.h>

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

MvtxPrototype2Align::MvtxPrototype2Align(const string &name ) :
  SubsysReco(name),
  clusters_(NULL),
  fdir_(""),
  runnumber_(0),
  apff_(true),
  m_afname(""),
  m_is_global(true),
  _timer(PHTimeServer::get()->insert_new(name))
{

}

int MvtxPrototype2Align::InitRun(PHCompositeNode* topNode)
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

int MvtxPrototype2Align::process_event(PHCompositeNode *topNode)
{

  _timer.get()->restart();

  //------
  //--- get cluster container
  //------
  clusters_ = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusters_)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //------
  //-- shift clusters
  //------
  TrkrClusterContainer::ConstRange clusrange = clusters_->getClusters();
  for ( TrkrClusterContainer::ConstIterator iter = clusrange.first;
        iter != clusrange.second;
        ++iter)
  {

    // get the key and check if we need to realign this cluster
    TrkrDefs::cluskey ckey = (iter->second)->getClusKey();
    TrkrDefs::hitsetkey hkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
    auto aligniter = alignmap_.find(hkey);

    if ( aligniter != alignmap_.end() )
    {
      // get a non-const pointer to the cluster
      TrkrClusterv1* clus = static_cast<TrkrClusterv1*>(clusters_->findCluster(ckey));

      if ( Verbosity() > 1 )
      {
        cout << " applying alignment to " << endl;
        clus->identify();
      }

      // apply the alignment
      clus->setX(clus->getX() - (aligniter->second).dx);//row and xaxis are oppositive
      clus->setY(clus->getY() - (aligniter->second).dy);
      clus->setZ(clus->getZ() - (aligniter->second).dz);

      if ( Verbosity() > 1 )
        clus->identify();
    }

  }


  //------
  //--- done
  //------
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}


void MvtxPrototype2Align::AddAlignmentPar(TrkrDefs::hitsetkey key, double dx, double dy, double dz)
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

    if ( Verbosity() > 0 )
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



void MvtxPrototype2Align::PrintAlignmentPars(std::ostream& os) const
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


int MvtxPrototype2Align::ReadAlignmentParFile()
{
  string flname(fdir_);
  if ( m_afname.empty() )
  {
    char fname[100];
    sprintf(fname, "beamcenter_%08d.txt", runnumber_);
    flname += fname;
  }
  else
  {
    flname += m_afname;
  }
  ifstream fin;
  fin.open(flname);
  if ( fin.is_open() )
  {
    cout << PHWHERE << " Reading alignment parameters from " << flname << endl;
    //MvtxDefUtil util;

    string line;
    // static bool is_first = true;
    // float bc0_x = 0;
    // float bc0_y = 0;
    // float bc0_z = 0;
    while ( getline(fin, line) )
    {
      int lyr, stave;
      float bcx, bcy, bcz;

      sscanf(line.c_str(), "%d %d %f %f %f", &lyr, &stave, &bcx, &bcy, &bcz);

      // if ( is_first )
      // {
      //   bc0_x = bcx;
      //   bc0_y = bcy;
      //   bc0_z = bcz;
      //   is_first = false;
      // }

      // double dx = bc0_x - bcx;
      // double dy = bc0_y - bcy;
      // double dz = bc0_z - bcz;

      for ( int chip = 0; chip < 9; ++chip )
      {
        AddAlignmentPar(
          MvtxDefs::genHitSetKey(lyr, 0, chip),
          double(bcx), double(bcy), double(bcz));
      }
    }
  }
  else
  {
    cout << PHWHERE << "ERROR: Unable to open file " << flname << endl;
    return 1;
  }

  if ( Verbosity() > 0 )
    PrintAlignmentPars();

  return 0;
}

// $Id: $

/*!
 * \file TpcPrototypeTrack.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "TpcPrototypeTrack.h"

#include <cmath>
#include <limits>

using namespace std;

TpcPrototypeTrack::TpcPrototypeTrack()
  : event(0)
  , trackID(-1)
  , chisq(NAN)
  , ndf(0)
  , px(NAN)
  , py(NAN)
  , pz(NAN)
  , x(NAN)
  , y(NAN)
  , z(NAN)
  , nCluster(0)
{
  for (int i = 0; i < nLayer; ++i)
  {
    clusterKey[i] = numeric_limits<uint64_t>::max();
    clusterlayer[i] = -1;
    clusterid[i] = -1;
    clusterX[i] = numeric_limits<float>::signaling_NaN();
    clusterY[i] = numeric_limits<float>::signaling_NaN();
    clusterZ[i] = numeric_limits<float>::signaling_NaN();
    clusterE[i] = numeric_limits<float>::signaling_NaN();
    clusterSizePhi[i] = numeric_limits<float>::signaling_NaN();
    clusterProjectionPhi[i] = numeric_limits<float>::signaling_NaN();
    clusterResidualPhi[i] = clusterResidualZ[i] = numeric_limits<float>::signaling_NaN();
  }
}

TpcPrototypeTrack::~TpcPrototypeTrack()
{
  // TODO Auto-generated destructor stub
}

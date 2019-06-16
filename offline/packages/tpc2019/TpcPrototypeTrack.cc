// $Id: $

/*!
 * \file TpcPrototypeTrack.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "TpcPrototypeTrack.h"

#include <limits>
#include <iostream>

using namespace std;

TpcPrototypeTrack::TpcPrototypeTrack()
  : trackID(-1)
  , chisq(NAN)
  , ndf(0)
  , nCluster(0)
{
  for (int i = 0; i < nLayer; ++i)
  {
    clusterX[i] = numeric_limits<float>::signaling_NaN();
    clusterY[i] = numeric_limits<float>::signaling_NaN();
    clusterZ[i] = numeric_limits<float>::signaling_NaN();
    clusterE[i] = numeric_limits<float>::signaling_NaN();
    clusterSizePhi[i] = 0;
    clusterResidualIsolated[i] = numeric_limits<float>::signaling_NaN();
  }
}

TpcPrototypeTrack::~TpcPrototypeTrack()
{
  // TODO Auto-generated destructor stub
}

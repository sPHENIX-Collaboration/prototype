#include "MvtxPrototype2Geom.h"

#include <iostream>


std::unique_ptr<MvtxPrototype2Geom> MvtxPrototype2Geom::sInstance;

void MvtxPrototype2Geom::Build()
{
  if (isBuilt()) {
    std::cout << "Already built" << std::endl;
    return; // already initialized
  }

  mLastChipIndex.resize(mNumberOfLayers);
  int numberOfChips = 0;
  for (int iLay = 0; iLay < mNumberOfLayers; ++iLay){
    numberOfChips += mNumberOfStaves * mNumberOfChips;
    mLastChipIndex[iLay] = numberOfChips - 1;
  }
  setSize(numberOfChips);
  fillMatrixCache();
}


void MvtxPrototype2Geom::fillMatrixCache()
{
  if (mSize < 1) {
    std::cout << "Build was not called yet. Calling now..." << std::endl;
    Build();
  }
  mL2G.setSize(mSize);
  for (int i = 0; i < mSize; i++) {
    TGeoHMatrix* hm = extractMatrixSensor(i);
    mL2G.setMatrix(*hm, i);
  }

  return;
}


TGeoHMatrix* MvtxPrototype2Geom::extractMatrixSensor(int index)
{
  int lay        = index / getNumberOfChipsInLay();
  int indexInLay = index % getNumberOfChipsInLay();
  int stv        = indexInLay / mNumberOfChips;
  int indexInStv = indexInLay % mNumberOfChips;

  float shift_dz = (2 * SegmentationAlpide::PassiveEdgeSide) + \
                   SegmentationAlpide::ActiveMatrixSizeCols + PitchChip_IB;

  float dz = (indexInStv - 4) * shift_dz;
  float dy = -lay * GapStave_TB;

  if (Verbose()>0) {
    std::cout << "Filling matrix for sensor in chip " << indexInStv << "  in stave " << stv;
    std::cout << "  of layer " << lay << "\n";
    std::cout << "shift in z " << dz << " and in dy " << dy << std::endl;
  }

  static TGeoHMatrix matTmp;
  matTmp = TGeoTranslation(0., dy, dz);
  static TGeoTranslation tra(0., -.5 * Segmentation::SensLayerThickness, 0.);
  matTmp *= tra;

  return &matTmp;
}

#include "MvtxPrototype2Geom.h"

#include <iostream>


std::unique_ptr<MvtxPrototype2Geom> MvtxPrototype2Geom::s_instance;

void MvtxPrototype2Geom::Build()
{
  if (isBuilt()) {
    std::cout << "Already built" << std::endl;
    return; // already initialized
  }

  m_lastChipIndex.resize(m_numOfLayers);
  int num_of_chips = 0;
  for (int iLay = 0; iLay < m_numOfLayers; ++iLay){
    num_of_chips += m_numOfStaves * m_numOfChips;
    m_lastChipIndex[iLay] = num_of_chips - 1;
  }
  setSize(num_of_chips);
  fillMatrixCache();
}


void MvtxPrototype2Geom::fillMatrixCache()
{
  if (m_size < 1) {
    std::cout << "Build was not called yet. Calling now..." << std::endl;
    Build();
  }
  m_l2G.setSize(m_size);
  for (int i = 0; i < m_size; i++) {
    TGeoHMatrix* hm = extractMatrixSensor(i);
    m_l2G.setMatrix(*hm, i);
  }

  return;
}


TGeoHMatrix* MvtxPrototype2Geom::extractMatrixSensor(int index)
{
  int lay        = index / getNumberOfChipsInLay();
  int indexInLay = index % getNumberOfChipsInLay();
  int stv        = indexInLay / m_numOfChips;
  int indexInStv = indexInLay % m_numOfChips;

  float shift_dz = (2 * SegmentationAlpide::PassiveEdgeSide) + \
                   SegmentationAlpide::ActiveMatrixSizeCols + s_pitchChip_IB;

  float dz = (indexInStv - 4) * shift_dz;
  float dx = lay * s_gapLayers_TB;

  if (Verbose()>0) {
    std::cout << "Filling matrix for sensor in chip " << indexInStv << "  in stave " << stv;
    std::cout << "  of layer " << lay << "\n";
    std::cout << "shift in z " << dz << " and in dy " << dx << std::endl;
  }

  static TGeoHMatrix matTmp;
  matTmp = TGeoTranslation(dx, 0., dz);
  static TGeoTranslation tra(.5 * Segmentation::SensLayerThickness, 0., 0.);
  matTmp *= tra;

  return &matTmp;
}

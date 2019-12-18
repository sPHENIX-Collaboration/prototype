// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTX_TB_SETUP_STAVE_GEOM_H
#define MVTX_TB_SETUP_STAVE_GEOM_H

#define NMAXRU 2
#define NMAXRUCHN 28
#define NLAYER 4
#define NCHIP 9
#define NROW 512
#define NCOL 1024
#define NREGION 32
#define NCOL_PER_REGION 32

#include <trackbase/TrkrDefs.h>
#include <mvtx/MvtxDefs.h>
#include <mvtx/SegmentationAlpide.h>

#include <TGeoMatrix.h>

typedef SegmentationAlpide Segmentation;

template <typename T = TGeoHMatrix>
class MatrixCache
{
 public:
  MatrixCache()  = default;
  ~MatrixCache() = default;
  MatrixCache(const MatrixCache& src) = delete;
  MatrixCache& operator=(const MatrixCache src) = delete;

  ///set size of matrix cache
  void setSize(int s)
  {
    if (!m_cache.size())
      m_cache.resize(s);
  }

  /// get the size of the cache
  int getSize() const { return m_cache.size(); }
  /// assign matrix to a slot
  void setMatrix(const T& mat, int sensID)
  {
    // assign matrix for given sensor. The cache must be booked in advance
    //if ((unsigned int)sensID >= mCache.size()) {
    //  LOG(FATAL) << "SensID " << sensID << " exceeds cache size of " << mCache.size();
    //}
    m_cache[sensID] = mat;
  }

  const T& getMatrix(int sensID) const { return m_cache[sensID]; }
  bool isFilled() const { return !m_cache.empty(); }

 private:
  std::vector<T> m_cache;
};

/** @Brief Class to create MVTX testbeam setup geometry
 *
 */
class MvtxPrototype2Geom
{
 public:
  MvtxPrototype2Geom():m_verbose(0)
  {
    if (s_instance) {
      exit(1);
    }
    Build();
  }

  ~MvtxPrototype2Geom() = default;
  MvtxPrototype2Geom(const MvtxPrototype2Geom& src) = delete;
  MvtxPrototype2Geom& operator=(const MvtxPrototype2Geom& geo) = delete;

  const MatrixCache<TGeoHMatrix>& getCacheL2G() const { return m_l2G; }
  const TGeoHMatrix& getMatrixL2G(int sensID) const { return m_l2G.getMatrix(sensID); }
  bool isBuilt() const { return m_size != 0; }
  int getSize() const { return m_size; }
  // before calling fillMatrixCache, detector implementation should set the size of the matrix cache
  void setSize(int s)
  {
  // set the size of the matrix cache, can be done only once
  //if (mSize != 0) {
  //  LOG(FATAL) << "Cache size (N sensors) was already set to " << mSize;
  //}
    m_size = s;
  }
  int getLastChipIndex(int lay) const { return m_lastChipIndex[lay]; }
  int getFirstChipIndex(int lay) const { return (lay == 0) ? 0 : m_lastChipIndex[lay - 1] + 1; }
  int getChipIndex(int lay, int chipInStave)
  {
    return getFirstChipIndex(lay) + chipInStave;
  }
  int getChipIndex(TrkrDefs::hitsetkey _key) {
    return getChipIndex(TrkrDefs::getLayer(_key), MvtxDefs::getChipId(_key));
  }
  int getLayer(int index) const { return index / m_numOfChips; }
  int getChipInLay(int index) const { return index % m_numOfChips; }

  static MvtxPrototype2Geom* Instance()
  {
    // get (create if needed) a unique instance of the object
    if (!s_instance)
      s_instance = std::unique_ptr<MvtxPrototype2Geom>(new MvtxPrototype2Geom());
    return s_instance.get();
  }

  void Build();

  void setVerbose(int v) { m_verbose = v; }
  int Verbose() { return m_verbose; }

  public:
    static constexpr float s_pitchChip_IB = 100.e-4; ///< nominal chip pitch in IB stave = 100 um.
    static constexpr float s_gapLayers_TB = 4.;       ///< nominal gap between staves in 2019 TB setup = 4 cm.

  //  static bool detectorToStave(int chip, int iRow, int iCol, int& sRow, int& sCol);
  //  static bool staveToDetector(int sRow, int sCol, int& chip, int& iRow, int& iCol);
    bool detectorToGlobal(int index, int iCol, int iRow, double* glo);
    bool detectorToGlobal(TrkrDefs::hitsetkey, TrkrDefs::hitkey, double*);

 private:
  static int getNumberOfChipsInLay()
  {
    return m_numOfStaves * m_numOfChips;
  }
  TGeoHMatrix* extractMatrixSensor(int index);
  void fillMatrixCache();

 private:
  int m_verbose;                     ///< verbosity level
  static std::unique_ptr<MvtxPrototype2Geom> s_instance;
  int m_size = 0;                    ///< prebooked number of sensors
  MatrixCache<TGeoHMatrix> m_l2G;    ///< Local to Global matrices

  static constexpr int m_numOfLayers = 4; ///< Number of layers
  static constexpr int m_numOfStaves = 1; ///< Number of staves/layer
  static constexpr int m_numOfChips  = 9; ///< Number of Chips/stave

  std::vector<int> m_lastChipIndex;

};

//_________________________________________________________________________________________________
/* inline bool MvtxPrototype2Geom::detectorToStave(int chip, int iRow, int iCol, int& sRow, int& sCol)
{
  //convert to concecutive row/col in a stave from row/col in chip
  if (iRow<0 || iRow>=Segmentation::NRows || iCol<0 || iCol>=Segmentation::NCols)
    return false;
  sRow = iRow;
  sCol = chip * Segmentation::NCols + iCol;
  return true;
}

//_________________________________________________________________________________________________
inline bool MvtxPrototype2Geom::staveToDetector(int sRow, int sCol, int& chip, int& iRow, int& iCol)
{
  //convert to chip row/col from consecutive stave from row/col in chip
  if (sRow<0 || sRow>=Segmentation::NRows || sCol<0 || sCol>=(m_numOfChips * Segmentation::NCols))
    return false;
  chip = sCol / Segmentation::NCols;
  iCol = sCol % Segmentation::NCols;
  iRow = sRow;
  return true;
}
*/
//_________________________________________________________________________________________________
inline bool MvtxPrototype2Geom::detectorToGlobal(int index, int iCol, int iRow, double* glo)
{
  if (iRow<0 || iRow>=Segmentation::NRows || iCol<0 || iCol>=Segmentation::NCols)
    return false;
  Point3D<float> locPoint;
  Segmentation::detectorToLocal(iRow, iCol, locPoint);
  const TGeoHMatrix& matL2G = Instance()->getMatrixL2G(index);
  //matL2G.Print();
  double loc[3];
  locPoint.GetCoordinates(loc);
  matL2G.LocalToMaster(loc, glo);
  return true;
}

inline bool MvtxPrototype2Geom::detectorToGlobal(TrkrDefs::hitsetkey _hitsetkey, TrkrDefs::hitkey _hitkey, double* glo)
{
  return detectorToGlobal(getChipIndex(_hitsetkey), MvtxDefs::getCol(_hitkey), MvtxDefs::getRow(_hitkey), glo);
}

#endif

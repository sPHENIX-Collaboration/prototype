#ifndef MVTX_TB_SETUP_STAVE_GEOM_H_
#define MVTX_TB_SETUP_STAVE_GEOM_H_

#define NMAXRU 2
#define NMAXRUCHN 28
#define NSTAVE 4
#define NCHIP 9
#define NROW 512
#define NCOL 1024
#define NREGION 32
#define NCOL_PER_REGION 32

#include "SegmentationAlpide.h"

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

  //set size of matrix cache
  void setSize(int s)
  {
    if (!mCache.size())
      mCache.resize(s);
  }

  /// get the size of the cache
  int getSize() const { return mCache.size(); }
  /// assign matrix to a slot
  void setMatrix(const T& mat, int sensID)
  {
    // assign matrix for given sensor. The cache must be booked in advance
    //if ((unsigned int)sensID >= mCache.size()) {
    //  LOG(FATAL) << "SensID " << sensID << " exceeds cache size of " << mCache.size();
    //}
    mCache[sensID] = mat;
  }

  const T& getMatrix(int sensID) const { return mCache[sensID]; }
  bool isFilled() const { return !mCache.empty(); }

 private:
  std::vector<T> mCache;
};

class MvtxPrototype2Geom
{
 public:
  MvtxPrototype2Geom():mVerbose(0)
  {
    if (sInstance) {
      exit(1);
    }
    Build();
  }

  ~MvtxPrototype2Geom() = default;
  MvtxPrototype2Geom(const MvtxPrototype2Geom& src) = delete;
  MvtxPrototype2Geom& operator=(const MvtxPrototype2Geom& geo) = delete;

  const MatrixCache<TGeoHMatrix>& getCacheL2G() const { return mL2G; }
  const TGeoHMatrix& getMatrixL2G(int sensID) const { return mL2G.getMatrix(sensID); }
  bool isBuilt() const { return mSize != 0; }
  int getSize() const { return mSize; }
  // before calling fillMatrixCache, detector implementation should set the size of the matrix cache
  void setSize(int s)
  {
  // set the size of the matrix cache, can be done only once
  //if (mSize != 0) {
  //  LOG(FATAL) << "Cache size (N sensors) was already set to " << mSize;
  //}
    mSize = s;
  }
  int getLastChipIndex(int lay) const { return mLastChipIndex[lay]; }
  int getFirstChipIndex(int lay) const { return (lay == 0) ? 0 : mLastChipIndex[lay - 1] + 1; }
  int getChipIndex(int lay, int chipInStv){
    return getFirstChipIndex(lay) + chipInStv;
  }
  int getLayer(int index) const { return index / mNumberOfChips; }
  int getChipInLay(int index) const { return index % mNumberOfChips; }

  static MvtxPrototype2Geom* Instance()
  {
    // get (create if needed) a unique instance of the object
    if (!sInstance)
      sInstance = std::unique_ptr<MvtxPrototype2Geom>(new MvtxPrototype2Geom());
    return sInstance.get();
  }
  void Build();
  void setVerbose(int v) { mVerbose = v; }
  int Verbose() { return mVerbose; }

  public:
    static constexpr float PitchChip_IB = 100.e-4; //nominal chip pitch in IB stave = 100 um.
    static constexpr float GapStave_TB = 4.; // nominal gap between staves in 2019 TB setup = 4 cm.

    static bool detectorToStave(int chip, int iRow, int iCol, int& sRow, int& sCol);
    static bool staveToDetector(int sRow, int sCol, int& chip, int& iRow, int& iCol);
    static bool detectorToGlobal(int index, int iCol, int iRow, double* glo);

 private:
  static int getNumberOfChipsInLay(){
    return mNumberOfStaves * mNumberOfChips;
  }
  TGeoHMatrix* extractMatrixSensor(int index);
  void fillMatrixCache();

 private:
  int mVerbose;
  static std::unique_ptr<MvtxPrototype2Geom> sInstance;
  int mSize = 0;               ///< prebooked number of sensors
  MatrixCache<TGeoHMatrix> mL2G;    ///< Local to Global matrices

  static constexpr int mNumberOfLayers = 4; //< Number of layers
  static constexpr int mNumberOfStaves = 1; //< Number of staves/layer
  static constexpr int mNumberOfChips  = 9; //< Number of Chips/stave

  std::vector<int> mLastChipIndex;

};

//_________________________________________________________________________________________________
inline bool MvtxPrototype2Geom::detectorToStave(int chip, int iRow, int iCol, int& sRow, int& sCol)
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
  if (sRow<0 || sRow>=Segmentation::NRows || sCol<0 || sCol>=(mNumberOfChips * Segmentation::NCols))
    return false;
  chip = sCol / Segmentation::NCols;
  iCol = sCol % Segmentation::NCols;
  iRow = sRow;
  return true;
}

//_________________________________________________________________________________________________
inline bool MvtxPrototype2Geom::detectorToGlobal(int index, int iCol, int iRow, double* glo)
{
  if (iRow<0 || iRow>=Segmentation::NRows || iCol<0 || iCol>=Segmentation::NCols)
    return false;
  Point3D<float> locPoint;
  Segmentation::detectorToLocal(iCol, iRow, locPoint);
  const TGeoHMatrix& matL2G = Instance()->getMatrixL2G(index);
  //matL2G.Print();
  double loc[3];
  locPoint.GetCoordinates(loc);
  matL2G.LocalToMaster(loc, glo);
  return true;
}

#endif

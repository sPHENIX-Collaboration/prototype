#ifndef __MvtxStandaloneTracking_h__
#define __MvtxStandaloneTracking_h__

#include <vector>
#include <set>

#include <TMatrixD.h>
#include <TVectorD.h>

class TrkrCluster;
class TrkrClusterContainer;
class PHCompositeNode;

class MvtxStandaloneTracking
{
public:

  struct MvtxTrack
  {
    std::vector<TrkrCluster*> ClusterList;
    std::vector<double> dx;
    std::vector<double> dz;
    double m_xy;
    double b_xy;
    double chi2_xy;
    double m_zy;
    double b_zy;
    double chi2_zy;
  };
  typedef std::vector<MvtxTrack> MvtxTrackList;

  //! ctor
  MvtxStandaloneTracking();

  //! dtor
  ~MvtxStandaloneTracking();

  //! run tracking
  void RunTracking(PHCompositeNode* topNode, MvtxTrackList &trklst, std::vector<int> &lyrs);

  //! reject ghosts
  void RunGhostRejection( MvtxTrackList &trklst);

  //! set association window
  void SetWindowX(float w) { window_x_ = w; }
  void SetWindowZ(float w) { window_z_ = w; }

  //! set ghost rejection
  void SetGhostRejection( bool yn ) {}

  //! set verbosity
  void Verbosity(int v ) { verbosity_ = v; }

private:

  //! Associate clusters into track candidates
  void AssociateClusters(MvtxTrackList &trklst, std::vector<int> &lyrs);

  //! Fit track in x vs y
  void TrackFitXY(MvtxTrack &trk);

  //! Fit track in z vs y
  void TrackFitZY(MvtxTrack &trk);

  //! Generalized least squares fitter
  TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L);


  double CalcSlope(double x0, double y0, double x1, double y1);
  double CalcIntecept(double x0, double y0, double m);
  double CalcProjection(double x, double m, double b);

  //! Print out track candidate information
  void PrintTrackCandidates(MvtxTrackList &trklst);

  //! Cluster node
  TrkrClusterContainer* clusters_;

  //! window size in x & z
  float window_x_;
  float window_z_;

  //! ghost rejection
  bool ghostrejection_;
  
  // verbosity
  int verbosity_;

  static const int NLYR = 4;

};

#endif //__MvtxStandaloneTracking_h__#
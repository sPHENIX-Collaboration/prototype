#ifndef __MvtxQAHisto_H__
#define __MvtxQAHisto_H__

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TFile;
class TH1F;
class TH2F;

//Brief: basic ntuple and histogram creation for sim evaluation
class MvtxQAHisto: public SubsysReco
{
 public: 
  //Default constructor
  MvtxQAHisto(const std::string &name="MvtxQAHisto");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile = file; }

 private:
  //output filename
  std::string _outfile;
   
  //Event counter
  int _event;

  //Get all the nodes
  void GetNodes(PHCompositeNode *);
  
  //flags
  unsigned int _flags;

	float beam_x[4];
	float beam_y[4];

	TH2F *h2d_hit[4];
	TH2F *h2d_hit_beam[4];
	TH2F *h2d_hit_trk[4];
	TH1F *h1d_hit_per_evt[4];

	TH2F *h2d_clus[4];
	TH2F *h2d_clus_beam[4];
	TH1F *h1d_clus_per_evt[4]; //cluster in a chip per event 
	TH1F *h1d_clus_size_x[4]; 
	TH1F *h1d_clus_size_z[4]; 

	TrkrClusterContainer *cluscon;
	TrkrHitSetContainer *hitsetcon;

	TH1F *h1d_trk_finder_x;
	TH1F *h1d_trk_finder_z;
	TH2F *h2d_trk_finder;

	TH1F *h1d_clus_associated;
	TH1F *h1d_clus_eff;


};

#endif //* __MvtxQAHisto_H__ *//

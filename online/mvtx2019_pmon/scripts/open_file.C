#if !defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__)

#include "mvtxOM.h"

R__LOAD_LIBRARY(libmvtxOM)

#endif

void open_file(int runnum, int nevents = 10001, const char* ftype = "longrun")
{
  gSystem->Load("libmvtxOM.so");
  //set_verbose(1);
  char filein[500];
  //int runnum=atoi(filename.Data());
  //sprintf(filein,"../beamtest2019/%s/%s_%08d-0000.prdf",ftype,ftype,runnum);
  //sprintf(filein,"/mnt/databkup_1/MVTX_testbeam2019/%s/%s_%08d-0000.prdf",ftype,ftype,runnum);
  sprintf(filein,"/sphenix/data/data01/mvtx/%s/%s_%08d-0000.prdf",ftype,ftype,runnum);
  pfileopen(filein);
  prun(nevents);
  OM();
  return;
}

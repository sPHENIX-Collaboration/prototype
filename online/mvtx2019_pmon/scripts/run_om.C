#if !defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__)

#include "mvtxOM.h"

R__LOAD_LIBRARY(libmvtxOM)

#endif

void run_om()
{
  gSystem->Load("libmvtxOM.so");
  //set_verbose(1);
  set_refresh(5000);
  rcdaqopen();
  OM();
  pstart();
  return;
}

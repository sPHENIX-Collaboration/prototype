//#ifdef __cline__
#include "mvtx.h"
//R__LOAD_LIBRARY(libpmonitor)
R__LOAD_LIBRARY(libmvtx)
//#endif

void mvtxBeamTestRunScan(const char *inputFile)
{
  pfileopen(inputFile);
  //OM();
  //rcdaqopen();
  //pstatus();
  //pstart();
  prun();
  OM();
  //sleep(15);
  print_canvas();
  get_alignment();
  exit(0);
}

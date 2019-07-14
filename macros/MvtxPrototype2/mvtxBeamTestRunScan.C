//#ifdef __cline__
#include "mvtxOM.h"
//R__LOAD_LIBRARY(libpmonitor)
R__LOAD_LIBRARY(libmvtxOM)
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

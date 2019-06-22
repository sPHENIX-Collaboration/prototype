#ifndef __MVTX_H__
#define __MVTX_H__

#include <pmonitor/pmonitor.h>   // Added by X. He following Martin's trick
#include <Event/Event.h>
#include <Event/EventTypes.h>

int process_event (Event *e); //++CINT
int process_histos(float thresh = 10);
int mask_pixels(float thresh = 5);
int analysis();
int OM();
int print_canvas();
int print_status();
void reset_histos();
void set_verbose(int v);
void set_refresh(int r);
void get_alignment();

#endif /* __MVTX_H__ */

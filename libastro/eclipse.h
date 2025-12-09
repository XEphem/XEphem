
#ifndef ECLIPSE_H
#define ECLIPSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "astro.h"

#define ADAY (1.0)
#define ANHOUR (ADAY / 24.0)
#define AMINUTE (ADAY / 24.0 / 60.0)
#define ASECOND (ADAY / 24.0 / 60.0 / 60.0)
#define AYEAR (ADAY * 365.2425)


void init_eclipse_now( Now * pn, double tick );
void init_eclipse_obj( Obj * po, PLCode c );
void init_eclipse_now_sun_moon( Now * pn, Obj * psun, Obj * pmoon, double tick );
void update_eclipse_now_sun_moon( Now * pn, Obj * psun, Obj * pmoon, double tick );
bool calculate_decD( Obj * psun, Obj * pmoon, double * pdecD, double * pr0, double * pr1 );
bool calculate_skyD( Obj * psun, Obj * pmoon, double decD, double * pskyD, double r0, double r1 );
int is_eclipsing( Obj * psun, Obj * pmoon );
void increment_eclipse_time( Now * pn, Obj * psun, Obj * pmoon, double offset );
bool find_eclipse_start_mid_stop( Now * pn, Obj * psun, Obj * pmoon, double * pstart, double * pmid, double * pstop );
bool get_eclipse_path_location( Now * pn, Obj * psun, Obj * pmoon, double * plt, double * plg );


#ifdef __cplusplus
}
#endif

#endif

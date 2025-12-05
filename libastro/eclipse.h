
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

void init_modifiable_instances( Now * np, Obj * po0, Obj * po1, Now * pdn, Obj *pdo0, Obj * pdo1 );
bool calculate_decD( Obj * po0, Obj * po1, double * pdecD, double * pr0, double * pr1 );
bool calculate_skyD( Obj * po0, Obj * po1, double decD, double * pskyD, double r0, double r1 );
int is_eclipsing( Obj * pobj1, Obj * pobj2 );
void increment_eclipse_time( Now * np, Obj * po0, Obj * po1, double offset );
bool find_eclipse_start_mid_stop( Now * np, Obj * po0, Obj * po1, double * pstart, double * pmid, double * pstop );
bool get_eclipse_path_location( Now * np, Obj * po0, Obj * po1, double * plt, double * plg );

#ifdef __cplusplus
}
#endif

#endif

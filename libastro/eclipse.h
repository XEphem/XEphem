
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
#define ASAROS (ADAY * 6585.3211)
#define ANINEX (ADAY * 10571.95)

#define ECLIPSESCANTIMERANGE (ANHOUR * 8)
#define ECLIPSESCANTIMEINCREMENT (ASECOND * 30)

#define MAXSAROSECLIPSES 100

typedef struct {
	double tickstart;
	double tickmid;
	double tickstop;
	double latstart;
	double lonstart;
	double latmid;
	double lonmid;
	double latstop;
	double lonstop;
} Eclipse;

typedef struct {
	int sarosnumber;
	Eclipse * pe[MAXSAROSECLIPSES];
	size_t eclipsecount;
} Saros;

typedef struct {
	int sarosnumbermin;
	int sarosnumbermax;
	size_t maxinexsaros;
	Saros ** ps;
	size_t saroscount;
} Inex;

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

Eclipse * malloc_eclipse_data( double tickstart, double tickmid, double tickstop, double latstart, double lonstart, double latmid, double lonmid, double latstop, double lonstop );
void free_eclipse_data( Eclipse * pe );
Saros * malloc_saros_data( int sarosnumber );
void free_saros_data( Saros * ps );
Inex * malloc_inex_data( int sarosnumbermin, int sarosnumbermax );
void free_inex_data( Inex * pi );

int compare_saros_eclipses( const void * pa, const void * pb );
int compare_saros_eclipses_reverse( const void * pa, const void * pb );
void sort_saros_data( Saros * ps );
void sort_saros_data_reverse( Saros * ps );
int compare_inex_saros( const void * pa, const void * pb );
int compare_inex_saros_reverse( const void * pa, const void * pb );
void sort_inex_data( Inex * pi );
void sort_inex_data_reverse( Inex * pi );

bool scan_for_eclipse( double * ptick );
bool scan_for_inex_eclipse( double * ptick, Saros * ps, double offset );
bool add_eclipse_to_saros_data( Saros * ps, double tick );
bool add_saros_to_inex_data( Inex * pi, int sarosnumber, double tick );
Inex * create_inex( double tickseed, int sarosnumberseed, int sarosnumbermin, int sarosnumbermax );


#ifdef __cplusplus
}
#endif

#endif

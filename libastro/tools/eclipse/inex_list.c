
#include <stdlib.h>
#include <math.h>

#include "libastro/astro.h"
#include "libastro/eclipse.h"


bool format_mjd_date( char * pd, size_t dmax, double tick ) {
	bool retval = false;
	if( pd ) {
		int m = 0;
		double d = 0.0;
		int y = 0;
		mjd_cal( tick, &m, &d, &y );
		if( snprintf( pd, dmax, "%5i.%02i.%02i", y, m, (int)d ) > 0 ) {
			retval = true;
		}
	}
	return retval;
}

bool format_mjd_time( char * pt, size_t tmax, double tick ) {
	bool retval = false;
	if( pt ) {
		int m = 0;
		double d = 0.0;
		int y = 0;
		int hr = 0;
		int mn = 0;
		int sc = 0;
		double t = 0.0;
		mjd_cal( tick, &m, &d, &y );
		t = d - floor( d );
		hr = t / ANHOUR;
		mn = (t - (hr * ANHOUR)) / AMINUTE;
		sc = (t - (hr * ANHOUR) - (mn * AMINUTE )) / ASECOND;
		if( snprintf( pt, tmax, "%02i:%02i:%02i", hr, mn, sc ) > 0 ) {
			retval = true;
		}
	}
	return retval;
}

bool format_lat( char * plt, size_t ltmax, double lt ) {
	bool retval = false;
	if( plt && snprintf( plt, ltmax, "%f", raddeg( lt ) ) > 0 ) {
		retval = true;
	}
	return retval;
}

bool format_lon( char * plg, size_t lgmax, double lg ) {
	bool retval = false;
	if( plg && snprintf( plg, lgmax, "%f", raddeg( lg ) ) > 0 ) {
		retval = true;
	}
	return retval;
}

bool format_mjd_date_time( char * pstate, size_t statemax, double tick ) {
	bool retval = false;
	if( pstate && statemax > 30 ) {
		char datebuf[16] = "";
		char timebuf[12] = "";
		if( format_mjd_date( datebuf, sizeof( datebuf ) - 1, tick ) ) {
			if( format_mjd_time( timebuf, sizeof( timebuf ) - 1, tick ) ) {
				if( snprintf( pstate, statemax, "%s %s", datebuf, timebuf ) > 0 ) {
					retval = true;
				}
			}
		}
	}
	return retval;
}

bool format_mjd_date_time_lat_lon( char * pstate, size_t statemax, double tick, double lt, double lg ) {
	bool retval = false;
	if( pstate && statemax > 60 ) {
		char datebuf[16] = "";
		char timebuf[12] = "";
		char latbuf[16] = "";
		char lonbuf[16] = "";
		if( format_mjd_date( datebuf, sizeof( datebuf ) - 1, tick ) ) {
			if( format_mjd_time( timebuf, sizeof( timebuf ) - 1, tick ) ) {
				if( format_lat( latbuf, sizeof( latbuf ) - 1, lt ) ) {
					if( format_lat( lonbuf, sizeof( lonbuf ) - 1, lg ) ) {
						if( snprintf( pstate, statemax, "%s %s %s %s", datebuf, timebuf, latbuf, lonbuf ) > 0 ) {
							retval = true;
						}
					}
				}
			}
		}
	}
	return retval;
}

void print_eclipse( Eclipse * pe ) {
	if( pe ) {
		char bufstart[32] = "";
		char bufmid[32] = "";
		char bufstop[32] = "";
		format_mjd_date_time( bufstart, sizeof( bufstart ), pe->tickstart );
		format_mjd_date_time( bufmid, sizeof( bufmid ), pe->tickmid );
		format_mjd_date_time( bufstop, sizeof( bufstop ), pe->tickstop );
		printf( "ECLIPSE %20s  %20s  %20s  %9.5f %10.5f  %9.5f %10.5f  %9.5f %10.5f\n", bufstart, bufmid, bufstop, raddeg( pe->latstart ), raddeg( pe->lonstart ), raddeg( pe->latmid ), raddeg( pe->lonmid ), raddeg( pe->latstop ), raddeg( pe->lonstop ) );
	}
}

void print_saros( Saros * ps ) {
	if( ps ) {
		for( size_t loop = 0; loop < ps->eclipsecount; ++loop ) {
			printf( "SAROS %03i ", ps->sarosnumber );
			print_eclipse( ps->pe[loop] );
		}
	}
}

void print_inex( Inex * pi ) {
	if( pi ) {
		for( size_t loop = 0; loop < pi->saroscount; ++loop ) {
			if( loop != 0 ) printf( "\n" );
			print_saros( pi->ps[loop] );
		}
	}
}

int main( int argc, char * argv[] ) {
	double tick = 36382; /* 8/11/1999 */
	int sarosnumber = 145; /* The 8/11/1999 eclipse part of Saros series 145. */
	int sarosnumbermin = 0;
	int sarosnumbermax = 180;
	if( argc == 3 || argc == 5 ) {
		sarosnumbermin = atoi( argv[1] );
		sarosnumbermax = atoi( argv[2] );
		if( argc > 3 ) tick = atof( argv[3] );
		if( argc > 4 ) sarosnumber = atoi( argv[4] );
		Inex * pi = create_inex( tick, sarosnumber, sarosnumbermin, sarosnumbermax );
		if( pi ) {
			print_inex( pi );
			free_inex_data( pi );
		}
	} else {
		printf( "Usage: %s sarosnumbermin sarosnumbermax [eclipsingMJD sarosnumberofeclipsingMJD]\n", argv[0] );
		printf( "\t Example: %s 0 180 (lists Saros series 0 thru 180)\n", argv[0] );
		printf( "\t Example: %s 144 146 (lists Saros series 144, 145, 146)\n", argv[0] );
		printf( "\t Example: %s 144 146 36382.0 145 (lists Saros series 144, 145, 146)\n", argv[0] );
		printf( "\t Example: %s 144 146 36382.0 42 (unconventional Saros series numbering)\n", argv[0] );
	}
	return 0;
}

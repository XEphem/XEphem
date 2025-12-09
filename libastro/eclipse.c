
#include <math.h>
#include <string.h>

#include "eclipse.h"


/* See also: GUI/xephem/tools/libastro_sample.c */

void init_eclipse_now( Now * pn, double tick ) {
	if( pn ) {
		memset( pn, 0, sizeof( Now ) );
		pn->n_mjd = tick;
		/* BUG:
		 * Changes to pn->n_lat and pn->n_lng alter eclipse path results.
		 * Why do LAT/LON changes impact eclipse path calculations?
		 * Do other changes impact path calculations?
		*/
	} else {
		fprintf( stderr, "ERROR: NULL passed to init_eclipse_now().\n" );
	}
}

void init_eclipse_obj( Obj * po, PLCode c ) {
	if( po ) {
		memset( po, 0, sizeof( Obj ) );
		po->o_type = PLANET;
		po->pl_code = c;
		// if( c == SUN ) strncpy( po->o_name, "Sun", MAXNM );
		// if( c == MOON ) strncpy( po->o_name, "Moon", MAXNM );
	} else {
		fprintf( stderr, "ERROR: NULL passed to init_eclipse_obj().\n" );
	}
}

void init_eclipse_now_sun_moon( Now * pn, Obj * psun, Obj * pmoon, double tick ) {
	if( pn ) init_eclipse_now( pn, tick );
	if( psun ) {
		init_eclipse_obj( psun, SUN );
		if( pn ) obj_cir( pn, psun );
	}
	if( pmoon ) {
		init_eclipse_obj( pmoon, MOON );
		if( pn ) obj_cir( pn, pmoon );
	}
}

void update_eclipse_now_sun_moon( Now * pn, Obj * psun, Obj * pmoon, double tick ) {
	if( pn ) {
		pn->n_mjd = tick;
		if( psun ) obj_cir( pn, psun );
		if( pmoon ) obj_cir( pn, pmoon );
	}
}

/* Extracted from XEphem e_soleclipse(), originally by Downey. */
bool calculate_decD( Obj * psun, Obj * pmoon, double * pdecD, double * pr0, double * pr1 ) {
	bool retval = false;
	if( psun && pmoon && pdecD && is_ssobj( psun ) && is_ssobj( pmoon ) ) {
		double r0 = psun->s_edist * (MAU / ERAD);
		double r1 = pmoon->s_edist * (MAU / ERAD);
		double decA = pmoon->s_gaedec - psun->s_gaedec;
		*pdecD = r0 * r1 * sin( decA ) / (r0 - r1);
		if( pr0 ) *pr0 = r0;
		if( pr1 ) *pr1 = r1;
		if( fabs( *pdecD ) < 1.0 ) {
			retval = true;
		}
	}
	return retval;
}

/* Extracted from XEphem e_soleclipse(), originally by Downey. */
bool calculate_skyD( Obj * psun, Obj * pmoon, double decD, double * pskyD, double r0, double r1 ) {
	bool retval = false;
	if( psun && pmoon && pskyD && is_ssobj( psun ) && is_ssobj( pmoon ) ) {
		double skyA = acos(
			sin( psun->s_gaedec ) * sin( pmoon->s_gaedec )
			+ cos( psun->s_gaedec ) * cos( pmoon->s_gaedec )
			* cos( psun->s_gaera - pmoon->s_gaera )
		);
		*pskyD = r0 * r1 * sin( skyA ) / (r0 - r1);
		if( fabs( *pskyD ) < 1.0 ) {
			retval = true;
		}
	}
	return retval;
}

int is_eclipsing( Obj * psun, Obj * pmoon ) {
	double decD = 0.0;
	double skyD = 0.0;
	double r0 = 0.0;
	double r1 = 0.0;
	bool validdecd = calculate_decD( psun, pmoon, &decD, &r0, &r1 );
	bool validskyd = calculate_skyD( psun, pmoon, decD, &skyD, r0, r1 );
	return (validdecd && validskyd);
}

void increment_eclipse_time( Now * pn, Obj * psun, Obj * pmoon, double offset ) {
	if( pn && psun && pmoon ) {
		double tick = pn->n_mjd;
		while( is_eclipsing( psun, pmoon ) ) {
			tick = pn->n_mjd;
			update_eclipse_now_sun_moon( pn, psun, pmoon, pn->n_mjd + offset );
		}
		update_eclipse_now_sun_moon( pn, psun, pmoon, tick );
	}
}

bool find_eclipse_start_mid_stop( Now * pn, Obj * psun, Obj * pmoon, double * pstart, double * pmid, double * pstop ) {
	bool retval = false;
	if( pn && psun && pmoon ) {
		init_eclipse_now_sun_moon( pn, psun, pmoon, pn->n_mjd );
		if( pstart ) {
			increment_eclipse_time( pn, psun, pmoon, -ANHOUR );
			increment_eclipse_time( pn, psun, pmoon, -AMINUTE );
			increment_eclipse_time( pn, psun, pmoon, -ASECOND );
			*pstart = pn->n_mjd + ASECOND;
			retval = true;
		}
		if( pstop ) {
			update_eclipse_now_sun_moon( pn, psun, pmoon, pn->n_mjd );
			increment_eclipse_time( pn, psun, pmoon, ANHOUR );
			increment_eclipse_time( pn, psun, pmoon, AMINUTE );
			increment_eclipse_time( pn, psun, pmoon, ASECOND );
			*pstop = pn->n_mjd - ASECOND;
			retval = true;
		}
		if( pmid && retval && pstart && pstop ) {
			*pmid = (*pstart + ((*pstop - *pstart) / 2.0));
		}
	}
	return retval;
}

/* Extracted from XEphem e_soleclipse(), originally by Downey. */
bool get_eclipse_path_location( Now * pn, Obj * psun, Obj * pmoon, double * plt, double * plg ) {
	bool retval = false;
	if( pn && psun && pmoon ) {
		double decD = 0.0;
		double skyD = 0.0;
		double r0 = 0.0;
		double r1 = 0.0;
		if( calculate_decD( psun, pmoon, &decD, &r0, &r1 ) ) {
			if( calculate_skyD( psun, pmoon, decD, &skyD, r0, r1 ) ) {
				double theta = acos( decD / skyD );
				double skyP = atan( skyD / r0 );

				/* FIX: when asin() goes nan, kludge with arbitrary 1.3 ? */
				/*
				double skyT = asin( skyD * r0 / (sqrt( r0 * r0 ) + sqrt( skyD * skyD )) ) - skyP;
				if( isnan( skyT ) ) {
					skyT = 1.3;
				}
				*/
				double skyT = asin( skyD * r0 / sqrt( r0 * r0 + skyD * skyD ) ) - skyP;

				if( ! isnan( theta ) && ! isnan( skyP ) && ! isnan( skyT ) ) {
					double sD = 0.0;
					double dRA = 0.0;
					double lt = 0.0;
					double lg = 0.0;
					double lst = 0.0;
					double gst = 0.0;
					double diffgaera = fabs( pmoon->s_gaera - psun->s_gaera );
					solve_sphere( theta, skyT, sin( psun->s_gaedec ), cos( psun->s_gaedec ), &sD, &dRA );

					/*
					* FIX: slightly arbitrary 4.0 kludge, to avoid eclipse path LAT/LON bugs
					* where the eclipse path location is occassionally incorrect
					* maybe because the Moon crosses the equator ?
					* e.g. 3/20/2034 8:40 UTC good, 11:57 UTC bad
					* https://eclipse.gsfc.nasa.gov/5MCSEmap/2001-2100/2034-03-20.gif
					*/
/*
					if (pmoon->s_gaera > psun->s_gaera)
					dRA = -dRA;
*/
					if (pmoon->s_gaera > psun->s_gaera) {
						if (diffgaera < 4.0) {
							dRA = -dRA;	/* eastward */
						}
					} else if( diffgaera > 4.0 ) {
						dRA = -dRA;	/* eastward */
					}

					lt = asin( sD );
					lst = psun->s_gaera - dRA;
					utc_gst( mjd_day( pn->n_mjd ), mjd_hr( pn->n_mjd ), &gst );
					lg = lst - hrrad( gst );
					while( lg < -PI ) lg += 2 * PI;
					while( lg >  PI ) lg -= 2 * PI;
					if( ! isnan( lt ) && ! isnan( lg ) ) {
						if( plt ) *plt = lt;
						if( plg ) *plg = lg;
						retval = true;
					}
				}
			}
		}
	}
	return retval;
}

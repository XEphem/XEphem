
#include <math.h>
#include <string.h>

#include "eclipse.h"

/*
 * Is there a better way to initialize/acquire
 * _modifiable_ instances (without disrupting XEphem state) ?
 */
void init_modifiable_instances( Now * np, Obj * po0, Obj * po1, Now * pdn, Obj *pdo0, Obj * pdo1 ) {
	if( np && pdn ) memcpy( (void *)pdn, (void *)np, sizeof( Now ) );
	if( po0 && pdo0 ) {
		memcpy( (void *)pdo0, (void *)po0, sizeof( Obj ) );
		if( np ) obj_cir( np, pdo0 );
	}
	if( po1 && pdo1 ) {
		memcpy( (void *)pdo1, (void *)po1, sizeof( Obj ) );
		if( np ) obj_cir( np, pdo1 );
	}
}

/* Extracted from XEphem e_soleclipse(), originally by Downey. */
bool calculate_decD( Obj * po0, Obj * po1, double * pdecD, double * pr0, double * pr1 ) {
	bool retval = false;
	if( po0 && po1 && pdecD && is_ssobj( po0 ) && is_ssobj( po1 ) ) {
		double r0 = po0->s_edist * (MAU / ERAD);
		double r1 = po1->s_edist * (MAU / ERAD);
		double decA = po1->s_gaedec - po0->s_gaedec;
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
bool calculate_skyD( Obj * po0, Obj * po1, double decD, double * pskyD, double r0, double r1 ) {
	bool retval = false;
	if( po0 && po1 && pskyD && is_ssobj( po0 ) && is_ssobj( po1 ) ) {
		double skyA = acos(
			sin( po0->s_gaedec ) * sin( po1->s_gaedec )
			+ cos( po0->s_gaedec ) * cos( po1->s_gaedec )
			* cos( po0->s_gaera - po1->s_gaera )
		);
		*pskyD = r0 * r1 * sin( skyA ) / (r0 - r1);
		if( fabs( *pskyD ) < 1.0 ) {
			retval = true;
		}
	}
	return retval;
}

int is_eclipsing( Obj * po0, Obj * po1 ) {
	double decD = 0.0;
	double skyD = 0.0;
	double r0 = 0.0;
	double r1 = 0.0;
	bool validdecd = calculate_decD( po0, po1, &decD, &r0, &r1 );
	bool validskyd = calculate_skyD( po0, po1, decD, &skyD, r0, r1 );
	return (validdecd && validskyd);
}

void increment_eclipse_time( Now * np, Obj * po0, Obj * po1, double offset ) {
	if( np && po0 && po1 ) {
		double tick = np->n_mjd;
		while( is_eclipsing( po0, po1 ) ) {
			tick = np->n_mjd;
			np->n_mjd += offset;
			obj_cir( np, po0 );
			obj_cir( np, po1 );
		}
		np->n_mjd = tick;
		obj_cir( np, po0 );
		obj_cir( np, po1 );
	}
}

bool find_eclipse_start_mid_stop( Now * np, Obj * po0, Obj * po1, double * pstart, double * pmid, double * pstop ) {
	bool retval = false;
	Now enow;
	Obj eo0;
	Obj eo1;
	init_modifiable_instances( np, po0, po1, &enow, &eo0, &eo1 );
	if( pstart ) {
		enow.n_mjd = np->n_mjd;
		obj_cir( &enow, &eo0 );
		obj_cir( &enow, &eo1 );
		increment_eclipse_time( &enow, &eo0, &eo1, -ANHOUR );
		increment_eclipse_time( &enow, &eo0, &eo1, -AMINUTE );
		increment_eclipse_time( &enow, &eo0, &eo1, -ASECOND );
		*pstart = enow.n_mjd + ASECOND;
		retval = true;
	}
	if( pstop ) {
		enow.n_mjd = np->n_mjd;
		obj_cir( &enow, &eo0 );
		obj_cir( &enow, &eo1 );
		increment_eclipse_time( &enow, &eo0, &eo1, ANHOUR );
		increment_eclipse_time( &enow, &eo0, &eo1, AMINUTE );
		increment_eclipse_time( &enow, &eo0, &eo1, ASECOND );
		*pstop = enow.n_mjd - ASECOND;
		retval = true;
	}
	if( pmid && retval && pstart && pstop ) {
		*pmid = (*pstart + ((*pstop - *pstart) / 2.0));
	}
	return retval;
}

/* Extracted from XEphem e_soleclipse(), originally by Downey. */
bool get_eclipse_path_location( Now * np, Obj * po0, Obj * po1, double * plt, double * plg ) {
	bool retval = false;
	if( np && po0 && po1 ) {
		double decD = 0.0;
		double skyD = 0.0;
		double r0 = 0.0;
		double r1 = 0.0;
		if( calculate_decD( po0, po1, &decD, &r0, &r1 ) ) {
			if( calculate_skyD( po0, po1, decD, &skyD, r0, r1 ) ) {
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
					double diffgaera = fabs( po1->s_gaera - po0->s_gaera );
					solve_sphere( theta, skyT, sin( po0->s_gaedec ), cos( po0->s_gaedec ), &sD, &dRA );

					/*
					* FIX: slightly arbitrary 4.0 kludge, to avoid eclipse path LAT/LON bugs
					* where the eclipse path location is occassionally incorrect
					* maybe because the Moon crosses the equator ?
					* e.g. 3/20/2034 8:40 UTC good, 11:57 UTC bad
					* https://eclipse.gsfc.nasa.gov/5MCSEmap/2001-2100/2034-03-20.gif
					*/
					/*
					* if (eo1.s_gaera > eo0.s_gaera)
					* dRA = -dRA;
					*/
					if (po1->s_gaera > po0->s_gaera) {
						if (diffgaera < 4.0) {
							dRA = -dRA;	/* eastward */
						}
					} else if( diffgaera > 4.0 ) {
						dRA = -dRA;	/* eastward */
					}

					lt = asin( sD );
					lst = po0->s_gaera - dRA;
					utc_gst( mjd_day( np->n_mjd ), mjd_hr( np->n_mjd ), &gst );
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

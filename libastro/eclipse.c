
#include <stdlib.h>
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
/*
		if( c == SUN ) strncpy( po->o_name, "Sun", MAXNM );
		if( c == MOON ) strncpy( po->o_name, "Moon", MAXNM );
*/
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

Eclipse * malloc_eclipse_data( double tickstart, double tickmid, double tickstop, double latstart, double lonstart, double latmid, double lonmid, double latstop, double lonstop ) {
	Eclipse * retval = malloc( sizeof( Eclipse ) );
	if( retval ) {
		retval->tickstart = tickstart;
		retval->tickmid = tickmid;
		retval->tickstop = tickstop;
		retval->latstart = latstart;
		retval->lonstart = lonstart;
		retval->latmid = latmid;
		retval->lonmid = lonmid;
		retval->latstop = latstop;
		retval->lonstop = lonstop;
	}
	return retval;
}

void free_eclipse_data( Eclipse * pe ) {
	if( pe ) {
		free( pe );
		pe = NULL;
	}
}

Saros * malloc_saros_data( int sarosnumber ) {
	Saros * retval = malloc( sizeof( Saros ) );
	if( retval ) {
		retval->sarosnumber = sarosnumber;
		for( size_t loop = 0; loop < MAXSAROSECLIPSES; ++loop ) {
			retval->pe[loop] = NULL;
		}
		retval->eclipsecount = 0;
	}
	return retval;
}

void free_saros_data( Saros * ps ) {
	if( ps ) {
		for( size_t loop = 0; loop < ps->eclipsecount; ++loop ) {
			free_eclipse_data( ps->pe[loop] );
			ps->pe[loop] = NULL;
		}
		free( ps );
		ps = NULL;
	}
}

Inex * malloc_inex_data( int sarosnumbermin, int sarosnumbermax ) {
	Inex * retval = NULL;
	Inex * pi = malloc( sizeof( Inex ) );
	if( pi ) {
		pi->sarosnumbermin = sarosnumbermin;
		pi->sarosnumbermax = sarosnumbermax;
		/* The + 1 is to allow for negative saros numbers. */
		pi->maxinexsaros = sarosnumbermax - sarosnumbermin + 1;
		pi->ps = malloc( sizeof( Saros * ) * pi->maxinexsaros );
		if( pi->ps ) {
			for( size_t loop = 0; loop < pi->maxinexsaros; ++loop ) {
				pi->ps[loop] = NULL;
			}
			pi->saroscount = 0;
			retval = pi;
		} else {
			free( pi );
			pi = NULL;
		}
	}
	return retval;
}

void free_inex_data( Inex * pi ) {
	if( pi ) {
		if( pi->ps ) {
			for( size_t loop = 0; loop < pi->saroscount; ++loop ) {
				free_saros_data( pi->ps[loop] );
				pi->ps[loop] = NULL;
			}
			free( pi->ps );
			pi->ps = NULL;
		}
		free( pi );
		pi = NULL;
	}
}

int compare_saros_eclipses( const void * pa, const void * pb ) {
	const Eclipse * pea = *(const Eclipse **)pa;
	const Eclipse * peb = *(const Eclipse **)pb;
	return pea->tickmid < peb->tickmid;
}

int compare_saros_eclipses_reverse( const void * pa, const void * pb ) {
	return compare_saros_eclipses( pb, pa );
}

void sort_saros_data( Saros * ps ) {
	if( ps ) {
		qsort( ps->pe, ps->eclipsecount, sizeof( Eclipse * ), compare_saros_eclipses );
	}
}

void sort_saros_data_reverse( Saros * ps ) {
	if( ps ) {
		qsort( ps->pe, ps->eclipsecount, sizeof( Eclipse * ), compare_saros_eclipses_reverse );
	}
}

int compare_inex_saros( const void * pa, const void * pb ) {
	const Saros * psa = *(const Saros **)pa;
	const Saros * psb = *(const Saros **)pb;
	return psa->sarosnumber < psb->sarosnumber;
}

int compare_inex_saros_reverse( const void * pa, const void * pb ) {
	return compare_inex_saros( pb, pa );
}

void sort_inex_data( Inex * pi ) {
	if( pi ) {
		qsort( pi->ps, pi->saroscount, sizeof( Saros * ), compare_inex_saros );
	}
}

void sort_inex_data_reverse( Inex * pi ) {
	if( pi ) {
		qsort( pi->ps, pi->saroscount, sizeof( Saros * ), compare_inex_saros_reverse );
	}
}

/*
 * ptick is a pointer to a MJD of interest.
 * If an eclipse is found *ptick is updated to an eclipsing MJD.
 * If generally useful, RANGE and INCREMENT could become function arguments.
 * Returns true if a total solar eclipse is found near the MJD.
*/
bool scan_for_eclipse( double * ptick ) {
	bool retval = false;
	if( ptick ) {
		Now enow;
		Obj esun;
		Obj emoon;
		init_eclipse_now_sun_moon( &enow, &esun, &emoon, *ptick );
		if( is_eclipsing( &esun, &emoon ) ) {
			retval = true;
		} else {
			double tickminusmin = *ptick - ECLIPSESCANTIMERANGE;
			double tickplusmax = *ptick + ECLIPSESCANTIMERANGE;
			double tickminus = *ptick;
			double tickplus = *ptick;
			double increment = ECLIPSESCANTIMEINCREMENT;
			while( ! retval ) {
				if( tickminus >= tickminusmin ) {
					update_eclipse_now_sun_moon( &enow, &esun, &emoon, tickminus );
					if( is_eclipsing( &esun, &emoon ) ) {
						*ptick = tickminus;
						retval = true;
						break;
					}
					tickminus -= increment;
				}
				if( tickplus <= tickplusmax ) {
					update_eclipse_now_sun_moon( &enow, &esun, &emoon, tickplus );
					if( is_eclipsing( &esun, &emoon ) ) {
						*ptick = tickplus;
						retval = true;
						break;
					}
					tickplus += increment;
				}
				if( tickminus < tickminusmin && tickplus > tickplusmax ) {
					break;
				}
			}
		}
	}
	return retval;
}

bool scan_for_inex_eclipse( double * ptick, Saros * ps, double offset ) {
	bool retval = false;
	if( ptick && ps && ps->eclipsecount > 0 ) {
		/*
		 * Middle out search.
		 * Because Saros middle eclipses are good Inex eclipse candidates.
		*/
		size_t indexminusmin = 0;
		size_t indexplusmax = ps->eclipsecount;
		size_t indexminus = indexplusmax / 2;
		size_t indexplus = indexplusmax / 2;
		double tick = 0.0;
		while( ! retval ) {
			if( indexminus >= indexminusmin ) {
				tick = ps->pe[indexminus]->tickmid + offset;
				if( scan_for_eclipse( &tick ) ) {
					*ptick = tick;
					retval = true;
					break;
				}
				indexminus -= 1;
			}
			if( indexplus < indexplusmax ) {
				tick = ps->pe[indexplus]->tickmid + offset;
				if( scan_for_eclipse( &tick ) ) {
					*ptick = tick;
					retval = true;
					break;
				}
				indexplus += 1;
			}
			if( indexminus < indexminusmin && indexplus >= indexplusmax ) {
				break;
			}
		}
	}
	return retval;
}

bool add_eclipse_to_saros_data( Saros * ps, double tick ) {
	bool retval = false;
	if( ps ) {
		if( ps->eclipsecount < MAXSAROSECLIPSES ) {
			Now enow;
			Obj esun;
			Obj emoon;
			double tickstart = 0.0;
			double tickmid = 0.0;
			double tickstop = 0.0;
			double latstart = 0.0;
			double lonstart = 0.0;
			double latmid = 0.0;
			double lonmid = 0.0;
			double latstop = 0.0;
			double lonstop = 0.0;
			init_eclipse_now_sun_moon( &enow, &esun, &emoon, tick );
			if( find_eclipse_start_mid_stop( &enow, &esun, &emoon, &tickstart, &tickmid, &tickstop ) ) {
				update_eclipse_now_sun_moon( &enow, &esun, &emoon, tickstart );
				get_eclipse_path_location( &enow, &esun, &emoon, &latstart, &lonstart );
				update_eclipse_now_sun_moon( &enow, &esun, &emoon, tickmid );
				get_eclipse_path_location( &enow, &esun, &emoon, &latmid, &lonmid );
				update_eclipse_now_sun_moon( &enow, &esun, &emoon, tickstop );
				get_eclipse_path_location( &enow, &esun, &emoon, &latstop, &lonstop );
				Eclipse * pe = malloc_eclipse_data( tickstart, tickmid, tickstop, latstart, lonstart, latmid, lonmid, latstop, lonstop );
				if( pe ) {
					ps->pe[ps->eclipsecount] = pe;
					ps->eclipsecount += 1;
					retval = true;
				}
			}
		} else {
			fprintf( stderr, "ERROR: exceeded default (%i) Saros Eclipse capacity.\n", MAXSAROSECLIPSES );
		}
	}
	return retval;
}

bool add_saros_to_inex_data( Inex * pi, int sarosnumber, double tick ) {
	bool retval = false;
	if( pi ) {
		if( pi->saroscount < pi->maxinexsaros ) {
			Saros * ps = malloc_saros_data( sarosnumber );
			if( ps ) {
				double tock = tick;
				if( scan_for_eclipse( &tock ) ) {
					if( add_eclipse_to_saros_data( ps, tock ) ) {
						tock = tick - ASAROS;
						while( scan_for_eclipse( &tock ) ) {
							if( add_eclipse_to_saros_data( ps, tock ) ) {
								tock -= ASAROS;
							} else {
								break;
							}
						}
						tock = tick + ASAROS;
						while( scan_for_eclipse( &tock ) ) {
							if( add_eclipse_to_saros_data( ps, tock ) ) {
								tock += ASAROS;
							} else {
								break;
							}
						}
						sort_saros_data_reverse( ps );
						pi->ps[pi->saroscount] = ps;
						pi->saroscount += 1;
						retval = true;
					} else {
						free_saros_data( ps );
					}
				} else {
					free_saros_data( ps );
				}
			}
		} else {
			fprintf( stderr, "ERROR: exceeded default (%zu) Inex Saros capacity.\n", pi->maxinexsaros );
		}
	}
	return retval;
}

/*
 * Allocates and returns a Inex->Saros->Eclipse data structure.
 * Call free_inex_data() to free the allocated memory.
 * tickseed is a MJD which should be during an eclipse.
 * sarosnumberseed is a Saros number which should correlate with the tickseed eclipse.
 * sarosnumbermin is the first Saros series of interest.
 * sarosnumbermax is the last Saros series of interest.
 * For example, the 8/11/1999 eclipse (36382) is part of Saros series 145.
 * Returns an allocated Inex->Saros->Eclipse data structure.
*/
Inex * create_inex( double tickseed, int sarosnumberseed, int sarosnumbermin, int sarosnumbermax ) {
	Inex * retval = NULL;
	if( sarosnumberseed >= sarosnumbermin && sarosnumberseed <= sarosnumbermax ) {
		Inex * pi = malloc_inex_data( sarosnumbermin, sarosnumbermax );
		if( pi ) {
			if( add_saros_to_inex_data( pi, sarosnumberseed, tickseed ) ) {
				for( int loop = sarosnumberseed - 1; loop >= pi->sarosnumbermin; --loop ) {
					double tick = 0.0;
					if( scan_for_inex_eclipse( &tick, pi->ps[pi->saroscount - 1], -ANINEX ) ) {
						if( ! add_saros_to_inex_data( pi, loop, tick ) ) {
							break;
						}
					}
				}
				sort_inex_data_reverse( pi );
				for( int loop = sarosnumberseed + 1; loop <= pi->sarosnumbermax; ++loop ) {
					double tick = 0.0;
					if( scan_for_inex_eclipse( &tick, pi->ps[pi->saroscount - 1], ANINEX ) ) {
						if( ! add_saros_to_inex_data( pi, loop, tick ) ) {
							break;
						}
					}
				}
				sort_inex_data_reverse( pi );
				retval = pi;
			}
		}
	}
	return retval;
}

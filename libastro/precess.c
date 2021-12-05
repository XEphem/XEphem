#include <stdio.h>
#include <math.h>

#include "astro.h"

static void precess_hiprec (double mjd1, double mjd2, double *ra, double *dec);

#define	DCOS(x)		cos(degrad(x))
#define	DSIN(x)		sin(degrad(x))
#define	DASIN(x)	raddeg(asin(x))
#define	DATAN2(y,x)	raddeg(atan2((y),(x)))

/* corrects ra and dec, both in radians, for precession from epoch 1 to epoch 2.
 * the epochs are given by their modified JDs, mjd1 and mjd2, respectively.
 * N.B. ra and dec are modifed IN PLACE.
 */
void
precess (
double mjd1, double mjd2,	/* initial and final epoch modified JDs */
double *ra, double *dec)	/* ra/dec for mjd1 in, for mjd2 out */
{
	precess_hiprec (mjd1, mjd2, ra, dec);
}
/*
 *
 * Precession formulas zeta_A, z_A and theta_A
 *	NB : all 3 functions get_zeta_A(), get_z_A(), get_theta_A()
 *	return arcseconds
 */
static double
do_precess (
double equinox, double *alpha, double *delta, int dir
)
{
	double T;
	double zeta_A, z_A, theta_A;
	double A, B, C;

	/* precession progresses about 1 arc second in .047 years */
	if (fabs (equinox - 2000.0) > .02) {
	    T = (equinox - 2000.0)/100.0;

	    /* Compute coefficients in arcseconds from */
	    /*	Astronomical Ephemeris 2020, p. B25 */
	    /* and convert them to degrees */

	    zeta_A= (2.650545 + T*(2306.083227+T*(0.2988499 + T*(0.01801828 +T*(-5.971e-6-3.173e-7*T)))))/3600.00 ;
	    z_A=(-2.650545+T*(2306.077181+T*(1.0927348+T*(0.01826837 +T*(-28.596e-6-2.904e-7*T)))))/3600.00 ;
	    theta_A=(T*(2004.191903-T*(0.4294934+T*(0.04182264+T*(7.089e-6+1.274e-7*T)))))/3600.00 ;

	    if (dir<0){
	    	/* If dir is negative, go from equinox to 2000.0 */
	    	A = DSIN(*alpha - z_A) * DCOS(*delta);
	    	B = DCOS(*alpha - z_A) * DCOS(theta_A) * DCOS(*delta)
	    	  + DSIN(theta_A) * DSIN(*delta);
	    	C = -DCOS(*alpha - z_A) * DSIN(theta_A) * DCOS(*delta)
	    	  + DCOS(theta_A) * DSIN(*delta);
	    	*alpha = DATAN2(A,B) - zeta_A;
	    }

	    if (dir>0){
	    	/* If dir is positive, go from 2000.0 to equinox */
	    	A = DSIN(*alpha + zeta_A) * DCOS(*delta);
	    	B = DCOS(*alpha + zeta_A) * DCOS(theta_A) * DCOS(*delta)
	    	    - DSIN(theta_A) * DSIN(*delta);
	    	C = DCOS(*alpha + zeta_A) * DSIN(theta_A) * DCOS(*delta)
	    	    + DCOS(theta_A) * DSIN(*delta);

	    	*alpha = DATAN2(A,B) + z_A;
	    }
	    range (alpha, 360.0);
	    *delta = DASIN(C);
	}
}

/*
 * Copyright (c) 1990 by Craig Counterman. All rights reserved.
 *
 * This software may be redistributed freely, not sold.
 * This copyright notice and disclaimer of warranty must remain
 *    unchanged. 
 *
 * No representation is made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty, to the extent permitted by applicable law.
 *
 */
static void
precess_hiprec (
double mjd1, double mjd2,	/* initial and final epoch modified JDs */
double *ra, double *dec)	/* ra/dec for mjd1 in, for mjd2 out */
{
	static double last_mjd1 = -213.432, last_from;
	static double last_mjd2 = -213.432, last_to;
	double alpha, delta;
	double from_equinox, to_equinox;

	/* convert mjds to years;
	 * avoid the remarkably expensive calls to mjd_year()
	 */
	if (last_mjd1 == mjd1)
	    from_equinox = last_from;
	else {
	    mjd_year (mjd1, &from_equinox);
	    last_mjd1 = mjd1;
	    last_from = from_equinox;
	}
	if (last_mjd2 == mjd2)
	    to_equinox = last_to;
	else {
	    mjd_year (mjd2, &to_equinox);
	    last_mjd2 = mjd2;
	    last_to = to_equinox;
	}
	/* convert coords from rads to degs */
	alpha = raddeg(*ra);
	delta = raddeg(*dec);

	/* From from_equinox to 2000.0, "backwards" inplace transformation : */
	do_precess(from_equinox, &alpha, &delta, -1);
	/* From 2000.0 to to_equinox, "forward" inplace transformation : */
	do_precess(to_equinox, &alpha, &delta, 1);

	/* convert result coords from degs to rads */
	*ra = degrad(alpha);
	*dec = degrad(delta);
}

#if 0
static void
precess_fast (
double mjd1, double mjd2,	/* initial and final epoch modified JDs */
double *ra, double *dec)	/* ra/dec for mjd1 in, for mjd2 out */
{
#define	N	degrad (20.0468/3600.0)
#define	M	hrrad (3.07234/3600.0)
	double nyrs;

	nyrs = (mjd2 - mjd1)/365.2425;
	*dec += N * cos(*ra) * nyrs;
	*ra += (M + (N * sin(*ra) * tan(*dec))) * nyrs;
	range (ra, 2.0*PI);
}
#endif


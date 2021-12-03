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
get_zeta_A (
double T)
{
	return 2.650545 + T*(2306.083227+T*(0.2988499 + T*(0.01801828 +T*(-5.971e-6-3.173e-7*T)))) ;
}

static double
get_z_A (
double T)
{
	return -2.650545+T*(2306.077181+T*(1.0927348+T*(0.01826837 +T*(-28.596e-6-2.904e-7*T)))) ;
}

static double
get_theta_A (
double T)
{
	return T*(2004.191903-T*(0.4294934+T*(0.04182264+T*(7.089e-6+1.274e-7*T)))) ;
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
 * Rigorous precession. From Astronomical Ephemeris 1989, p. B18
 *
 * 96-06-20 Hayo Hase <hase@wettzell.ifag.de>: theta_a corrected
 */
static void
precess_hiprec (
double mjd1, double mjd2,	/* initial and final epoch modified JDs */
double *ra, double *dec)	/* ra/dec for mjd1 in, for mjd2 out */
{
	static double last_mjd1 = -213.432, last_from;
	static double last_mjd2 = -213.432, last_to;
	double zeta_A, z_A, theta_A;
	double T;
	double A, B, C;
	double alpha, delta;
	double alpha_in, delta_in;
	double from_equinox, to_equinox;
	double alpha2000, delta2000;

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

	/* convert coords in rads to degs */
	alpha_in = raddeg(*ra);
	delta_in = raddeg(*dec);

	/* precession progresses about 1 arc second in .047 years */
	/* From from_equinox to 2000.0 */
	if (fabs (from_equinox-2000.0) > .02) {
	    T = (from_equinox - 2000.0)/100.0;
	    /* Get coefficients in arcseconds and convert them to degrees */
	    zeta_A  = get_zeta_A(T)/3600.0;
	    z_A  = get_z_A(T)/3600.0;
	    theta_A  = get_theta_A(T)/3600.0;

	    A = DSIN(alpha_in - z_A) * DCOS(delta_in);
	    B = DCOS(alpha_in - z_A) * DCOS(theta_A) * DCOS(delta_in)
	      + DSIN(theta_A) * DSIN(delta_in);
	    C = -DCOS(alpha_in - z_A) * DSIN(theta_A) * DCOS(delta_in)
	      + DCOS(theta_A) * DSIN(delta_in);

	    alpha2000 = DATAN2(A,B) - zeta_A;
	    range (&alpha2000, 360.0);
	    delta2000 = DASIN(C);
	} else {
	    /* should get the same answer, but this could improve accruacy */
	    alpha2000 = alpha_in;
	    delta2000 = delta_in;
	};


	/* From 2000.0 to to_equinox */
	if (fabs (to_equinox - 2000.0) > .02) {
	    T = (to_equinox - 2000.0)/100.0;
	    /* Get coefficients in arcseconds and convert them to degrees */
	    zeta_A  = get_zeta_A(T)/3600.0;
	    z_A  = get_z_A(T)/3600.0;
	    theta_A  = get_theta_A(T)/3600.0;

	    A = DSIN(alpha2000 + zeta_A) * DCOS(delta2000);
	    B = DCOS(alpha2000 + zeta_A) * DCOS(theta_A) * DCOS(delta2000)
	      - DSIN(theta_A) * DSIN(delta2000);
	    C = DCOS(alpha2000 + zeta_A) * DSIN(theta_A) * DCOS(delta2000)
	      + DCOS(theta_A) * DSIN(delta2000);

	    alpha = DATAN2(A,B) + z_A;
	    range(&alpha, 360.0);
	    delta = DASIN(C);
	} else {
	    /* should get the same answer, but this could improve accruacy */
	    alpha = alpha2000;
	    delta = delta2000;
	};

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


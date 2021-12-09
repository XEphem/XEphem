#include <stdio.h>
#include <math.h>

#include "astro.h"

/*
 *
 * Precession computation to or from J2000 depending on parameter dir :
 *	if dir < 0 go from equinox to J2000
 *	if dir > 0 go from J2000 to equinox
 */
static double
do_precess (
double equinox, double *alpha, double *delta, int dir
)
{
	double T;
	double zeta_A, z_A, theta_A;
	double A, B, C;
	/* precession progresses about 1 arc second in .047 years
	 * that is 17 days
	 * threshold must be at most 1.7 days for 0.1 arcsec accuracy */

	double threshhold=1;

	if (fabs (equinox - J2000) > threshhold) {
	    T = (equinox - J2000)/36525.0;

	    /* Compute coefficients in arcseconds from */
	    /*	Astronomical Ephemeris 2020, p. B25 */
	    /* and convert them to radians */

	    zeta_A= (2.650545 + T*(2306.083227+T*(0.2988499 + T*(0.01801828 +T*(-5.971e-6-3.173e-7*T)))))*arcsecrad;
	    z_A= (-2.650545+T*(2306.077181+T*(1.0927348+T*(0.01826837 +T*(-28.596e-6-2.904e-7*T)))))*arcsecrad;
	    theta_A=(T*(2004.191903-T*(0.4294934+T*(0.04182264+T*(7.089e-6+1.274e-7*T)))))*arcsecrad;

	    if (dir<0){
	    	/* If dir is negative, go from equinox to 2000.0 */
	    	A = sin(*alpha - z_A) * cos(*delta);
	    	B = cos(*alpha - z_A) * cos(theta_A) * cos(*delta)
	    	  + sin(theta_A) * sin(*delta);
	    	C = -cos(*alpha - z_A) * sin(theta_A) * cos(*delta)
	    	  + cos(theta_A) * sin(*delta);
	    	*alpha = atan2(A,B) - zeta_A;
	    }

	    if (dir>0){
	    	/* If dir is positive, go from 2000.0 to equinox */
	    	A = sin(*alpha + zeta_A) * cos(*delta);
	    	B = cos(*alpha + zeta_A) * cos(theta_A) * cos(*delta)
	    	    - sin(theta_A) * sin(*delta);
	    	C = cos(*alpha + zeta_A) * sin(theta_A) * cos(*delta)
	    	    + cos(theta_A) * sin(*delta);

	    	*alpha = atan2(A,B) + z_A;
	    }
	    range (alpha, TWOPI);
	    *delta = asin(C);
	}
}

/* corrects ra and dec, both in radians, for precession from epoch 1 to epoch 2.
 * the epochs are given by their modified JDs, from_mjd and to_mjd, respectively.
 * N.B. ra and dec are modifed IN PLACE.
 */
void
precess (
double from_mjd, double to_mjd,	/* initial and final epoch modified JDs */
double *ra, double *dec)	/* ra/dec for mjd1 in, for mjd2 out */
{
	/* From mjd1  to J2000, "backwards" inplace transformation : */
	do_precess(from_mjd, ra, dec, -1);
	/* From J2000 to mjd2, "forward" inplace transformation : */
	do_precess(to_mjd, ra, dec, 1);
}


/* VSOP87 planetary theory
 *
 * currently uses version VSOP87D:
 * heliocentric spherical, mean ecliptic of date.
 *
 * calculation of rates (daily changes) is optional;
 * see header file for the necessary #define's
 *
 * rough orientation on calculation time, miliseconds
 * on an HP 715/75, all planets Mercury to Neptune, prec=0.0:
 *
 *      terms	with rates	without rates
 * 	3598	11		7.1
 *      31577	51		44
 *
 * with secular terms for JD 2232395.0  19/12/1399 0h TDB:
 *
 *	FULL PRECISION code (31577 terms), milliseconds
 *	prec	terms	rates	no rates
 *	1e-8	15086	62	36
 *	1e-7	10105	44	25
 *	1e-6	3725	20	13
 *	1e-5	1324	11	7.8
 *	1e-4	443	7.0	6.0
 *	1e-3	139	6.0	5.0
 *
 *	REDUCED PRECISION code (3598 terms), milliseconds
 *	prec	terms	rates	no rates
 *	1e-7	2463	9.9	5.5
 *	1e-6	1939	8.0	4.5
 *	1e-5	1131	4.9	2.9
 *	1e-4	443	2.2	1.5
 *	1e-3	139	1.0	0.9
 */

#include <math.h>

#include "astro.h"
#include "vsop87.h"

#define VSOP_A1000	365250.0	/* days per millenium */
#define VSOP_MAXALPHA	5		/* max degree of time */

/******************************************************************
 * adapted from BdL FORTRAN Code; stern
 *
 *    Reference : Bureau des Longitudes - PBGF9502
 *
 *    Object :  calculate a VSOP87 position for a given time.
 *
 *    Input :
 *
 *    mj       modified julian date, counted from J1900.0
 *             time scale : dynamical time TDB.
 *
 *    obj	object number as in astro.h, NB: not for pluto
 *
 *    prec     relative precision
 *
 *             if prec is equal to 0 then the precision is the precision
 *                p0 of the complete solution VSOP87.
 *                Mercury    p0 =  0.6 10**-8
 *                Venus      p0 =  2.5 10**-8
 *                Earth      p0 =  2.5 10**-8
 *                Mars       p0 = 10.0 10**-8
 *                Jupiter    p0 = 35.0 10**-8
 *                Saturn     p0 = 70.0 10**-8
 *                Uranus     p0 =  8.0 10**-8
 *                Neptune    p0 = 42.0 10**-8
 *
 *             if prec is not equal to 0, let us say in between p0 and
 *             10**-3, the precision is :
 *                for the positions :
 *                - prec*a0 au for the distances.
 *                - prec rad for the other variables.
 *                for the velocities :
 *                - prec*a0 au/day for the distances.
 *                - prec rad/day for the other variables.
 *                  a0 is the semi-major axis of the body.
 *
 *    Output :
 *
 *    ret[6]     array of the results (double).
 *
 *             for spherical coordinates :
 *                 1: longitude (rd)
 *                 2: latitude (rd)
 *                 3: radius (au)
 *		#if VSOP_GETRATE:
 *                 4: longitude velocity (rad/day)
 *                 5: latitude velocity (rad/day)
 *                 6: radius velocity (au/day)
 *
 *    return:     error index (int)
 *                 0: no error.
 *		   2: object out of range [MERCURY .. NEPTUNE, SUN]
 ******************************************************************/
/******************************************************************
 * added entire vsop87 data set (31577 terms)
 * change use of indexes for pointers
 * delete VSOP_SCALE macro
 * delete precision control as always called with 0 (full precision)
 * Gustavo A. Corradi - Dec - 2022 (gcorrad@gmail.com)
 ******************************************************************/

extern t_vsop87v 
    *mercury_vs[], *venus_vs[], *mars_vs[], *jupiter_vs[],
    *saturn_vs[], *uranus_vs[], *neptune_vs[], *earth_vs[];

int
vsop87 (double mj, int obj, double *ret)
{
    static t_vsop87v **map_obj[] = {	
        mercury_vs, venus_vs, mars_vs, jupiter_vs,
	    saturn_vs, uranus_vs, neptune_vs, 0, earth_vs,
    };
    double t[VSOP_MAXALPHA+1];			/* powers of time */
    double term, arg, mj0;
#if VSOP_GETRATE
    double termdot;
#endif
    t_vsop87 *p; /* VSOP87 term pointer */
    t_vsop87v *pp, **vsp; /* VSOP87 variable and planet pointer */
    int i, cooidx, alpha, maxt; /* misc indexes */

    if (obj == PLUTO || obj > SUN)
	    return (2);

    vsp = map_obj[obj];
    /* zero result array */
    for (i = 0; i < 6; ++i) ret[i] = 0.0;

    /* time and its powers */
    mj0 = (mj - J2000)/VSOP_A1000;
    t[0] = 1.0;
    t[1] = mj0;
    for (i = 2; i <= VSOP_MAXALPHA; ++i) t[i] = t[i-1] * mj0;

    /* do the term summation; first the spatial dimensions */
    for (cooidx = 0; cooidx < 3; ++cooidx) {

        /* then the powers of time */
        pp = vsp[cooidx];
	    for (alpha = 0; pp->vars; ++alpha, ++pp) {

            maxt = pp->maxt;
            p = pp->vars;
            term = 0.0;
#if VSOP_GETRATE
            termdot = 0.0;  
#endif
            for (i = 0; i < maxt; ++i, ++p) {

		        arg = p->B + p->r * mj0;
		        term += p->L * cos(arg);
#if VSOP_GETRATE
		        termdot += -p->r * p->L * sin(arg);
#endif
	    }

	    ret[cooidx] += t[alpha] * term;
#if VSOP_GETRATE
	    ret[cooidx + 3] += t[alpha] * termdot +
		    ((alpha > 0) ? alpha * t[alpha - 1] * term : 0.0);
#endif
	} /* alpha */
    } /* cooidx */

#if VSOP_SPHERICAL
    /* reduce longitude to 0..2pi */
    ret[0] -= floor(ret[0]/(2.*PI)) * (2.*PI);
#endif

#if VSOP_GETRATE
    /* convert millenium rate to day rate */
    for (i = 3; i < 6; ++i) ret[i] /= VSOP_A1000;
#endif

#if VSOP_SPHERICAL
    /* reduction from dynamical equinox of VSOP87 to FK5;
     */
    {
	double L1, c1, s1;
	L1 = ret[0] - degrad(13.97 * t[1] - 0.031 * t[2]);
	c1 = cos(L1); s1 = sin(L1);
	ret[0] += degrad(-0.09033 + 0.03916 * (c1 + s1) * tan(ret[1]))/3600.0;
	ret[1] += degrad(0.03916 * (c1 - s1))/3600.0;
    }
#endif

    return (0);
}


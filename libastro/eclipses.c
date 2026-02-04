#include <math.h>
#include "astro.h"

/* compute the center of totality of a solar eclipse, if any.
 * return 1 if an eclipse was detected, 0 otherwise.
 *
 * N.B. this code is geometrically correct for all objects, but in practice
 * only the sun and moon are computed accurately enough by xephem to make this
 * worth while for now. in particular, I have tried tests against several S&T
 * tables of minor planet occultations and while the asteroids are computed
 * well enough for visual identification remember that at even 1 AU the earth
 * only subtends 18 arc seconds and the asteroids are not computed *that*
 * accurately (especially since we do not yet include perturbations).
 *
 * I will try to describe the algorithm: the executive summary is that we
 * are striving for the spherical arc subtended by the intersection of a line
 * from one object and the earth's center and the line from the other object
 * to the earth's center.
 *
 * N.B. I tried just computing the intersection of a line connecting the two
 * objects and a unit sphere but it suffered from terrible numerical
 * instability.
 *
 * start in a plane defined by the center of the earth, the north pole and
 * object obj0. label the center of the earth as O, the location of object 0
 * as P0, and place object 1 someplace off the line O-P0 and label it P1.
 * what you have actually placed is the location of P1 as it projected onto
 * this place; ie, we are only working with dec here.  define decA as the
 * angle P1-O-P0; it is so named because it is the amount of declination
 * subtended from P0 to P1. Project the line P0-P1 back to a line
 * perpendicular to the line O-P0 at O. decD is the distance from O to the
 * point where P0-P1 intersects the line. if it is less than the earth radius
 * we have an occultation! now do all this again only this time place
 * yourself in a plane defined by the real locations of O, P0 and P1. and
 * repeat everything except this time use the real angled subtended in the
 * sky between P0 and P1 (not just the dec difference). this angle we define
 * as skyA (and its projection back onto a plane perpendicular to P0-P1
 * through O we call skyD). what we want next is the spherical angle
 * subtended between the point at which O-P0 intersects the earth's surface
 * (which is just the geocentric coords of P0) and the point where a line
 * from the tip of skyD to P1 intersects the earth's surface. we call this
 * skyT (I used tau in my original sketch) and I will let you work out the
 * trig (it's just planar trig since you are working in the O-P0-P1 plane).
 * this gives us the spherical angle between the two lines and the earth
 * surface; now all we need is the angle.image yourself at P0 now looking
 * right at O. we see decD as a vertical line and SkyD as a line going off
 * from O at an angle somewhere. the angle between these lines we define as
 * theta. knowing decD and skyD and knowing that there is a right angle at
 * the tip of decD between O and the tip of skyD we can compute the angle
 * between them. theta.  now just use a little spherical trig to find where
 * our arc ends up, compute the new RA, compute longitude by subtracting
 * gst, set latitude to dec, and return 1 for success!
 */
int
solar_eclipse_center(sunp, moonp, np, lt_out, lg_out)
Obj *sunp, *moonp;
Now *np;
double *lt_out, *lg_out;
{
	Obj obj0, obj1;			/* use copies */
	double r0, r1;			/* dist to objects, in earth radii */
	double theta;			/* angle between projections */
	double decD, decA;		/* dec-only proj dist and angle */
	double skyD, skyA, skyP, skyT;	/* full sky projection */
	Now now;			/* local copy to compute EOD info */
	double lst, gst;		/* local and UTC time */
	double lt, lg;			/* lat/long */
	double sD, dRA;

	obj0 = *sunp;
	obj1 = *moonp;
	now = *np;

	now.n_epoch = EOD;
	(void) obj_cir (&now, &obj0);
	if (is_ssobj(&obj0))
	    r0 = obj0.s_edist*(MAU/ERAD);	/* au to earth radii */
	else
	    r0 = 1e7;				/* way past pluto */

	(void) obj_cir (&now, &obj1);
	if (is_ssobj(&obj1))
	    r1 = obj1.s_edist*(MAU/ERAD);	/* au to earth radii */
	else
	    r1 = 1e7;				/* way past pluto */

	decA = obj1.s_gaedec - obj0.s_gaedec;
	decD = r0*r1*sin(decA)/(r0 - r1);	/* similar triangles */
	if (fabs(decD) >= 1.0)
	    return 0;

	skyA = acos (sin(obj0.s_gaedec)*sin(obj1.s_gaedec) +
				    cos(obj0.s_gaedec)*cos(obj1.s_gaedec) *
					cos(obj0.s_gaera-obj1.s_gaera));
	skyD = r0*r1*sin(skyA)/(r0 - r1);	/* similar triangles */
	if (fabs(skyD) >= 1.0)
	    return 0;

	/* skyP is angle subtended by skyD as seen from obj0 (I called it psi).
	 * skyT is angle subtended by line from earth center to obj0 to a
	 *   point where the line from obj0 to the tip of skyD intersects the
	 *   earth surface (I called it tau).
	 */
	skyP = atan(skyD/r0);
	skyT = asin(skyD*r0/sqrt(r0*r0+skyD*skyD)) - skyP;

	theta = acos(decD/skyD);
	solve_sphere (theta, skyT, sin(obj0.s_gaedec), cos(obj0.s_gaedec),
								    &sD, &dRA);

	lt = asin(sD);

	if (obj1.s_gaera > obj0.s_gaera)
	    dRA = -dRA;	/* eastward */

	lst = obj0.s_gaera - dRA;
	utc_gst (mjd_day(mjd), mjd_hr(mjd), &gst);
	lg = lst - hrrad(gst);
	while (lg < -PI) lg += 2*PI;
	while (lg >  PI) lg -= 2*PI;

        *lt_out = lt;
        *lg_out = lg;
        return 1;
}

#include "astro.h"

int main() {
     char out[1024];

     printf("Precession J2000 -> J2015\n");
     double ra, pra;
     double dec, pdec;
     for (double i=0.0; i<7.0; i+=1.0) {
          ra = pra = i;
          dec = pdec = -i;
          precess(J2000, J2015, &pra, &pdec);
          printf("RA  %+.16f -> %+.14f\n", ra, pra);
          printf("Dec %+.16f -> %+.14f\n", dec, pdec);
     }

     printf("\nPluto discontinuity when Chapront (1995) ends on 2247/10/1\n");
     double mjd0, mjd1;
     cal_mjd(5, 16, 2247, &mjd0);
     cal_mjd(9, 27, 2248, &mjd1);
     Now n;
     n.n_lat = n.n_lng = n.n_tz = n.n_temp = n.n_pressure = n.n_elev
          = n.n_dip = 0.0;
     n.n_epoch = J2000;
     Obj o;
     o.o_type = PLANET;
     o.pl_code = PLUTO;
     o.pl_moon = X_PLANET;

     n.n_mjd = mjd0;
     obj_cir(&n, &o);
     fs_sexa(out, radhr(o.s_gaera), 2, 360000);
     printf("RA = %s\n", out);

     n.n_mjd = mjd1;
     obj_cir(&n, &o);
     fs_sexa(out, radhr(o.s_gaera), 2, 360000);
     printf("RA = %s\n", out);

/* From Elwood's email, this is how bad the discontinuity was before the
   commit 139b2ea2 'Revert "Astronomical Almanac 2020 Pluto elements...':
5/16/2247 00:00 UTC    4.1.0 RA 16:52:53.90  Dec -10:36:26.5
                       3.6.7 RA 16:52:53.90  Dec -10:36:26.5
9/27/2248 00:00 UTC    4.1.0 RA 11:52:36.30  Dec 16:13:59.6
                       3.6.7 RA 16:46:16.01  Dec -11:05:43.7
*/

     printf("\nEclipse track for 2034 March 20 at 20-minute intervals:\n");

     Obj sun, moon;
     Now now;
     double lt, lg;

     sun.o_type = PLANET;
     sun.pl_code = SUN;
     sun.pl_moon = X_PLANET;

     moon.o_type = PLANET;
     moon.pl_code = MOON;
     moon.pl_moon = X_PLANET;

     cal_mjd(3, 20, 2034, &(now.n_mjd));
     now.n_mjd += 9.0 / 24.0;  /* advance 9 more hours */

     int r;
     while (1) {
          r = solar_eclipse_center(&sun, &moon, &now, &lt, &lg);
          if (!r)
               break;
          printf("Lat %.1f Lon %.1f\n", raddeg(lt), raddeg(lg));
          now.n_mjd += 1.0 / 24.0 / 3.0;  /* 1/3 of an hour is 20 minutes */
     }

     return 0;
}

#include "astro.h"

int main() {
     printf("Precession J2000 -> J2015\n");
     double ra, pra;
     double dec, pdec;
     for (double i=0.0; i<7.0; i+=1.0) {
          ra = pra = i;
          dec = pdec = -i;
          precess(J2000, J2015, &pra, &pdec);
          printf("RA  %+.16f -> %+.16f\n", ra, pra);
          printf("Dec %+.16f -> %+.16f\n", dec, pdec);
     }
     return 0;
}

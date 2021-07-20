#include <R.h>
#include <Rmath.h>

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }

void F77_SUB(normrnd)(double r[], int *n) {
   int i;
   for (i = 0; i<*n; i++) {
      r[i] = rnorm(0.0, 1.0); //norm_rand();
   }
}

#include <R.h>
#include <Rdefines.h>

void F77_SUB(prtrac)(int *neval, double *fx, double *nmg, int *n, double x[]) {
   
   SEXP pv, call;
   
   Rprintf("Iteration %d, F(x) = %21.16e, max|g(x)| = %11.4e\n", *neval, *fx, *nmg);
   
   Rprintf("Parameter:\n");
   
   PROTECT( pv = NEW_NUMERIC(*n));
   Memcpy(REAL(pv), x, *n);
   
   PROTECT( call = lang2( install("print"), pv) );
   PROTECT( eval( call, R_GlobalEnv) );
   
   UNPROTECT(3);
   
   Rprintf("\n");
}

void F77_SUB(prterr)(double *f) {
   Rprintf("DEBUG ::: f = %21.16e\n", *f);
}

void F77_SUB(prteri)(int *f) {
   Rprintf("DEBUG ::: f = %i\n", *f);
}

/*
 void F77_SUB(prline)(double *a, double sl[]) {
 Rprintf(" Line search: alpha =%11.4e, dphi(0) =%11.4e, dphi(1) =%11.4e\n",
 *a, sl[0], sl[1]);
 }
 
 void F77_SUB(prconv)() {
 Rprintf(" Optimization has converged\n");
 }
 
 void F77_SUB(prfail)(int *neval) {
 Rprintf(" Optimization stopped after %d function evaluations\n", *neval);
 }
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <Rdefines.h>

#ifdef SUPPORT_OPENMP
   #include <omp.h>
#endif

/*
      SUBROUTINE MAIN(INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,NIMAT,OMAT,
     $                NOMAT,NOBS,NSET,INFO)
*/
extern void F77_NAME(main)(int[], int*, double[], int*, double[], int*, double[], int*, double[], int*, int*, int*, double[], int*);

// Define a "global" variable for the current environment
extern SEXP R_envir;
int cnt;

SEXP sxp,sxm,suo,st;
SEXP R_envir;
SEXP R_amat;
SEXP R_bmat;
SEXP R_cmat;
SEXP R_dmat;
SEXP R_smat;
SEXP R_sigmat;

void F77_SUB(evlmat)(int *nmat, double xm[], int *nxm, double xp[], int *nxp, double uo[], int *nuo, double *t, double matrix[], int *N, int *M) {
   /*
   
      nmat is indicating which matrix should be evaluated and returned

         A      : nmat =  1
         B      : nmat =  2
         C      : nmat =  3
         D      : nmat =  4

         fvecxj : nmat =  5
         fvecuj : nmat =  6
         hvecxj : nmat =  7
         hvecuj : nmat =  8

         SIGMAT : nmat =  9
         SMAT   : nmat = 10

   */

   cnt++;

   int i;
   for (i=0;i<*nxp;i++) REAL(sxp)[i] = xp[i];
   for (i=0;i<*nxm;i++) REAL(sxm)[i] = xm[i];
   for (i=0;i<*nuo;i++) REAL(suo)[i] = uo[i];
                        REAL(st)[0] = *t;
   SEXP R_fcall,ans;

   switch (*nmat) {
      case 1:
      case 5:
         PROTECT(R_fcall = lang5(R_amat,sxp,sxm,suo,st));
         break;
      case 2:
      case 6:
         PROTECT(R_fcall = lang5(R_bmat,sxp,sxm,suo,st));
         break;
      case 3:
      case 7:
         PROTECT(R_fcall = lang5(R_cmat,sxp,sxm,suo,st));
         break;
      case 4:
      case 8:
         PROTECT(R_fcall = lang5(R_dmat,sxp,sxm,suo,st));
         break;
      case 9:
         PROTECT(R_fcall = lang5(R_sigmat,sxp,sxm,suo,st));
         break;
      case 10:
         PROTECT(R_fcall = lang5(R_smat,sxp,sxm,suo,st));
         break;
   }
   PROTECT(ans = eval(R_fcall, R_envir));

   for (i=0;i<((*M)*(*N));i++) matrix[i] = REAL(ans)[i];

   UNPROTECT(2);
}


SEXP ctsmdriver(SEXP amat, SEXP bmat, SEXP cmat, SEXP dmat, SEXP smat, SEXP sigmat, SEXP rho) {

   R_envir = rho;

   int  nint  = asInteger(findVar(install(".nint"), rho)),
        ndoub = asInteger(findVar(install(".ndoub"), rho)),
        ntmat = asInteger(findVar(install(".ntmat"), rho)),
        nimat = asInteger(findVar(install(".nimat"), rho)),
        nomat = asInteger(findVar(install(".nomat"), rho)),
        nobs  = asInteger(findVar(install(".nobs"), rho)),
        nset  = asInteger(findVar(install(".nset"), rho),
        threads = asInteger(findVar(install(".threads"),rho)));

   SEXP ints   = findVar(install(".ints"), rho),
        doubls = findVar(install(".doubls"), rho),
        tmat   = findVar(install(".tmat"), rho),
        imat   = findVar(install(".imat"), rho),
        omat   = findVar(install(".omat"), rho),
        trace  = findVar(install(".trace"), rho),
        info   = findVar(install(".info"), rho);


   PROTECT(sxp = allocVector(REALSXP,3));
   PROTECT(sxm = allocVector(REALSXP,14));
   PROTECT(suo = allocVector(REALSXP,2));
   PROTECT(st = NEW_NUMERIC(1));

   R_amat = amat;
   R_bmat = bmat;
   R_cmat = cmat;
   R_dmat = dmat;
   R_smat = smat;
   R_sigmat = sigmat;

   int cpus = 1;

   cnt = 0;

#ifdef SUPPORT_OPENMP
   omp_set_num_threads(threads);
   cpus = omp_get_num_procs();
#endif

// extern void F77_NAME(main)(int int[], int nint, double doubls[], int ndoub, double tmat[], int ntmat, double imat[], int nimat, double omat[], int nomat, int nobs, int nset, int info);
   F77_CALL(main)(INTEGER(ints),&nint,
                  REAL(doubls),&ndoub,
                  REAL(tmat),&ntmat,
                  REAL(imat),&nimat,
                  REAL(omat),&nomat,
                  &nobs,&nset,REAL(trace),INTEGER(info));


   UNPROTECT(4);

   return R_NilValue;
}

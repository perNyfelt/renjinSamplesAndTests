#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <Rdefines.h>

/*
#ifdef _OPENMP
   #include <omp.h>
   
   #if !defined(WIN32)
      #define CSTACK_DEFNS 7
      #include "Rinterface.h"
   #endif
#endif
*/
   
/*
      SUBROUTINE MAIN(INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,NIMAT,OMAT,
     $                NOMAT,NOBS,NSET,INFO)
 */
extern void F77_NAME(main)(int[], int*, double[], int*, double[], int*, double[], int*, double[], int*, int*, int*, double[], int*);

/*
      SUBROUTINE LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,
     $                 NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,FPEN,FPRIOR,
     $                 EPSM,JOB,FUNK,INFO)
 */
extern void F77_NAME(llike)(
		int*, // NX 
		double[], // X
		int[], int*, // INTS, NINT
		double[], int*, // DOUBLS, NDOUB
		double[], int*, // TMAT, NTMAT
		double[], int*, // IMAT, NIMAT
		double[], int*, // OMAT, NOMAT
		int[], int*, // NOBS, NSET
		int*, // NMISST
		double*, // FPEN
		double*, // FPRIOR
		double*, // EPSM
		int*, // JOB
		double*, // FUNK
		int*, // INFO
      double[] // VX0
		);

SEXP ctsmdriver(SEXP rho) {

   int  nint  = asInteger(findVar(install(".nint"), rho)),
        ndoub = asInteger(findVar(install(".ndoub"), rho)),
        ntmat = asInteger(findVar(install(".ntmat"), rho)),
        nimat = asInteger(findVar(install(".nimat"), rho)),
        nomat = asInteger(findVar(install(".nomat"), rho)),
        nset  = asInteger(findVar(install(".nset"), rho)),
        threads = asInteger(findVar(install(".threads"),rho));

   SEXP ints   = findVar(install(".ints"), rho),
        doubls = findVar(install(".doubls"), rho),
        tmat   = findVar(install(".tmat"), rho),
        imat   = findVar(install(".imat"), rho),
        nobs   = findVar(install(".nobs"), rho),
        omat   = findVar(install(".omat"), rho),
        trace  = findVar(install(".trace"), rho),
        cpus   = findVar(install(".cpus"), rho),
        info   = findVar(install(".info"), rho);
/*
#ifdef _OPENMP
   omp_set_dynamic(0);
   omp_set_num_threads(threads);
   
   #if !defined(WIN32)
      R_CStackLimit=(uintptr_t)-1;
   #endif
   
   INTEGER(cpus)[0] = omp_get_num_procs();
#else
   threads = 1;
   INTEGER(cpus)[0] = threads;
#endif
 */

// extern void F77_NAME(main)(int int[], int nint, double doubls[], int ndoub, double tmat[], int ntmat, double imat[], int nimat, double omat[], int nomat, int nobs, int nset, int info);
   F77_CALL(main)(INTEGER(ints),&nint,
                  REAL(doubls),&ndoub,
                  REAL(tmat),&ntmat,
                  REAL(imat),&nimat,
                  REAL(omat),&nomat,
                  INTEGER(nobs),&nset,REAL(trace),INTEGER(info));

   return R_NilValue;
}

SEXP ctsmllike(SEXP rho) {
	
	int  nint  = asInteger(findVar(install(".nint") , rho)),
	     ndoub = asInteger(findVar(install(".ndoub"), rho)),
	     ntmat = asInteger(findVar(install(".ntmat"), rho)),
	     nimat = asInteger(findVar(install(".nimat"), rho)),
	     nomat = asInteger(findVar(install(".nomat"), rho)),
	     nset  = asInteger(findVar(install(".nset") , rho)),
	     job   = asInteger(findVar(install(".job")  , rho));

    SEXP ints   = findVar(install(".ints"), rho),
	     doubls = findVar(install(".doubls"), rho),
	     tmat   = findVar(install(".tmat"), rho),
	     imat   = findVar(install(".imat"), rho),
	     omat   = findVar(install(".omat"), rho),
	     nobs   = findVar(install(".nobs"), rho),
	     info   = findVar(install(".info"), rho),
        vx0   = findVar(install(".vx0"), rho);
    
    int nx,nmisst;
    double x,funk,fpen,fprior,epsm;

    x = 0;
    nx = 1;
    epsm = -1;
    nmisst = 0;
    
    /*
    extern void F77_NAME(llike)(
    		int*, // NX 
    		double[], // X
    		int[], int*, // INTS, NINT
    		double[], int*, // DOUBLS, NDOUB
    		double[], int*, // TMAT, NTMAT
    		double[], int*, // IMAT, NIMAT
    		double[], int*, // OMAT, NOMAT
    		double[], int*, // NOBS, NSET
    		int*, // NMISST
    		double*, // FPEN
    		double*, // FPRIOR
    		double*, // EPSM
    		int*, // JOB
    		double*, // FUNK
    		int* // INFO
    		);
    */
    
	F77_CALL(llike)(&nx, &x, INTEGER(ints), &nint,
			REAL(doubls), &ndoub,
            REAL(tmat),&ntmat,
            REAL(imat),&nimat,
            REAL(omat),&nomat,
			INTEGER(nobs), &nset,
			&nmisst,
			&fpen,
			&fprior,
			&epsm,
			&job,
			&funk,
			INTEGER(info),
         REAL(vx0));
	   
	return R_NilValue;
}

SEXP ctsmllike2(SEXP rho) {
	
	int  nint  = asInteger(findVar(install(".nint") , rho)),
	     ndoub = asInteger(findVar(install(".ndoub"), rho)),
	     ntmat = asInteger(findVar(install(".ntmat"), rho)),
	     nimat = asInteger(findVar(install(".nimat"), rho)),
	     nomat = asInteger(findVar(install(".nomat"), rho)),
	     nset  = asInteger(findVar(install(".nset") , rho)),
        nx    = asInteger(findVar(install(".nx") , rho)),
	     job   = asInteger(findVar(install(".job")  , rho));

    SEXP x   = findVar(install(".x"), rho),
        ints   = findVar(install(".ints"), rho),
	     doubls = findVar(install(".doubls"), rho),
	     tmat   = findVar(install(".tmat"), rho),
	     imat   = findVar(install(".imat"), rho),
	     omat   = findVar(install(".omat"), rho),
	     nobs   = findVar(install(".nobs"), rho),
	     nmisst = findVar(install(".nmisst"), rho),
	     f      = findVar(install(".f"), rho),
	     fpen   = findVar(install(".fpen"), rho),
	     fprior = findVar(install(".fprior"), rho),
	     info   = findVar(install(".info"), rho),
        vx0   = findVar(install(".vx0"), rho);

    double epsm = -1;
    
    /*
    extern void F77_NAME(llike)(
    		int*, // NX 
    		double[], // X
    		int[], int*, // INTS, NINT
    		double[], int*, // DOUBLS, NDOUB
    		double[], int*, // TMAT, NTMAT
    		double[], int*, // IMAT, NIMAT
    		double[], int*, // OMAT, NOMAT
    		double[], int*, // NOBS, NSET
    		int*, // NMISST
    		double*, // FPEN
    		double*, // FPRIOR
    		double*, // EPSM
    		int*, // JOB
    		double*, // FUNK
    		int* // INFO
    		);
    */
    
	F77_CALL(llike)(&nx, REAL(x), INTEGER(ints), &nint,
			REAL(doubls), &ndoub,
            REAL(tmat),&ntmat,
            REAL(imat),&nimat,
            REAL(omat),&nomat,
			INTEGER(nobs), &nset,
			INTEGER(nmisst),
			REAL(fpen),
			REAL(fprior),
			&epsm,
			&job,
			REAL(f),
			INTEGER(info),
         REAL(vx0));
	   
	return R_NilValue;
}

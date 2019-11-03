C-----------------------------------------------------------------------C
C                                                                       C
C     V A R I O U S    M A T H E M A T I C A L    U T I L I T I E S     C
C                                                                       C
C-----------------------------------------------------------------------C
      SUBROUTINE OPTIM(N,X,XMINS,XMAXS,EPS,ETA,MAXFUN,ITR,W,IW,ICONTR,
     $          INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,
     $               NOBS,NSET,NMISST,FPEN,FPRIOR,EPSM,TRACE,VINFO,INFO)
C-----------------------------------------------------------------------C
C     UNCONSTRAINED MINIMIZATION OF A SCALAR FUNCTION.                  C
C-----------------------------------------------------------------------C
      INTEGER          N,MAXFUN,ITR,IW,ICONTR,NINT,INTS(NINT),NDOUB,
     $                 NTMAT,NIMAT,NOMAT,NSET,NOBS(NSET),NMISST,
     $                 VINFO(N),INFO
      DOUBLE PRECISION X(N),EPS,ETA,W(IW),DOUBLS(NDOUB),TMAT(NTMAT),
     $                 IMAT(NIMAT),OMAT(NOMAT),FPEN,FPRIOR,EPSM
      INTEGER          IL,NNH,NH,ND,NW,NXA,NGA,NXB,NGB,NG,IPRINT,MAXF1
      DOUBLE PRECISION OD2,OD3,F
      DOUBLE PRECISION TRACE(*),XMINS(N),XMAXS(N)
C
      IL=(N*(N+15))/2+1
  10  CONTINUE
      NNH=(N*(N+1))/2
      NH=2
      ND=NH+NNH
      NW=ND+N
      NXA=NW+N
      NGA=NXA+N
      NXB=NGA+N
      NGB=NXB+N
      NG=NGB+N
      IPRINT=1
      MAXF1=MAXFUN+1
      OD3=ETA**(1.0D0/3.0D0)
      OD2=DSQRT(ETA)
      F=0.0D0
      CALL VA13CD(N,X,XMINS,XMAXS,F,W(NG),EPS,OD2,OD3,MAXF1,ITR,W(NH),
     $ NNH,W(ND),W(NW),W(NXA),W(NGA),W(NXB),W(NGB),IPRINT,INTS,NINT,
     $ DOUBLS,NDOUB,TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,
     $ NMISST,FPEN,FPRIOR,EPSM,TRACE,VINFO,INFO)
      IF (INFO.NE.0) RETURN
      W(1)=F
      IF (MAXF1.GT.MAXFUN) ICONTR=2
      MAXFUN=MAXF1
      RETURN
      END
C
      SUBROUTINE VA13CD (N,X,XMINS,XMAXS,F,G,ACC,OD2,OD3,MAXFUN,ITR,H,
     $                 IH,D,W,XA,GA,XB,GB,IPRINT,INTS,NINT,DOUBLS,NDOUB,
     $                   TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,
     $                   NMISST,FPEN,FPRIOR,EPSM,TRACE,VINFO,INFO)
C-----------------------------------------------------------------------C
C     THIS IS A MODIFIED VERSION OF THE SUBROUTINE WITH THE SAME NAME   C
C     FROM THE HARWELL LIBRARY.                                         C
C-----------------------------------------------------------------------C
      INTEGER          N,MAXFUN,ITR,IH,NINT,INTS(NINT),NDOUB,NTMAT,
     $                 NIMAT,NOMAT,NSET,NOBS(NSET),NMISST,VINFO(N),INFO
      DOUBLE PRECISION X(N),F,G(N),ACC,OD2,OD3,H(IH),D(N),W(N),XA(N),
     $                 GA(N),XB(N),GB(N),DOUBLS(NDOUB),TMAT(NTMAT),
     $                 IMAT(NIMAT),OMAT(NOMAT),FPEN,FPRIOR,EPSM
      INTEGER          I,IP,IPRA,IPRINT,IR,ISFV,K,MD,MODE,NFUN,NP
      DOUBLE PRECISION C,DFF,FA,DGA,STMIN,STEPBD,STEPLB,FMIN,GMIN,STEP,
     $                 FB,GL1,GL2,DGB
      INTEGER          TRPOS
      DOUBLE PRECISION TRACE(*),MAXNGA,XTMP(N),XMINS(N),XMAXS(N)
C
      MD=1
      MODE=1
C
C     CALCULATE THE INITIAL FUNCTION VALUE
C
   10 CALL FDF(N,X,G,F,OD2,OD3,MD,INTS,NINT,DOUBLS,NDOUB,
     $         TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,
     $         NMISST,FPEN,FPRIOR,EPSM,VINFO,INFO)
      IF (INFO.NE.0) RETURN
      NFUN=1
      ITR=0
      NP=N+1
C
C     SET THE HESSIAN TO A DIAGONAL MATRIX DEPENDING ON SCALE(.)
C
      IF (MODE.GE.2) GOTO 60
   20 C=0.0D0
      DO 30 I=1,N
   30   C=DMAX1(C,DABS(G(I)))
      IF (C.LE.0.0D0) C=1.0D0
      K=(N*NP)/2
      C=1D-2*C
      DO 40 I=1,K
   40   H(I)=0.0D0
      K=1
      DO 50 I=1,N
        H(K)=C
   50   K=K+NP-I
      GOTO 100
C
C     FACTORIZE THE GIVEN HESSIAN MATRIX
C
   60 IF (MODE.GE.3) GOTO 80
      CALL MC11BD(H,IH,N,K)
      IF (K.GE.N) GOTO 100
   70 GOTO 20
C
C     CHECK THAT THE GIVEN DIAGONAL ELEMENTS ARE POSITIVE
C
   80 K=1
      DO 90 I=1,N
        IF (H(K).LE.0.0D0) GOTO 70
   90   K=K+NP-I
C
C     SET SOME VARIABLES FOR THE FIRST ITERATION
C
  100 DFF=0.0D0
      IPRA=IABS(IPRINT)
      IP=IABS(IPRA-1)
  110 FA=F
      ISFV=1
      DO 120 I=1,N
        XA(I)=X(I)
  120   GA(I)=G(I)
C
C     BEGIN THE ITERATION BY GIVING THE REQUIRED PRINTING
C
  130 IP=IP+1
      IF (IP.NE.IPRA) GOTO 140
      IP=0
  140 ITR=ITR+1
C
C     SAVE THE CURRENT H (DIAGONAL), GA, C, XA, FA
C              SIZE    N +  N + 1 + N + 1 = 3*N+2
C
      TRPOS=(ITR-1)*(3*N+2)
      K=1
      DO I=1,N
        TRACE(TRPOS+I) = H(K)
        K=K+NP-I
      ENDDO
      MAXNGA=0.0D0
      DO I=1,N
        TRACE(TRPOS+N+I)=GA(I)
        MAXNGA=DMAX1(MAXNGA,DABS(GA(I)))
      ENDDO
      TRACE(TRPOS+2*N+1)=C
C     Copy and backtransform the current parameter estimate      
      DO I=1,N
         XTMP(I) = XA(I)
         CALL ATF(XTMP(I),XMINS(I),XMAXS(I))
         TRACE(TRPOS+2*N+1+I) = XTMP(I)
      ENDDO
      TRACE(TRPOS+3*N+2) = FA
C     Print current info      
      CALL PRTRAC(ITR,FA,MAXNGA,N,XTMP)
C
C     CALCULATE THE SEARCH DIRECTION OF THE ITERATION
C
      DO 150 I=1,N
  150   D(I)=-GA(I)
      CALL MC11ED(H,IH,N,D,W,N)
C
C     CALCULATE A LOWER BOUND ON THE STEP-LENGTH
C     AND THE INITIAL DIRECTIONAL DERIVATIVE
C
      C=0.0D0
      DGA=0.0D0
      DO 160 I=1,N
        C=DMAX1(C,DABS(D(I)))
  160   DGA=DGA+GA(I)*D(I)
C
C     TEST IF THE SEARCH DIRECTION IS DOWNHILL
C
      IF (DGA .GE. 0.0D0) GOTO 240
C
C     SET THE INITIAL STEP-LENGTH OF THE LINE SEARCH
C
      STMIN=0.0D0
      STEPBD=0.0D0
      STEPLB=ACC/C
      FMIN=FA
      GMIN=DGA
      STEP=1.0D0
      IF (DFF.LE.0.0D0) STEP=DMIN1(STEP,1.0D0/C)
      IF (DFF.GT.0.0D0) STEP=DMIN1(STEP,(DFF+DFF)/(-DGA))
  170 C=STMIN+STEP
C
C     TEST WHETHER FDF HAS BEEN CALLED MAXFUN TIMES
C
      IF (NFUN.EQ.MAXFUN) GOTO 250
      NFUN=NFUN+1
C
C     CALCULATE ANOTHER FUNCTION VALUE AND GRADIENT
C
      DO 180 I=1,N
  180   XB(I)=XA(I)+C*D(I)
      CALL FDF(N,XB,GB,FB,OD2,OD3,MD,INTS,NINT,DOUBLS,NDOUB,
     $         TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,
     $         NMISST,FPEN,FPRIOR,EPSM,VINFO,INFO)
      IF (INFO.NE.0) RETURN
C
C     STORE THIS FUNCTION VALUE IF IT IS THE SMALLEST SO FAR
C
      ISFV=MIN0(2,ISFV)
      IF (FB.GT.F) THEN
        IF (ITR.GT.2*N) MD=2
        GOTO 220
      END IF
      IF (FB.LT.F) GOTO 200
      GL1=0.0D0
      GL2=0.0D0
      DO 190 I=1,N
        GL1=GL1+G(I)**2
  190   GL2=GL2+GB(I)**2
      IF (GL2.GE.GL1) GOTO 220
  200 ISFV=3
      F=FB
      DO 210 I=1,N
        X(I)=XB(I)
  210   G(I)=GB(I)
C
C     CALCULATE THE DIRECTIONAL DERIVATIVE AT THE NEW POINT
C
  220 DGB=0.0D0
      DO 230 I=1,N
  230   DGB=DGB+GB(I)*D(I)
C
C     BRANCH IF WE HAVE FOUND A NEW LOWER BOUND ON THE STEP-LENGTH
C
      IF (FB-FA.LE.0.1D0*C*DGA) GOTO 280
C
C     FINISH THE ITERATION IF THE CURRENT STEP IS STEPLB
C
      IF (STEP.GT.STEPLB) GOTO 270
  240 IF (ISFV.GE.2) GOTO 110
C
C     AT THIS STAGE THE WHOLE CALCULATION IS COMPLETE
C
  250 CONTINUE
  260 MAXFUN=NFUN
      RETURN
C
C     CALCULATE A NEW STEP-LENGTH BY CUBIC INTERPOLATION
C
  270 STEPBD=STEP
      C=GMIN+DGB-3D0*(FB-FMIN)/STEP
      C=GMIN/(C+GMIN-DSQRT(C*C-GMIN*DGB))
      STEP=STEP*DMAX1(0.1D0,C)
      GOTO 170
C
C     SET THE NEW BOUNDS ON THE STEP-LENGTH
C
  280 STEPBD=STEPBD-STEP
      STMIN=C
      FMIN=FB
      GMIN=DGB
C
C     CALCULATE A NEW STEP-LENGTH BY EXTRAPOLATION
C
      STEP=9D0*STMIN
      IF (STEPBD.GT.0.0D0) STEP=.5D0*STEPBD
      C=DGA+3D0*DGB-4D0*(FB-FA)/STMIN
      IF (C.GT.0.0D0) STEP=DMIN1(STEP,STMIN*DMAX1(1.0D0,-DGB/C))
      IF (DGB.LT.0.7D0*DGA) GOTO 170
C
C     TEST FOR CONVERGENCE OF THE ITERATIONS
C
      ISFV=4-ISFV
      IF (STMIN+STEP.LE.STEPLB) GOTO 240
C
C     REVISE THE SECOND DERIVATIVE MATRIX
C
      IR=-N
      DO 290 I=1,N
        XA(I)=XB(I)
        XB(I)=GA(I)
        D(I)=GB(I)-GA(I)
  290   GA(I)=GB(I)
      CALL MC11AD(H,IH,N,XB,1.0D0/DGA,W,IR,1,0.0D0)
      IR=-IR
      CALL MC11AD(H,IH,N,D,1.0D0/(STMIN*(DGB-DGA)),D,IR,0,0.0D0)
C
C     BRANCH IF THE RANK OF THE NEW MATRIX IS DEFICIENT
C
      IF (IR.LT.N) GOTO 250
C
C     BEGIN ANOTHER ITERATION
C
      DFF=FA-FB
      FA=FB
      GOTO 130
      END
C
*----------------------------------------------------------------------|
      SUBROUTINE DGPADM( IDEG,M,T,H,LDH,WSP,LWSP,IPIV,IEXPH,NS,IFLAG )

      INTEGER IDEG, M, LDH, LWSP, IEXPH, NS, IFLAG, IPIV(M)
      DOUBLE PRECISION T, H(LDH,M), WSP(LWSP)

*-----PURPOSE----------------------------------------------------------|
*
*     COMPUTES EXP(T*H), THE MATRIX EXPONENTIAL OF A GENERAL MATRIX IN
*     FULL, USING THE IRREDUCIBLE RATIONAL PADE APPROXIMATION TO THE
*     EXPONENTIAL FUNCTION EXP(X) = R(X) = (+/-)( I + 2*(Q(X)/P(X)) ),
*     COMBINED WITH SCALING-AND-SQUARING.
*
*-----ARGUMENTS--------------------------------------------------------|
*
*     IDEG      : (INPUT) THE DEGREE OF THE DIAGONAL PADE TO BE USED.
*                 A VALUE OF 6 IS GENERALLY SATISFACTORY.
*
*     M         : (INPUT) ORDER OF H.
*
*     H(LDH,M)  : (INPUT) ARGUMENT MATRIX.
*
*     T         : (INPUT) TIME-SCALE (CAN BE < 0).
*
*     WSP(LWSP) : (WORKSPACE/OUTPUT) LWSP .GE. 4*M*M+IDEG+1.
*
*     IPIV(M)   : (WORKSPACE)
*
*     IEXPH     : (OUTPUT) NUMBER SUCH THAT WSP(IEXPH) POINTS TO EXP(TH)
*                 I.E., EXP(TH) IS LOCATED AT WSP(IEXPH ... IEXPH+M*M-1)
*
*                 NOTE: IF THE ROUTINE WAS CALLED WITH WSP(IPTR),
*                       THEN EXP(TH) WILL START AT WSP(IPTR+IEXPH-1).
*
*     NS        : (OUTPUT) NUMBER OF SCALING-SQUARING USED.
*
*     IFLAG     : (OUTPUT) EXIT FLAG.
*                      0 - NO PROBLEM
*                     <0 - PROBLEM
*
*----------------------------------------------------------------------|
*     ROGER B. SIDJE (RBS@MATHS.UQ.EDU.AU)
*     EXPOKIT: SOFTWARE PACKAGE FOR COMPUTING MATRIX EXPONENTIALS.
*     ACM - TRANSACTIONS ON MATHEMATICAL SOFTWARE, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*     MODIFIED TO RETURN AN ERROR FLAG INSTEAD OF TERMINATING, WHEN
*     TRYING TO COMPUTE THE EXPONENTIAL OF A MATRIX WITH TOO LARGE
*     LARGE ELEMENTS.
*
*     N.R.KRISTENSEN, KT, TECHNICAL UNIVERSITY OF DENMARK, 2000
*----------------------------------------------------------------------|
      INTEGER MM,I,J,K,IH2,IP,IQ,IUSED,IFREE,IODD,ICOEF,IPUT,IGET
      DOUBLE PRECISION HNORM,SCALE,SCALE2,CP,CQ
      INTRINSIC INT,ABS,DBLE,LOG,MAX
*---  CHECK RESTRICTIONS ON INPUT PARAMETERS ...
      MM = M*M
      IFLAG = 0
      IF ( LDH.LT.M ) IFLAG = -1
      IF ( LWSP.LT.4*MM+IDEG+1 ) IFLAG = -2
      IF (IFLAG.NE.0) RETURN
*
*---  INITIALISE POINTERS ...
*
      ICOEF = 1
      IH2 = ICOEF + (IDEG+1)
      IP  = IH2 + MM
      IQ  = IP + MM
      IFREE = IQ + MM
*
*---  SCALING: SEEK NS SUCH THAT ||T*H/2^NS|| < 1/2;
*     AND SET SCALE = T/2^NS ...
*
      DO I = 1,M
         WSP(I) = 0.0D0
      ENDDO
      DO J = 1,M
         DO I = 1,M
            WSP(I) = WSP(I) + ABS( H(I,J) )
         ENDDO
      ENDDO
      HNORM = 0.0D0
      DO I = 1,M
         HNORM = MAX( HNORM,WSP(I) )
      ENDDO
      HNORM = ABS( T*HNORM )
      IF (HNORM.EQ.0.0D0) THEN
         IFLAG=-1
         RETURN
        ENDIF
      NS = MAX( 0,INT(LOG(HNORM)/LOG(2.0D0))+2 )
      SCALE = T / DBLE(2**NS)
      SCALE2 = SCALE*SCALE
*
*---  COMPUTE PADE COEFFICIENTS ...
*
      I = IDEG+1
      J = 2*IDEG+1
      WSP(ICOEF) = 1.0D0
      DO K = 1,IDEG
         WSP(ICOEF+K) = (WSP(ICOEF+K-1)*DBLE( I-K ))/DBLE( K*(J-K) )
      ENDDO
*
*---  H2 = SCALE2*H*H ...
*
      CALL DGEMM( 'N','N',M,M,M,SCALE2,H,LDH,H,LDH,0.0D0,WSP(IH2),M )
*
*---  INITIALIZE P (NUMERATOR) AND Q (DENOMINATOR) ...
*
      CP = WSP(ICOEF+IDEG-1)
      CQ = WSP(ICOEF+IDEG)
      DO J = 1,M
         DO I = 1,M
            WSP(IP + (J-1)*M + I-1) = 0.0D0
            WSP(IQ + (J-1)*M + I-1) = 0.0D0
         ENDDO
         WSP(IP + (J-1)*(M+1)) = CP
         WSP(IQ + (J-1)*(M+1)) = CQ
      ENDDO
*
*---  APPLY HORNER RULE ...
*
      IODD = 1
      K = IDEG - 1
 100  CONTINUE
      IUSED = IODD*IQ + (1-IODD)*IP
      CALL DGEMM( 'N','N',M,M,M, 1.0D0,WSP(IUSED),M,
     $             WSP(IH2),M, 0.0D0,WSP(IFREE),M )
      DO J = 1,M
         WSP(IFREE+(J-1)*(M+1)) = WSP(IFREE+(J-1)*(M+1))+WSP(ICOEF+K-1)
      ENDDO
      IP = (1-IODD)*IFREE + IODD*IP
      IQ = IODD*IFREE + (1-IODD)*IQ
      IFREE = IUSED
      IODD = 1-IODD
      K = K-1
      IF ( K.GT.0 )  GOTO 100
*
*---  OBTAIN (+/-)(I + 2*(P\Q)) ...
*
      IF ( IODD .EQ. 1 ) THEN
         CALL DGEMM( 'N','N',M,M,M, SCALE,WSP(IQ),M,
     $                H,LDH, 0.0D0,WSP(IFREE),M )
         IQ = IFREE
      ELSE
         CALL DGEMM( 'N','N',M,M,M, SCALE,WSP(IP),M,
     $                H,LDH, 0.0D0,WSP(IFREE),M )
         IP = IFREE
      END IF
      CALL DAXPY( MM, -1.0D0,WSP(IP),1, WSP(IQ),1 )
      CALL DGESV( M,M, WSP(IQ),M, IPIV, WSP(IP),M, IFLAG )
      IF (IFLAG.NE.0) RETURN
      CALL DSCAL( MM, 2.0D0, WSP(IP), 1 )
      DO J = 1,M
         WSP(IP+(J-1)*(M+1)) = WSP(IP+(J-1)*(M+1)) + 1.0D0
      ENDDO
      IPUT = IP
      IF ( NS.EQ.0 .AND. IODD.EQ.1 ) THEN
         CALL DSCAL( MM, -1.0D0, WSP(IP), 1 )
         GOTO 200
      END IF
*
*--   SQUARING : EXP(T*H) = (EXP(T*H))^(2^NS) ...
*
      IODD = 1
      DO K = 1,NS
         IGET = IODD*IP + (1-IODD)*IQ
         IPUT = (1-IODD)*IP + IODD*IQ
         CALL DGEMM( 'N','N',M,M,M, 1.0D0,WSP(IGET),M, WSP(IGET),M,
     $                0.0D0,WSP(IPUT),M )
         IODD = 1-IODD
      ENDDO
 200  CONTINUE
      IEXPH = IPUT
      END
*----------------------------------------------------------------------|
      SUBROUTINE MATADD(A, B, C, N, M)
C----------------------------------------------------------------------C
C     ADD MATRICES A(N,M) + B(N,M) => C(N,M)                           C
C----------------------------------------------------------------------C
      INTEGER          N, M
      DOUBLE PRECISION A(*), B(*), C(*)
      INTEGER          I
      DO 10 I = 1, N*M
 10   C(I) = A(I) + B(I)
      RETURN
      END
C
      SUBROUTINE MATSUB(A, B, C, N, M)
C----------------------------------------------------------------------C
C     SUBTRACT MATRICES A(N,M) - B(N,M) => C(N,M)                      C
C----------------------------------------------------------------------C
      INTEGER          N, M
      DOUBLE PRECISION A(*), B(*), C(*)
      INTEGER          I
      DO 10 I = 1, N*M
 10   C(I) = A(I) - B(I)
      RETURN
      END
C
      SUBROUTINE TF(X, XMIN, XMAX)
C----------------------------------------------------------------------C
C     TRANSFORMATION OF THE INTERVAL (XMIN,XMAX) -> R                  C
C----------------------------------------------------------------------C
      DOUBLE PRECISION X, XMIN, XMAX
      X = DLOG((X - XMIN)/(XMAX - X))
      RETURN
      END
C
      SUBROUTINE ATF(X, XMIN, XMAX)
C----------------------------------------------------------------------C
C     INVERSE TRANSFORMATION OF TF                                     C
C----------------------------------------------------------------------C
      DOUBLE PRECISION X, XMIN, XMAX, C
      C = DEXP(X)
      X = (XMIN + C*XMAX)/(1.0D0 + C)
      RETURN
      END
C
      SUBROUTINE DTF(X, XMIN, XMAX, Y)
C----------------------------------------------------------------------C
C     DERIVATIVE OF TF                                                 C
C----------------------------------------------------------------------C
      DOUBLE PRECISION X, XMIN, XMAX, Y
      Y = (XMAX - XMIN)/((XMAX - X)*(X - XMIN))
      RETURN
      END
C
      SUBROUTINE EXPAND(A, NA, N, L, D)
C----------------------------------------------------------------------C
C     EXPANDS AN L*D*L' FACTORED MATRIX STORED THE WAY MC11 FROM       C
C     THE HARWELL LIBRARY 1989 DOES, INTO AN L(N,N) MATRIX AND         C
C     A D(N) VECTOR.                                                   C
C----------------------------------------------------------------------C
      INTEGER          NA, N
      DOUBLE PRECISION A(NA), L(N, N), D(N)
      INTEGER          I, J, K
      K = 1
      DO 5 I = 1, N
         DO 6 J = I, N
            L(J, I) = A(K)
            K = K + 1
            IF (J .EQ. I) GOTO 6
            L(I, J) = 0.0D0
 6       CONTINUE
         D(I) = L(I, I)
         L(I, I) = 1.0D0
 5    CONTINUE
      RETURN
      END
C
      SUBROUTINE UPDATE(A, B, C, D, NA, NNA, NC, NNC, W, DW, LC, EPSM)
C----------------------------------------------------------------------C
C     UPDATES THE MATRIX  A = B*C*B' + D  , WHERE A(NNA) AND D(NNA)    C
C     OF ORDER NA ARE FACTORIZED AND C(NNC) OF ORDER NC IS FACTORIZED  C
C     AND B(NA,NC) IS A GENERAL MATRIX.                                C
C----------------------------------------------------------------------C
      INTEGER          NA, NNA, NC, NNC
      DOUBLE PRECISION A(NNA),B(NA,NC),C(NNC),D(NNA),W(NA*(NA+NC)),
     $                 DW(NA+NC),LC(NC*NC),EPSM
      CALL EXPAND(C,NNC,NC,LC,DW(1))
      CALL DGEMM('N','N',NA,NC,NC,1.0D0,B,NA,LC,NC,0.0D0,W,NA)
      CALL EXPAND(D,NNA,NA,W(NC*NA+1),DW(NC+1))
      CALL MWGS(W,NA,NA+NC,DW,A,EPSM)
      RETURN
      END
C
      SUBROUTINE MWGS(W, RW, CW, D, A, EPSM)
C----------------------------------------------------------------------C
C     PERFORM THE MODIFIED WEIGHTED GRAM-SCHMIDT ORTHOGONALIZATION,    C
C     AND MATRIX FACTORIZATION. A MODIFIED VERSION OF THEOREM          C
C     VI.4.1 IN BIERMAN, G.J.: FACTORIZATION METHODS FOR DISCRETE      C
C     SEQUENTIAL ESTIMATION, 1977.                                     C
C                                                                      C
C     A SYSTEM A=W*D*W', WHERE W(RW,CW) IS A GENERAL MATRIX, AND D(CW) C
C     IS A VECTOR REPRESENTING A DIAGONAL MATRIX, IS TRANSFORMED TO    C
C     A=LS*DS*LS', WHERE LS IS A UNIT LOWER TRIANGULAR MATRIX OF       C
C     DIMENSION RW, AND DS IS A DIAGONAL MATRIX.                       C
C                                                                      C
C     ON OUTPUT LS AND DS ARE STORED IN A IN THE SAME ORDER AS IN      C
C     MC11 FROM THE HARWELL LIBRARY, 1989.                             C
C----------------------------------------------------------------------C
      INTEGER          RW, CW
      DOUBLE PRECISION W(RW,CW), D(CW), A(RW*(RW+1)/2), EPSM
      INTEGER          I, J, K, NI, NJ, IV
      DOUBLE PRECISION EPS, AD, DELT, AMIN
      EPS = 1D1*EPSM
      AMIN = 1D-50
      NI = 1
      DO 10  J = 2, RW
 5       A(NI) = 0.0D0
         DO 20  I = 1, CW
            IV = J-1
            A(NI) = A(NI) + D(I)*W(IV, I)**2
 20      CONTINUE
         IF(A(NI).LE.0.0D0) THEN
            AD = 0.0D0
            DO 50 I = 1,CW
               AD = AD + DABS(D(I))*W(IV, I)**2
 50         CONTINUE
            DELT = (AMIN - A(NI))/AD
            DELT = DMAX1(DELT,EPS)
            DO 55 I = 1, CW
               D(I) = D(I) + DELT*DABS(D(I))
 55         CONTINUE
            GOTO 5
         END IF
         NJ = NI + 1
         DO 25  K = J, RW
            A(NJ) = 0.0D0
            DO 30  I = 1, CW
               A(NJ) = A(NJ) + W(K, I)*D(I)*W(J-1, I)
 30         CONTINUE
            A(NJ) = A(NJ)/A(NI)
            DO 35  I = 1, CW
               W(K,  I) = W(K, I) - A(NJ )*W(J-1, I)
 35         CONTINUE
            NJ = NJ + 1
 25      CONTINUE
         NI = NI + RW - J + 2
 10   CONTINUE
      A(NI) = 0.0D0
      DO 40  I = 1, CW
         IV = RW
         A(NI) = A(NI) + D(I)*W(IV, I)**2
 40   CONTINUE
      IF(A(NI).LE.0.0D0) A(NI) = AMIN
      RETURN
      END
C
      SUBROUTINE TDIST(T, DF, PT)
C----------------------------------------------------------------------C
C     CALCULATES CUMULATIVE PROBABILITIES 'PT' IN STUDENT'S            C
C     T-DISTRIBUTION WITH 'DF' DEGREES OF FREEDOM FOR THE VALUE 'T'.   C
C                                                                      C
C     THIS APPROXIMATION IS BEST FOR LARGE DF.                         C
C----------------------------------------------------------------------C
      INTEGER          DF
      DOUBLE PRECISION T, PT
      DOUBLE PRECISION X
      X = T*(1.0D0 - 1.0D0/(4D0*DF))/DSQRT(1.0D0 + T*T/(2D0*DF))
      CALL NDIST(X, PT)
      RETURN
      END
C
      SUBROUTINE NDIST(X, PX)
C----------------------------------------------------------------------C
C     CALCULATES CUMULATIVE PROBABILITIES 'PX' IN THE STANDARD         C
C     GAUSSIAN DISTRIBUTION FOR THE VALUE 'X'.                         C
C                                                                      C
C     THE APPROXIMATION HAS THE ERROR: <7.5D-8.                        C
C----------------------------------------------------------------------C
      DOUBLE PRECISION X, PX
      DOUBLE PRECISION PI, P, B1, B2, B3, B4, B5, ZX, T
      PI =  DACOS(-1.0D0)
      P  =  0.2316419D0
      B1 =  0.319381530D0
      B2 = -0.356563782D0
      B3 =  1.781477937D0
      B4 = -1.821255978D0
      B5 =  1.330274429D0
      ZX = (1.0D0/DSQRT(2D0*PI))*DEXP(-X*X/2D0)
      T = 1.0D0/(1.0D0 + P*X)
      IF(X .LT. 0.0D0) T = 1.0D0/(1.0D0 - P*X)
      PX = 1.0D0 - ZX*T*(B1 + T*(B2 + T*(B3 + T*(B4 + T*B5))))
      IF(X .LT. 0.0D0) PX = 1.0D0 - PX
      RETURN
      END
C
      SUBROUTINE MC11AD(A,NA,N,Z,SIG,W,IR,MK,EPS)
C----------------------------------------------------------------------C
C     THIS IS A MODIFIED VERSION OF THE SUBROUTINE WITH THE SAME NAME  C
C     FROM THE HARWELL LIBRARY.                                        C
C----------------------------------------------------------------------C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NA),Z(N),W(N)
C     UPDATE FACTORS GIVEN IN A BY SIG*Z*ZTRANSPOSE
      IF (N .GT. 1) GOTO 1
      A(1)=A(1)+SIG *Z(1)**2
      IR=1
      IF (A(1) .GT. 0.0D0) RETURN
      A(1)=0.0D0
      IR=0
      RETURN
    1 CONTINUE
      NP=N+1
      IF (SIG .GT. 0.0D0) GOTO 40
      IF (SIG .EQ. 0.0D0  .OR.  IR .EQ. 0) RETURN
      TI=1.0D0/SIG
      IJ=1
      IF (MK .EQ. 0) GOTO 10
      DO 7 I=1,N
        IF (A(IJ) .NE. 0.0D0) TI=TI+W(I)**2/A(IJ)
    7   IJ=IJ+NP-I
      GOTO 20
   10 CONTINUE
      DO 11 I=1,N
   11   W(I)=Z(I)
      DO 15 I=1,N
        IP=I+1
        V=W(I)
        IF (A(IJ) .GT. 0.0D0) GOTO 12
        W(I)=0.0D0
        IJ=IJ+NP-I
        GOTO 15
   12   CONTINUE
        TI=TI+V**2/A(IJ)
        IF (I .EQ. N) GOTO 14
        DO 13 J=IP,N
          IJ=IJ+1
   13     W(J)=W(J)-V*A(IJ)
   14   IJ=IJ+1
   15 CONTINUE
   20 CONTINUE
      IF (IR .LE. 0) GOTO 21
      IF (TI .GT. 0.0D0) GOTO 22
      IF (MK-1) 40,40,23
   21 TI=0.0D0
      IR=-IR-1
      GOTO 23
   22 TI=EPS/SIG
      IF (EPS .EQ. 0.0D0) IR=IR-1
   23 CONTINUE
      MM=1
      TIM=TI
      DO 30 I=1,N
        J=NP-I
        IJ=IJ-I
        IF (A(IJ) .NE. 0.0D0) TIM=TI-W(J)**2/A(IJ)
        W(J)=TI
   30   TI=TIM
      GOTO 41
   40 CONTINUE
      MM=0
      TIM=1.0D0/SIG
   41 CONTINUE
      IJ=1
      DO 66 I=1,N
        IP=I+1
        V=Z(I)
        IF (A(IJ) .GT. 0.0D0) GOTO 53
        IF (IR .GT. 0  .OR.  SIG .LT. 0.0D0  .OR.  V .EQ. 0.0D0) GOTO 52
        IR=1-IR
        A(IJ)=V**2/TIM
        IF (I .EQ. N) RETURN
        DO 51 J=IP,N
          IJ=IJ+1
   51     A(IJ)=Z(J)/V
        RETURN
   52   CONTINUE
        TI=TIM
        IJ=IJ+NP-I
        GOTO 66
   53   CONTINUE
        AL=V/A(IJ)
        IF (MM) 54,54,55
   54   TI=TIM+V*AL
        GOTO 56
   55   TI=W(I)
   56   CONTINUE
        R=TI/TIM
        A(IJ)=A(IJ)*R
        IF (R .EQ. 0.0D0) GOTO 70
        IF (I .EQ. N) GOTO 70
        B=AL/TI
        IF (R .GT. 4D0) GOTO 62
        DO 61 J=IP,N
          IJ=IJ+1
          Z(J)=Z(J)-V*A(IJ)
   61     A(IJ)=A(IJ)+B*Z(J)
        GOTO 64
   62   GM=TIM/TI
        DO 63 J=IP,N
          IJ=IJ+1
          Y=A(IJ)
          A(IJ)=B*Z(J)+Y*GM
   63     Z(J)=Z(J)-V*Y
   64   CONTINUE
        TIM=TI
        IJ=IJ+1
   66 CONTINUE
   70 CONTINUE
      IF (IR .LT. 0) IR=-IR
      RETURN
C     FACTORIZE A MATRIX GIVEN IN A
      ENTRY MC11BD(A,NA,N,IR)
      IR=N
      IF (N .GT. 1) GOTO 100
      IF (A(1) .GT. 0.0D0) RETURN
      A(1)=0.0D0
      IR=0
      RETURN
  100 CONTINUE
      NP=N+1
      II=1
      DO 104 I=2,N
        AA=A(II)
        NI=II+NP-I
        IF (AA .GT. 0.0D0) GOTO 101
        A(II)=0.0D0
        IR=IR-1
        II=NI+1
        GOTO 104
  101   CONTINUE
        IP=II+1
        II=NI+1
        JK=II
        DO 103 IJ=IP,NI
          V=A(IJ)/AA
          DO 102 IK=IJ,NI
            A(JK)=A(JK)-A(IK)*V
  102       JK=JK+1
  103     A(IJ)=V
  104 CONTINUE
      IF (A(II) .GT. 0.0D0) RETURN
      A(II)=0.0D0
      IR=IR-1
      RETURN
C     MULTIPLY OUT THE FACTORS GIVEN IN A
      ENTRY MC11CD(A,NA,N)
      IF (N .EQ. 1) RETURN
      NP=N+1
      II=N*NP/2
      DO 202 NIP=2,N
        JK=II
        NI=II-1
        II=II-NIP
        AA=A(II)
        IP=II+1
        IF (AA .GT. 0.0D0) GOTO 203
        DO 204 IJ=IP,NI
  204     A(IJ)=0.0D0
        GOTO 202
  203   CONTINUE
        DO 201 IJ=IP,NI
          V=A(IJ)*AA
          DO 200 IK=IJ,NI
            A(JK)=A(JK)+A(IK)*V
  200     JK=JK+1
  201     A(IJ)=V
  202 CONTINUE
      RETURN
C     MULTIPLY A VECTOR Z BY THE FACTORS GIVEN IN A
      ENTRY MC11DD(A,NA,N,Z,W)
      IF (N .GT. 1) GOTO 300
      Z(1)=Z(1)*A(1)
      W(1)=Z(1)
      RETURN
  300 CONTINUE
      NP=N+1
      II=1
      N1=N-1
      DO 303 I=1,N1
        Y=Z(I)
        IF (A(II) .EQ. 0.0D0) GOTO 302
        IJ=II
        IP=I+1
        DO 301 J=IP,N
          IJ=IJ+1
  301     Y=Y+Z(J)*A(IJ)
  302   Z(I)=Y*A(II)
        W(I)=Z(I)
  303   II=II+NP-I
      Z(N)=Z(N)*A(II)
      W(N)=Z(N)
      DO 311 K=1,N1
        I=N-K
        II=II-NP+I
        IF (Z(I) .EQ. 0.0D0) GOTO 311
        IP=I+1
        IJ=II
        Y=Z(I)
        DO 310 J=IP,N
          IJ=IJ+1
  310     Z(J)=Z(J)+A(IJ)*Z(I)
  311 CONTINUE
      RETURN
C     MULTIPLY A VECTOR Z BY THE INVERSE OF THE FACTORS GIVEN IN A
      ENTRY MC11ED(A,NA,N,Z,W,IR)
      IF (IR .LT. N) RETURN
      W(1)=Z(1)
      IF (N .GT. 1) GOTO 400
      Z(1)=Z(1)/A(1)
      RETURN
  400 CONTINUE
      DO 402 I=2,N
        IJ=I
        I1=I-1
        V=Z(I)
        DO 401 J=1,I1
          V=V-A(IJ)*Z(J)
  401     IJ=IJ+N-J
        W(I)=V
  402   Z(I)=V
      Z(N)=Z(N)/A(IJ)
      NP=N+1
      DO 411 NIP=2,N
        I=NP-NIP
        II=IJ-NIP
        V=Z(I)/A(II)
        IP=I+1
        IJ=II
        DO 410 J=IP,N
          II=II+1
  410     V=V-A(II)*Z(J)
  411   Z(I)=V
      RETURN
C     COMPUTE THE INVERSE MATRIX FROM FACTORS GIVEN IN A
      ENTRY MC11FD(A,NA,N,IR)
      IF (IR .LT. N) RETURN
      A(1)=1.0D0/A(1)
      IF (N .EQ. 1) RETURN
      NP=N+1
      N1=N-1
      II=2
      DO 511 I=2,N
        A(II)=-A(II)
        IJ=II+1
        IF (I .EQ. N) GOTO 502
        DO 501 J=I,N1
          IK=II
          JK=IJ
          V=A(IJ)
          DO 500 K=I,J
            JK=JK+NP-K
            V=V+A(IK)*A(JK)
  500       IK=IK+1
          A(IJ)=-V
  501     IJ=IJ+1
  502   CONTINUE
        A(IJ)=1.0D0/A(IJ)
        II=IJ+1
        AA=A(IJ)
        IJ=I
        IP=I+1
        NI=N-I
        DO 511 J=2,I
          V=A(IJ)*AA
          IK=IJ
          K=IJ-IP+J
          I1=IJ-1
          NIP=NI+IJ
          DO 510 JK=K,I1
            A(JK)=A(JK)+V*A(IK)
  510       IK=IK+NIP-JK
          A(IJ)=V
  511   IJ=IJ+NP-J
      RETURN
      END

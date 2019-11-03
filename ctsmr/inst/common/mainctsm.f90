      SUBROUTINE MAIN(INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,NIMAT,OMAT,&
     &                NOMAT,NOBS,NSET,TRACE,INFO)
!----------------------------------------------------------------------C
!     ESTIMATION OF PARAMETERS IN STOCHASTIC DIFFERENTIAL EQUATIONS.   C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,NSET,&
     &                 NOBS(NSET),INFO
      DOUBLE PRECISION DOUBLS(NDOUB),TMAT(NTMAT),IMAT(NIMAT),OMAT(NOMAT)
      INTEGER          DF,I,ICONTR,ID(NPARAM),IR,ITR,J,K,MAXFUN,NA,&
     &                 NEVAL,NG,NMISST,NPR,NPRM,NX,NDATA,VINFO(NPARAM)
      DOUBLE PRECISION C1,C2,D(NPARAM,NPARAM),DET,DL(NPARAM),&
     &                 DPEN(NPARAM),EPS,EPSM,ETA,F,FPEN,FPRIOR,&
     &                 H(NPARAM,NPARAM),HK((NPARAM*(NPARAM+1))/2),&
     &                 LAMBDA,PRDET,PRCOVI(NPARAM*NPARAM),&
     &                 PRMEAN(NPARAM),PT(NPARAM),SD(NPARAM),T(NPARAM),&
     &                 TMP(NPARAM,NPARAM),TMP1(NPARAM),&
     &                 TMP2(NPARAM*(NPARAM+1)/2),&
     &                 WK(NPARAM*(NPARAM+15)/2+1),XMM,X(NPARAM),&
     &                 XM(NPARAM),XMAX(NPARAM),XMAXS(NPARAM),&
     &                 XMIN(NPARAM),XMINS(NPARAM),TRACE(*)
      DOUBLE PRECISION DLAMCH
      EXTERNAL         DLAMCH
!     ..................................................................
!                                           UNPACKING INTEGER ARGUMENTS.
!
      MAXFUN=INTS(1)
      ITR=INTS(2)
      NEVAL=INTS(3)
      NPR=INTS(9)
      DO 1 I=1,NPARAM
         ID(I)=INTS(10+I)
   1  CONTINUE
!     ..................................................................
!                                            UNPACKING DOUBLE ARGUMENTS.
!
      EPS=DOUBLS(1)
      ETA=DOUBLS(2)
      DET=DOUBLS(3)
      F=DOUBLS(4)
      DO 3 I=1,NPARAM
         DL(I)=DOUBLS(4+I)
   3  CONTINUE
      DO 4 I=1,NPARAM
         DPEN(I)=DOUBLS(NPARAM+4+I)
   4  CONTINUE
      DO 5 I=1,NPARAM
      DO 5 J=1,NPARAM
         H(J,I)=DOUBLS(2*NPARAM+4+(I-1)*NPARAM+J)
   5  CONTINUE
      DO 6 I=1,NPARAM
         PT(I)=DOUBLS((NPARAM+2)*NPARAM+4+I)
   6  CONTINUE
      DO 7 I=1,NPARAM
         SD(I)=DOUBLS((NPARAM+3)*NPARAM+4+I)
   7  CONTINUE
      DO 8 I=1,NPARAM
         T(I)=DOUBLS((NPARAM+4)*NPARAM+4+I)
   8  CONTINUE
      FPEN=DOUBLS((NPARAM+5)*NPARAM+5)
      FPRIOR=DOUBLS((NPARAM+5)*NPARAM+6)
      LAMBDA=DOUBLS((NPARAM+5)*NPARAM+8)
      XMM=DOUBLS((NPARAM+5)*NPARAM+10)
      DO 11 I=1,NPARAM*NPARAM
         PRCOVI(I)=DOUBLS((NPARAM+5)*NPARAM+11+I)
  11  CONTINUE
      PRDET=DOUBLS((2*NPARAM+5)*NPARAM+12)
      DO 12 I=1,NPARAM
         PRMEAN(I)=DOUBLS((2*NPARAM+5)*NPARAM+12+I)
  12  CONTINUE
      DO 13 I=1,NPARAM
         XM(I)=DOUBLS((2*NPARAM+6)*NPARAM+12+I)
  13  CONTINUE
      DO 14 I=1,NPARAM
         XMAX(I)=DOUBLS((2*NPARAM+7)*NPARAM+12+I)
  14  CONTINUE
      DO 15 I=1,NPARAM
         XMIN(I)=DOUBLS((2*NPARAM+8)*NPARAM+12+I)
  15  CONTINUE
!     ..................................................................
!                 SETTING INITIAL VALUE OF INFO AND WORKING ARRAY VINFO.
!
      INFO=0
      DO 17 I=1,NPARAM
         VINFO(I)=0
  17  CONTINUE
!     ..................................................................
!                                 CALCULATION OF X, XMAXS, XMINS AND NX.
!
      DO 18 I=1,NPARAM
         X(I)=0.0D0
         XMINS(I)=0.0D0
         XMAXS(I)=0.0D0
  18  CONTINUE
      NX=0
      DO 20 I=1,NPARAM
         IF (ID(I).NE.0) THEN
            NX=NX+1
            X(NX)=XM(I)
            XMINS(NX)=XMIN(I)
            XMAXS(NX)=XMAX(I)
         END IF
  20  CONTINUE
      NA=(NX*(NX+1))/2
!     ..................................................................
!           CALCULATION OF THE INVERTED PRIOR COVARIANCE AND PRIOR MEAN.
!
      DO 30 I=1,NPR
      DO 30 J=1,NPR
         PRCOVI((I-1)*NPR+J)=PRCOVI((I-1)*NPR+J)*PRMEAN(I)*PRMEAN(J)
  30  CONTINUE
      IF (NPR.GT.0) THEN
         NPRM=NPR*(NPR+1)/2
         DO 31 I=1,(NPARAM*(NPARAM+1)/2)
            TMP2(I)=0.0D0
  31     CONTINUE
         K=1
         DO 32 I=1,NPR
         DO 32 J=I,NPR
            TMP2(K)=PRCOVI((I-1)*NPR+J)
            K=K+1
  32     CONTINUE
         CALL MC11BD(TMP2,NPRM,NPR,IR)
         IF (IR.NE.NPR) THEN
            INFO=5
            RETURN
         END IF
         CALL EXPAND(TMP2,NPRM,NPR,PRCOVI,PRMEAN)
         PRDET=1.0D0
         DO 33 I=1,NPR
            PRDET=PRDET*PRMEAN(I)
  33     CONTINUE
         CALL MC11FD(TMP2,NPRM,NPR,NPR)
         K=1
         DO 34 I=1,NPR
         DO 34 J=I,NPR
            PRCOVI((I-1)*NPR+J)=TMP2(K)
            PRCOVI((J-1)*NPR+I)=PRCOVI((I-1)*NPR+J)
            K=K+1
  34     CONTINUE
         J=0
         DO 35 I=1, NPARAM
            IF (ID(I).EQ.2) THEN
               J=J+1
               PRMEAN(J)=XM(I)
            END IF
  35     CONTINUE
      END IF
!     ..................................................................
!                        PACKING PRCOVI, PRDET, PRMEAN, XMAXS AND XMINS.
!
      DO 36 I=1,NPARAM*NPARAM
         DOUBLS((NPARAM+5)*NPARAM+11+I)=PRCOVI(I)
  36  CONTINUE
      DOUBLS((2*NPARAM+5)*NPARAM+12)=PRDET
      DO 37 I=1,NPARAM
         DOUBLS((2*NPARAM+5)*NPARAM+12+I)=PRMEAN(I)
  37  CONTINUE
      DO 38 I=1,NPARAM
         DOUBLS((2*NPARAM+7)*NPARAM+12+I)=XMAXS(I)
  38  CONTINUE
      DO 39 I=1,NPARAM
         DOUBLS((2*NPARAM+8)*NPARAM+12+I)=XMINS(I)
  39  CONTINUE
!     ..................................................................
!                                  TRANSFORMATION OF INITIAL PARAMETERS.
!
      DO 40 I=1,NX
         CALL TF(X(I),XMINS(I),XMAXS(I))
  40  CONTINUE
!     ..................................................................
!                                                   INITIALIZING ARRAYS.
!
      DO 41 I=1,(NPARAM*(NPARAM+15)/2+1)
         WK(I)=0.0D0
  41  CONTINUE
      DO 42 I=1,(NPARAM*(NPARAM+1)/2)
         HK(I)=0.0D0
  42  CONTINUE
      DO 43 I=1,NPARAM
      DO 43 J=1,NPARAM
         D(J,I)=0.0D0
  43  CONTINUE
!     ..................................................................
!                                                          OPTIMISATION.
!
      ICONTR=0
      NMISST=0
      EPSM=DLAMCH('E')
      CALL OPTIM(NX,X,XMINS,XMAXS,EPS,ETA,MAXFUN,ITR,WK,&
     &           NPARAM*(NPARAM+15)/2+1,ICONTR,INTS,NINT,DOUBLS,NDOUB,&
     &           TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,FPEN,&
     &           FPRIOR,EPSM,TRACE,VINFO,INFO)
      IF (INFO.NE.0) RETURN
      IF (ICONTR.NE.0) THEN
         INFO=ICONTR
         RETURN
      END IF
!     ..................................................................
!                                                UNPACKING THE SOLUTION.
!
      F=WK(1)
      DO 45 I=1,NA
         HK(I)=WK(I+1)
  45  CONTINUE
      CALL MC11FD(WK(2),NA,NX,NX)
      K = 2
      DO 50 I=1,NX
      DO 50 J=I,NX
         D(I,J)=WK(K)
         D(J,I)=D(I,J)
         K=K+1
  50  CONTINUE
!     ..................................................................
!                                   BACK-TRANSFORMATION OF THE SOLUTION.
!
      DO 55 I=1,NX
         CALL ATF(X(I),XMINS(I),XMAXS(I))
  55  CONTINUE
!     ..................................................................
!                                                     CALCULATION OF XM.
!
      NX=0
      DO 60 I=1,NPARAM
         IF (ID(I).NE.0) THEN
            NX=NX+1
            XM(I)=X(NX)
         END IF
  60  CONTINUE
!     ..................................................................
!                            BACK-TRANSFORMATION OF THE VARIANCE MATRIX.
!
      C1=1.0D0
      C2=1.0D0
      DO 65 J=1,NX
         CALL DTF(X(J),XMINS(J),XMAXS(J),C2)
      DO 65 I=1,NX
         CALL DTF(X(I),XMINS(I),XMAXS(I),C1)
         D(I,J)=D(I,J)/(C1*C2)
  65  CONTINUE
!     ..................................................................
!                         CALCULATION OF THE DETERMINANT OF THE HESSIAN.
!
      K=1
      DET=1.0D0
      DO 70 I=1,NX
         CALL DTF(X(I),XMINS(I),XMAXS(I),C1)
         DET=DET*(C1**2)*HK(K)
         K=K+NX-I+1
  70  CONTINUE
      IF (DET.LE.0.0D0) THEN
         DET=0.0D0
      ELSE
         DET=-DLOG(DET)
      END IF
!     ..................................................................
!                        CALCULATION OF THE STANDARD DEVIATIONS AND T'S.
!
      NDATA=0
      DO 75 I=1,NSET
         NDATA=NDATA+NOBS(I)*S
  75  CONTINUE
      DF=NDATA-NMISST-NX
      NG=1+(NX*(NX+13))/2
      DO 80 I=1,NX
         DL(I)=0.0D0
         IF (D(I,I).LE.0.0D0) D(I,I)=1D-100
         SD(I)=DSQRT(D(I,I))
         T(I)=X(I)/SD(I)
         CALL TDIST(DABS(T(I)),DF,PT(I))
         PT(I)=2.0D0*(1.0D0-PT(I))
         CALL DTF(X(I),XMINS(I),XMAXS(I),DL(I))
         DL(I)=X(I)*WK(I+NG)*DL(I)
         DPEN(I)=0.0D0
         C1=DABS(XMINS(I))
         C2=DABS(XMAXS(I))
         IF (C1.LT.XMM) C1=XMM
         IF (C2.LT.XMM) C2=XMM
         DPEN(I)=-C1/(X(I)-XMINS(I))**2+C2/(XMAXS(I)-X(I))**2
         DPEN(I)=LAMBDA*DPEN(I)*X(I)
  80  CONTINUE
!     ..................................................................
!                                 CALCULATION OF THE CORRELATION MATRIX.
!
      DO 90 I=1,NX
      DO 90 J=1,NX
         H(J,I)=D(J,I)/(DSQRT(D(I,I)*D(J,J)))
  90  CONTINUE
!     ..................................................................
!                                               PACKING INTEGER RESULTS.
!
      INTS(2)=ITR
      NEVAL=MAXFUN
      INTS(3)=NEVAL
!     ..................................................................
!                                                PACKING DOUBLE RESULTS.
!
      DOUBLS(3)=DET
      DOUBLS(4)=F
      DO 903 I=1,NPARAM
         DOUBLS(4+I)=DL(I)
 903  CONTINUE
      DO 904 I=1,NPARAM
         DOUBLS(NPARAM+4+I)=DPEN(I)
 904  CONTINUE
      DO 905 I=1,NPARAM
      DO 905 J=1,NPARAM
         DOUBLS(2*NPARAM+4+(I-1)*NPARAM+J)=H(J,I)
 905  CONTINUE
      DO 906 I=1,NPARAM
         DOUBLS((NPARAM+2)*NPARAM+4+I)=PT(I)
 906  CONTINUE
      DO 907 I=1,NPARAM
         DOUBLS((NPARAM+3)*NPARAM+4+I)=SD(I)
 907  CONTINUE
      DO 908 I=1,NPARAM
         DOUBLS((NPARAM+4)*NPARAM+4+I)=T(I)
 908  CONTINUE
      DOUBLS((NPARAM+5)*NPARAM+5)=FPEN
      DOUBLS((NPARAM+5)*NPARAM+6)=FPRIOR
      DO 913 I=1,NPARAM
         DOUBLS((2*NPARAM+6)*NPARAM+12+I)=XM(I)
 913  CONTINUE
      DO 914 I=1,NPARAM
         DOUBLS((2*NPARAM+7)*NPARAM+12+I)=XMAX(I)
 914  CONTINUE
      DO 915 I=1,NPARAM
         DOUBLS((2*NPARAM+8)*NPARAM+12+I)=XMIN(I)
 915  CONTINUE
      RETURN
      END
!
      SUBROUTINE LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,&
     &                 NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,FPEN,FPRIOR,&
     &                 EPSM,JOB,FUNK,INFO,VX0)
!----------------------------------------------------------------------C
!     CALCULATION OF LIKELIHOOD FUNCTION.                              C
!----------------------------------------------------------------------C

      use covariance


      INCLUDE          'global.h'
      INTEGER          NX,NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,NSET,&
     &                 NOBS(NSET),NMISST,JOB,INFO
      DOUBLE PRECISION X(NX),DOUBLS(NDOUB),TMAT(NTMAT),IMAT(NIMAT),&
     &                 OMAT(NOMAT),EPSM,FUNK
      INTEGER          CTSFLG,HN,I,ID(NPARAM),IDEG,II,IR,J,JJ,K,KK,LL,&
     &                 MISS,MV,NDATA,NIEKF,NN,NONLIN,NPR,NSEKF,PRED,&
     &                 SVDFLG
      DOUBLE PRECISION A(N,N),ALPHA(M+1),ALPHAS(M+1),B(N,M+1),C(S,N),C1,&
     &                 C2,D(S,M+1),DELT,DET,DSIGMA(N),EKFEPS,EPS(S),&
     &                 F(N),FPEN,FPRIOR,GAM(N,N),GAMH(N,N),H(S),HC,&
     &                 LAMBDA,LOGLIK,ODEEPS,PHI(N,N),PPDF,&
     &                 PRCOVI(NPARAM*NPARAM),PRDET,PRMEAN(NPARAM),&
     &                 R1M(N,N),R1S(N*(N+1)/2),R1SM(N,N),R2(S*(S+1)/2),&
     &                 R2M(S,S),SVDEPS,T,TF,TMP(NPARAM),TMP1(NPARAM),&
     &                 TMP2(S,S),TMP3(S),TMP4,TMP5(N),TMP6(S),TS,TSAMP,&
     &                 TSAMPS,U(N,N),UO(M+1),UOS(M+1),UP(M+1),UPS(M+1),&
     &                 VXO(N*(N+1)/2),VXP(N*(N+1)/2),VXP0S,&
     &                 VXPS(N*(N+1)/2),VYP(S*(S+1)/2),VYPS(S*(S+1)/2),&
     &                 XCUR(N),XM(NPARAM),XMAX(NPARAM),XMIN(NPARAM),XMM,&
     &                 XO(N),XP(N),XP0(N),XPREV(N),XPS(N),XT(NPARAM),&
     &                 Y(S),YP(S)
      DOUBLE PRECISION DLAMCH
      EXTERNAL         DLAMCH
      COMMON           /COM/XM,ALPHA,UO,T,HN
!$OMP THREADPRIVATE(/COM/)
      COMMON           /BWCOM/ODEEPS,TF,TSAMP,XCUR,XPREV,NONLIN
!$OMP THREADPRIVATE(/BWCOM/)
      INTEGER          OUTOPT,SCOV,OCOV,CURPOS,ESTVAR

      double precision, optional, dimension(N*(N+1)/2), intent(in) :: VX0
!
!  RJ 2011 : Check if the user requested an interrupt 
      CALL RCHKUSR
      
      OUTOPT = 0
      SCOV = 0
      OCOV = 0

      call clear_cov_counter()
      
      IF (JOB.NE.0) THEN
!heck output requirements
         IF (OMAT(1).EQ.0D0) THEN
!           Output states, outputs and standard deviances
         ELSE IF (OMAT(1).EQ.1D0) THEN
!           Output states + full covariance and outputs + sd
            OUTOPT = 1
!            SCOV = 1
            SCOV = N*(N-1)/2
         ELSE IF (OMAT(1).EQ.2D0) THEN
!           Output states + sd and outputs + full covariance
            OUTOPT = 2
!            OCOV = 1
            OCOV = S*(S-1)/2
         ELSE IF (OMAT(1).EQ.3D0) THEN
!           Output states, outputs and full covariances
            OUTOPT = 3
!            SCOV = 1
            SCOV = N*(N-1)/2
!            OCOV = 1
            OCOV = S*(S-1)/2
         END IF
         OMAT(1) = 0D0
      END IF
      
!
!     ..................................................................
!                                                SETTING INITIAL VALUES.
!
      INFO=0
      FUNK=0.0D0
      NN=N*(N+1)/2
      IF (EPSM.LT.0.0D0) THEN
         EPSM=DLAMCH('E')
      END IF
!     ..................................................................
!                                           UNPACKING INTEGER ARGUMENTS.
!
      CTSFLG=INTS(4)
      HN=INTS(5)
      NONLIN=INTS(6)
      IDEG=INTS(7)
      NSEKF=INTS(8)
      NPR=INTS(9)

!      if (PRESENT(VX0)) then
!         ESTVAR = 2
!      else
         ESTVAR=INTS(10)
!      endif

      DO 1 I=1,NPARAM
         ID(I)=INTS(10+I)
   1  CONTINUE
      NIEKF=INTS(NINT)
!     ..................................................................
!                                            UNPACKING DOUBLE ARGUMENTS.
!
      HC=DOUBLS((NPARAM+5)*NPARAM+7)
      LAMBDA=DOUBLS((NPARAM+5)*NPARAM+8)
      SVDEPS=DOUBLS((NPARAM+5)*NPARAM+9)
      XMM=DOUBLS((NPARAM+5)*NPARAM+10)
      VXP0S=DOUBLS((NPARAM+5)*NPARAM+11)
      DO 11 I=1,NPARAM*NPARAM
         PRCOVI(I)=DOUBLS((NPARAM+5)*NPARAM+11+I)
  11  CONTINUE
      PRDET=DOUBLS((2*NPARAM+5)*NPARAM+12)
      DO 12 I=1,NPARAM
         PRMEAN(I)=DOUBLS((2*NPARAM+5)*NPARAM+12+I)
  12  CONTINUE
      DO 13 I=1,NPARAM
         XM(I)=DOUBLS((2*NPARAM+6)*NPARAM+12+I)
  13  CONTINUE
      DO 14 I=1,NPARAM
         XMAX(I)=DOUBLS((2*NPARAM+7)*NPARAM+12+I)
  14  CONTINUE
      DO 15 I=1,NPARAM
         XMIN(I)=DOUBLS((2*NPARAM+8)*NPARAM+12+I)
  15  CONTINUE
      EKFEPS=DOUBLS(NDOUB-1)
      ODEEPS=DOUBLS(NDOUB)
!     ..................................................................
!                          SPECIFICATION OF PARAMETERS OF THE ALGORITHM.
!
      SVDFLG=-1
!     ..................................................................
!                                     BACK-TRANSFORMATION OF PARAMETERS.
!
      IF (JOB.EQ.0) THEN
         DO 20 I=1,NPARAM
            XT(I)=0.0D0
  20     CONTINUE
         DO 30 I=1,NX
            XT(I)=X(I)
            CALL ATF(X(I),XMIN(I),XMAX(I))
  30     CONTINUE
         J=1
         DO 31 I=1,NPARAM
            IF (ID(I).NE.0) THEN
               XM(I)=X(J)
               J=J+1
            END IF
  31     CONTINUE
      END IF
!     ..................................................................
!                                                LOOP FOR SEVERAL FILES.
!
      PPDF=0.0D0
      NDATA=0
      NMISST=0
      DO 32 K=1,NSET
!     ..................................................................
!                         INITIALIZATIONS BEFORE LIKELIHOOD CALCULATION.
!
      LOGLIK=0.0D0
      MV=0
      T=0.0D0
      DO 35 II=1,M+1
         UO(II)=0.0D0
         UP(II)=0.0D0
         ALPHA(II)=0.0D0
  35  CONTINUE
      DO 36 II=1,N
         XO(II)=0.0D0
         XP(II)=0.0D0
         XP0(II)=0.0D0
  36  CONTINUE
!     ..................................................................
!                            SPECIFICATION OF THE ACTUAL DYNAMIC SYSTEM.
!
      TSAMP=(TMAT(NDATA+2)-TMAT(NDATA+1))/DBLE(NSEKF)
      IF (NONLIN.EQ.0) THEN
         CALL DEFSYS(XM,XP0,XP0,UO,T,A,B,C,D,R1M,R2M,F,H,NONLIN,1)
      ELSE
         T=TMAT(NDATA+1)
         DO 38 II=1,M
            UO(II)=IMAT(NDATA*(M+S)+II)
  38     CONTINUE
         CALL DEFSYS(XM,XP0,XP0,UO,T,A,B,C,D,R1M,R2M,F,H,NONLIN,0)
         CALL DEFSYS(XM,XP0,XP0,UO,T,A,B,C,D,R1M,R2M,F,H,NONLIN,3)
      END IF
      DO 40 II=1,N
         XP(II)=XP0(II)
  40  CONTINUE
!     ..................................................................
!                                   CALCULATION OF THE INITIAL VARIANCE.
!
      IF (ESTVAR.EQ.0) THEN
!        Use the integral solution.
         IF (NONLIN.LE.2) THEN
            CALL EXPA(A,R1M,TSAMP,PHI,R1SM,IDEG,INFO)
            IF (INFO.NE.0) RETURN
            DELT=0.0D0
            LL=1
  45        KK=1
            DO 50 JJ=1,N
            DO 50 II=JJ,N
               R1S(KK)=R1SM(II,JJ)
               IF (II.EQ.JJ) R1S(KK)=DABS(R1S(KK))*(1.0D0+DELT)
               VXP(KK)=VXP0S*R1S(KK)
               KK=KK+1
  50        CONTINUE
            CALL MC11BD(R1S,NN,N,IR)
            IF (IR.NE.N) THEN
               IF (LL.LE.300) THEN
                  DELT=(1.0D1**LL)*EPSM
                  LL=LL+1
                  GOTO 45
               ELSE
                  INFO=30
                  RETURN
               END IF
            END IF
         ELSE
            DO 55 II=1,NN
               VXO(II)=0.0D0
               VXP(II)=0.0D0
  55        CONTINUE
         CALL PROPA2(XP0,XP0,VXO,VXP,T,TSAMP,NONLIN,ODEEPS,EPSM,0,INFO)
            DO 60 II=1,NN
               VXP(II)=VXP(II)*VXP0S
  60        CONTINUE
         END IF
      ELSE if (ESTVAR == 2) THEN
         VXO = 0.0D0
         VXP = VX0
      ELSE
!        Set the covariance matrix to k*I
         KK=1
         DO JJ=1,N
            DO II=JJ,N
               VXO(KK) = 0.0D0
               IF (II.EQ.JJ) THEN
                  VXP(KK) = DEXP(XM(NPARAM))
               ELSE
                  VXP(KK) = 0.0D0
               END IF
               KK=KK+1
            END DO
         END DO
      END IF
!     ..................................................................
!                                 FACTORIZATION OF THE INITIAL VARIANCE.
!
      CALL MC11BD(VXP,NN,N,IR)
      IF (IR.NE.N) THEN
         INFO=30
         RETURN
      END IF
!     ..................................................................
!                                       LOOP FOR LIKELIHOOD CALCULATION.
!
      DO 80 I=1,NOBS(K)
         IF ((CTSFLG.EQ.1).AND.(I.NE.NOBS(K))) THEN
            TSAMP=(TMAT(NDATA+I+1)-TMAT(NDATA+I))/DBLE(NSEKF)
         END IF
         DO 81 II=1,M
            UO(II)=IMAT((NDATA+I-1)*(M+S)+II)
            IF ((HN.EQ.1).AND.(I.NE.NOBS(K))) THEN
               UP(II)=IMAT((NDATA+I)*(M+S)+II)
               ALPHA(II)=(UP(II)-UO(II))/(TSAMP*DBLE(NSEKF))
               UP(II)=UO(II)+ALPHA(II)*TSAMP
            END IF
  81     CONTINUE
         DO 82 II=1,S
            Y(II)=IMAT((NDATA+I-1)*(M+S)+M+II)
  82     CONTINUE
         IF (NONLIN.GE.1) THEN
            T=TMAT(NDATA+I)
         END IF
!     ..................................................................
!                   KALMAN FILTER CALL FOR REGULAR PARAMETER ESTIMATION.
!
         IF (JOB.EQ.0) THEN
            PRED=0
            CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,R1M,R2M,&
     &                  EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,SVDFLG,U,&
     &                  DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,MISS,EPSM,&
     &                  NIEKF,EKFEPS,ODEEPS,INFO)
            IF (INFO.NE.0) RETURN
!     ..................................................................
!             KALMAN FILTER CALL FOR STATE ESTIMATION (PURE SIMULATION).
!
         ELSEIF (JOB.EQ.-1) THEN
            CALL MC11CD(VXP,NN,N)
            CURPOS=(NDATA+I-1)*(2*N+SCOV+2*S+OCOV)
            DO J=1,N
               OMAT(CURPOS+J)=XP(J)
            END DO
            CURPOS=CURPOS+N
            IF (outopt.EQ.0) THEN
               II=1
               DO J=1,N
                  OMAT(CURPOS+J)=DSQRT(VXP(II))
                  II=II+N-J+1
               END DO
               CURPOS=CURPOS+N
            ELSE
               DO J=1,(N*(N+1)/2)
                  OMAT(CURPOS+J)=VXP(J)
               END DO
               CURPOS=CURPOS+N*(N+1)/2
            END IF
!           Update the current pos with J from the DO's
            
            CALL MC11BD(VXP,NN,N,IR)
            PRED=3
            CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,R1M,R2M,&
     &                  EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,SVDFLG,U,&
     &                  DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,MISS,EPSM,&
     &                  NIEKF,EKFEPS,ODEEPS,INFO)
            IF (INFO.NE.0) RETURN
            CALL MC11CD(VYPS,S*(S+1)/2,S)
            
            DO J=1,S
               OMAT(CURPOS+J)=YP(J)
            END DO
            CURPOS=CURPOS+S
            IF (OCOV.EQ.0) THEN
               II=1
               DO J=1,S
                  OMAT(CURPOS+J)=DSQRT(VYPS(II))
                  II=II+S-J+1
               END DO
            ELSE
               DO J=1,(S*(S+1)/2)
                  OMAT(CURPOS+J)=VYPS(J)
               END DO
            END IF
!     ..................................................................
!           KALMAN FILTER CALL FOR STATE ESTIMATION (K-STEP PREDICTION).
!
         ELSEIF (JOB.GE.1) THEN
            IF (I.EQ.1) THEN
               CALL MC11CD(VXP,NN,N)
!               II=1
!               DO 87 J=1,N
!                  OMAT((NDATA+I-1)*(2*N+2*S)+J)=XP(J)
!                  OMAT((NDATA+I-1)*(2*N+2*S)+N+J)=DSQRT(VXP(II))
!                  II=II+N-J+1
!  87           CONTINUE

               call store_covariance(omat,nomat,xp,vxp,N,outopt)

               CALL MC11BD(VXP,NN,N,IR)
               PRED=1
               CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,R1M,&
     &                     R2M,EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,&
     &                     SVDFLG,U,DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,&
     &                     MISS,EPSM,NIEKF,EKFEPS,ODEEPS,INFO)
               CALL MC11CD(VYPS,S*(S+1)/2,S)
!               II=1
!               DO 88 J=1,S
!                  OMAT((NDATA+I-1)*(2*N+2*S)+2*N+J)=YP(J)
!                  OMAT((NDATA+I-1)*(2*N+2*S)+2*N+S+J)=DSQRT(VYPS(II))
!                  II=II+S-J+1
!  88           CONTINUE

               call store_covariance(omat,nomat,yp,vyps,S,outopt)

            END IF
            PRED=0
            CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,R1M,R2M,&
     &                  EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,SVDFLG,U,&
     &                  DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,MISS,EPSM,&
     &                  NIEKF,EKFEPS,ODEEPS,INFO)
            IF (INFO.NE.0) RETURN
            DO 89 II=1,N
               XPS(II)=XP(II)
  89        CONTINUE
            DO 90 II=1,NN
               VXPS(II)=VXP(II)
  90        CONTINUE
            TSAMPS=TSAMP
            TS=T
            DO 91 II=1,M+1
               UOS(II)=UO(II)
               UPS(II)=UP(II)
               ALPHAS(II)=ALPHA(II)
  91        CONTINUE
            DO 92 J=1,MIN(JOB,NOBS(K)-I)
               IF ((I.EQ.1).AND.(J.NE.1)) THEN
                CALL MC11CD(VYPS,S*(S+1)/2,S)
!                II=1
!                DO 93 JJ=1,S
!                 OMAT((NDATA+I+J-2)*(2*N+2*S)+2*N+JJ)=YP(JJ)
!                 OMAT((NDATA+I+J-2)*(2*N+2*S)+2*N+S+JJ)=DSQRT(VYPS(II))
!                 II=II+S-JJ+1
!  93            CONTINUE
                 call store_covariance(omat,nomat,yp,vyps,S,outopt)
               END IF
               DO 94 JJ=1,(NSEKF-1)
                  IF (HN.EQ.1) THEN
                     DO 95 II=1,M
                        UO(II)=UP(II)
                        UP(II)=UO(II)+ALPHA(II)*TSAMP
  95                 CONTINUE
                  END IF
                  T=T+JJ*TSAMP
                  PRED=2
                  CALL KALMAN(XPS,VXPS,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,&
     &                        R1M,R2M,EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,&
     &                        R2,SVDFLG,U,DSIGMA,GAM,GAMH,IDEG,SVDEPS,&
     &                        CTSFLG,MISS,EPSM,NIEKF,EKFEPS,ODEEPS,INFO)
                  IF (INFO.NE.0) RETURN
  94           CONTINUE
               IF ((I.EQ.1).AND.(J.NE.JOB)) THEN
                  CALL MC11CD(VXPS,NN,N)
!                  II=1
!                  DO 96 JJ=1,N
!                     OMAT((NDATA+I+J-1)*(2*N+2*S)+JJ)=XPS(JJ)
!                     OMAT((NDATA+I+J-1)*(2*N+2*S)+N+JJ)=DSQRT(VXPS(II))
!                     II=II+N-JJ+1
!  96              CONTINUE
                  call store_covariance(omat,nomat,xps,vxps,N,outopt)
                  CALL MC11BD(VXPS,NN,N,IR)
               END IF
               IF ((CTSFLG.EQ.1).AND.((I+J).NE.NOBS(K))) THEN
                  TSAMP=(TMAT(NDATA+I+J+1)-TMAT(NDATA+I+J))/DBLE(NSEKF)
               END IF
               DO 97 II=1,M
                  UO(II)=IMAT((NDATA+I+J-1)*(M+S)+II)
                  IF ((HN.EQ.1).AND.((I+J).NE.NOBS(K))) THEN
                     UP(II)=IMAT((NDATA+I+J)*(M+S)+II)
                     ALPHA(II)=(UP(II)-UO(II))/(TSAMP*DBLE(NSEKF))
                     UP(II)=UO(II)+ALPHA(II)*TSAMP
                  END IF
  97           CONTINUE
               IF (NONLIN.GE.1) THEN
                  T=TMAT(NDATA+I+J)
               END IF
               IF (J.EQ.JOB) THEN
                  PRED=1
               ELSE
                  PRED=3
               END IF
               CALL KALMAN(XPS,VXPS,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,R1M,&
     &                     R2M,EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,&
     &                     SVDFLG,U,DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,&
     &                     MISS,EPSM,NIEKF,EKFEPS,ODEEPS,INFO)
               IF (INFO.NE.0) RETURN
  92        CONTINUE
            IF (I.LE.(NOBS(K)-JOB)) THEN
               CALL MC11CD(VXPS,NN,N)
!               II=1
!               DO 98 J=1,N
!                  OMAT((NDATA+I+JOB-1)*(2*N+2*S)+J)=XPS(J)
!                  OMAT((NDATA+I+JOB-1)*(2*N+2*S)+N+J)=DSQRT(VXPS(II))
!                  II=II+N-J+1
!  98           CONTINUE
               call store_covariance(omat,nomat,xps,vxps,N,outopt)
               CALL MC11BD(VXPS,NN,N,IR)
               CALL MC11CD(VYPS,S*(S+1)/2,S)
!               II=1
!               DO 99 J=1,S
!                 OMAT((NDATA+I+JOB-1)*(2*N+2*S)+2*N+J)=YP(J)
!                 OMAT((NDATA+I+JOB-1)*(2*N+2*S)+2*N+S+J)=DSQRT(VYPS(II))
!                 II=II+S-J+1
!  99           CONTINUE
               call store_covariance(omat,nomat,yp,vyps,S,outopt)
            END IF
            TSAMP=TSAMPS
            T=TS
            DO 100 II=1,M+1
               UO(II)=UOS(II)
               UP(II)=UPS(II)
               ALPHA(II)=ALPHAS(II)
 100        CONTINUE
!     ..................................................................
!                   KALMAN FILTER CALL FOR STATE ESTIMATION (FILTERING).
!
         ELSEIF (JOB.EQ.-2) THEN
            PRED=0
            CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,R1M,R2M,&
     &                  EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,SVDFLG,U,&
     &                  DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,MISS,EPSM,&
     &                  NIEKF,EKFEPS,ODEEPS,INFO)
            IF (INFO.NE.0) RETURN

            CALL MC11CD(VXO,NN,N)

            CURPOS=(NDATA+I-1)*(2*N+SCOV)

            DO J=1,N
               OMAT(CURPOS+J)=XO(J)
            END DO
            CURPOS=CURPOS+N
            IF (outopt.EQ.0) THEN
               II=1
               DO J=1,N
                  OMAT(CURPOS+J)=DSQRT(VXO(II))
                  II=II+N-J+1
               END DO
            ELSE
               DO J=1,(N*(N+1)/2)
                  OMAT(CURPOS+J)=VXO(J)
               END DO
            END IF
            CALL MC11BD(VXO,NN,N,IR)
!     ..................................................................
!                   KALMAN FILTER CALL FOR STATE ESTIMATION (SMOOTHING).
!
         ELSE
            IF (I.EQ.1) THEN
               DO 102 II=1,N
                  OMAT(NDATA*(N+NN)+II)=XP(II)
 102           CONTINUE
               DO 103 II=1,NN
                  OMAT(NDATA*(N+NN)+N+II)=VXP(II)
 103           CONTINUE
               PRED=0
               CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,1,R1M,&
     &                     R2M,EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,&
     &                     SVDFLG,U,DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,&
     &                     MISS,EPSM,NIEKF,EKFEPS,ODEEPS,INFO)
               IF (INFO.NE.0) RETURN
               DO 104 J=1,NOBS(K)-1
                  IF ((CTSFLG.EQ.1).AND.(J.NE.(NOBS(K)-1))) THEN
                     TSAMP=(TMAT(NDATA+J+2)-TMAT(NDATA+J+1))
                  END IF
                  DO 105 II=1,M
                     UO(II)=IMAT((NDATA+J)*(M+S)+II)
                     IF ((HN.EQ.1).AND.(J.NE.(NOBS(K)-1))) THEN
                        UP(II)=IMAT((NDATA+J+1)*(M+S)+II)
                        ALPHA(II)=(UP(II)-UO(II))/TSAMP
                     END IF
 105              CONTINUE
                  DO 106 II=1,S
                     Y(II)=IMAT((NDATA+J)*(M+S)+M+II)
 106              CONTINUE
                  T=TMAT(NDATA+J+1)
                  DO 107 II=1,N
                     OMAT((NDATA+J)*(N+NN)+II)=XP(II)
 107              CONTINUE
                  DO 108 II=1,NN
                     OMAT((NDATA+J)*(N+NN)+N+II)=VXP(II)
 108              CONTINUE
                  PRED=0
                  CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,1,R1M,&
     &                        R2M,EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,&
     &                        SVDFLG,U,DSIGMA,GAM,GAMH,IDEG,SVDEPS,&
     &                        CTSFLG,MISS,EPSM,NIEKF,EKFEPS,ODEEPS,INFO)
                  IF (INFO.NE.0) RETURN
 104           CONTINUE
               DO 109 II=1,N
                  XP(II)=0.0D0
 109           CONTINUE
               DO 110 II=1,NN
                  VXP(II)=0.0D0
 110           CONTINUE
               TF=TMAT(NDATA+NOBS(K))
               DO 111 J=NOBS(K),1,-1
                  IF ((CTSFLG.EQ.1).AND.(J.NE.1)) THEN
                     TSAMP=(TMAT(NDATA+J)-TMAT(NDATA+J-1))
                  END IF
                  DO 112 II=1,M
                     UO(II)=IMAT((NDATA+J-1)*(M+S)+II)
                     IF (J.NE.1) THEN
                        UP(II)=IMAT((NDATA+J-2)*(M+S)+II)
                     END IF
                     ALPHA(II)=(UP(II)-UO(II))/TSAMP
 112              CONTINUE
                  DO 113 II=1,S
                     Y(II)=IMAT((NDATA+J-1)*(M+S)+M+II)
 113              CONTINUE
                  T=TMAT(NDATA+J)
                  DO 114 II=1,N
                     XCUR(II)=OMAT((NDATA+J-1)*(N+NN)+II)
                     IF (J.NE.1) THEN
                        XPREV(II)=OMAT((NDATA+J-2)*(N+NN)+II)
                     END IF
 114              CONTINUE
                  CALL BWFLTR(XP,VXP,XO,VXO,Y,INFO)
                  IF (INFO.NE.0) RETURN
                  DO 115 II=1,N
                     OMAT((NDATA+NOBS(K)+J-1)*(N+NN)+II)=XO(II)
 115              CONTINUE
                  DO 116 II=1,NN
                     OMAT((NDATA+NOBS(K)+J-1)*(N+NN)+N+II)=VXO(II)
 116              CONTINUE
 111           CONTINUE
            END IF
            DO 117 II=1,N
               XO(II)=OMAT((NDATA+I-1)*(N+NN)+II)
               XP(II)=OMAT((NDATA+NOBS(K)+I-1)*(N+NN)+II)
 117        CONTINUE
            DO 118 II=1,NN
               VXO(II)=OMAT((NDATA+I-1)*(N+NN)+N+II)
               VXP(II)=OMAT((NDATA+NOBS(K)+I-1)*(N+NN)+N+II)
 118        CONTINUE
            CALL MC11ED(VXO,NN,N,XO,TMP5,N)
            CALL MATADD(XO,XP,XPS,N,1)
            CALL MC11FD(VXO,NN,N,N)
            CALL MATADD(VXO,VXP,VXPS,NN,1)
            CALL MC11BD(VXPS,NN,N,IR)
            IF (IR.NE.N) THEN
               INFO=30
               RETURN
            END IF
            CALL MC11ED(VXPS,NN,N,XPS,TMP5,N)
            CALL MC11FD(VXPS,NN,N,N)
            II=1
            DO 119 J=1,N
               OMAT((NDATA+I-1)*2*N+J)=XPS(J)
               OMAT((NDATA+I-1)*2*N+N+J)=DSQRT(VXPS(II))
               II=II+N-J+1
 119        CONTINUE
            IF (I.EQ.NOBS(K)) THEN
               DO 120 J=((NDATA+NOBS(K))*2*N+1),((NDATA+NOBS(K))*(N+NN))
                  OMAT(J)=0.0D0
 120           CONTINUE
            END IF
         END IF
!     ..................................................................
!                           NEGATIVE OF CONDITIONAL LIKELIHOOD FUNCTION.
!
         IF (JOB.EQ.0) THEN
           IF (MISS.NE.S) THEN
              DO 125 II=1,S
                 TMP3(II)=EPS(II)
                 DO 125 JJ=1,S
                    TMP2(JJ,II)=0.0D0
 125          CONTINUE
              CALL MC11ED(VYP,(S-MISS)*(S-MISS+1)/2,S-MISS,TMP3,TMP6,&
     &                    S-MISS)
              TMP4=0.0D0
              DO 130 II=1,S
                 TMP4=TMP4+EPS(II)*TMP3(II)
 130          CONTINUE
              CALL EXPAND(VYP,(S-MISS)*(S-MISS+1)/2,S-MISS,TMP2,TMP3)
              DET=1.0D0
              DO 135 II=1,(S-MISS)
                 DET=DET*TMP3(II)
 135          CONTINUE
!     ..................................................................
!                                                   OUTLIER ELIMINATION.
!
              IF (TMP4.GE.HC**2) TMP4=HC*(2*DSQRT(TMP4)-HC)
              TMP4=0.5D0*(DLOG(DET)+TMP4)
              LOGLIK=LOGLIK+TMP4
           END IF
           MV=MV+MISS
         END IF
!     ..................................................................
!                                    EXTENDED KALMAN FILTER SUBSAMPLING.
!
         IF ((NSEKF.GE.2).AND.(I.NE.NOBS(K))) THEN
            DO 140 JJ=1,(NSEKF-1)
               IF (HN.EQ.1) THEN
                  DO 145 II=1,M
                     UO(II)=UP(II)
                     UP(II)=UO(II)+ALPHA(II)*TSAMP
 145              CONTINUE
               END IF
               T=T+JJ*TSAMP
               PRED=2
               CALL KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,I,R1M,&
     &                     R2M,EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,&
     &                     SVDFLG,U,DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,&
     &                     MISS,EPSM,NIEKF,EKFEPS,ODEEPS,INFO)
               IF (INFO.NE.0) RETURN
 140        CONTINUE
         END IF
  80  CONTINUE
!     ..................................................................
!                                                CONSTANT TERM ADDITION.
!
      IF (JOB.EQ.0) THEN
         TMP4=6.283185308D0
         LOGLIK=LOGLIK+0.5D0*DBLE(NOBS(K)*S-MV)*DLOG(TMP4)
         PPDF=PPDF+LOGLIK
         NMISST=NMISST+MV
      END IF
      NDATA=NDATA+NOBS(K)
  32  CONTINUE
!     ..................................................................
!                                   CHECK FOR SUFFICIENT AMOUNT OF DATA.
!
      IF (JOB.EQ.0) THEN
         IF ((NDATA*S-NMISST).LE.NX) THEN
            INFO=10
            RETURN
         END IF
      END IF
!     ..................................................................
!                                       PENALTY AND PRIOR TERM ADDITION.
!
      FPEN=0.0D0
      FPRIOR=0.0D0
      IF (JOB.EQ.0) THEN
         DO 165 I=1,NX
            C1=DABS(XMIN(I))
            C2=DABS(XMAX(I))
            IF (C1.LT.XMM) C1=XMM
            IF (C2.LT.XMM) C2=XMM
            if (x(i) == xmin(i)) then
               fpen = fpen + c1 / EPSM
               call dblepr('ctsmr :: WARNING :: The lower boundary has been hit!', -1, 0.0d0, 1)
            else
               fpen = fpen + c1 / (x(i) - xmin(i))
            end if
            if (x(i) == xmax(i)) then
               fpen = fpen + c2 / EPSM
               call dblepr('ctsmr :: WARNING :: The upper boundary has been hit!', -1, 0.0d0, 1)
            else
               fpen = fpen + c2 / (xmax(i) - x(i))
            end if
!            FPEN=FPEN+C1/(X(I)-XMIN(I))+C2/(XMAX(I)-X(I))
 165     CONTINUE
         FPEN=LAMBDA*(FPEN-DBLE(NX))
         IF (NPR.GT.0) THEN
            DO 166 I=1,NPARAM
               TMP(I)=0.0D0
               TMP1(I)=0.0D0
 166        CONTINUE
            J=0
            DO 167 I=1,NPARAM
               IF (ID(I).EQ.2) THEN
                  J=J+1
                  TMP(J)=XM(I)-PRMEAN(J)
               END IF
 167        CONTINUE
            CALL DGEMV('T',NPR,NPR,1.0D0,PRCOVI,NPR,TMP,1,0.0D0,TMP1,1)
            DO 168 I=1,NPR
               FPRIOR=FPRIOR+TMP1(I)*TMP(I)
 168        CONTINUE
            FPRIOR=0.5D0*(FPRIOR+DBLE(NPR)*DLOG(TMP4)+DLOG(PRDET))
         END IF
         PPDF=PPDF+FPEN+FPRIOR
      END IF
      IF (DABS(PPDF).LT.1D300) THEN
         FUNK=PPDF
      ELSE
         INFO=20
         RETURN
      END IF
!     ..................................................................
!                                          TRANSFORMATION OF PARAMETERS.
!
      IF (JOB.EQ.0) THEN
         DO 340 I=1,NX
            X(I)=XT(I)
 340     CONTINUE
      END IF
      RETURN
      END
!
      SUBROUTINE DEFSYS(XM,XP0,XP,UO,T,A,B,C,D,R1M,R2M,F,H,NONLIN,&
     &                  SWITCH)
!----------------------------------------------------------------------C
!     SPECIFICATION OF THE ACTUAL DYNAMIC SYSTEM.                      C
!----------------------------------------------------------------------C
!     IF NONLIN=0/1 :                                                  C
!                                                                      C
!     DX = A(T)*X*DT + B(T)*U*DT + DW  ;    E[DW*(DW)'] = R1M          C
!     Y  = C(T)*X + D(T)*U + E         ;    E[E*(E)']   = R2M          C
!                                                                      C
!     X:   THE PARAMETER VECTOR.                                       C
!     T:   THE CURRENT TIME INSTANT.                                   C
!     A:   THE MATRIX A(N,N).                                          C
!     B:   THE MATRIX B(N,M).                                          C
!     C:   THE MATRIX C(S,N).                                          C
!     D:   THE MATRIX D(S,M).                                          C
!     R1M: THE MATRIX R1(N,N).                                         C
!     R2M: THE MATRIX R2(S,S).                                         C
!----------------------------------------------------------------------C
!     IF NONLIN=2 :                                                    C
!                                                                      C
!     DX = F(X,U,T) + DW         ;    E[DW*(DW)'] = R1M                C
!     Y  = H(X,U,T) + E          ;    E[E*(E)']   = R2M                C
!                                                                      C
!     X:   THE PARAMETER VECTOR.                                       C
!     T:   THE CURRENT TIME INSTANT.                                   C
!     F:   THE VECTOR F(N).                                            C
!     A:   THE MATRIX A(N,N).        =DF/DX                            C
!     B:   THE MATRIX B(N,M).        =DF/DU                            C
!     C:   THE MATRIX C(S,N).        =DH/DX                            C
!     D:   THE MATRIX D(S,M).        =DH/DU                            C
!     R1M: THE MATRIX R1(N,N).                                         C
!     R2M: THE MATRIX R2(S,S).                                         C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NONLIN,SWITCH
      DOUBLE PRECISION XM(NPARAM),XP0(N),XP(N),UO(M+1),T,A(N,N),&
     &                 B(N,M+1),C(S,N),D(S,M+1),R1M(N,N),R2M(S,S),F(N),&
     &                 H(S)
      INTEGER          I,J
      DOUBLE PRECISION JACOB1(N,N),JACOB2(N,M+1),JACOB3(S,N),&
     &                 JACOB4(S,M+1),SIGMA(N,N)
!
      IF (SWITCH.EQ.0) THEN
         DO 10 I=1,N
            XP0(I)=XM(I)
  10     CONTINUE
         RETURN
      END IF
      IF (SWITCH.EQ.1) THEN
         DO 20 I=1,N
            XP0(I)=XM(I)
  20     CONTINUE
         CALL AMAT(XM,XP,UO,T,A,NPARAM,N,M+1)
         CALL BMAT(XM,XP,UO,T,B,NPARAM,N,M+1)
         CALL CMAT(XM,XP,UO,T,C,NPARAM,N,M+1,S)
         CALL DMAT(XM,XP,UO,T,D,NPARAM,N,M+1,S)
         CALL SIGMAT(XM,UO,T,SIGMA,NPARAM,N,M+1)
         CALL DGEMM('N','T',N,N,N,1.0D0,SIGMA,N,SIGMA,N,0.0D0,R1M,N)
         CALL SMAT(XM,UO,T,R2M,NPARAM,N,M+1,S)
!         CALL SMAT(XM,UO,T,R2M,NPARAM,S,M+1)
         RETURN
      END IF
      IF (SWITCH.EQ.2) THEN
         IF (NONLIN.LE.1) THEN
            CALL CMAT(XM,XP,UO,T,C,NPARAM,N,M+1,S)
            CALL DMAT(XM,XP,UO,T,D,NPARAM,N,M+1,S)
         ELSE
            CALL HVECXJ(XM,XP,UO,T,JACOB3,H,NPARAM,N,M+1,S)
!            CALL HVECXJ(XM,XP,UO,T,H,NPARAM,N,M+1,S,JACOB3,S,YJACO3,
!     $                  LYJAC3,RFS3,IFS3,LFS3)
            DO 30 I=1,N
            DO 30 J=1,S
               C(J,I)=JACOB3(J,I)
  30        CONTINUE
            CALL HVECUJ(XM,XP,UO,T,JACOB4,NPARAM,N,M+1,S)
!            CALL HVECUJ(XM,XP,UO,T,H,NPARAM,N,M+1,S,JACOB4,S,YJACO4,
!     $                  LYJAC4,RFS4,IFS4,LFS4)
            DO 40 I=1,M
            DO 40 J=1,S
               D(J,I)=JACOB4(J,I)
  40        CONTINUE
         END IF
         CALL SMAT(XM,UO,T,R2M,NPARAM,N,M+1,S)
!         CALL SMAT(XM,UO,T,R2M,NPARAM,S,M+1)
         RETURN
      END IF
      IF (SWITCH.EQ.3) THEN
         IF (NONLIN.LE.1) THEN
            CALL AMAT(XM,XP,UO,T,A,NPARAM,N,M+1)
            CALL BMAT(XM,XP,UO,T,B,NPARAM,N,M+1)
         ELSE
            CALL FVECXJ(XM,XP,UO,T,JACOB1,F,NPARAM,N,M+1)
!            CALL FVECXJ(XM,XP,UO,T,F,NPARAM,N,M+1,JACOB1,N,YJACO1,
!     $                  LYJAC1,RFS1,IFS1,LFS1)
            DO 50 I=1,N
            DO 50 J=1,N
               A(J,I)=JACOB1(J,I)
  50        CONTINUE
            CALL FVECUJ(XM,XP,UO,T,JACOB2,NPARAM,N,M+1)
!            CALL FVECUJ(XM,XP,UO,T,F,NPARAM,N,M+1,JACOB2,N,YJACO2,
!     $                  LYJAC2,RFS2,IFS2,LFS2)
            DO 60 I=1,M
            DO 60 J=1,N
               B(J,I)=JACOB2(J,I)
  60        CONTINUE
         END IF
         CALL SIGMAT(XM,UO,T,SIGMA,NPARAM,N,M+1)
         CALL DGEMM('N','T',N,N,N,1.0D0,SIGMA,N,SIGMA,N,0.0D0,R1M,N)
         RETURN
      END IF
      RETURN
      END
!
      SUBROUTINE KALMAN(XP,VXP,XO,VXO,YP,VYP,VYPS,Y,A,B,C,D,OBSNO,R1M,&
     &           R2M,EPS,PRED,F,H,NONLIN,TSAMP,PHI,R1S,R2,SVDFLG,U,&
     &           DSIGMA,GAM,GAMH,IDEG,SVDEPS,CTSFLG,MISS,EPSM,NIEKF,&
     &           EKFEPS,ODEEPS,INFO)
!----------------------------------------------------------------------C
!     FORWARD KALMAN FILTER UPDATE OF STATES AND COVARIANCES.          C
!----------------------------------------------------------------------C
!     PRED=0 = (XP,VXP,XO,VXO,YP,VYP).                                 C
!     PRED=1 = (YP,VYPS).                                              C
!     PRED=2 = (XP,VXP).                                               C
!     PRED=3 = (XP,VXP,YP,VYPS).                                       C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          OBSNO,PRED,NONLIN,SVDFLG,IDEG,CTSFLG,MISS,NIEKF,&
     &                 INFO
      DOUBLE PRECISION XP(N),VXP(N*(N+1)/2),YP(S),VYP(S*(S+1)/2),&
     &                 VYPS(S*(S+1)/2),Y(S),A(N,N),B(N,M+1),C(S,N),&
     &                 D(S,M+1),R1M(N,N),R2M(S,S),EPS(S),XO(N),&
     &                 VXO(N*(N+1)/2),F(N),H(S),TSAMP,U(N,N),DSIGMA(N),&
     &                 PHI(N,N),R1S(N*(N+1)/2),R2(S*(S+1)/2),GAM(N,N),&
     &                 GAMH(N,N),SVDEPS,EPSM,EKFEPS,ODEEPS
      INTEGER          HN,I,IEKF,II,IR,J,JJ,K,KK,LL,NN,NS
      DOUBLE PRECISION ALPHA(M+1),CMISS(S*N),DELT,DW1(S+N),DW2(N+S),&
     &                 DW3(N+N),E(S,S),EPSTMP(S),ERR,ETANEW(N),&
     &                 ETAOLD(N),KF(N*S),LC1(N*N),LC2(S*S),LC3(N*N),&
     &                 NVYP(S*(S+1)/2),R1SM(N,N),R2MISS(S*(S+1)/2),T,&
     &                 TMP(S),TMP1(S),TMP2(S*S),TMP3(N*(N+1)/2),&
     &                 TMP4(N*N),TMP5(N*S),TMP6(N),TMP7(S,S),TMP8(S,N),&
     &                 TMP9(N),TMP10(S),UO(M+1),W1(S*(S+N)),W2(N*(N+S)),&
     &                 W3(N*(N+N)),XM(NPARAM),XP0(N),YPTMP(S)
      COMMON           /COM/XM,ALPHA,UO,T,HN
!$OMP THREADPRIVATE(/COM/)
      IEKF=1
      MISS=0
      NN=N*(N+1)/2
      NS=S*(S+1)/2
      IF ((NONLIN.EQ.0).AND.(OBSNO.NE.1)) THEN
         GOTO 2
      ELSE
         IF (NONLIN.GE.1) THEN
            CALL DEFSYS(XM,XP0,XP,UO,T,A,B,C,D,R1M,R2M,F,H,NONLIN,2)
         END IF
         KK=1
         DO 1 JJ=1,S
         DO 1 II=JJ,S
            R2(KK)=R2M(II,JJ)
            KK=KK+1
   1     CONTINUE
         CALL MC11BD(R2,NS,S,IR)
         IF (IR.NE.S) THEN
            INFO=40
            RETURN
         END IF
      END IF
   2  IF (PRED.GE.2) THEN
         MISS=S
         IF (PRED.EQ.2) GOTO 100
      END IF
!     ..................................................................
!                      CALCULATE THE VARIANCE OF THE FULL OUTPUT VECTOR.
!
      IF (PRED.NE.0) THEN
         CALL UPDATE(VYPS,C,VXP,R2,S,NS,N,NN,W1,DW1,LC1,EPSM)
      END IF
!     ..................................................................
!                          CALCULATE THE MEAN OF THE FULL OUTPUT VECTOR.
!
   3  IF (NONLIN.LE.1) THEN
         CALL DGEMV('N',S,N,1.0D0,C,S,XP,1,0.0D0,TMP,1)
         CALL DGEMV('N',S,M,1.0D0,D,S,UO,1,0.0D0,TMP1,1)
         CALL MATADD(TMP,TMP1,YP,S,1)
      ELSE
         DO 4 I=1,S
            YP(I)=H(I)
   4     CONTINUE
      END IF
      IF (PRED.EQ.3) GOTO 100
      IF (PRED.EQ.1) RETURN
!     ..................................................................
!                    CHECK FOR MISSING DATA IN THE ACTUAL OUTPUT VECTOR.
!
      IF (IEKF.EQ.1) THEN
         MISS=0
         DO 5 I=1,S
            IF (Y(I).GE.1D300) THEN
               IF (MISS.EQ.0) THEN
                  DO 6 K=1,S
                  DO 6 J=1,S
                     IF (K.EQ.J) THEN
                        E(J,K)=1.0D0
                     ELSE
                        E(J,K)=0.0D0
                     END IF
   6              CONTINUE
               END IF
               MISS=MISS+1
               Y(I)=YP(I)
               DO 7 K=1,S
                  DO 8 J=I,(S-1)
                     E(J,K)=E(J+1,K)
   8              CONTINUE
                  E(S,K)=0.0D0
   7           CONTINUE
            END IF
   5     CONTINUE
      END IF
!     ..................................................................
!                          CALCULATE THE ACTUAL PREDICTION ERROR VECTOR.
!
      CALL MATSUB(Y,YP,EPS,S,1)
!     ..................................................................
!                       HANDLE MISSING DATA IN THE ACTUAL OUTPUT VECTOR.
!
      IF (MISS.EQ.S) THEN
         GOTO 100
      ELSEIF (MISS.GE.1) THEN
         CALL DGEMV('N',S,S,1.0D0,E,S,EPS,1,0.0D0,TMP,1)
         CALL DCOPY(S,TMP,1,EPS,1)
         IF (IEKF.EQ.1) THEN
            CALL DGEMM('N','N',S,S,S,1.0D0,E,S,R2M,S,0.0D0,TMP2,S)
            CALL DGEMM('N','T',S,S,S,1.0D0,TMP2,S,E,S,0.0D0,TMP7,S)
            DO 9 I=1,NS
               R2MISS(I)=0.0D0
   9        CONTINUE
            K=1
            DO 10 J=1,(S-MISS)
            DO 10 I=J,(S-MISS)
               R2MISS(K)=TMP7(I,J)
               K=K+1
  10        CONTINUE
            CALL MC11BD(R2MISS,(S-MISS)*(S-MISS+1)/2,S-MISS,IR)
            IF (IR.NE.(S-MISS)) THEN
               INFO=40
               RETURN
            END IF
         END IF
         CALL DGEMM('N','N',S,N,S,1.0D0,E,S,C,S,0.0D0,TMP8,S)
         DO 11 I=1,S*N
            CMISS(I)=0.0D0
  11     CONTINUE
         K=1
         DO 12 I=1,N
         DO 12 J=1,(S-MISS)
            CMISS(K)=TMP8(J,I)
            K=K+1
  12     CONTINUE
      END IF
!     ..................................................................
!                    CALCULATE THE VARIANCE OF THE ACTUAL OUTPUT VECTOR.
!
      IF (MISS.EQ.0) THEN
         CALL UPDATE(VYP,C,VXP,R2,S,NS,N,NN,W1,DW1,LC1,EPSM)
      ELSE
         DO 14 I=1,NS
            VYP(I)=0.0D0
  14     CONTINUE
         CALL UPDATE(VYP,CMISS,VXP,R2MISS,S-MISS,(S-MISS)*(S-MISS+1)/2,&
     &               N,NN,W1,DW1,LC1,EPSM)
      END IF
!     ..................................................................
!                                             CALCULATE THE KALMAN GAIN.
!
      DO 15 I=1,NS
         NVYP(I)=VYP(I)
  15  CONTINUE
      IF (MISS.EQ.0) THEN
         CALL MC11FD(NVYP,NS,S,S)
         K=1
         DO 20 I=1,S
         DO 20 J=I,S
            TMP2((J-1)*S+I)=NVYP(K)
            TMP2((I-1)*S+J)=TMP2((J-1)*S+I)
            K=K+1
  20     CONTINUE
         DO 30 I=1,NN
            TMP3(I)=VXP(I)
  30     CONTINUE
         CALL MC11CD(TMP3,NN,N)
         K=1
         DO 40 I=1,N
         DO 40 J=I,N
            TMP4((J-1)*N+I)=TMP3(K)
            TMP4((I-1)*N+J)=TMP4((J-1)*N+I)
            K=K+1
  40     CONTINUE
         CALL DGEMM('N','T',N,S,N,1.0D0,TMP4,N,C,S,0.0D0,TMP5,N)
         CALL DGEMM('N','N',N,S,S,1.0D0,TMP5,N,TMP2,S,0.0D0,KF,N)
      ELSE
         CALL MC11FD(NVYP,(S-MISS)*(S-MISS+1)/2,S-MISS,S-MISS)
         DO 50 I=1,S*S
            TMP2(I)=0.0D0
  50     CONTINUE
         K=1
         DO 60 I=1,(S-MISS)
         DO 60 J=I,(S-MISS)
            TMP2((J-1)*(S-MISS)+I)=NVYP(K)
            TMP2((I-1)*(S-MISS)+J)=TMP2((J-1)*(S-MISS)+I)
            K=K+1
  60     CONTINUE
         DO 70 I=1,NN
            TMP3(I)=VXP(I)
  70     CONTINUE
         CALL MC11CD(TMP3,NN,N)
         K=1
         DO 80 I=1,N
         DO 80 J=I,N
            TMP4((J-1)*N+I)=TMP3(K)
            TMP4((I-1)*N+J)=TMP4((J-1)*N+I)
            K=K+1
  80     CONTINUE
         CALL DGEMM('N','T',N,S,N,1.0D0,TMP4,N,TMP8,S,0.0D0,TMP5,N)
         DO 90 I=1,(N*S)
            KF(I)=0.0D0
  90     CONTINUE
         CALL DGEMM('N','N',N,S-MISS,S-MISS,1.0D0,TMP5,N,TMP2,S-MISS,&
     &              0.0D0,KF,N)
      END IF
!     ..................................................................
!                                       CALCULATE THE A POSTERIORI MEAN.
!
 100  IF (MISS.EQ.S) THEN
         DO 110 I=1,N
            XO(I)=XP(I)
 110     CONTINUE
      ELSE
         CALL DGEMV('N',N,S,1.0D0,KF,N,EPS,1,0.0D0,TMP6,1)
         CALL MATADD(TMP6,XP,XO,N,1)
         IF (NONLIN.GE.2) THEN
            IF (IEKF.LT.NIEKF) THEN
               IF (IEKF.EQ.1) THEN
                  DO 112 I=1,N
                     ETAOLD(I)=XP(I)
 112              CONTINUE
                  DO 113 I=1,S
                     EPSTMP(I)=EPS(I)
                     YPTMP(I)=YP(I)
 113              CONTINUE
               END IF
               CALL MATSUB(XP,ETAOLD,TMP9,N,1)
               CALL DGEMV('N',S,N,1.0D0,C,S,TMP9,1,0.0D0,TMP10,1)
               CALL DGEMV('N',N,S,1.0D0,KF,N,TMP10,1,0.0D0,TMP9,1)
               CALL MATSUB(XO,TMP9,ETANEW,N,1)
               ERR=0.0D0
               DO 114 I=1,N
                  ERR=ERR+DABS(1.0D0-ETAOLD(I)/ETANEW(I))
 114           CONTINUE
               IF (ERR.GE.EKFEPS) THEN
                  CALL DEFSYS(XM,XP0,ETANEW,UO,T,A,B,C,D,R1M,R2M,F,H,&
     &                        NONLIN,2)
                  DO 116 I=1,N
                     ETAOLD(I)=ETANEW(I)
 116              CONTINUE
                  IEKF=IEKF+1
                  GOTO 3
               ELSE
                  DO 117 I=1,S
                     EPS(I)=EPSTMP(I)
                     YP(I)=YPTMP(I)
 117              CONTINUE
                  DO 118 I=1,N
                    XO(I)=ETANEW(I)
 118              CONTINUE
               END IF
            END IF
         END IF
      END IF
!     ..................................................................
!                                   CALCULATE THE A POSTERIORI VARIANCE.
!
      IF (MISS.EQ.S) THEN
         DO 120 I=1,NN
            VXO(I)=VXP(I)
 120     CONTINUE
      ELSE
         DO 130 I=1,NS
            NVYP(I)=VYP(I)
 130     CONTINUE
         K=1
         DO 140 I=1,(S-MISS)
            NVYP(K)=-VYP(K)
            K=K+S-MISS-I+1
 140     CONTINUE
         CALL UPDATE(VXO,KF,NVYP,VXP,N,NN,S-MISS,&
     &               (S-MISS)*(S-MISS+1)/2,W2,DW2,LC2,EPSM)
      END IF
!     ..................................................................
!                                       CALCULATE THE A PRIORI VARIANCE.
!
      IF ((NONLIN.EQ.0).AND.(CTSFLG.EQ.0).AND.(OBSNO.NE.1)) THEN
         GOTO 150
      ELSE
         IF (NONLIN.GE.1) THEN
            CALL DEFSYS(XM,XP0,XO,UO,T,A,B,C,D,R1M,R2M,F,H,NONLIN,3)
         END IF
         IF (NONLIN.LE.2) THEN
            CALL EXPA(A,R1M,TSAMP,PHI,R1SM,IDEG,INFO)
            IF (INFO.NE.0) RETURN
            DELT=0.0D0
            LL=1
 141        KK=1
            DO 142 JJ=1,N
            DO 142 II=JJ,N
               R1S(KK)=R1SM(II,JJ)
               IF (II.EQ.JJ) R1S(KK)=DABS(R1S(KK))*(1.0D0+DELT)
               KK=KK+1
 142        CONTINUE
            CALL MC11BD(R1S,NN,N,IR)
            IF (IR.NE.N) THEN
               IF (LL.LE.300) THEN
                  DELT=(1.0D1**LL)*EPSM
                  LL=LL+1
                  GOTO 141
               ELSE
                  INFO=30
                  RETURN
               END IF
            END IF
         END IF
      END IF
 150  IF (NONLIN.LE.2) THEN
         CALL UPDATE(VXP,PHI,VXO,R1S,N,NN,N,NN,W3,DW3,LC3,EPSM)
      END IF
!     ..................................................................
!                                           CALCULATE THE A PRIORI MEAN.
!
      IF (NONLIN.LE.2) THEN
         CALL PROPA1(A,B,F,PHI,GAM,GAMH,U,DSIGMA,UO,ALPHA,TSAMP,XO,&
     &               XP,NONLIN,HN,SVDFLG,SVDEPS,CTSFLG,INFO)
      ELSE
         CALL PROPA2(XO,XP,VXO,VXP,T,TSAMP,NONLIN,ODEEPS,EPSM,1,INFO)
      END IF
      RETURN
      END
!
      SUBROUTINE EXPA(A,R1M,T,PHI,R1SM,IDEG,INFO)
!----------------------------------------------------------------------C
!     CALCULATION OF MATRIX EXPONENTIAL.                               C
!----------------------------------------------------------------------C
      INCLUDE 'global.h'
      INTEGER IDEG,INFO
      DOUBLE PRECISION A(N,N),R1M(N,N),T,PHI(N,N),R1SM(N,N)
      INTEGER IPIV(2*N),IEXPH,NS,I,J
      DOUBLE PRECISION H(2*N,2*N),WSP(4*2*N*2*N+IDEG+1),F3(N,N),G2(N,N)
!
         DO 10 I=1,N
         DO 10 J=1,N
            H(J,I)=-A(J,I)
  10     CONTINUE
         DO 11 I=1,N
         DO 11 J=1,N
            H(J+N,I)=0.0D0
  11     CONTINUE
         DO 12 I=1,N
         DO 12 J=1,N
            H(J+N,I+N)=A(I,J)
  12     CONTINUE
         DO 13 I=1,N
         DO 13 J=1,N
            H(J,I+N)=R1M(J,I)
  13     CONTINUE
         CALL DGPADM(IDEG,2*N,T,H,2*N,WSP,4*2*N*2*N+IDEG+1,IPIV,IEXPH,&
     &            NS,INFO)
         IF (INFO.NE.0) THEN
            INFO=50
            RETURN
         END IF
         DO 20 I=1,N
         DO 30 J=1,N
            G2(J,I)=WSP(IEXPH+2*N*N)
            IEXPH=IEXPH+1
  30     CONTINUE
         DO 40 J=1,N
            F3(J,I)=WSP(IEXPH+2*N*N)
            PHI(I,J)=F3(J,I)
            IEXPH=IEXPH+1
  40     CONTINUE
  20     CONTINUE
         CALL DGEMM('T','N',N,N,N,1.0D0,F3,N,G2,N,0.0D0,R1SM,N)
         DO 45 I=1,N
         DO 45 J=I,N
            R1SM(J,I)=5D-1*(R1SM(J,I)+R1SM(I,J))
            R1SM(I,J)=R1SM(J,I)
  45     CONTINUE
      RETURN
      END
!
      SUBROUTINE PROPA1(A,B,F,PHI,GAM,GAMH,U,DSIGMA,UO,ALPHA,TSAMP,XO,&
     &                  XP,NONLIN,HN,SVDFLG,SVDEPS,CTSFLG,INFO)
!----------------------------------------------------------------------C
!     SVD-BASED SOLUTION OF FORWARD PROPAGATION EQUATION FOR MEAN.     C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NONLIN,HN,SVDFLG,CTSFLG,INFO
      DOUBLE PRECISION A(N,N),B(N,M+1),F(N),PHI(N,N),GAM(N,N),GAMH(N,N),&
     &                 U(N,N),DSIGMA(N),UO(M+1),ALPHA(M+1),TSAMP,XO(N),&
     &                 XP(N),SVDEPS
      CHARACTER*1      EQUED
      INTEGER          GAMFLG,I,J,K,IPIV(N),IWORK(N)
      DOUBLE PRECISION ZO(N),ZP(N),ATILDE(N,N),BTILDE(N,M+1),FTILDE(N),&
     &                 PHITIL(N,N),A1(N,N),A2(N,N),B1(N,M+1),B2(N,M+1),&
     &                 F1(N),F2(N),PHI1(N,N),AINV(N,N),PHII(N,N),&
     &                 AIPHII(N,N),TMP(N,N),TMP1(N),TMP2(N),TMP3(N),&
     &                 TMP4(N),TMP5(N,N),TMP6(N,N),TMP7(N,N),TMP8(N,N),&
     &                 TMP9(N),VT(N,N),UNIT(N,N),AF(N,N),R(N),C(N),&
     &                 ANORM,RCOND,FERR(N),BERR(N),WORK(5*N),WORK1(4*N)
      DOUBLE PRECISION DLANGE
      EXTERNAL         DLANGE
!
      GAMFLG=SVDFLG
      IF ((NONLIN.EQ.0).AND.(SVDFLG.NE.-1)) THEN
!     ..................................................................
!                                          SVD ALREADY CALCULATED (LTI).
!
         I=SVDFLG
      ELSE
         IF (N.EQ.1) THEN
!     ..................................................................
!                                       CHECKING SINGULARITY (SCALAR A).
!
            IF (A(1,1).EQ.0.0D0) THEN
               I=0
               GOTO 10
            END IF
         ELSE
!     ..................................................................
!                                 COMPUTING SVD IF NEEDED (NONSCALAR A).
!
            ANORM=DLANGE('1',N,N,A,N,WORK1)
            CALL DGECON('1',N,A,N,ANORM,RCOND,WORK1,IWORK,INFO)
            IF (INFO.NE.0) THEN
               INFO=60
               RETURN
            END IF
            IF (RCOND.LE.SVDEPS) THEN
               CALL DCOPY(N*N,A,1,TMP,1)
               CALL DGESVD('A','A',N,N,TMP,N,DSIGMA,U,N,VT,N,WORK,5*N,&
     &                     INFO)
               IF (INFO.NE.0) THEN
                  INFO=70
                  RETURN
               END IF
               IF (DSIGMA(1).EQ.0.0D0) THEN
                  I=0
                  GOTO 10
               ELSE
                  DO 5 I=1,N-1
                     IF ((DSIGMA(I+1)/DSIGMA(1)).LE.SVDEPS) THEN
                        GOTO 10
                     END IF
   5              CONTINUE
               END IF
            END IF
         END IF
         I=N
      END IF
  10  SVDFLG=I
      IF ((I.GT.0).AND.(I.LT.N)) THEN
!     ..................................................................
!                                   INITIALIZING VARIABLES (SINGULAR A).
!
         DO 20 J=1,N
            F1(J)=0.0D0
            F2(J)=0.0D0
         DO 15 K=1,(M+1)
            B1(J,K)=0.0D0
            B2(J,K)=0.0D0
            BTILDE(J,K)=0.0D0
  15     CONTINUE
         DO 20 K=1,N
            A1(K,J)=0.0D0
            AINV(K,J)=0.0D0
            A2(K,J)=0.0D0
            PHI1(K,J)=0.0D0
            UNIT(K,J)=0.0D0
  20     CONTINUE
!     ..................................................................
!                                                      COMPUTING ATILDE.
!
         CALL DGEMM('N','N',N,N,N,1.0D0,A,N,U,N,0.0D0,TMP,N)
         CALL DGEMM('T','N',N,N,N,1.0D0,U,N,TMP,N,0.0D0,ATILDE,N)
!     ..................................................................
!                                                      COMPUTING BTILDE.
!
         CALL DGEMM('T','N',N,M+1,N,1.0D0,U,N,B,N,0.0D0,BTILDE,N)
!     ..................................................................
!                                                      COMPUTING FTILDE.
!
         CALL DGEMV('T',N,N,1.0D0,U,N,F,1,0.0D0,FTILDE,1)
!     ..................................................................
!                                                    COMPUTING PHITILDE.
!
         CALL DGEMM('N','N',N,N,N,1.0D0,PHI,N,U,N,0.0D0,TMP,N)
         CALL DGEMM('T','N',N,N,N,1.0D0,U,N,TMP,N,0.0D0,PHITIL,N)
!     ..................................................................
!                   COMPUTING I, A1, A1INV, A2, B1, B2, PHI1, F1 AND F2.
!
         DO 30 J=1,N
            IF (J.LE.I) THEN
               UNIT(J,J)=1.0D0
               F1(J)=FTILDE(J)
               DO 40 K=1,(M+1)
                  B1(J,K)=BTILDE(J,K)
  40           CONTINUE
               DO 50 K=1,I
                  A1(J,K)=ATILDE(J,K)
                  PHI1(J,K)=PHITIL(J,K)
  50           CONTINUE
            ELSE
               F2(J-I)=FTILDE(J)
               DO 60 K=1,(M+1)
                  B2(J-I,K)=BTILDE(J,K)
  60           CONTINUE
               DO 70 K=1,I
                  A2(K,J-I)=ATILDE(K,J)
  70           CONTINUE
            END IF
  30     CONTINUE
         CALL DCOPY(N*N,PHI1,1,PHII,1)
         DO 80 J=1,I
            PHII(J,J)=PHII(J,J)-1.0D0
  80     CONTINUE
         CALL DCOPY(N*N,A1,1,TMP,1)
         CALL DGESVX('E','N',I,I,TMP,N,AF,N,IPIV,EQUED,R,C,UNIT,N,&
     &                AINV,N,RCOND,FERR,BERR,WORK1,IWORK,INFO)
         IF (INFO.NE.0) THEN
            INFO=80
            RETURN
         END IF
         CALL DGEMM('N','N',N,N,N,1.0D0,AINV,N,PHII,N,0.0D0,AIPHII,N)
      ELSEIF (I.EQ.N) THEN
!     ..................................................................
!                                INITIALIZING VARIABLES (NONSINGULAR A).
!
         DO 85 J=1,N
         DO 85 K=1,N
            UNIT(K,J)=0.0D0
  85     CONTINUE
         CALL DCOPY(N*N,PHI,1,PHII,1)
         DO 90 J=1,N
            UNIT(J,J)=1.0D0
            PHII(J,J)=PHII(J,J)-1.0D0
  90     CONTINUE
         CALL DCOPY(N*N,A,1,TMP,1)
         CALL DGESVX('E','N',N,N,TMP,N,AF,N,IPIV,EQUED,R,C,UNIT,N,&
     &                AINV,N,RCOND,FERR,BERR,WORK1,IWORK,INFO)
         IF (INFO.NE.0) THEN
            INFO=80
            RETURN
         END IF
         CALL DGEMM('N','N',N,N,N,1.0D0,AINV,N,PHII,N,0.0D0,AIPHII,N)
      END IF
!     ..................................................................
!                                                COMPUTING THE SOLUTION.
!
      IF (NONLIN.LE.1) THEN
!     ..................................................................
!                                                         LTI/LTV CASES.
!
         IF (((NONLIN.EQ.0).AND.(CTSFLG.EQ.0)).AND.(GAMFLG.NE.-1)) THEN
!     ..................................................................
!                                       GAM AND GAMH ALREADY CALCULATED.
!
            GOTO 300
         END IF
         IF (I.EQ.0) THEN
!     ..................................................................
!                                                                ZERO A.
!
!                                                          COMPUTING GAM
            DO 92 J=1,N
            DO 92 K=1,N
               GAM(K,J)=0.0D0
  92        CONTINUE
            DO 94 J=1,N
               GAM(J,J)=TSAMP
  94        CONTINUE
            IF (HN.EQ.1) THEN
!
!                                                         COMPUTING GAMH
!
               DO 96 J=1,N
               DO 96 K=1,N
                  GAMH(K,J)=0.0D0
  96           CONTINUE
               DO 98 J=1,N
                  GAMH(J,J)=5D-1*(TSAMP**2)
  98           CONTINUE
            END IF
         ELSEIF ((I.GT.0).AND.(I.LT.N)) THEN
!     ..................................................................
!                                                            SINGULAR A.
!
!                                                          COMPUTING GAM
            CALL DCOPY(N*N,AIPHII,1,TMP6,1)
            DO 100 J=1,I
               TMP6(J,J)=TMP6(J,J)-TSAMP
 100        CONTINUE
            CALL DGEMM('N','N',N,N,N,1.0D0,AINV,N,TMP6,N,0.0D0,TMP7,N)
            CALL DGEMM('N','N',N,N,N,1.0D0,TMP7,N,A2,N,0.0D0,TMP6,N)
            DO 110 J=1,I
               DO 120 K=1,I
                  GAM(K,J)=AIPHII(K,J)
 120           CONTINUE
               DO 130 K=1,(N-I)
                  GAM(J,I+K)=TMP6(J,K)
 130           CONTINUE
 110        CONTINUE
            DO 140 J=1,(N-I)
               DO 150 K=1,N
                  GAM(I+J,K)=0.0D0
 150           CONTINUE
               GAM(I+J,I+J)=TSAMP
 140        CONTINUE
            CALL DGEMM('N','N',N,N,N,1.0D0,U,N,GAM,N,0.0D0,TMP5,N)
            CALL DGEMM('N','T',N,N,N,1.0D0,TMP5,N,U,N,0.0D0,GAM,N)
            IF (HN.EQ.1) THEN
!
!                                                         COMPUTING GAMH
!
               CALL DCOPY(N*N,PHI1,1,TMP6,1)
               DO 160 J=1,I
               DO 160 K=1,I
                  TMP6(K,J)=TMP6(K,J)*TSAMP
 160           CONTINUE
               CALL MATSUB(TMP6,AIPHII,TMP7,N,N)
               CALL DGEMM('N','N',N,N,N,1D0,AINV,N,TMP7,N,0.0D0,TMP8,N)
               CALL DCOPY(N*N,TMP8,1,TMP6,1)
               DO 170 J=1,I
                  TMP6(J,J)=TMP6(J,J)-5D-1*(TSAMP**2)
 170           CONTINUE
               CALL DGEMM('N','N',N,N,N,1D0,AINV,N,TMP6,N,0.0D0,TMP7,N)
               CALL DGEMM('N','N',N,N,N,1.0D0,TMP7,N,A2,N,0.0D0,TMP6,N)
               DO 180 J=1,I
               DO 190 K=1,I
                  GAMH(K,J)=TMP8(K,J)
 190           CONTINUE
               DO 200 K=1,(N-I)
                  GAMH(J,I+K)=TMP6(J,K)
 200           CONTINUE
 180           CONTINUE
               DO 210 J=1,(N-I)
               DO 220 K=1,N
                  GAMH(I+J,K)=0.0D0
 220           CONTINUE
               GAMH(I+J,I+J)=5D-1*(TSAMP**2)
 210           CONTINUE
               CALL DGEMM('N','N',N,N,N,1.0D0,U,N,GAMH,N,0.0D0,TMP5,N)
               CALL DGEMM('N','T',N,N,N,1.0D0,TMP5,N,U,N,0.0D0,GAMH,N)
            END IF
         ELSE
!     ..................................................................
!                                                         NONSINGULAR A.
!
!                                                          COMPUTING GAM
            DO 230 J=1,I
            DO 230 K=1,I
               GAM(K,J)=AIPHII(K,J)
 230        CONTINUE
            IF (HN.EQ.1) THEN
!
!                                                         COMPUTING GAMH
!
               CALL DCOPY(N*N,PHI,1,TMP6,1)
               DO 240 J=1,N
               DO 240 K=1,N
                  TMP6(K,J)=TMP6(K,J)*TSAMP
 240           CONTINUE
               CALL MATSUB(TMP6,AIPHII,TMP7,N,N)
               CALL DGEMM('N','N',N,N,N,1D0,AINV,N,TMP7,N,0.0D0,TMP8,N)
               DO 250 J=1,N
               DO 250 K=1,N
                  GAMH(K,J)=TMP8(K,J)
 250           CONTINUE
            END IF
         END IF
!     ..................................................................
!                                                          COMPUTING XP.
!
 300     CALL DGEMV('N',N,N,1.0D0,PHI,N,XO,1,0.0D0,TMP1,1)
         CALL DGEMV('N',N,M+1,1.0D0,B,N,UO,1,0.0D0,TMP2,1)
         IF (HN.EQ.1) THEN
            CALL DGEMV('N',N,M+1,1.0D0,B,N,ALPHA,1,0.0D0,TMP3,1)
            DO 310 J=1,N
               TMP4(J)=TMP3(J)*TSAMP+TMP2(J)
 310        CONTINUE
            CALL DGEMV('N',N,N,1.0D0,GAM,N,TMP4,1,0.0D0,TMP2,1)
            CALL DGEMV('N',N,N,1.0D0,GAMH,N,TMP3,1,0.0D0,TMP4,1)
            CALL MATADD(TMP1,TMP2,TMP3,N,1)
            CALL MATSUB(TMP3,TMP4,XP,N,1)
         ELSE
            CALL DGEMV('N',N,N,1.0D0,GAM,N,TMP2,1,0.0D0,TMP3,1)
            CALL MATADD(TMP1,TMP3,XP,N,1)
         END IF
      ELSE
!     ..................................................................
!                                                        NONLINEAR CASE.
!
         IF (I.EQ.0) THEN
!     ..................................................................
!                                                                ZERO A.
!
            IF (HN.EQ.1) THEN
               CALL DGEMV('N',N,M+1,1.0D0,B,N,ALPHA,1,0.0D0,TMP1,1)
               DO 400 J=1,N
                  XP(J)=XO(J)+F(J)*TSAMP+5D-1*TMP1(J)*(TSAMP**2)
 400           CONTINUE
            ELSE
               DO 410 J=1,N
                  XP(J)=XO(J)+F(J)*TSAMP
 410           CONTINUE
            END IF
         ELSEIF ((I.GT.0).AND.(I.LT.N)) THEN
!     ..................................................................
!                                                            SINGULAR A.
!
            CALL DGEMV('T',N,N,1.0D0,U,N,XO,1,0.0D0,ZO,1)
            IF (HN.EQ.1) THEN
               CALL DGEMV('N',N,M+1,1.0D0,B1,N,ALPHA,1,0.0D0,TMP1,1)
               CALL DGEMV('N',N,N,1.0D0,A2,N,F2,1,0.0D0,TMP2,1)
               CALL MATADD(TMP1,TMP2,TMP3,N,1)
               CALL DGEMV('N',N,M+1,1.0D0,B2,N,ALPHA,1,0.0D0,TMP1,1)
               CALL DGEMV('N',N,N,1.0D0,A2,N,TMP1,1,0.0D0,TMP2,1)
               CALL DGEMV('N',N,N,1.0D0,AINV,N,TMP2,1,0.0D0,TMP4,1)
               CALL MATADD(TMP4,TMP3,TMP2,N,1)
               CALL DGEMV('N',N,N,1.0D0,AINV,N,TMP2,1,0.0D0,TMP3,1)
               CALL MATADD(TMP3,F1,TMP2,N,1)
               CALL DGEMV('N',N,N,1.0D0,AIPHII,N,TMP2,1,0.0D0,TMP9,1)
               DO 420 J=1,I
                  ZP(J)=ZO(J)-5D-1*TMP4(J)*(TSAMP**2)-TMP3(J)*TSAMP&
     &                  +TMP9(J)
 420           CONTINUE
               DO 430 J=1,(N-I)
                  ZP(J+I)=ZO(J+I)+F2(J)*TSAMP+5D-1*TMP1(J)&
     &                    *(TSAMP**2)
 430           CONTINUE
            ELSE
               CALL DGEMV('N',N,N,1.0D0,A2,N,F2,1,0.0D0,TMP1,1)
               CALL DGEMV('N',N,N,1.0D0,AINV,N,TMP1,1,0.0D0,TMP2,1)
               CALL MATADD(TMP2,F1,TMP1,N,1)
               CALL DGEMV('N',N,N,1.0D0,AIPHII,N,TMP1,1,0.0D0,TMP3,1)
               DO 440 J=1,I
                  ZP(J)=ZO(J)-TMP2(J)*TSAMP+TMP3(J)
 440           CONTINUE
               DO 450 J=1,(N-I)
                  ZP(J+I)=ZO(J+I)+F2(J)*TSAMP
 450           CONTINUE
            END IF
            CALL DGEMV('N',N,N,1.0D0,U,N,ZP,1,0.0D0,XP,1)
         ELSE
!     ..................................................................
!                                                         NONSINGULAR A.
!
            IF (HN.EQ.1) THEN
               CALL DGEMV('N',N,M+1,1.0D0,B,N,ALPHA,1,0.0D0,TMP1,1)
               CALL DGEMV('N',N,N,1.0D0,AINV,N,TMP1,1,0.0D0,TMP2,1)
               CALL MATADD(TMP2,F,TMP1,N,1)
               CALL DGEMV('N',N,N,1.0D0,AIPHII,N,TMP1,1,0.0D0,TMP3,1)
               DO 460 J=1,N
                  XP(J)=XO(J)-TMP2(J)*TSAMP+TMP3(J)
 460           CONTINUE
            ELSE
               CALL DGEMV('N',N,N,1.0D0,AIPHII,N,F,1,0.0D0,TMP1,1)
               DO 470 J=1,N
                  XP(J)=XO(J)+TMP1(J)
 470           CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END
!
      SUBROUTINE PROPA2(XO,XP,VXO,VXP,T,TSAMP,NONLIN,ODEEPS,EPSM,CHOL,&
     &                  INFO)
!----------------------------------------------------------------------C
!     ODE SOLUTION OF FORWARD PROPAGATION EQUATIONS.                   C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NONLIN,CHOL,INFO,LIW,LRW
      DOUBLE PRECISION XO(N),XP(N),VXO(N*(N+1)/2),VXP(N*(N+1)/2),T,&
     &                 TSAMP,ODEEPS,EPSM
      PARAMETER        (LIW=20+N*(N+1)/2+N)
      PARAMETER        (LRW=32+9*(N*(N+1)/2+N)+(N*(N+1)/2+N)**2)
      INTEGER          COUNT,I,II,IR,ISTATE,IWORK(LIW),JJ,KK,LL,MF
      DOUBLE PRECISION DELT,RWORK(LRW),TO,TOUT,Y(N*(N+1)/2+N)
      EXTERNAL         RES,JAC
      TO=T
      TOUT=T+TSAMP
!     ..................................................................
!                                            SPECIFY INITIAL CONDITIONS.
!
      DO 10 I=1,N
         Y(I)=XO(I)
  10  CONTINUE
      IF (CHOL.EQ.1) THEN
         CALL MC11CD(VXO,N*(N+1)/2,N)
      END IF
      DO 20 I=1,(N*(N+1)/2)
         Y(N+I)=VXO(I)
  20  CONTINUE
!     ..................................................................
!                                               SPECIFY SOLUTION METHOD.
!
      IF (NONLIN.EQ.3) THEN
         MF=10
      ELSE
         MF=21
      END IF
!     ..................................................................
!                                                    CALCULATE SOLUTION.
!
      ISTATE=1
      COUNT=1
  30  CALL DLSODE(RES,N*(N+1)/2+N,Y,TO,TOUT,1,ODEEPS,ODEEPS,1,ISTATE,0,&
     &            RWORK,LRW,IWORK,LIW,JAC,MF)
      IF (ISTATE.NE.2) THEN
         IF ((ISTATE.EQ.-1).AND.(COUNT.LE.1000)) THEN
            ISTATE=2
            COUNT=COUNT+1
            GOTO 30
         ELSE
            INFO=90
            RETURN
         END IF
      END IF
!     ..................................................................
!                                                        STORE SOLUTION.
!
      DO 40 I=1,N
         XP(I)=Y(I)
  40  CONTINUE
      IF (CHOL.EQ.0) THEN
         DO 50 I=1,(N*(N+1)/2)
            VXP(I)=Y(N+I)
  50     CONTINUE
      ELSE
         DELT=0.0D0
         LL=1
  60     KK=1
         DO 70 JJ=1,N
         DO 70 II=JJ,N
            VXP(KK)=Y(N+KK)
            IF (II.EQ.JJ) VXP(KK)=DABS(VXP(KK))*(1.0D0+DELT)
            KK=KK+1
  70     CONTINUE
         CALL MC11BD(VXP,N*(N+1)/2,N,IR)
         IF (IR.NE.N) THEN
            IF (LL.LE.300) THEN
               DELT=(1.0D1**LL)*EPSM
               LL=LL+1
               GOTO 60
            ELSE
               INFO=30
               RETURN
            END IF
         END IF
      END IF
      RETURN
      END
!
      SUBROUTINE PROPA3(XO,XP,T,TSAMP,NONLIN,ODEEPS,INFO)
!----------------------------------------------------------------------C
!     ODE SOLUTION OF FORWARD PROPAGATION EQUATION FOR MEAN.           C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NONLIN,INFO,LIW,LRW
      DOUBLE PRECISION XO(N),XP(N),T,TSAMP,ODEEPS
      PARAMETER        (LIW=20+N)
      PARAMETER        (LRW=32+9*N+N**2)
      INTEGER          COUNT,I,ISTATE,IWORK(LIW),MF
      DOUBLE PRECISION RWORK(LRW),TO,TOUT,Y(N)
      EXTERNAL         DYNRES,DYNJAC
      TO=T
      TOUT=T+TSAMP
!     ..................................................................
!                                            SPECIFY INITIAL CONDITIONS.
!
      DO 10 I=1,N
         Y(I)=XO(I)
  10  CONTINUE
!     ..................................................................
!                                               SPECIFY SOLUTION METHOD.
!
      IF (NONLIN.EQ.3) THEN
         MF=10
      ELSE
         MF=21
      END IF
!     ..................................................................
!                                                    CALCULATE SOLUTION.
!
      ISTATE=1
      COUNT=1
  20  CALL DLSODA(DYNRES,N,Y,TO,TOUT,1,ODEEPS,ODEEPS,1,ISTATE,0,RWORK,&
     &            LRW,IWORK,LIW,DYNJAC,2)
      IF (ISTATE.NE.2) THEN
         IF ((ISTATE.EQ.-1).AND.(COUNT.LE.1000)) THEN
            ISTATE=2
            COUNT=COUNT+1
            GOTO 20
         ELSE
            INFO=90
            RETURN
         END IF
      END IF
!     ..................................................................
!                                                        STORE SOLUTION.
!
      DO 30 I=1,N
         XP(I)=Y(I)
  30  CONTINUE
      RETURN
      END
!
      SUBROUTINE RES(NEQ,T,Y,YDOT)
!----------------------------------------------------------------------C
!     RHS VECTOR FOR ODE SOLUTION OF FW PROPAGATION EQUATIONS.         C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NEQ
      DOUBLE PRECISION T,Y(NEQ),YDOT(NEQ)
      INTEGER          HN,I,II,J,JJ,K,KK
      DOUBLE PRECISION ALPHA(M+1),F(N),JACOB(N,N),SIGMA(N,N),&
     &                 TO,U(M+1),UO(M+1),XM(NPARAM)
      COMMON           /COM/XM,ALPHA,UO,TO,HN
!$OMP THREADPRIVATE(/COM/)
      IF (HN.EQ.0) THEN
         DO 5 I=1,M+1
            U(I)=UO(I)
   5     CONTINUE
      ELSE
         DO 10 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*(T-TO)
  10     CONTINUE
      END IF
      CALL FVECXJ(XM,Y,U,T,JACOB,F,NPARAM,N,M+1)
!      CALL FVECXJ(XM,Y,U,T,F,NPARAM,N,M+1,JACOB,N,YJACO,LYJAC,
!     $            RFS,IFS,LFS)
      CALL SIGMAT(XM,U,T,SIGMA,NPARAM,N,M+1)
      DO 20 I=1,N
         YDOT(I)=F(I)
  20  CONTINUE
      K=N+1
      DO 30 I=1,N
      DO 30 J=I,N
         YDOT(K)=0.0D0
         DO 40 II=1,N
            JJ=MAX(II,I)
            KK=MIN(II,I)
            YDOT(K)=YDOT(K)+JACOB(J,II)*Y(KK*N+(KK-KK*KK)/2+JJ)
            JJ=MAX(J,II)
            KK=MIN(J,II)
            YDOT(K)=YDOT(K)+JACOB(I,II)*Y(KK*N+(KK-KK*KK)/2+JJ)+&
     &              SIGMA(J,II)*SIGMA(I,II)
  40     CONTINUE
         K=K+1
  30  CONTINUE
      RETURN
      END
!
      SUBROUTINE JAC(NEQ,T,Y,ML,MU,PD,NROWPD)
!----------------------------------------------------------------------C
!     JACOBIAN MATRIX FOR ODE SOLUTION OF FW PROPAGATION EQUATIONS.    C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          ML,MU,NEQ,NROWPD
      DOUBLE PRECISION PD(NROWPD,NEQ),T,Y(NEQ)
      INTEGER          HN,I,II,III,J,JJ,K,KK
      DOUBLE PRECISION ALPHA(M+1),F(N),JACOB(N,N),TO,U(M+1),&
     &                 UO(M+1),XM(NPARAM)
      COMMON           /COM/XM,ALPHA,UO,TO,HN
!$OMP THREADPRIVATE(/COM/)
      IF (HN.EQ.0) THEN
         DO 5 I=1,M+1
            U(I)=UO(I)
   5     CONTINUE
      ELSE
         DO 10 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*(T-TO)
  10     CONTINUE
      END IF
      CALL FVECXJ(XM,Y,U,T,JACOB,F,NPARAM,N,M+1)
!      CALL FVECXJ(XM,Y,U,T,F,NPARAM,N,M+1,JACOB,N,YJACO,LYJAC,&
!     &            RFS,IFS,LFS)
      DO 20 I=1,N
      DO 20 J=1,N
         PD(J,I)=JACOB(J,I)
  20  CONTINUE
      K=N+1
      DO 30 I=1,N
      DO 30 J=I,N
         DO 40 II=1,N
            JJ=MAX(II,I)
            KK=MIN(II,I)
            DO 50 III=N+1,N+N*(N+1)/2
               IF (III.EQ.(KK*N+(KK-KK*KK)/2+JJ)) THEN
                  PD(K,III)=PD(K,III)+JACOB(J,II)
               END IF
  50        CONTINUE
            JJ=MAX(J,II)
            KK=MIN(J,II)
            DO 60 III=N+1,N+N*(N+1)/2
               IF (III.EQ.(KK*N+(KK-KK*KK)/2+JJ)) THEN
                  PD(K,III)=PD(K,III)+JACOB(I,II)
               END IF
  60        CONTINUE
  40     CONTINUE
         K=K+1
  30  CONTINUE
      RETURN
      END
!
      SUBROUTINE DYNRES(NEQ,T,Y,YDOT)
!----------------------------------------------------------------------C
!     RHS VECTOR FOR ODE SOLUTION OF MEAN PROPAGATION EQUATION.        C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NEQ
      DOUBLE PRECISION T,Y(NEQ),YDOT(NEQ)
      INTEGER          HN,I,NONLIN
      DOUBLE PRECISION ALPHA(M+1),F(N),JACOB(N,N),ODEEPS,TF,TO,&
     &                 TSAMP,U(M+1),UO(M+1),XCUR(N),XM(NPARAM),XPREV(N)
      COMMON           /COM/XM,ALPHA,UO,TO,HN
!$OMP THREADPRIVATE(/COM/)
      COMMON           /BWCOM/ODEEPS,TF,TSAMP,XCUR,XPREV,NONLIN
!$OMP THREADPRIVATE(/BWCOM/)
      IF (HN.EQ.0) THEN
         DO 5 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*TSAMP
   5     CONTINUE
      ELSE
         DO 10 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*(T-(TO-TSAMP))
  10     CONTINUE
      END IF
      CALL FVECXJ(XM,Y,U,T,JACOB,F,NPARAM,N,M+1)
      DO 20 I=1,N
         YDOT(I)=F(I)
  20  CONTINUE
      RETURN
      END
!
      SUBROUTINE DYNJAC(NEQ,T,Y,ML,MU,PD,NROWPD)
!----------------------------------------------------------------------C
!     JACOBIAN MATRIX FOR ODE SOLUTION OF MEAN PROPAGATION EQUATION.   C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          ML,MU,NEQ,NROWPD
      DOUBLE PRECISION PD(NROWPD,NEQ),T,Y(NEQ)
      INTEGER          HN,I,J,NONLIN
      DOUBLE PRECISION ALPHA(M+1),F(N),JACOB(N,N),ODEEPS,TF,TO,&
     &                 TSAMP,U(M+1),UO(M+1),XCUR(N),XM(NPARAM),XPREV(N)
      COMMON           /COM/XM,ALPHA,UO,TO,HN
!$OMP THREADPRIVATE(/COM/)
      COMMON           /BWCOM/ODEEPS,TF,TSAMP,XCUR,XPREV,NONLIN
!$OMP THREADPRIVATE(/BWCOM/)
      IF (HN.EQ.0) THEN
         DO 5 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*TSAMP
   5     CONTINUE
      ELSE
         DO 10 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*(T-(TO-TSAMP))
  10     CONTINUE
      END IF
      CALL FVECXJ(XM,Y,U,T,JACOB,F,NPARAM,N,M+1)
      DO 20 I=1,N
      DO 20 J=1,N
         PD(J,I)=JACOB(J,I)
  20  CONTINUE
      RETURN
      END
!
      SUBROUTINE BWFLTR(SP,VSP,SO,VSO,Y,INFO)
!----------------------------------------------------------------------C
!     BACKWARD KALMAN FILTER UPDATE OF STATES AND COVARIANCES.         C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          INFO
      DOUBLE PRECISION SP(N),VSP(N*(N+1)/2),SO(N),VSO(N*(N+1)/2),Y(S)
      INTEGER          HN,I,IR,J,K,MISS,NN,NONLIN,NS
      DOUBLE PRECISION A(N,N),ALPHA(M+1),B(N,M+1),C(S,N),CMISS(S,N),&
     &                 D(S,M+1),E(S,S),EPS(S),F(N),H(S),ODEEPS,R1M(N,N),&
     &                 R2(S*(S+1)/2),R2M(S,S),R2MISS(S*(S+1)/2),T,&
     &                 TF,TMP(S),TMP1(S,S),TMP2(S,S),TMP3(N),TMP4(S,N),&
     &                 TMP5(N,N),TMP6(S),TSAMP,UO(M+1),XCUR(N),&
     &                 XM(NPARAM),XP0(N),XPREV(N)
      COMMON           /COM/XM,ALPHA,UO,T,HN
!$OMP THREADPRIVATE(/COM/)
      COMMON           /BWCOM/ODEEPS,TF,TSAMP,XCUR,XPREV,NONLIN
!$OMP THREADPRIVATE(/BWCOM/)
      MISS=0
      NN=N*(N+1)/2
      NS=S*(S+1)/2
      CALL DEFSYS(XM,XP0,XCUR,UO,T,A,B,C,D,R1M,R2M,F,H,NONLIN,2)
!     ..................................................................
!                          CALCULATE THE MEAN OF THE FULL OUTPUT VECTOR.
!
      CALL DGEMV('N',S,N,-1.0D0,C,S,XCUR,1,1.0D0,H,1)
!     ..................................................................
!                    CHECK FOR MISSING DATA IN THE ACTUAL OUTPUT VECTOR.
!
      MISS=0
      DO 1 I=1,S
         IF (Y(I).GE.1D300) THEN
            IF (MISS.EQ.0) THEN
               DO 2 K=1,S
               DO 2 J=1,S
                  IF (K.EQ.J) THEN
                     E(J,K)=1.0D0
                  ELSE
                     E(J,K)=0.0D0
                  END IF
   2           CONTINUE
            END IF
            MISS=MISS+1
            Y(I)=H(I)
            DO 3 K=1,S
               DO 4 J=I,(S-1)
                  E(J,K)=E(J+1,K)
   4           CONTINUE
               E(S,K)=0.0D0
   3        CONTINUE
         END IF
   1  CONTINUE
!     ..................................................................
!                          CALCULATE THE ACTUAL PREDICTION ERROR VECTOR.
!
      CALL MATSUB(Y,H,EPS,S,1)
!     ..................................................................
!                       HANDLE MISSING DATA IN THE ACTUAL OUTPUT VECTOR.
!
      IF (MISS.EQ.S) THEN
         GOTO 8
      ELSEIF (MISS.GE.1) THEN
         CALL DGEMV('N',S,S,1.0D0,E,S,EPS,1,0.0D0,TMP,1)
         CALL DCOPY(S,TMP,1,EPS,1)
         CALL DGEMM('N','N',S,S,S,1.0D0,E,S,R2M,S,0.0D0,TMP1,S)
         CALL DGEMM('N','T',S,S,S,1.0D0,TMP1,S,E,S,0.0D0,TMP2,S)
         DO 5 I=1,NS
            R2MISS(I)=0.0D0
   5     CONTINUE
         K=1
         DO 6 J=1,(S-MISS)
         DO 6 I=J,(S-MISS)
            R2MISS(K)=TMP2(I,J)
            K=K+1
   6     CONTINUE
         CALL MC11BD(R2MISS,(S-MISS)*(S-MISS+1)/2,S-MISS,IR)
         IF (IR.NE.(S-MISS)) THEN
            INFO=40
            RETURN
         END IF
         CALL DGEMM('N','N',S,N,S,1.0D0,E,S,C,S,0.0D0,CMISS,S)
      ELSE
         K=1
         DO 7 J=1,S
         DO 7 I=J,S
            R2(K)=R2M(I,J)
            K=K+1
   7     CONTINUE
         CALL MC11BD(R2,NS,S,IR)
         IF (IR.NE.S) THEN
            INFO=40
            RETURN
         END IF
      END IF
!     ..................................................................
!                                   CALCULATE THE A POSTERIORI VARIANCE.
!
   8  IF (MISS.EQ.S) THEN
         DO 9 I=1,NN
            VSO(I)=VSP(I)
   9     CONTINUE
      ELSEIF (MISS.GE.1) THEN
         CALL MC11FD(R2MISS,(S-MISS)*(S-MISS+1)/2,S-MISS,S-MISS)
         K=1
         DO 10 J=1,(S-MISS)
         DO 10 I=J,(S-MISS)
            TMP2(I,J)=R2MISS(K)
            TMP2(J,I)=R2MISS(K)
            K=K+1
  10     CONTINUE
         CALL DGEMM('N','N',S-MISS,N,S-MISS,1.0D0,TMP2,S,CMISS,S,0.0D0,&
     &              TMP4,S)
         CALL DGEMM('T','N',N,N,S-MISS,1.0D0,CMISS,S,TMP4,S,0.0D0,TMP5,&
     &              N)
         K=1
         DO 11 J=1,N
         DO 11 I=J,N
            VSO(K)=VSP(K)+TMP5(I,J)
            K=K+1
  11     CONTINUE
         CALL MC11BD(R2MISS,(S-MISS)*(S-MISS+1)/2,S-MISS,IR)
      ELSE
         CALL MC11FD(R2,NS,S,S)
         K=1
         DO 12 J=1,S
         DO 12 I=J,S
            TMP2(I,J)=R2(K)
            TMP2(J,I)=R2(K)
            K=K+1
  12     CONTINUE
         CALL DGEMM('N','N',S,N,S,1.0D0,TMP2,S,C,S,0.0D0,TMP4,S)
         CALL DGEMM('T','N',N,N,S,1.0D0,C,S,TMP4,S,0.0D0,TMP5,N)
         K=1
         DO 13 J=1,N
         DO 13 I=J,N
            VSO(K)=VSP(K)+TMP5(I,J)
            K=K+1
  13     CONTINUE
         CALL MC11BD(R2,NS,S,IR)
      END IF
!     ..................................................................
!                                       CALCULATE THE A POSTERIORI MEAN.
!
      IF (MISS.EQ.S) THEN
         DO 14 I=1,N
            SO(I)=SP(I)
  14     CONTINUE
      ELSE
         IF (MISS.GE.1) THEN
            CALL MC11DD(R2MISS,(S-MISS)*(S-MISS+1)/2,S-MISS,EPS,TMP6)
            CALL DGEMV('T',S,N,1.0D0,CMISS,S,EPS,1,0.0D0,TMP3,1)
         ELSE
            CALL MC11DD(R2,NS,S,EPS,TMP6)
            CALL DGEMV('T',S,N,1.0D0,C,S,EPS,1,0.0D0,TMP3,1)
         END IF
         CALL MATADD(TMP3,SP,SO,N,1)
      END IF
!     ..................................................................
!                              CALCULATE THE A PRIORI MEAN AND VARIANCE.
!
      CALL BWPROP(SO,SP,VSO,VSP,TF-T,TSAMP,NONLIN,ODEEPS,INFO)
      RETURN
      END
!
      SUBROUTINE BWPROP(SO,SP,VSO,VSP,T,TSAMP,NONLIN,ODEEPS,INFO)
!----------------------------------------------------------------------C
!     ODE SOLUTION OF BACKWARD PROPAGATION EQUATIONS.                  C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NONLIN,INFO,LIW,LRW
      DOUBLE PRECISION SO(N),SP(N),VSO(N*(N+1)/2),VSP(N*(N+1)/2),T,&
     &                 TSAMP,ODEEPS
      PARAMETER        (LIW=20+N*(N+1)/2+N)
      PARAMETER        (LRW=32+9*(N*(N+1)/2+N)+(N*(N+1)/2+N)**2)
      INTEGER          COUNT,I,ISTATE,IWORK(LIW),MF
      DOUBLE PRECISION RWORK(LRW),TO,TOUT,Y(N*(N+1)/2+N)
      EXTERNAL         BWRES,BWJAC
      TO=T
      TOUT=T+TSAMP
!     ..................................................................
!                                            SPECIFY INITIAL CONDITIONS.
!
      DO 10 I=1,N
         Y(I)=SO(I)
  10  CONTINUE
      DO 20 I=1,(N*(N+1)/2)
         Y(N+I)=VSO(I)
  20  CONTINUE
!     ..................................................................
!                                               SPECIFY SOLUTION METHOD.
!
      IF (NONLIN.EQ.3) THEN
         MF=10
      ELSE
         MF=21
      END IF
!     ..................................................................
!                                                    CALCULATE SOLUTION.
!
      ISTATE=1
      COUNT=1
  30  CALL DLSODE(BWRES,N*(N+1)/2+N,Y,TO,TOUT,1,ODEEPS,ODEEPS,1,ISTATE,&
     &            0,RWORK,LRW,IWORK,LIW,BWJAC,MF)
      IF (ISTATE.NE.2) THEN
         IF ((ISTATE.EQ.-1).AND.(COUNT.LE.1000)) THEN
            ISTATE=2
            COUNT=COUNT+1
            GOTO 30
         ELSE
            INFO=90
            RETURN
         END IF
      END IF
!     ..................................................................
!                                                        STORE SOLUTION.
!
      DO 40 I=1,N
         SP(I)=Y(I)
  40  CONTINUE
      DO 50 I=1,(N*(N+1)/2)
         VSP(I)=Y(N+I)
  50  CONTINUE
      RETURN
      END
!
      SUBROUTINE BWRES(NEQ,T,Y,YDOT)
!----------------------------------------------------------------------C
!     RHS VECTOR FOR ODE SOLUTION OF BW PROPAGATION EQUATIONS.         C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          NEQ
      DOUBLE PRECISION T,Y(NEQ),YDOT(NEQ)
      INTEGER          HN,I,II,INFO,ISAV(37),J,JJ,JJJ,K,KK,KKK,&
     &                 LL,NONLIN
      DOUBLE PRECISION ALPHA(M+1),F(N),JACOB(N,N),ODEEPS,R1M(N,N),&
     &                 RSAV(218),SIGMA(N,N),TF,TO,TSAMP,U(M+1),&
     &                 UO(M+1),XCUR(N),XF(N),XM(NPARAM),XPREV(N)
      COMMON           /COM/XM,ALPHA,UO,TO,HN
!$OMP THREADPRIVATE(/COM/)
      COMMON           /BWCOM/ODEEPS,TF,TSAMP,XCUR,XPREV,NONLIN
!$OMP THREADPRIVATE(/BWCOM/)
      IF (HN.EQ.0) THEN
         DO 5 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*TSAMP
   5     CONTINUE
      ELSE
         DO 10 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*(T-(TF-TO))
  10     CONTINUE
      END IF
      CALL DSRCOM(RSAV,ISAV,1)
      CALL PROPA3(XPREV,XF,TO-TSAMP,TF-T-TO+TSAMP,NONLIN,ODEEPS,INFO)
      CALL DSRCOM(RSAV,ISAV,2)
      CALL FVECXJ(XM,XF,U,TF-T,JACOB,F,NPARAM,N,M+1)
      CALL SIGMAT(XM,U,TF-T,SIGMA,NPARAM,N,M+1)
      CALL DGEMM('N','T',N,N,N,1.0D0,SIGMA,N,SIGMA,N,0.0D0,R1M,N)
      CALL DGEMV('N',N,N,-1.0D0,JACOB,N,XF,1,1.0D0,F,1)
      CALL DGEMV('N',N,N,1.0D0,R1M,N,Y,1,0.0D0,XF,1)
      DO 15 I=1,N
         YDOT(I)=0.0D0
         DO 20 J=1,N
            JJ=MAX(I,J)
            KK=MIN(I,J)
            YDOT(I)=YDOT(I)+JACOB(J,I)*Y(J)&
     &              -Y(KK*N+(KK-KK*KK)/2+JJ)*XF(J)&
     &              -Y(KK*N+(KK-KK*KK)/2+JJ)*F(J)
  20     CONTINUE
  15  CONTINUE
      K=N+1
      DO 30 I=1,N
      DO 30 J=I,N
         YDOT(K)=0.0D0
         DO 40 II=1,N
            JJ=MAX(J,II)
            KK=MIN(J,II)
            YDOT(K)=YDOT(K)+JACOB(II,I)*Y(KK*N+(KK-KK*KK)/2+JJ)
            DO 50 LL=1,N
               JJJ=MAX(LL,I)
               KKK=MIN(LL,I)
               YDOT(K)=YDOT(K)-Y(KK*N+(KK-KK*KK)/2+JJ)*R1M(II,LL)&
     &                 *Y(KKK*N+(KKK-KKK*KKK)/2+JJJ)
  50        CONTINUE
            JJ=MAX(II,I)
            KK=MIN(II,I)
            YDOT(K)=YDOT(K)+JACOB(II,J)*Y(KK*N+(KK-KK*KK)/2+JJ)
  40     CONTINUE
         K=K+1
  30  CONTINUE
      RETURN
      END
!
      SUBROUTINE BWJAC(NEQ,T,Y,ML,MU,PD,NROWPD)
!----------------------------------------------------------------------C
!     JACOBIAN MATRIX FOR ODE SOLUTION OF BW PROPAGATION EQUATIONS.    C
!----------------------------------------------------------------------C
      INCLUDE          'global.h'
      INTEGER          ML,MU,NEQ,NROWPD
      DOUBLE PRECISION PD(NROWPD,NEQ),T,Y(NEQ)
      INTEGER          HN,I,II,III,INFO,ISAV(37),J,JJ,JJJ,K,KK,&
     &                 KKK,LL,NONLIN
      DOUBLE PRECISION ALPHA(M+1),F(N),JACOB(N,N),ODEEPS,R1M(N,N),&
     &                 RSAV(218),SIGMA(N,N),TF,TO,TSAMP,U(M+1),&
     &                 UO(M+1),XCUR(N),XF(N),XM(NPARAM),XPREV(N)
      COMMON           /COM/XM,ALPHA,UO,TO,HN
!$OMP THREADPRIVATE(/COM/)
      COMMON           /BWCOM/ODEEPS,TF,TSAMP,XCUR,XPREV,NONLIN
!$OMP THREADPRIVATE(/BWCOM/)
      IF (HN.EQ.0) THEN
         DO 5 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*TSAMP
   5     CONTINUE
      ELSE
         DO 10 I=1,M+1
            U(I)=UO(I)+ALPHA(I)*(T-(TF-TO))
  10     CONTINUE
      END IF
      CALL DSRCOM(RSAV,ISAV,1)
      CALL PROPA3(XPREV,XF,TO-TSAMP,TF-T-TO+TSAMP,NONLIN,ODEEPS,INFO)
      CALL DSRCOM(RSAV,ISAV,2)
      CALL FVECXJ(XM,XF,U,TF-T,JACOB,F,NPARAM,N,M+1)
      CALL SIGMAT(XM,U,TF-T,SIGMA,NPARAM,N,M+1)
      CALL DGEMM('N','T',N,N,N,1.0D0,SIGMA,N,SIGMA,N,0.0D0,R1M,N)
      CALL DGEMV('N',N,N,-1.0D0,JACOB,N,XF,1,1.0D0,F,1)
      CALL DGEMV('N',N,N,1.0D0,R1M,N,Y,1,0.0D0,XF,1)
      DO 15 I=1,N
         DO 20 J=1,N
            PD(I,J)=JACOB(J,I)
            DO 25 K=1,N
               JJ=MAX(I,K)
               KK=MIN(I,K)
               PD(I,J)=PD(I,J)-Y(KK*N+(KK-KK*KK)/2+JJ)*R1M(K,J)
  25        CONTINUE
  20     CONTINUE
         II=N+1
         DO 30 K=1,N
         DO 30 J=K,N
            IF (I.EQ.J) THEN
               PD(I,II)=-XF(K)-F(K)
            END IF
            II=II+1
  30     CONTINUE
  15  CONTINUE
      K=N+1
      DO 35 I=1,N
      DO 35 J=I,N
         DO 40 II=1,N
            JJ=MAX(J,II)
            KK=MIN(J,II)
            DO 45 III=N+1,N+N*(N+1)/2
               IF (III.EQ.(KK*N+(KK-KK*KK)/2+JJ)) THEN
                  PD(K,III)=PD(K,III)+JACOB(II,I)
               END IF
  45        CONTINUE
            DO 50 LL=1,N
               JJJ=MAX(LL,I)
               KKK=MIN(LL,I)
               DO 55 III=N+1,N+N*(N+1)/2
                  IF (III.EQ.(KK*N+(KK-KK*KK)/2+JJ)) THEN
                     IF (III.EQ.(KKK*N+(KKK-KKK*KKK)/2+JJJ)) THEN
                        PD(K,III)=PD(K,III)-2.0D0*R1M(II,LL)&
     &                            *Y(KKK*N+(KKK-KKK*KKK)/2+JJJ)
                     ELSE
                        PD(K,III)=PD(K,III)-R1M(II,LL)&
     &                            *Y(KKK*N+(KKK-KKK*KKK)/2+JJJ)
                     END IF
                  ELSE
                     IF (III.EQ.(KKK*N+(KKK-KKK*KKK)/2+JJJ)) THEN
                        PD(K,III)=PD(K,III)-Y(KK*N+(KK-KK*KK)/2+JJ)&
     &                            *R1M(II,LL)
                     END IF
                  END IF
  55           CONTINUE
  50        CONTINUE
            JJ=MAX(II,I)
            KK=MIN(II,I)
            DO 60 III=N+1,N+N*(N+1)/2
               IF (III.EQ.(KK*N+(KK-KK*KK)/2+JJ)) THEN
                  PD(K,III)=PD(K,III)+JACOB(II,J)
               END IF
  60        CONTINUE
  40     CONTINUE
         K=K+1
  35  CONTINUE
      RETURN
      END

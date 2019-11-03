      SUBROUTINE FDF(NX,X,DF,F,OD2,OD3,MD,INTS,NINT,DOUBLS,NDOUB,&
     &               TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,&
     &               NMISST,FPEN,FPRIOR,EPSM,VINFO,INFO)
!-----------------------------------------------------------------------C
!     CALCULATION OF THE FUNCTION VALUE F AND THE GRADIENT DF OF F      C
!     AT THE POINT X(NX).                                               C
!                                                                       C
!     MD = 1  =>   FORWARD DIFFERENCE GRADIENT APPROXIMATION            C
!     MD = 2  =>   CENTRAL DIFFERENCE GRADIENT APPROXIMATION            C
!-----------------------------------------------------------------------C

      implicit none

      INTEGER          NX,MD,NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,&
     &                 NSET,NOBS(NSET),NMISST,VINFO(NX),INFO
      DOUBLE PRECISION X(NX),DF(NX),F,OD2,OD3,DOUBLS(NDOUB),TMAT(NTMAT),&
     &                 IMAT(NIMAT),OMAT(NOMAT),FPEN,FPRIOR,EPSM
      LOGICAL          FIRST
      INTEGER          I,ISTOP
      DOUBLE PRECISION XMIN,EPS,SFMIN,BASE,T,RND,EMIN,RMIN,EMAX,RMAX,&
     &                 PREC
      COMMON           /SPECS1/FIRST
      SAVE             /SPECS1/
!$OMP THREADPRIVATE(/SPECS1/)
      COMMON           /SPECS2/EPS,SFMIN,BASE,T,RND,EMIN,RMIN,&
     &                         EMAX,RMAX,PREC
      SAVE             /SPECS2/
!$OMP THREADPRIVATE(/SPECS2/)

      interface

      SUBROUTINE LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,&
     &                 NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,FPEN,FPRIOR,&
     &                 EPSM,JOB,FUNK,INFO,VX0)


      INTEGER          NX,NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,NSET,&
     &                 NOBS(NSET),NMISST,JOB,INFO
      DOUBLE PRECISION X(NX),DOUBLS(NDOUB),TMAT(NTMAT),IMAT(NIMAT),&
     &                 OMAT(NOMAT),EPSM,FUNK
      DOUBLE PRECISION FPEN,FPRIOR
      double precision, optional, dimension(:), intent(in) :: vx0

      end SUBROUTINE LLIKE
      end interface
!
      XMIN = 1D-1
      CALL LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,&
     &           TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,&
     &           NMISST,FPEN,FPRIOR,EPSM,0,F,INFO)
      IF (INFO.NE.0) RETURN
      IF (MD.EQ.1) THEN
!
!     FORWARD DIFFERENCE APPROXIMATION TO GRADIENT.
!
!$OMP PARALLEL DO PRIVATE(I,NMISST) FIRSTPRIVATE(X) SHARED(DF,VINFO)
      DO 4 I=1,NX
         CALL FWDIFF(NX,X,I,XMIN,OD2,F,DF,INTS,NINT,DOUBLS,NDOUB,&
     &               TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,&
     &               EPSM,VINFO)
   4  CONTINUE
!$OMP END PARALLEL DO
      ELSE
!
!     CENTRAL DIFFERENCE APPROXIMATION TO GRADIENT.
!
!$OMP PARALLEL DO PRIVATE(I) FIRSTPRIVATE(X) SHARED(DF,VINFO)
      DO 5 I=1,NX
         CALL CTDIFF(NX,X,I,XMIN,OD3,DF,INTS,NINT,DOUBLS,NDOUB,&
     &               TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,&
     &               EPSM,VINFO)
   5  CONTINUE
!$OMP END PARALLEL DO
      END IF
!
!     CHECK TO SEE IF EVERYTHING IS OK - REPORT IF NOT
!
      DO 10 I=1,NX
         IF (VINFO(I).NE.0) THEN
            INFO=I*100+VINFO(I)
            RETURN
         END IF
  10  CONTINUE
      RETURN
      END
!
      SUBROUTINE CTDIFF(NX,X,I,XMIN,OD3,DF,INTS,NINT,DOUBLS,NDOUB,&
     &                  TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,&
     &                  EPSM,VINFO)
!-----------------------------------------------------------------------C
!     CALCULATION OF THE GRADIENT DF OF F AT THE POINT X(NX).           C
!     CENTRAL DIFFERENCE APPROXIMATION.                                 C
!-----------------------------------------------------------------------C
      INTEGER NX,I,NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,NSET,&
     &        NOBS(NSET),NMISST,VINFO(NX)
      DOUBLE PRECISION X(NX),XMIN,OD3,DF(NX),DOUBLS(NDOUB),TMAT(NTMAT),&
     &                 IMAT(NIMAT),OMAT(NOMAT),EPSM
      INTEGER JOB,INFO
      DOUBLE PRECISION C,H,FF,FB,FPEN,FPRIOR

      interface

      SUBROUTINE LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,&
     &                 NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,FPEN,FPRIOR,&
     &                 EPSM,JOB,FUNK,INFO,VX0)


      INTEGER          NX,NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,NSET,&
     &                 NOBS(NSET),NMISST,JOB,INFO
      DOUBLE PRECISION X(NX),DOUBLS(NDOUB),TMAT(NTMAT),IMAT(NIMAT),&
     &                 OMAT(NOMAT),EPSM,FUNK
      DOUBLE PRECISION FPEN,FPRIOR
      double precision, optional, dimension(:), intent(in) :: vx0

      end SUBROUTINE LLIKE
      end interface

         JOB=0
         C=X(I)
         IF (DABS(C).GT.XMIN) THEN
            H=C*OD3
         ELSE
            H=XMIN*OD3
         END IF
         X(I)=C+H
         CALL LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,&
     &              TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,&
     &              FPEN,FPRIOR,EPSM,JOB,FF,INFO)
         VINFO(I)=INFO
         IF (INFO.NE.0) RETURN
         X(I)=C-H
         CALL LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,&
     &              TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,&
     &              FPEN,FPRIOR,EPSM,JOB,FB,INFO)
         VINFO(I)=INFO
         IF (INFO.NE.0) RETURN
         DF(I)=(FF-FB)/(2D0*H)
         X(I)=C
      RETURN
      END
!
      SUBROUTINE FWDIFF(NX,X,I,XMIN,OD2,F,DF,INTS,NINT,DOUBLS,NDOUB,&
     &                  TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,&
     &                  EPSM,VINFO)
!-----------------------------------------------------------------------C
!     CALCULATION OF THE GRADIENT DF OF F AT THE POINT X(NX).           C
!     FORWARD DIFFERENCE APPROXIMATION.                                 C
!-----------------------------------------------------------------------C
      INTEGER NX,I,NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,NSET,&
     &        NOBS(NSET),NMISST,VINFO(NX)
      DOUBLE PRECISION X(NX),XMIN,OD2,DF(NX),F,DOUBLS(NDOUB),&
     &                 TMAT(NTMAT),IMAT(NIMAT),OMAT(NOMAT),EPSM
      INTEGER JOB,INFO
      DOUBLE PRECISION C,H,FF,FPEN,FPRIOR

      interface

      SUBROUTINE LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,TMAT,NTMAT,IMAT,&
     &                 NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,FPEN,FPRIOR,&
     &                 EPSM,JOB,FUNK,INFO,VX0)


      INTEGER          NX,NINT,INTS(NINT),NDOUB,NTMAT,NIMAT,NOMAT,NSET,&
     &                 NOBS(NSET),NMISST,JOB,INFO
      DOUBLE PRECISION X(NX),DOUBLS(NDOUB),TMAT(NTMAT),IMAT(NIMAT),&
     &                 OMAT(NOMAT),EPSM,FUNK
      DOUBLE PRECISION FPEN,FPRIOR
      double precision, optional, dimension(:), intent(in) :: vx0

      end SUBROUTINE LLIKE
      end interface

         JOB=0
         C=X(I)
         IF (DABS(C).GT.XMIN) THEN
            H=C*OD2
         ELSE
            H=DSIGN(XMIN,C)*OD2
         END IF
         X(I)=C-H
         CALL LLIKE(NX,X,INTS,NINT,DOUBLS,NDOUB,&
     &              TMAT,NTMAT,IMAT,NIMAT,OMAT,NOMAT,NOBS,NSET,NMISST,&
     &              FPEN,FPRIOR,EPSM,JOB,FF,INFO)
         VINFO(I)=INFO
         IF (INFO.NE.0) RETURN
         DF(I)=(F-FF)/H
         X(I)=C
      RETURN
      END

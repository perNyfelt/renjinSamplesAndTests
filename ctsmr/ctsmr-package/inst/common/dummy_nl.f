c
c This file contains dummy functions for the system matrices used
c in for the linear Kalman filter.
c 
      SUBROUTINE AMAT(XM,XP,UO,T,A,NPARAM,N,M)
      INTEGER NPARAM,N,M
      DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,A(N,N)
      INTEGER I,J
      DO 10 I=1,N
      DO 10 J=1,N
         A(J,I)=0.0D0
  10  CONTINUE
      RETURN
      END
c
      SUBROUTINE BMAT(XM,XP,UO,T,B,NPARAM,N,M)
      INTEGER NPARAM,N,M
      DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,B(N,M)
      INTEGER I,J
      DO 10 I=1,M
      DO 10 J=1,N
         B(J,I)=0.0D0
  10  CONTINUE
      RETURN
      END
c
      SUBROUTINE CMAT(XM,XP,UO,T,C,NPARAM,N,M,S)
      INTEGER NPARAM,N,M,S
      DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,C(S,N)
      INTEGER I,J
      DO 10 I=1,N
      DO 10 J=1,S
         C(J,I)=0.0D0
  10  CONTINUE
      RETURN
      END
c
      SUBROUTINE DMAT(XM,XP,UO,T,D,NPARAM,N,M,S)
      INTEGER NPARAM,N,M,S
      DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,D(S,M)
      INTEGER I,J
      DO 10 I=1,M
      DO 10 J=1,S
         D(J,I)=0.0D0
  10  CONTINUE
      RETURN
      END

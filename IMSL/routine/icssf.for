C   IMSL ROUTINE NAME   - ICSSF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ICSSCV
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSSF  (N,RO,H,Z,D,A)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               RO,H,Z(1),D(1),A(N,2,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,N1,ITEMP
      REAL               A11,A12,A1,A21,A22,A2,ALFA,B1,BETA,D1,DET,
     *                   Z1,XB,YB,XV,YV,X,AUX,AO,BO,TOL
C                                  CALCULATE THE CUBIC SPLINE TO
C                                    SMOOTH THE GIVEN DATA Z, WITH
C                                    SMOOTHING PARAMETER RO.
C                                  FIRST EXECUTABLE STATEMENT
      DATA               TOL/Z3C100000/
      ALFA = H*H*H*RO/6.0
      IF (ALFA.LT.10.*TOL) GO TO 20
      DET = 1.0/(1.0+ALFA+ALFA)
      Z1 = ALFA*Z(1)*DET
      BETA = 4.0+ALFA
      Z(1) = 2.0*Z1
      D(1) = -3.0*Z1
      A(1,1,1) = -DET
      A(1,1,2) = DET
      A(1,2,1) = -3.0*ALFA*DET
      A(1,2,2) = (ALFA-1.0)*DET
      N1 = N-1
      DO 5 I=2,N1
         ITEMP = I-1
         A1 = 2.0*A(ITEMP,1,1)+A(ITEMP,2,1)
         A2 = 2.0*A(ITEMP,1,2)+A(ITEMP,2,2)
         B1 = 2.0*Z(ITEMP)+D(ITEMP)
         A22 = BETA+A1
         A12 = A2
         A21 = -A1-A(ITEMP,1,1)
         A11 = 4.0-A2-A(ITEMP,1,2)
         Z1 = ALFA*Z(I)+B1
         D1 = -(B1+Z(ITEMP))
         DET = 1.0/(A11*A22-A12*A21)
         A11 = A11*DET
         A12 = -A12*DET
         A21 = -A21*DET
         A22 = A22*DET
         A(I,1,1) = -2.0*A11-3.0*A12
         A(I,1,2) = A11+A12
         A(I,2,1) = -2.0*A21-3.0*A22
         A(I,2,2) = A21+A22
         Z(I) = A11*Z1+A12*D1
         D(I) = A21*Z1+A22*D1
    5 CONTINUE
      A1 = 2.0*A(N1,1,1)+A(N1,2,1)
      A2 = 2.0*A(N1,1,2)+A(N1,2,2)
      B1 = 2.0*Z(N1)+D(N1)
      Z1 = ALFA*Z(N)+B1
      A22 = ALFA+2.0+A1
      A12 = -1.0+A2
      A21 = -3.0-A1-A(N1,1,1)
      A11 = 2.0-A2-A(N1,1,2)
      D1 = -(B1+Z(N1))
      DET = 1.0/(A11*A22-A12*A21)
      A11 = A11*DET
      A12 = -A12*DET
      A21 = -A21*DET
      A22 = A22*DET
      Z(N) = A11*Z1+A12*D1
      D(N) = A21*Z1+A22*D1
      DO 10 J=1,N1
         I = N-J
         ITEMP = I+1
         Z(I) = Z(I)-(A(I,1,1)*Z(ITEMP)+A(I,1,2)*D(ITEMP))
         D(I) = D(I)-(A(I,2,1)*Z(ITEMP)+A(I,2,2)*D(ITEMP))
   10 CONTINUE
      DET = 1.0/H
      DO 15 I=1,N
         D(I) = D(I)*DET
   15 CONTINUE
      RETURN
C                                   STRAIGHT LINE CASE
   20 XB = (N+1)/2.*H
      YB = 0.
      YV = 0.
      XV = 0.
      X = 0.
      DO 25 I=1,N
         X = X+H
         AUX = X-XB
         XV = XV+AUX*AUX
         YB = YB+Z(I)
         YV = YV+Z(I)*AUX
   25 CONTINUE
      XV = XV/N
      YV = YV/N
      YB = YB/N
      AO = YV/XV
      BO = YB-AO*XB
      X = 0.
      DO 30 I=1,N
         D(I) = AO
         X = X+H
         Z(I) = AO*X+BO
   30 CONTINUE
      RETURN
      END

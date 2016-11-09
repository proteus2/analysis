C   IMSL ROUTINE NAME   - ZRPQLE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   REQD. IMSL ROUTINES - ZRPQLH
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLE (ITYPE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITYPE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN
      REAL               ARE,ETA,RMRE
      REAL               P(101),QP(101),RK(101),QK(101),SVK(101)
      REAL               SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
C                                  THIS ROUTINE CALCULATES SCALAR
C                                    QUANTITIES USED TO COMPUTE THE
C                                    NEXT K POLYNOMIAL AND NEW
C                                    ESTIMATES OF THE QUADRATIC
C                                    COEFFICIENTS
C                                  ITYPE - INTEGER VARIABLE SET HERE
C                                    INDICATING HOW THE CALCULATIONS
C                                    ARE NORMALIZED TO AVOID OVERFLOW
C                                  SYNTHETIC DIVISION OF K BY THE
C                                    QUADRATIC 1,U,V
C                                  FIRST EXECUTABLE STATEMENT
      CALL ZRPQLH (N,U,V,RK,QK,C,D)
      IF (ABS(C).GT.ABS(RK(N))*100.*ETA) GO TO 5
      IF (ABS(D).GT.ABS(RK(N-1))*100.*ETA) GO TO 5
      ITYPE = 3
C                                  TYPE=3 INDICATES THE QUADRATIC IS
C                                    ALMOST A FACTOR OF K
      RETURN
    5 IF (ABS(D).LT.ABS(C)) GO TO 10
      ITYPE = 2
C                                  TYPE=2 INDICATES THAT ALL FORMULAS
C                                    ARE DIVIDED BY D
      E = RA/D
      F = C/D
      G = U*RB
      H = V*RB
      A3 = (RA+G)*E+H*(RB/D)
      A1 = RB*F-RA
      A7 = (F+U)*RA+H
      RETURN
   10 ITYPE = 1
C                                  TYPE=1 INDICATES THAT ALL FORMULAS
C                                    ARE DIVIDED BY C
      E = RA/C
      F = D/C
      G = U*E
      H = V*RB
      A3 = RA*E+(H/C+G)*RB
      A1 = RB-RA*(D/C)
      A7 = RA+G*D+H*F
      RETURN
      END

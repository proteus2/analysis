C   IMSL ROUTINE NAME   - ZRPQLG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLG (ITYPE,UU,VV)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITYPE
      REAL               UU,VV
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN
      REAL               ARE,ETA,RMRE
      REAL               P(101),QP(101),RK(101),QK(101),SVK(101)
      REAL               SR,SI,U,V,RA,RB,C,D,A1,A2,A3,ONE,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   A4,A5,B1,B2,C1,C2,C3,C4,TEMP,ZERO
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
      DATA               ZERO/0.0E0/, ONE/1.0E0/
C                                  COMPUTE NEW ESTIMATES OF THE
C                                    QUADRATIC COEFFICIENTS USING THE
C                                    SCALARS COMPUTED IN ZRPQLE
C                                  USE FORMULAS APPROPRIATE TO SETTING
C                                    OF TYPE.
C                                  FIRST EXECUTABLE STATEMENT
      IF (ITYPE.EQ.3) GO TO 15
      IF (ITYPE.EQ.2) GO TO 5
      A4 = RA+U*RB+H*F
      A5 = C+(U+V*F)*D
      GO TO 10
    5 A4 = (RA+G)*F+H
      A5 = (F+U)*C+V*D
C                                  EVALUATE NEW QUADRATIC COEFFICIENTS.
C
   10 B1 = -RK(N)/P(NN)
      B2 = -(RK(N-1)+B1*P(N))/P(NN)
      C1 = V*B2*A1
      C2 = B1*A7
      C3 = B1*B1*A3
      C4 = C1-C2-C3
      TEMP = A5+B1*A4-C4
      IF (ONE*TEMP.EQ.ZERO) GO TO 15
      UU = U-(U*(C3+C2)+V*(B1*A1+B2*A7))/TEMP
      VV = V*(1+C4/TEMP)
      RETURN
C                                  IF TYPE=3 THE QUADRATIC IS ZEROED
   15 UU = ZERO
      VV = ZERO
      RETURN
      END

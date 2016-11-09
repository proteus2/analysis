C   IMSL ROUTINE NAME   - ZRPQLI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
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
      SUBROUTINE ZRPQLI (RA,B1,C,SR,SI,RLR,RLI)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               RA,B1,C,SR,SI,RLR,RLI
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               RB,D,E,ZERO,ONE,TWO
      DATA               ZERO,ONE,TWO/0.0,1.0,2.0/
C                                  CALCULATE THE ZEROS OF THE QUADRATIC
C                                    A*Z**2 + B1*Z + C. THE QUADRATIC
C                                    FORMULA, MODIFIED TO AVOID
C                                    OVERFLOW, IS USED TO FIND THE
C                                    LARGER ZERO IF THE ZEROS ARE REAL
C                                    AND BOTH ZEROS ARE COMPLEX.
C                                  THE SMALLER REAL ZERO IS FOUND
C                                    DIRECTLY FROM THE PRODUCT OF THE
C                                    ZEROS C/A
C                                  FIRST EXECUTABLE STATEMENT
      IF (RA.NE.ZERO) GO TO 10
      SR = ZERO
      IF (B1.NE.ZERO) SR = -C/B1
      RLR = ZERO
    5 SI = ZERO
      RLI = ZERO
      RETURN
   10 IF (C.NE.ZERO) GO TO 15
      SR = ZERO
      RLR = -B1/RA
      GO TO 5
C                                  COMPUTE DISCRIMINANT AVOIDING
C                                    OVERFLOW
   15 RB = B1/TWO
      IF (ABS(RB).LT.ABS(C)) GO TO 20
      E = ONE-(RA/RB)*(C/RB)
      D = SQRT(ABS(E))*ABS(RB)
      GO TO 25
   20 E = RA
      IF (C.LT.ZERO) E = -RA
      E = RB*(RB/ABS(C))-E
      D = SQRT(ABS(E))*SQRT(ABS(C))
   25 IF (E.LT.ZERO) GO TO 30
C                                  REAL ZEROS
      IF (RB.GE.ZERO) D = -D
      RLR = (-RB+D)/RA
      SR = ZERO
      IF (RLR.NE.ZERO) SR = (C/RLR)/RA
      GO TO 5
C                                  COMPLEX CONJUGATE ZEROS
   30 SR = -RB/RA
      RLR = SR
      SI = ABS(D/RA)
      RLI = -SI
      RETURN
      END

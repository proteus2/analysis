C   IMSL ROUTINE NAME   - ZCPQLL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION ZCPQLL (CR,CI)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               CR,CI
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               AR,AI,ONE,RSQ2
      DATA               ONE/1.0/
      DATA               RSQ2/1.414214/
C                                  FIRST EXECUTABLE STATEMENT
      AR = ABS(CR)
      AI = ABS(CI)
C                                  MODULUS OF A COMPLEX NUMBER AVOIDING
C                                    OVERFLOW
      IF (AR.GE.AI) GO TO 5
      ZCPQLL = AI*SQRT(ONE+(AR/AI)**2)
      RETURN
    5 IF (AR.LE.AI) GO TO 10
      ZCPQLL = AR*SQRT(ONE+(AI/AR)**2)
      RETURN
   10 ZCPQLL = AR*RSQ2
      RETURN
      END

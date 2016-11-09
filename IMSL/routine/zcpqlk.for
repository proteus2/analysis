C   IMSL ROUTINE NAME   - ZCPQLK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - ZCPQLM
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPQLK (AR,AI,BR,BI,CR,CI)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               AR,AI,BR,BI,CR,CI
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               R,D,T,RINFP,ZERO,ONE
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (BR.NE.ZERO.OR.BI.NE.ZERO) GO TO 5
C                                  COMPLEX DIVISION C = A/B, AVOIDING
C                                    OVERFLOW
C                                  DIVISION BY ZERO, C = INFINITY
      CALL ZCPQLM (T,RINFP,T,T)
      CR = RINFP
      CI = RINFP
      RETURN
    5 IF (ABS(BR).GE.ABS(BI)) GO TO 10
      R = BR/BI
      D = BI+R*BR
      D = ONE/D
      CR = (AR*R+AI)*D
      CI = (AI*R-AR)*D
      RETURN
   10 R = BI/BR
      D = BR+R*BI
      D = ONE/D
      CR = (AR+AI*R)*D
      CI = (AI-AR*R)*D
      RETURN
      END

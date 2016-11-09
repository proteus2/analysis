C   IMSL ROUTINE NAME   - CTPR1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE CTPR
C
C   USAGE               - FUNCTION CTPR1(N)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION CTPR1 (N)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IFLAG,I
      DOUBLE PRECISION   HALF,ONE,X,XP1,RTPILG,TS,TWL
      DATA               RTPILG/.9189385332046727D0/
      DATA               HALF/0.5D0/,ONE/1.0D0/,TWL/12.0D0/,TS/360.0D0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  USE STIRLINGS APPROXIMATION
      X = N
      XP1 = X + ONE
      CTPR1= (X + HALF) * DLOG(XP1) - XP1 + RTPILG
     *       + ONE/(TWL*XP1) - ONE/(TS*XP1*XP1*XP1)
      RETURN
      END

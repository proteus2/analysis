C   IMSL ROUTINE NAME   - GGNQF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NORMAL RANDOM DEVIATE GENERATOR - FUNCTION
C                           FORM OF GGNML
C
C   USAGE               - FUNCTION GGNQF (DSEED)
C
C   ARGUMENTS    GGNQF  - RESULTANT NORMAL (0,1) DEVIATE.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION GGNQF (DSEED)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IER
      REAL               GGNQFX,XGGNQF
      DOUBLE PRECISION   D2P31M,D2PN31
C                                  D2P31M = (2**31) - 1
C                                  D2PN31 = (2**31)
      DATA               D2P31M/2147483647.D0/
      DATA               D2PN31/2147483648.D0/
C                                  GENERATE A RANDOM (0,1) DEVIATE.
C                                  16807 = (7**5)
C                                  FIRST EXECUTABLE STATEMENT
      DSEED = DMOD(16807.D0*DSEED,D2P31M)
      GGNQFX = DSEED / D2PN31
C                                  TRANSFORM TO NORMAL DEVIATE
      CALL MDNRIS(GGNQFX,XGGNQF,IER)
      GGNQF = XGGNQF
      RETURN
      END

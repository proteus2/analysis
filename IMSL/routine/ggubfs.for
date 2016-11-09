C   IMSL ROUTINE NAME   - GGUBFS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - BASIC UNIFORM (0,1) RANDOM NUMBER GENERATOR -
C                           FUNCTION FORM OF GGUBS
C
C   USAGE               - FUNCTION GGUBFS (DSEED)
C
C   ARGUMENTS    GGUBFS - RESULTANT DEVIATE.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION GGUBFS (DSEED)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   D2P31M,D2P31
C                                  D2P31M=(2**31) - 1
C                                  D2P31 =(2**31)(OR AN ADJUSTED VALUE)
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31 /2147483648.D0/
C                                  FIRST EXECUTABLE STATEMENT
      DSEED = DMOD(16807.D0*DSEED,D2P31M)
      GGUBFS = DSEED / D2P31
      RETURN
      END

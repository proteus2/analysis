C   IMSL ROUTINE NAME   - GGUBT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - UNIFORM (0,1) PSEUDO-RANDOM NUMBER GENERATOR
C                           USING ALTERNATE MULTIPLIER.
C
C   USAGE               - CALL GGUBT (DSEED,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0,2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF DEVIATES TO BE GENERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           PSEUDO-RANDOM UNIFORM (0,1) DEVIATES
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGUBT (DSEED,NR,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               R(NR)
      DOUBLE PRECISION   DSEED,DSEED1,DSEED2
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      DOUBLE PRECISION   D2P31M,D2P31
C                                  D2P31M=(2**31) - 1
C                                  D2P31 =(2**31)(OR AN ADJUSTED VALUE)
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31/2147483648.D0/
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I=1,NR
         DSEED1 = DMOD(32768.D0*DSEED,D2P31M)
         DSEED2 = DMOD(23166.D0*DSEED,D2P31M)
         DSEED  = DMOD(12121.D0*DSEED1+DSEED2,D2P31M)
         R(I)   = DSEED / D2P31
    5 CONTINUE
      RETURN
      END

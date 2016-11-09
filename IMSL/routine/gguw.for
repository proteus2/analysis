C   IMSL ROUTINE NAME   - GGUW
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - UNIFORM (0,1) RANDOM NUMBER GENERATOR WITH
C                           SHUFFLING
C
C   USAGE               - CALL GGUW (DSEED,NR,IOPT,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT.  NUMBER OF DEVIATES TO BE GENERATED
C                           ON THIS CALL.
C                IOPT   - INPUT.  OPTION SWITCH.
C                           ON INITIAL ENTRY, IOPT MUST BE ONE.
C                           ON SUBSEQUENT CALLS, IOPT MUST BE ZERO.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           SHUFFLED FLOATING POINT (0,1) DEVIATES.
C
C   REQD. IMSL ROUTINES - NONE
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGUW   (DSEED,NR,IOPT,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,IOPT
      REAL               R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J
      REAL               WK(128)
      DOUBLE PRECISION   D2P31M,D2P31
C                                  D2P31M=(2**31) - 1
C                                  D2P31 =(2**31)(OR AN ADJUSTED VALUE)
      DATA               WK/128*0.0/
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31 /2147483648.D0/
C                                  INITIAL ENTRY - FILL TABLE
C                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT .NE. 1)   GO TO 10
      DO 5  I = 1,128
         DSEED= DMOD( 16807.D0*DSEED,D2P31M )
    5 WK(I)   = DSEED / D2P31
   10 CONTINUE
      DO 15 I = 1,NR
C                                  GENERATE A NEW UNIFORM (0,1) DEVIATE
C                                  16807 = (7**5)
         DSEED= DMOD(16807.D0*DSEED,D2P31M)
C                                  NOW SHUFFLE
         J    = DMOD(DSEED,128.D0)+1.D0
         X    = DSEED / D2P31
         R(I) = WK(J)
         WK(J)= X
   15 CONTINUE
      RETURN
      END

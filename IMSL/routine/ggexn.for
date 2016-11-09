C   IMSL ROUTINE NAME   - GGEXN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - EXPONENTIAL RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGEXN (DSEED,XM,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                XM     - INPUT MEAN VALUE.
C                NR     - INPUT NUMBER OF DEVIATES TO BE GENERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           EXPONENTIAL DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGEXN   (DSEED,XM,NR,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               XM,R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
C                                  GET THE RANDOM NUMBERS
C                                  FIRST EXECUTABLE STATEMENT
      CALL GGUBS(DSEED,NR,R)
C                                  TRANSFORM TO EXPONENTIAL DEVIATES
      DO 5 I=1,NR
    5 R(I) = -XM*ALOG(R(I))
      RETURN
      END

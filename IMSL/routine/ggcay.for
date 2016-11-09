C   IMSL ROUTINE NAME   - GGCAY
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - CAUCHY RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGCAY (DSEED,NR,WK,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT. NUMBER OF DEVIATES TO BE GENERATED.
C                WK     - WORK AREA VECTOR OF LENGTH 3*NR.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           CAUCHY DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGCAY  (DSEED,NR,WK,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               WK(1),R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K,M,I
      REAL               PT5,TWO,Y1,Y2,Y
      DATA               PT5/0.5/,TWO/2.0/
C                                  FIRST EXECUTABLE STATEMENT
      K = -1
      M = 3*NR
      I = 1
C                                  OBTAIN M UNIFORM DEVIATES
    5 K = K+2
      IF (K .GE. M) K = 1
      IF (K .EQ. 1) CALL GGUBS(DSEED,M,WK)
      Y1 = TWO*(WK(K)-PT5)
      Y2 = WK(K+1)
      Y = Y1*Y1+Y2*Y2
C                                  REJECTION TEST
      IF (Y .GT. 1.) GO TO 5
C                                  TAKE RATIO
      R(I) = Y1/Y2
      I = I+1
      IF (I .LE. NR) GO TO 5
      RETURN
      END

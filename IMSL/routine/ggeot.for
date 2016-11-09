C   IMSL ROUTINE NAME   - GGEOT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GEOMETRIC RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGEOT (DSEED,NR,P,WK,IR)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF DEVIATES TO BE GENERATED.
C                P      - INPUT. PROBABILITY OF GETTING A SUCCESS ON
C                           ANY TRIAL.
C                WK     - VECTOR OF LENGTH NR USED AS WORK STORAGE.
C                IR     - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           GEOMETRIC DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGEOT  (DSEED,NR,P,WK,IR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,IR(NR)
      REAL               P,WK(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K,I,ICNT
C                                  FIRST EXECUTABLE STATEMENT
      K = 0
      I = 1
    5 ICNT = 1
   10 K = K+1
      IF (K .GT. NR) K = 1
      IF (K .EQ. 1) CALL GGUBS(DSEED,NR,WK)
      IF (WK(K) .LE. P) GO TO 15
      ICNT = ICNT+1
      GO TO 10
   15 IR(I) = ICNT
      I = I+1
      IF (I .LE. NR) GO TO 5
      RETURN
      END

C   IMSL ROUTINE NAME   - GGTRA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - TRIANGULAR DISTRIBUTION RANDOM DEVIATE
C                           GENERATOR
C
C   USAGE               - CALL GGTRA (DSEED,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF DEVIATES TO BE GENERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           TRIANGULAR DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGTRA  (DSEED,NR,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
C                                  GENERATE NR (0,1) UNIFORM DEVIATES
C                                  FIRST EXECUTABLE STATEMENT
      CALL GGUBS(DSEED,NR,R)
C                                  TRANSFORM UNIFORM DEVIATES TO
C                                  TRIANGULAR DEVIATES
      DO 10 I = 1,NR
         IF (R(I) .LE. .5) GO TO 5
         R(I) = 1.-SQRT(.5*(1.-R(I)))
         GO TO 10
    5    R(I) = SQRT(.5*R(I))
   10 CONTINUE
      RETURN
      END

C   IMSL ROUTINE NAME   - GGUD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - DISCRETE UNIFORM RANDOM NUMBER GENERATOR
C
C   USAGE               - CALL GGUD   (DSEED,K,NR,IR)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0,2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                K      - INPUT PARAMETER OF DISCRETE UNIFORM
C                          DISTRIBUTION.  THE INTEGERS 1,2,...,K
C                          OCCUR WITH EQUAL PROBABILITY.  K MUST
C                          BE POSITIVE.
C                NR     - INPUT NUMBER OF RANDOM NUMBERS TO BE
C                           GENERATED.
C                IR     - OUTPUT VECTOR OF LENGTH NR CONTAINING
C                           THE UNIFORMLY DISTRIBUTED INTEGERS.
C
C   REQD. IMSL ROUTINES - GGUBFS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGUD   (DSEED,K,NR,IR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,NR,IR(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               FK,ONE
      DATA               ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      FK = FLOAT(K)
      DO 5 I=1,NR
         IR(I) = FK*GGUBFS(DSEED)+ONE
    5 CONTINUE
      RETURN
      END

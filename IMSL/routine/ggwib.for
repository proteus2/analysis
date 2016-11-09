C   IMSL ROUTINE NAME   - GGWIB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - WEIBULL RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGWIB (DSEED,A,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                A      - INPUT. SHAPE PARAMETER FOR THE DESIRED WEIBULL
C                           FUNCTION. A MUST BE GREATER THAN 0. NO
C                           PROGRAM CHECK IS MADE FOR THIS CONDITION.
C                NR     - INPUT NUMBER OF DEVIATES TO BE GENERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           WEIBULL DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGWIB  (DSEED,A,NR,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               A,R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               A1
C                                  FIRST EXECUTABLE STATEMENT
      A1 = 1./A
C                                  GENERATE NR (0,1) UNIFORM DEVIATES
      CALL GGUBS(DSEED,NR,R)
C                                  TRANSFORM UNIFORM DEVIATES TO WEIBULL
      DO 5 I = 1,NR
         R(I) = (-ALOG(R(I)))**A1
    5 CONTINUE
      RETURN
      END

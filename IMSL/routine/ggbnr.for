C   IMSL ROUTINE NAME   - GGBNR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NEGATIVE BINOMIAL RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGBNR (DSEED,K,P,NR,WK,IR)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                K      - NEGATIVE BINOMIAL PARAMETER (NUMBER OF
C                           SUCCESSES) AND MUST BE A NON-NEGATIVE
C                           INTEGER. (INPUT)
C                P      - NEGATIVE BINOMIAL PARAMETER (PROBABILITY OF
C                           A SUCCESS ON EACH TRIAL) AND MUST BE IN
C                           THE INCLUSIVE (0,1) RANGE. (INPUT)
C                NR     - NUMBER OF NEGATIVE BINOMIAL DEVIATES TO BE
C                           GENERATED. (INPUT)
C                WK     - WORK VECTOR OF LENGTH NR
C                IR     - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           NEGATIVE BINOMIAL DEVIATES
C
C   REQD. IMSL ROUTINES - GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGBNR  (DSEED,K,P,NR,WK,IR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,NR,IR(NR)
      REAL               P,WK(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,KK
      REAL               P1,R1,R,U,SUM,XI,XD,EPS
      DATA               EPS/1.E-6/
C                                  FIRST EXECUTABLE STATEMENT
      P1 = 1.-P
      R1 = P**K
      CALL GGUBS(DSEED,NR,WK)
      DO 15 I = 1,NR
         IR(I) = K
         R = R1
         U = WK(I)
         IF (U .LE. R) GO TO 15
         SUM = R
         XI = K
         XD = 1.
         KK = K
    5    KK = KK+1
         R = XI*P1*R/XD
         IF (R .LT. EPS*SUM) GO TO 10
         SUM = SUM+R
         XI = XI+1.
         XD = XD+1.
         IF (U .GT. SUM) GO TO 5
         IR(I) = KK
         GO TO 15
   10    CALL GGUBS(DSEED,1,WK(I))
         IR(I) = KK+1+9.*WK(I)
   15 CONTINUE
      RETURN
      END

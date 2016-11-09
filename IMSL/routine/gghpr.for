C   IMSL ROUTINE NAME   - GGHPR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - HYPERGEOMETRIC RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGHPR (DSEED,N,L,M,NR,WK,IR)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                N      - NUMBER OF ITEMS IN THE SAMPLE. (INPUT)
C                L      - NUMBER OF ITEMS IN THE LOT. (INPUT)
C                M      - NUMBER OF DEFECTIVE ITEMS IN THE LOT. (INPUT)
C                NR     - NUMBER OF RANDOM NUMBERS DESIRED.(INPUT)
C                WK     - VECTOR OF LENGTH NR USED AS WORK STORAGE.
C                IR     - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           HYPERGEOMETRIC DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGHPR  (DSEED,N,L,M,NR,WK,IR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,L,M,NR,IR(NR)
      REAL               WK(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K,I,MM,ICNT,IBAD,KK
      REAL               TOT1,TOT
C                                  FIRST EXECUTABLE STATEMENT
      K = 0
      TOT1 = L
      I = 1
    5 MM = M
      TOT = TOT1
      ICNT = 1
      IBAD = 0
   10 K = K+1
      IF (K .GT. NR) K = 1
      IF (K .EQ. 1) CALL GGUBS(DSEED,NR,WK)
      KK = TOT*WK(K)+1.
      IF (KK .GT. MM) GO TO 15
      IBAD = IBAD+1
      MM = MM-1
   15 IF(IBAD .EQ. M .OR. ICNT .EQ. N) GO TO 20
      TOT = TOT-1.
      ICNT = ICNT+1
      GO TO 10
   20 IR(I) = IBAD
      I = I+1
      IF (I .LE. NR) GO TO 5
      RETURN
      END

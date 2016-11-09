C   IMSL ROUTINE NAME   - MDGAM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - GAMMA PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDGAM (X,P,PROB,IER)
C
C   ARGUMENTS    X      - INPUT VALUE TO WHICH GAMMA IS TO BE INTEGRATED
C                P      - INPUT GAMMA PARAMETER
C                PROB   - OUTPUT PROBABILITY THAT A RANDOM VARIABLE
C                           HAVING A GAMMA DISTRIBUTION WITH PARAMETER P
C                           WILL BE LESS THAN OR EQUAL TO X.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES  X IS LESS THAN ZERO.
C                           IER = 130 INDICATES  P IS LESS THAN OR EQUAL
C                             TO ZERO.
C
C   REQD. IMSL ROUTINES - MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDGAM   (X,P,PROB,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               X,P,PROB
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               V(6),V1(6),PNLG,CNT,YCNT,XMIN,AX,BIG,CUT,
     *                   Y,Z,RATIO,REDUC
      EQUIVALENCE        (V(3),V1(1))
      DATA               XMIN/-174.673/
C                                  TEST  X AND  P
C                                  FIRST EXECUTABLE STATEMENT
      PROB = 0.0
      IF (X  .GE. 0.0) GO TO 5
      IER = 129
      GO TO 9000
    5 IF (P .GT. 0.0) GO TO 10
      IER = 130
      GO TO 9000
   10 IER = 0
      IF (X  .EQ. 0.0) GO TO 9005
C                                  DEFINE LOG-GAMMA AND INITIALIZE
      PNLG = ALGAMA(P)
      CNT = P * ALOG(X)
      YCNT = X + PNLG
      IF ((CNT-YCNT) .GT. XMIN) GO TO 15
      AX = 0.0
      GO TO 20
   15 AX = EXP(CNT-YCNT)
   20 BIG = 1.E35
      CUT = 1.E-8
C                                  CHOOSE ALGORITHMIC METHOD
      IF ((X  .LE. 1.0) .OR. (X  .LT. P )) GO TO 40
C                                  CONTINUED FRACTION EXPANSION
      Y = 1.0 - P
      Z = X  + Y + 1.0
      CNT = 0.0
      V(1) = 1.0
      V(2) = X
      V(3) = X + 1.0
      V(4) = Z * X
      PROB = V(3)/V(4)
   25 CNT = CNT + 1.0
      Y = Y + 1.0
      Z = Z + 2.0
      YCNT = Y * CNT
      V(5) = V1(1) * Z - V(1) * YCNT
      V(6) = V1(2) * Z - V(2) * YCNT
      IF (V(6) .EQ. 0.0) GO TO 50
      RATIO = V(5)/V(6)
      REDUC = ABS(PROB-RATIO)
      IF (REDUC .GT. CUT) GO TO 30
      IF (REDUC .LE. RATIO*CUT) GO TO 35
   30 PROB = RATIO
      GO TO 50
   35 PROB = 1.0 - PROB * AX
      GO TO 9005
C                                  SERIES EXPANSION
   40 RATIO =  P
      CNT = 1.0
      PROB = 1.0
   45 RATIO = RATIO + 1.0
      CNT = CNT *  X/RATIO
      PROB = PROB + CNT
      IF (CNT .GT. CUT) GO TO 45
      PROB = PROB * AX/P
      GO TO 9005
   50 DO 55 I=1,4
         V(I) = V1(I)
   55 CONTINUE
      IF (ABS(V(5)) .LT. BIG) GO TO 25
C                                  SCALE TERMS DOWN TO PREVENT OVERFLOW
      DO 60 I=1,4
         V(I) = V(I)/BIG
   60 CONTINUE
      GO TO 25
 9000 CONTINUE
      CALL UERTST (IER,6HMDGAM )
 9005 RETURN
      END

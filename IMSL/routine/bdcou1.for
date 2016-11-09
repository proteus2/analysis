C   IMSL ROUTINE NAME   - BDCOU1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TALLY OF OBSERVATIONS INTO A ONE-WAY
C                           FREQUENCY TABLE
C
C   USAGE               - CALL BDCOU1 (X,N,K,DIV,BU,BL,TAB,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF DATA TO BE TALLIED INTO
C                           CATEGORIES. X HAS LENGTH N.
C                N      - NUMBER OF OBSERVATIONS. LENGTH OF X. (INPUT)
C                K      - INPUT. IF K IS LESS THAN ZERO, EQUIDISTANT
C                           CATEGORIES ARE USED. SEE PARAMETERS BU,BL.
C                           OTHERWISE, K IS THE LENGTH OF THE DIV VECTOR
C                DIV    - INPUT VECTOR OF LENGTH K. IF K IS LESS THAN
C                           ZERO, DIV MAY BE DECLARED WITH LENGTH ONE,
C                           AND IS UNUSED. DIV SHOULD BE ORDERED
C                           MONOTONIC INCREASING ON ENTRY, AND IF DIV(I)
C                           IS LESS THAN X(J) IS LESS THAN OR EQUAL TO
C                           DIV(I+1), A COUNT IS ADDED TO  TAB(I+1).
C                BU     - INPUT UPPER BOUND OF THE EQUIDISTANT
C                           CATEGORIZATION. IF K.GT.0, BU IS THE EPSILON
C                           USED TO DETERMINE IF X(J) IS EQUAL TO SOME
C                           DIV(I). SEE DOCUMENT.
C                BL     - INPUT LOWER BOUND OF THE EQUIDISTANT
C                           CATEGORIZATION. IF K.GT.0, BL IS USED AS THE
C                           OPTION SWITCH TO ALLOW MATCHES, WITHIN AN
C                           EPSILON DIFFERENCE. BL.GT.0 IMPLIES THAT BU
C                           WILL BE USED AS THE EPSILON.  BL.LE.0
C                           IMPLIES THAT MATCHING IS NOT DESIRED.
C                TAB    - OUTPUT TALLY VECTOR OF LENGTH K+1 (IF DIV IS
C                           USED), AND OF LENGTH ABS(K)+2 (IF BU,BL ARE
C                           USED).
C                         ON INPUT, TAB MUST BE A NULL VECTOR.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS LENGTH OF X IS LESS THAN OR
C                             EQUAL TO 1.
C                           IER=130 MEANS K EQUAL TO ZERO.
C                           IER=131 MEANS K IS LESS THAN 0 AND BU IS
C                             LESS THAN BL.
C                         WARNING ERROR
C                           IER=36 MEANS THE X VECTOR IS CONSTANT WHEN K
C                             IS LESS THAN ZERO AND BU EQUAL TO BL. IN
C                             THIS CASE  TAB(ABS(K)+2)=THE LENGTH OF X.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDCOU1 (X,N,K,DIV,BU,BL,TAB,IER)
C
      REAL               X(1),DIV(1),TAB(1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER=0
C                                  CHECK ERRORS
      IF (N.GT.1) GO TO 5
      IER=129
      GO TO 9000
    5 IF (K.NE.0) GO TO 10
      IER=130
      GO TO 9000
   10 IF (K.LT.0.AND.BU.LT.BL) GO TO 15
      GO TO 20
   15 IER=131
      GO TO 9000
   20 IABSK=IABS(K)
   25 IF (K.LT.0) GO TO 50
      IF (BL.GT.0.0) GO TO 40
C                                  CASE WHERE K.GT.0 AND BL.LE.0
      KP1=K+1
      DO 35 I=1,K
         DO 35 J=1,N
            IF (I.NE.1) GO TO 30
            IF (X(J).LE.DIV(1))  TAB(1)= TAB(1)+1.0
            IF (X(J).GT.DIV(K))  TAB(KP1)= TAB(KP1)+1.0
            GO TO 35
   30       IF (X(J).GT.DIV(I-1).AND.X(J).LE.DIV(I))  TAB(I)= TAB(I)+1.0
   35 CONTINUE
      GO TO 9005
C                                  CASE WHERE K.GT.0 AND BL.GT.0
   40 DO 45 I=1,K
         DO 45 J=1,N
            IF (ABS(X(J)-DIV(I)).GT.BU) GO TO 45
            TAB(I) = TAB(I)+1.0
   45 CONTINUE
      GO TO 9005
   50 IF (BU.EQ.BL) GO TO 65
C                                  CASE WHERE K.LT.0 AND BU.GT.BL
      BUL=BU-BL
      DO 60 J=1,N
         IF (X(J).LT.BL.OR.X(J).GE.BU) GO TO 55
         I=((X(J)-BL)*IABSK)/BUL+2.0
         TAB(I)=TAB(I)+1.0
         GO TO 60
   55    IF (X(J).LT.BL)  TAB(1)= TAB(1)+1.0
         IF (X(J).GE.BU)  TAB(IABSK+2)= TAB(IABSK+2)+1.0
   60 CONTINUE
      GO TO 9005
C                                  CASE WHERE K.LT.0 AND BU=BL, FIND
C                                  THE MINIMUM AND MAXIMUM VALUES OF
C                                  X VECTOR
   65 XMIN=X(1)
      XMAX=X(1)
      DO 70 J=2,N
         IF (X(J).LT.XMIN) XMIN=X(J)
         IF (X(J).GT.XMAX) XMAX=X(J)
   70 CONTINUE
      BL=XMIN
      BU=XMAX
      IF (BL.NE.BU) GO TO 75
      TAB (IABSK+2)=N
      IER=36
      GO TO 9000
   75 BUL=BU-BL
      DO 85 J=1,N
         IF (X(J).GE.BU) GO TO 80
         I=((X(J)-BL)*IABSK)/BUL+2.0
         TAB(I)=TAB(I)+1.0
         GO TO 85
   80    TAB( IABSK+2)= TAB(IABSK+2)+1.0
   85 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BDCOU1')
 9005 RETURN
      END

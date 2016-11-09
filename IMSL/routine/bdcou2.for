C   IMSL ROUTINE NAME   - BDCOU2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TALLY OF OBSERVATIONS INTO A TWO-WAY
C                           FREQUENCY TABLE
C
C   USAGE               - CALL BDCOU2 (X,Y,N,K1,K2,DIVX,DIVY,XU,XL,
C                           YU,YL,IT,TAB,IER)
C
C   ARGUMENTS    X      - 1ST INPUT VECTOR OF DATA TO BE TALLIED.
C                           X HAS LENGTH N.
C                Y      - 2ND INPUT VECTOR OF DATA TO BE TALLIED.
C                           Y HAS LENGTH N.
C                N      - INPUT.  NUMBER OF OBSERVATIONS TO TALLY.
C                K1     - INPUT.  IF K1 IS LESS THAN ZERO, ABS(K1)
C                           EQUIDISTANT CATEGORIES ARE USED.  SEE
C                           PARAMETERS XU AND XL.  OTHERWISE, K1 IS THE
C                           LENGTH OF THE DIVX VECTOR.
C                K2     - INPUT.  IF K2 IS LESS THAN ZERO, ABS(K2)
C                           EQUIDISTANT CATEGORIES ARE USED.  SEE
C                           PARAMETERS YU AND YL.  OTHERWISE, K2 IS THE
C                           LENGTH OF THE DIVY VECTOR.
C                DIVX   - INPUT VECTOR OF LENGTH K1.  IF K1 IS LESS THAN
C                           ZERO, DIVX MAY BE DECLARED OF LENGTH ONE,
C                           AND IS UNUSED.  DIVX SHOULD BE ORDERED
C                           MONOTONIC INCREASING ON ENTRY, AND IF
C                           DIVX(I) IS LESS THAN X(J) IS LESS THAN OR
C                           EQUAL TO DIVX(I+1), A COUNT IS ADDED TO
C                           THE I+1ST ROW OF TAB.  THE COLUMN SUBSCRIPT
C                           OF TAB IS CHOSEN BY Y(J).
C                DIVY   - INPUT VECTOR OF LENGTH K2.  IF K2 IS LESS THAN
C                           ZERO, DIVY MAY BE DECLARED OF LENGTH ONE,
C                           AND IS UNUSED.  DIVY SHOULD BE ORDERED
C                           MONOTONIC INCREASING ON ENTRY, AND IF
C                           DIVY(I) IS LESS THAN Y(K) IS LESS THAN OR
C                           EQUAL TO DIVY(I+1), A COUNT IS ADDED TO
C                           TAB(I+1,K+1) WHERE I IS CHOSEN BY X(J).
C                XU     - INPUT UPPER BOUND OF THE EQUIDISTANT
C                           CATEGORIZATION X.  IF K1 IS GREATER THAN
C                           ZERO, XU IS THE EPSILON USED TO DETERMINE
C                           THAT X(J) IS EQUAL TO SOME DIVX(I).
C                XL     - INPUT LOWER BOUND OF THE EQUIDISTANT
C                           CATEGORIZATION X.  IF K1 IS GREATER THAN
C                           ZERO, XL IS USED AS THE SWITCH TO ALLOW
C                           MATCHES WITHIN AN EPSILON DIFFERENCE.
C                           XL GREATER THAN ZERO IMPLIES THAT XU WILL
C                           BE USED AS THE EPSILON.
C                YU     - INPUT UPPER BOUND OF THE EQUIDISTANT
C                           CATEGORIZATION Y.  IF K2 IS GREATER THAN
C                           ZERO, YU IS THE EPSILON USED TO DETERMINE
C                           THAT Y(J) IS EQUAL TO SOME DIVY(I).
C                YL     - INPUT LOWER BOUND OF THE EQUIDISTANT
C                           CATEGORIZATION Y.  IF K2 IS GREATER THAN
C                           ZERO, YL IS USED AS THE SWITCH TO ALLOW
C                           MATCHES WITHIN AN EPSILON DIFFERENCE.
C                           YL GREATER THAN ZERO IMPLIES THAT YU WILL
C                           BE USED AS THE EPSILON.
C                IT     - ROW DIMENSION OF MATRIX TAB EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                TAB    - OUTPUT TALLY MATRIX.  TAB MUST BE AT LEAST
C                           (K1+1,K2+1) IF DIVX AND DIVY ARE USED, AND
C                           AT LEAST (ABS(K1)+2,ABS(K2)+2) IF XU, XL,
C                           YU, AND YL ARE USED.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS XU IS LESS THAN XL WHEN K1 IS
C                             LESS THAN ZERO.
C                           IER=130 MEANS YU IS LESS THAN YL WHEN K2 IS
C                             LESS THAN ZERO.
C                           IER=131 MEANS N IS LESS THAN ONE.
C                           IER=132 MEANS K1 OR K2 IS EQUAL TO ZERO.
C                         WARNING ERROR
C                           IER=37 MEANS THE X VECTOR IS CONSTANT WHEN
C                             K1 IS LESS THAN ZERO AND XU=XL.
C                           IER=38 MEANS THE Y VECTOR IS CONSTANT WHEN
C                             K2 IS LESS THAN ZERO AND YU=YL.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDCOU2 (X,Y,N,K1,K2,DIVX,DIVY,XU,XL,YU,YL,IT,TAB,IER)
C
      REAL               X(1),Y(1),DIVX(1),DIVY(1),TAB(IT,1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  TEST FOR INPUT ERROR
      IF(N .LT. 1) GO TO 105
      IF(K1.EQ.0 .OR. K2.EQ.0) GO TO 110
      IER = 0
      IF(K1.GT.0 .AND. K2.GT.0) GO TO 30
      IF(K2 .GT. 0) GO TO 15
      IF(YU .LT. YL) GO TO 100
      IF(YU .GT. YL) GO TO 10
      YL = Y(1)
      YU = Y(1)
C                                  COMPUTE Y MAX AND MIN ONLY
      DO 5 I=2,N
         IF(Y(I) .LT. YL) YL = Y(I)
         IF(Y(I) .GT. YU) YU = Y(I)
    5 CONTINUE
   10 YUMYL = YU-YL
      IABK2 = -K2
   15 IF(K1 .GT. 0) GO TO 30
      IF(XU .LT. XL) GO TO 95
      IF(XU .GT. XL) GO TO 25
C                                  COMPUTE X MAX AND MIN ONLY
      XL = X(1)
      XU = X(1)
      DO 20 I=2,N
         IF(X(I) .LT. XL) XL = X(I)
         IF(X(I) .GT. XU) XU = X(I)
   20 CONTINUE
   25 XUMXL = XU-XL
      IABK1 = -K1
C                                  DO THE TALLYING
   30 DO 85 I=1,N
         J2 = 1
C                                  COMPUTE SECOND SUBSCRIPT
         IF(K2 .GT. 0) GO TO 35
         J2 = IABK2+2
         IF(YUMYL .EQ. 0.0) GO TO 55
         J2 = ((Y(I)-YL)*IABK2)/YUMYL+2
         IF(J2 .LT. 1) J2 = 1
         IF(J2 .GT. IABK2+2) J2 = IABK2+2
         GO TO 55
   35    IF(YL .LE. 0.0) GO TO 45
C                                  YL IS GREATER THAN ZERO
         DO 40 I1=1,K2
            IF(ABS(Y(I)-DIVY(I1)) .LE. YU) GO TO 55
   40    J2 = J2+1
         GO TO 85
C                                  YL IS LESS THAN OR EQUAL TO ZERO
   45    IF(Y(I) .LE. DIVY(1)) GO TO 55
         DO 50 I1=2,K2
            J2 = J2+1
            IF(Y(I).GT.DIVY(I1-1) .AND. Y(I).LE.DIVY(I1)) GO TO 55
   50    CONTINUE
         J2 = K2+1
   55    J1 = 1
C                                  COMPUTE FIRST SUBSCRIPT
         IF(K1 .GT. 0) GO TO 60
         J1 = IABK1+2
         IF(XUMXL .EQ. 0.0) GO TO 80
         J1 = ((X(I)-XL)*IABK1)/XUMXL+2
         IF(J1 .LT. 1) J1 = 1
         IF(J1 .GT. IABK1+2) J1 = IABK1+2
         GO TO 80
   60    IF(XL .LE. 0.0) GO TO 70
C                                  XL IS GREATER THAN ZERO
         DO 65 I1=1,K1
            IF(ABS(X(I)-DIVX(I1)) .LE. XU) GO TO 80
   65    J1 = J1+1
         GO TO 85
C                                  XL IS LESS THAN OR EQUAL TO ZERO
   70    IF(X(I) .LE. DIVX(1)) GO TO 80
         DO 75 I1=2,K1
            J1 = J1+1
            IF(X(I).GT.DIVX(I1-1) .AND. X(I).LE.DIVX(I1)) GO TO 80
   75    CONTINUE
         J1 = K1+1
C                                  INCREMENT THE TALLY MATRIX
   80    TAB(J1,J2) = TAB(J1,J2)+1.0
   85 CONTINUE
      IF(K1.LT.0 .AND. XU.EQ.XL) GO TO 90
      IF(K2.LT.0 .AND. YU.EQ.YL) GO TO 115
      GO TO 9005
   90 IER = 37
      GO TO 9000
   95 IER = 129
      GO TO 9000
  100 IER = 130
      GO TO 9000
  105 IER = 131
      GO TO 9000
  110 IER = 132
      GO TO 9000
  115 IER = 38
 9000 CONTINUE
      CALL UERTST(IER,'BDCOU2')
 9005 RETURN
      END

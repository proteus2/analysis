C   IMSL ROUTINE NAME   - GGAMR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - ONE PARAMETER GAMMA RANDOM DEVIATE GENERATOR,
C                           AND USABLE AS THE BASIS FOR TWO PARAMETER
C                           GAMMA, EXPONENTIAL, CHI-SQUARED, CHI, BETA,
C                           T, AND F DEVIATE GENERATION
C
C   USAGE               - CALL GGAMR (DSEED,A,NR,WK,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT.  DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                A      - INPUT. SHAPE PARAMETER FOR THE DESIRED GAMMA
C                           FUNCTION. A MUST BE GREATER THAN ZERO. NO
C                           PROGRAM CHECK IS MADE FOR THIS CONDITION.
C                NR     - INPUT.  NUMBER OF DEVIATES TO BE GENERATED.
C                WK     - WORK VECTOR OF LENGTH 2*NR IF A .LT. 1.0.
C                           IF A .GE. 1.0, WK IS NOT USED AND MAY BE
C                           DIMENSIONED AS WK(1).
C                R      - OUTPUT. VECTOR OF LENGTH NR CONTAINING THE
C                           GAMMA DEVIATES.
C
C   REQD. IMSL ROUTINES - GGNML,GGUBS,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGAMR (DSEED,A,NR,WK,R)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               A,WK(1),R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,K,M
      REAL               A1,ALETA,D,E,EPS,EX,F1,F2,F4,F5,ONE,P,P1,P10,
     *                   P2,P3,P4,P5,P6,P7,P8,P9,TEMP,TEMP1,TEMP2,U,
     *                   U1(1),U2,V,W(1),W1,W2,X,X1,X2,X3,X4,X5,XETA,
     *                   XLL,XLR,XM,ZERO
      DATA               XETA /Z00100000/
      DATA               EPS/Z3C100000/
      DATA               E/2.718282/
      DATA               ONE /1.0/,XLL /1.0/,XLR,XM,A1,X1,X2,X3,X4,X5,
     *                   F1,F2,F4,F5,ZERO,P1,P2,P3,P4,P5,P6,P7,P8,P9,
     *                   P10 /23*0.0/
C                                  FIRST EXECUTABLE STATEMENT
      I = 1
      M = NR+NR
      XM = A-1.
      IF (A.GT.1.) GO TO 50
C                                  AD HOC PROCEDURES
      K = A
      IF (A.EQ.1.) GO TO 10
      IF (A.NE..5) GO TO 20
C                                  FOR A=.5 USE
C                                  NORMAL**2 DEVIATE DIVIDED BY 2
      CALL GGNML(DSEED,NR,R)
      DO 5 II=1,NR
         R(II) = R(II)*R(II)*.5
    5 CONTINUE
      GO TO 120
C                                  FOR A EQUAL TO ONE, GENERATE
C                                  EXPONENTIAL DEVIATES
   10 CALL GGUBS(DSEED,NR,R)
      DO 15 II=1,NR
         R(II) = -ALOG(R(II))
   15 CONTINUE
      GO TO 120
C                                  A LESS THAN ONE
   20 TEMP1 = (E+A)/E
      TEMP2 = 1./A
      ALETA = ALOG(XETA)
      K = -1
   25 K = K+2
      IF (K.GE.M) K = 1
C                                  GENERATE A UNIFORM (0,1) DEVIATE AND
C                                  A UNIT EXPONENTIAL PSEUDO-RANDOM
C                                  DEVIATE
      IF (K.EQ.1) CALL GGUBS(DSEED,M,WK)
      EX = -ALOG(WK(K))
C                                  REJECTION TEST
      P = TEMP1*WK(K+1)
      IF (P.GT.1.) GO TO 40
C                                  SMALL X    CHECK FOR UNDERFLOW
      IF (TEMP2*ALOG(P).GE.ALETA) GO TO 30
      TEMP = XETA
      GO TO 35
   30 TEMP = P**TEMP2
   35 X = TEMP
      GO TO 45
C                                  LARGE X
   40 X = -ALOG((TEMP1-P)*TEMP2)
      TEMP = -XM*ALOG(X)
   45 IF (TEMP.GT.EX) GO TO 25
      R(I) = X
      I = I+1
      IF (I.LE.NR) GO TO 25
      GO TO 120
C                                  A GREATER THAN ONE
   50 CONTINUE
      IF (A.EQ.A1) GO TO 60
C                                  INITIALIZATION
      X3 = XM
      D = SQRT(X3)
      X2 = ZERO
      X1 = ZERO
      F1 = ZERO
      F2 = ZERO
      XLL = ONE
      IF (D*(ONE+EPS).GE.X3) GO TO 55
      X2 = X3-D
      X1 = X2*(ONE-ONE/D)
      XLL = ONE-X3/X1
      F1 = EXP(X3*ALOG(X1/X3)+X3-X1)
      F2 = EXP(X3*ALOG(X2/X3)+X3-X2)
   55 X4 = X3+D
      X5 = X4*(ONE+ONE/D)
      XLR = ONE-X3/X5
      F4 = EXP(X3*ALOG(X4/X3)+X3-X4)
      F5 = EXP(X3*ALOG(X5/X3)+X3-X5)
C                                  CALCULATE PROBABILITY FOR
C                                  EACH OF THE TEN REGIONS
      P1 = F2*(X3-X2)
      P2 = F4*(X4-X3)+P1
      P3 = F1*(X2-X1)+P2
      P4 = F5*(X5-X4)+P3
      P5 = (ONE-F2)*(X3-X2)+P4
      P6 = (ONE-F4)*(X4-X3)+P5
      P7 = (F2-F1)*(X2-X1)*.5+P6
      P8 = (F4-F5)*(X5-X4)*.5+P7
      P9 = -F1/XLL+P8
      P10 = F5/XLR+P9
      A1 = A
C                                  GENERATE ONE UNIFORM(0,1) DEVIATE
   60 CALL GGUBS(DSEED,1,U1)
      U = U1(1)*P10
C                                  THE FOUR REGIONS WITH ZERO
C                                  PROBABILITY OF REJECTION
      IF (U.GT.P4) GO TO 85
      IF (U.GT.P1) GO TO 70
      X = X2+U/F2
   65 R(I) = X
      I = I+1
      IF (I.LE.NR) GO TO 60
      GO TO 120
   70 IF (U.GT.P2) GO TO 75
      X = X3+(U-P1)/F4
      GO TO 65
   75 IF (U.GT.P3) GO TO 80
      X = X1+(U-P2)/F1
      GO TO 65
   80 X = X4+(U-P3)/F5
      GO TO 65
C                                  THE TWO REGIONS USING
C                                  RECTANGULAR REJECTION
   85 CALL GGUBS(DSEED,1,W)
      W1 = W(1)
      IF (U.GT.P5) GO TO 90
      X = X2+(X3-X2)*W1
      IF ((U-P4)/(P5-P4).LE.W1) GO TO 65
      V = F2+(U-P4)/(X3-X2)
      GO TO 115
   90 IF (U.GT.P6) GO TO 95
      X = X3+(X4-X3)*W1
      IF ((P6-U)/(P6-P5).GE.W1) GO TO 65
      V = F4+(U-P5)/(X4-X3)
      GO TO 115
C                                  THE TWO TRIANGULAR REGIONS
   95 IF (U.GT.P8) GO TO 105
      CALL GGUBS(DSEED,1,W)
      W2 = W(1)
      IF (W2.GT.W1) W1 = W2
      IF (U.GT.P7) GO TO 100
      X = X1+(X2-X1)*W1
      V = F1+2.*W1*(U-P6)/(X2-X1)
      IF (V.LE.F2*W1) GO TO 65
      GO TO 115
  100 X = X5-W1*(X5-X4)
      V = F5+2.*W1*(U-P7)/(X5-X4)
      IF (V.LE.F4*W1) GO TO 65
      GO TO 115
C                                  THE TWO EXPONENTIAL REGIONS
  105 IF (U.GE.P9) GO TO 110
      U = (P9-U)/(P9-P8)
      X = X1-ALOG(U)/XLL
      IF (X.LE.0.) GO TO 60
      IF (W1.LT.(XLL*(X1-X)+ONE)/U) GO TO 65
      V = W1*F1*U
      GO TO 115
  110 IF (P10.EQ.U) GO TO 55
      U = (P10-U)/(P10-P9)
      X = X5-ALOG(U)/XLR
      IF (W1.LT.(XLR*(X5-X)+ONE)/U) GO TO 65
      V = W1*F5*U
C                                  PERFORM THE STANDARD REJECTION
  115 IF (ALOG(V).GT.X3*ALOG(X/X3)+X3-X) GO TO 60
      GO TO 65
  120 RETURN
      END

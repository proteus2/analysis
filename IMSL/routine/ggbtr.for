C   IMSL ROUTINE NAME   - GGBTR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - BETA RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGBTR (DSEED,P,Q,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                P      - INPUT FIRST BETA PARAMETER. P MUST BE GREATER
C                           THAN ZERO.
C                Q      - INPUT SECOND BETA PARAMETER. Q MUST BE GREATER
C                           THAN ZERO.
C                NR     - INPUT NUMBER OF RANDOM NUMBERS TO BE GEN-
C                           ERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           BETA(P,Q) DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBFS,GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGBTR  (DSEED,P,Q,NR,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               P,Q,R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               UA(2),P1,Q1,X,F1,F2,T,S1,S2,GAMA,RLAM,TEMPA,
     1                   TEMPAA,TEMPB,TEMPBB,TEMPAB,TEMPBA,TEMP
      LOGICAL            PMIN
      DATA               EPS/Z3C100000/
      DATA               CONS1/1.386294/
      DATA               CONS2/2.609438/
      DATA               ONE/1.0/,TWO/2.0/
      DATA               ZERO/0.0/,FIVE/5.0/
C                                  FIRST EXECUTABLE STATEMENT
      I = 0
      IDX = 0
      IF (P.EQ.ONE.OR.Q.EQ.ONE) GO TO 20
      IF (P.LT.ONE) GO TO 5
      IF (Q.GT.ONE) GO TO 45
      GO TO 160
    5 IF (Q.GT.ONE) GO TO 150
C                                  JOHNKS METHOD - BOTH P AND Q LESS
C                                    THAN 1
      P1 = 1.0/P
      Q1 = 1.0/Q
   10 I = I+1
C                                  OBTAIN UNIFORM DEVIATES
   15 CALL GGUBS (DSEED,2,UA)
      X = UA(1)**P1
      Y = UA(2)**Q1
      Y = X+Y
C                                  REJECT OR ACCEPT
      IF (Y.GT.ONE) GO TO 15
      R(I) = X/Y
      IF (I.LT.NR) GO TO 10
      GO TO 9005
C                                  USE THE INVERSE OF THE CUMULATIVE
C                                    DISTRIBUTION FUNCTION FOR P OR Q
C                                    EQUAL TO ONE
   20 CALL GGUBS (DSEED,NR,R)
      IF (P.NE.ONE) GO TO 25
      IF (Q.NE.ONE) GO TO 35
      GO TO 9005
C                                  Q IS ONE
   25 P1 = ONE/P
      DO 30 I=1,NR
         R(I) = R(I)**P1
   30 CONTINUE
      GO TO 9005
C                                  P IS ONE
   35 Q1 = ONE/Q
      DO 40 I=1,NR
         R(I) = ONE-R(I)**Q1
   40 CONTINUE
      GO TO 9005
   45 IF (NR.LT.4) GO TO 120
C                                  SCHMEISERS METHOD - FOR BOTH P AND Q
C                                    GREATER THAN 1 AND NR GREATER THAN
C                                    OR EQUAL TO 4
      I = 1
      PP = P-ONE
      QQ = Q-ONE
      RR = PP+QQ
      TRR = -RR-RR
      C = RR*ALOG(RR)-PP*ALOG(PP)-QQ*ALOG(QQ)
      X1 = ZERO
      X2 = ZERO
      F1 = ZERO
      F2 = ZERO
      F4 = ZERO
      F5 = ZERO
      X3 = PP/RR
      X4 = ONE
      X5 = ONE
C
      IF (RR.LE.(ONE+EPS)) GO TO 55
      D = SQRT(PP*QQ/(RR-ONE))/RR
      IF (D*(ONE+EPS).GE.X3) GO TO 50
      X2 = X3-D
      X1 = X2-(X2*(ONE-X2))/(PP-RR*X2)
      A1 = PP/X1-QQ/(ONE-X1)
      F1 = EXP(C+PP*ALOG(X1)+QQ*ALOG(ONE-X1))
      F2 = EXP(C+PP*ALOG(X2)+QQ*ALOG(ONE-X2))
   50 IF (D.GE.(ONE-X3-EPS)) GO TO 55
      X4 = X3+D
      X5 = X4-((X4*(ONE-X4))/(PP-RR*X4))
      A5 = QQ/(ONE-X5)-PP/X5
      F4 = EXP(C+PP*ALOG(X4)+QQ*ALOG(ONE-X4))
      F5 = EXP(C+PP*ALOG(X5)+QQ*ALOG(ONE-X5))
C                                  CALCULATE AREAS OF THE TEN REGIONS
   55 P1 = ZERO
      IF (F2.GT.ZERO) P1 = F2*(X3-X2)
      P2 = P1
      IF (F4.GT.ZERO) P2 = F4*(X4-X3)+P1
      P3 = P2
      IF (F1.GT.ZERO) P3 = F1*(X2-X1)+P2
      P4 = P3
      IF (F5.GT.ZERO) P4 = F5*(X5-X4)+P3
      P5 = (ONE-F2)*(X3-X2)+P4
      P6 = (ONE-F4)*(X4-X3)+P5
      P7 = P6
      IF (F2.GT.F1) P7 = (F2-F1)*(X2-X1)*.5+P6
      P8 = P7
      IF (F4.GT.F5) P8 = (F4-F5)*(X5-X4)*.5+P7
      P9 = P8
      IF (F1.GT.ZERO) P9 = F1/A1+P8
      P10 = P9
      IF (F5.GT.ZERO) P10 = F5/A5+P9
C                                  REJECTION PROCEDURE BEGINS HERE
   60 U = GGUBFS(DSEED)*P10
C                                  THE FOUR REGIONS WITH ZERO
C                                    PROBABILITY OF REJECTION
      IF (U.GT.P4) GO TO 80
      IF (U.GT.P1) GO TO 65
      X = X2+U/F2
      GO TO 115
   65 IF (U.GT.P2) GO TO 70
      X = X3+(U-P1)/F4
      GO TO 115
   70 IF (U.GT.P3) GO TO 75
      X = X1+(U-P2)/F1
      GO TO 115
   75 X = X4+(U-P3)/F5
      GO TO 115
C                                  THE TWO REGIONS USING RECTANGULAR
C                                    REJECTION
   80 W = GGUBFS(DSEED)
      IF (U.GT.P5) GO TO 85
      X = X2+(X3-X2)*W
      IF ((U-P4)/(P5-P4).LE.W) GO TO 115
      V = F2+(U-P4)/(X3-X2)
      GO TO 110
   85 IF (U.GT.P6) GO TO 90
      X = X3+(X4-X3)*W
      IF ((P6-U)/(P6-P5).GE.W) GO TO 115
      V = F4+(U-P5)/(X4-X3)
      GO TO 110
C                                  THE TWO TRIANGULAR REGIONS
   90 IF (U.GT.P8) GO TO 100
      W2 = GGUBFS(DSEED)
      IF (W2.GT.W) W = W2
      IF (U.GT.P7) GO TO 95
      X = X1+(X2-X1)*W
      V = F1+2.*W*(U-P6)/(X2-X1)
      IF (V.LT.F2*W) GO TO 115
      GO TO 110
   95 X = X5-W*(X5-X4)
      V = F5+2.*W*(U-P7)/(X5-X4)
      IF (V.LE.F4*W) GO TO 115
      GO TO 110
C                                  THE TWO EXPONENTIAL REGIONS
  100 IF (U.GT.P9) GO TO 105
      U = (P9-U)/(P9-P8)
      X = X1+ALOG(U)/A1
      IF (X.LE.ZERO) GO TO 60
      IF (W.LE.(A1*(X-X1)+ONE)/U) GO TO 115
      V = W*F1*U
      GO TO 110
  105 U = (P10-U)/(P10-P9)
      X = X5-ALOG(U)/A5
      IF (X.GE.ONE) GO TO 60
      IF (W.LE.(A5*(X5-X)+ONE)/U) GO TO 115
      V = W*F5*U
C                                  CHECK EASY REJECTION VIA COMPARISON
C                                    WITH NORMAL DENSITY FUNCTIO
  110 ALV = ALOG(V)
      IF (ALV.GT.(X-X3)*(X-X3)*TRR) GO TO 60
C                                  PERFORM THE STANDARD REJECTION
C
      IF (ALV.GT.(C+PP*ALOG(X)+QQ*ALOG(ONE-X))) GO TO 60
  115 R(I) = X
      I = I+1
      IF (I.LE.NR) GO TO 60
      GO TO 9005
C                                  BB METHOD - BOTH P AND Q GREATER
C                                    THAN 1 AND NR LESS THAN 4
  120 PMIN = .FALSE.
      A = Q
      B = P
      IF (Q.LT.P) GO TO 125
      PMIN = .TRUE.
      A = P
      B = Q
  125 ALPHA = A+B
      BETA = SQRT((ALPHA-TWO)/(TWO*A*B-ALPHA))
      GAMA = A+ONE/BETA
      I = 1
  130 CALL GGUBS (DSEED,2,UA)
      V = BETA*ALOG(UA(1)/(ONE-UA(1)))
      W = A*EXP(V)
      Z = UA(1)*UA(1)*UA(2)
      X = GAMA*V-CONS1
      S = A+X-W
      TEMP2 = B+W
      IF ((S+CONS2).GE.FIVE*Z) GO TO 135
      T = ALOG(Z)
      IF (S.GE.T) GO TO 135
      TEMP1 = X+ALPHA*ALOG(ALPHA/TEMP2)
      IF (TEMP1.LT.T) GO TO 130
  135 IF (PMIN) GO TO 140
      R(I) = B/TEMP2
      GO TO 145
  140 R(I) = W/TEMP2
  145 I = I+1
      IF (I.LE.NR) GO TO 130
      GO TO 9005
C                                  ATKINSONS METHOD - ONE OF P AND Q
C                                    LESS THAN 1 AND THE OTHER GREATER
C                                    THAN 1
  150 A = P
      B = Q
  155 TEMPA = A-ONE
      TEMPB = B-ONE
      T = -TEMPA/(B+ONE-A)
      S1 = ONE
      S2 = T**TEMPA
      TEMP = B*T
      GAMA = TEMP/(TEMP+A*((ONE-T)**B))
      GO TO 165
  160 A = Q
      B = P
      IDX = 1
      GO TO 155
  165 PMA = P-A
      QMB = Q-B
      OMT = ONE-T
      OMG = ONE-GAMA
      AIN = ONE/A
      BIN = ONE/B
  170 CALL GGUBS (DSEED,2,UA)
      IF (UA(1).GT.GAMA) GO TO 180
      X = T*(UA(1)/GAMA)**AIN
      F1 = (ONE-X)**TEMPB
      IF (IDX.LE.1) GO TO 175
      F1 = F1*X**PMA
  175 IF (S1*UA(2).GT.F1) GO TO 170
      GO TO 190
  180 X = ONE-OMT*((ONE-UA(1))/OMG)**BIN
      F2 = X**TEMPA
      IF (IDX.LE.1) GO TO 185
      F2 = F2*(ONE-X)**QMB
  185 IF (S2*UA(2).GT.F2) GO TO 170
  190 I = I+1
      R(I) = X
      IF (IDX.EQ.1) R(I) = ONE-X
      IF (I.LT.NR) GO TO 170
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - CTRBYC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - ANALYSIS OF A CONTINGENCY TABLE.
C
C   USAGE               - CALL CTRBYC (A,IRC,JRC,IR,IC,STAT,IER)
C
C   ARGUMENTS    A      - INPUT CONTINGENCY TABLE CONTAINING OBSERVED
C                           COUNTS. THIS (IR+1) BY (IC+1) BY 3 ARRAY
C                           WILL CONTAIN THE EXPECTED VALUES AND CHI-
C                           SQUARED CONTRIBUTIONS AS WELL AS THE ROW
C                           AND COLUMN SUMS, ON OUTPUT.
C                IRC    - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                JRC    - COLUMN DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IR     - NUMBER OF ROWS USED FOR CALCULATING
C                           STATISTICS. (INPUT)
C                IC     - NUMBER OF COLUMNS USED FOR CALCULATING
C                           STATISTICS. (INPUT)
C                STAT   - OUTPUT VECTOR OF LENGTH 6 CONTAINING,
C                           (1) RESULTING CHI-SQUARED STATISTIC
C                           (2) DEGREES OF FREEDOM OF CHI-SQUARED
C                           (3) EXACT MEAN
C                           (4) EXACT STANDARD DEVIATION
C                           (5) EXACT PROBABILITY OF ALL EXTREME
C                                 DISTRIBUTIONS OF A 2X2 TABLE WHEN
C                                 SOME CELL HAS AN EXPECTED VALUE LESS
C                                 THAN 5.
C                           (6) TOCHERS PROBABILITY SATISFYING THE
C                                 SAME CONDITIONS AS (5) ABOVE.
C                           IF CONDITIONS IN (5) ARE NOT SATISFIED
C                           STAT(5),STAT(6) ARE SET TO ZERO.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           TERMINAL ERROR
C                             IER=129 IMPLIES ALL ROW, COLUMN SUMS ARE
C                               COMPUTED BUT ONE IS ZERO.
C                           WARNING ERROR
C                             IER=34 IMPLIES 20 PERCENT OF EXPECTED
C                               VALUES LESS THAN 5.
C                             IER=35 IMPLIES DEGREES OF FREEDOM GREATER
C                               THAN 30, AND USER SHOULD CONSIDER USING
C                               EXACT MEAN, STANDARD DEVIATION, AND
C                               NORMAL DISTRIBUTION FUNCTION. (ALSO
C                               IMPLIES IER=34 CONDITION MET).
C                             IER=36 IMPLIES SOME EXPECTED VALUE LESS
C                               THAN 2 AND STATISTICS MAY NOT BE
C                               VALID. (ALSO IMPLIES IER=34
C                               CONDITION MET).
C                             IER=37 IMPLIES SOME EXPECTED VALUE LESS
C                               THAN 1 AND STATISTICS MAY NOT BE VALID.
C                             IER=38 IMPLIES THAT THE TOTAL FREQUENCY OF
C                               THE TABLE IS THREE OR LESS.  STAT(3) AND
C                               STAT(4) ARE SET TO -1.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ALTHOUGH IT MAY BE ASSUMED THAT THE OBSERVATIONS IN A
C                ARE INTEGERS, THE USER IS CAUTIONED THAT CTRBYC
C                ASSUMES POSITIVE REAL NUMBERS AND NO TEST IS MADE FOR
C                NEGATIVE INPUT DATA.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CTRBYC (A,IRC,JRC,IR,IC,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IRC,JRC,IR,IC,IER
      REAL               A(IRC,JRC,1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            COL,DEL,ELT1,ELT2,ELT5,I1,IA(9),IER1,IINFP,IP1,
     1                   IRJ,IRR(9),IRRJ,IRRJP1,IRR1,IRR9,ITEMP,
     2                   I,J,K,MINA,NC1,NR1,ROW
      REAL               CM1,CP1,C,FF,RC,RM1,RNM1,RN,RP1,RPC,R,S,
     1                   T1,T2,T3,T4,T5,T6,TEMP,T,X,Z1,Z2,Z3,Z4,ZAPPRX
      REAL               NUM,PF,PA,DENOM,F(9)
      DATA               ZAPPRX/0.001/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IER1 = 0
      STAT(1) = 0.0
      STAT(5) = 0.0
      STAT(6) = 0.0
      NR1 = IR + 1
      NC1 = IC + 1
      A(NR1,NC1,1) = 0.0
      DO 5 I = 1, IR
        A(I,NC1,1) = 0.0
    5 CONTINUE
C                                  COMPUTE N , N , N I. .J ..
      IF (IR.LE.IC) GO TO 20
      DO 15 J = 1, IC
        A(NR1,J,1) = 0.0
        DO 10 I = 1, IR
          A(NR1,J,1) = A(NR1,J,1) + A(I,J,1)
          A(I,NC1,1) = A(I,NC1,1) + A(I,J,1)
          IF (J.EQ.IC.AND.A(I,NC1,1).LE.ZAPPRX) IER = 129
   10   CONTINUE
        IF (A(NR1,J,1).LE.ZAPPRX) IER = 129
        A(NR1,NC1,1) = A(NR1,NC1,1) + A(NR1,J,1)
   15 CONTINUE
      GO TO 40
   20 DO 25 J = 1, IC
        A(NR1,J,1) = 0.0
   25 CONTINUE
      DO 35 I = 1, IR
        A(I,NC1,1) = 0.0
        DO 30 J = 1, IC
          A(I,NC1,1) = A(I,NC1,1) + A(I,J,1)
          A(NR1,J,1) = A(NR1,J,1) + A(I,J,1)
          IF (I.EQ.IR.AND.A(NR1,J,1).LE.ZAPPRX) IER = 129
   30   CONTINUE
        IF (A(I,NC1,1).LE.ZAPPRX) IER = 129
        A(NR1,NC1,1) = A(NR1,NC1,1) + A(I,NC1,1)
   35 CONTINUE
   40 IF (IER.EQ.129) GO TO 9000
      ELT1 = 0
      ELT2 = 0
      ELT5 = 0
C                                  COMPUTE E = (N * N )/N IJ I. .J ..
      T1 = 1.0/A(NR1,NC1,1)
      DO 55 I = 1, IR
        DO 50 J = 1, IC
          A(I,J,2) = A(NR1,J,1)*A(I,NC1,1)*T1
          IF (A(I,J,2).GE.5.0) GO TO 45
          ELT5 = ELT5 + 1
          IF (A(I,J,2).GE.2.0) GO TO 45
          ELT2 = 1
          IF (A(I,J,2).LT.1.0) ELT1 = 1
C                                  2 2 COMPUTE X = (A - E ) / E IJ IJ
C                                    IJ
C
   45     A(I,J,3) = (A(I,J,1)-A(I,J,2))**2/A(I,J,2)
          STAT(1) = STAT(1) + A(I,J,3)
   50   CONTINUE
   55 CONTINUE
C                                  2 2 2 COMPUTE E ,E ,E ,X ,X ,X I. .J
C                                    .. I. .J ..
      A(NR1,NC1,2) = 0.0
      A(NR1,NC1,3) = 0.0
      IF (IR.LE.IC) GO TO 75
      DO 60 I = 1, IR
        A(I,NC1,2) = 0.0
        A(I,NC1,3) = 0.0
   60 CONTINUE
      DO 70 J = 1, IC
        A(NR1,J,2) = 0.0
        A(NR1,J,3) = 0.0
        DO 65 I = 1, IR
          A(NR1,J,2) = A(NR1,J,2) + A(I,J,2)
          A(NR1,J,3) = A(NR1,J,3) + A(I,J,3)
          A(I,NC1,2) = A(I,NC1,2) + A(I,J,2)
          A(I,NC1,3) = A(I,NC1,3) + A(I,J,3)
   65   CONTINUE
        A(NR1,NC1,2) = A(NR1,NC1,2) + A(NR1,J,2)
        A(NR1,NC1,3) = A(NR1,NC1,3) + A(NR1,J,3)
   70 CONTINUE
      GO TO 95
   75 DO 80 J = 1, IC
        A(NR1,J,2) = 0.0
        A(NR1,J,3) = 0.0
   80 CONTINUE
      DO 90 I = 1, IR
        A(I,NC1,2) = 0.0
        A(I,NC1,3) = 0.0
        DO 85 J = 1, IC
          A(I,NC1,2) = A(I,NC1,2) + A(I,J,2)
          A(I,NC1,3) = A(I,NC1,3) + A(I,J,3)
          A(NR1,J,2) = A(NR1,J,2) + A(I,J,2)
          A(NR1,J,3) = A(NR1,J,3) + A(I,J,3)
   85   CONTINUE
        A(NR1,NC1,2) = A(NR1,NC1,2) + A(I,NC1,2)
        A(NR1,NC1,3) = A(NR1,NC1,3) + A(I,NC1,3)
   90 CONTINUE
   95 CONTINUE
      S = 0.0
C                                  COMPUTE S = SUM(OVER I) OF 1/N I.
      DO 100 I = 1, IR
        S = S + 1.0/A(I,NC1,1)
  100 CONTINUE
      T = 0.0
C                                  COMPUTE T = SUM(OVER J) OF 1/N .J
      DO 105 J = 1, IC
        T = T + 1.0/A(NR1,J,1)
  105 CONTINUE
      R = IR
      C = IC
      RM1 = R - 1.0
      CM1 = C - 1.0
C                                  COMPUTE DF = (R-1)(C-1)
      STAT(2) = RM1*CM1
C                                  20 PERCENT OF E LESS THAN 5 IJ
      IF (R*C*0.2.LT.ELT5) GO TO 110
C                                  SOME E LESS THAN 1 IJ
      IF (ELT1.NE.0) IER = 37
      GO TO 120
  110 IER1 = 34
      IF (ELT2.EQ.0) GO TO 115
C                                  SOME E LESS THAN 2 IJ
      IER = 36
      GO TO 120
C                                  DF GREATER THAN 30
  115 IF (STAT(2).GT.30.0) IER = 35
C                                  BEGIN EXACT MEAN,VARIANCE
  120 RPC = R + C
      RP1 = NR1
      CP1 = NC1
      RC = R*C
      RN = A(NR1,NC1,1)
      RNM1 = RN - 1.0
C                                  EXACT MEAN
      STAT(3) = -1.0
      STAT(4) = -1.0
      IF (RN.LE.3.0) IER = 38
      IF (RN.LE.3.0) GO TO 125
      STAT(3) = STAT(2)*RN/RNM1
      T5 = RN/(RN-3.0)
      T2 = RM1 + RM1
      T1 = T2*CM1
      T3 = CM1 + CM1
      T2 = RC*RC + T2 + T3
      T3 = -(RC*STAT(2)+RPC*(RPC-2.0))
      T3 = T3 + T3
      T4 = RC*(R-2.0)*(C-2.0)
      T1 = ((T1*RN+T2)*RN+T3)*RN - T4
      STAT(4) = T1/(RNM1*(RN-2.0)*RNM1)
      T1 = -(RP1*RP1-3.0)*T - (CP1*CP1-3.0-T)*S
      T2 = R*(RM1-1.0)*T + C*(CM1-1.0)*S
      T3 = S*T
      T6 = (((T3*RN+T1)*RN+T2)*RN)/(RNM1*(RN-2.0))
C                                  STD DEV EQUALS SQRT(VARIANCE)
      TEMP = STAT(4) + T6
      IF (TEMP.LT.0.0) TEMP = 0.0
      STAT(4) = SQRT(T5*TEMP)
  125 IF (STAT(2).NE.1.0) GO TO 9000
C                                  MUST BE A 2X2
      IF (ELT5.LE.0) GO TO 9000
      DO 130 I = 5, 9
        IA(I) = 0
  130 CONTINUE
      MINA = A(1,1,1)
C                                  FIND SMALLEST ELEMENT
      IRJ = -2
      DO 140 J = 1, 2
        IRJ = IRJ + 2
        DO 135 I = 1, 2
C                                  STORE THE 2 BY 2 TABLE IN AN INTEGER
C                                    VECTOR IA
          K = IRJ + I
          IA(K) = A(I,J,1) + 0.5
          IA(4+I) = IA(4+I) + IA(K)
          IA(6+J) = IA(6+J) + IA(K)
          IF (MINA.GT.IA(K)) GO TO 135
          MINA = IA(K)
          ROW = I
          COL = J
  135   CONTINUE
        IA(9) = IA(9) + IA(6+J)
  140 CONTINUE
      DO 145 I = 1, 9
        IRR(I) = I
  145 CONTINUE
C                                  OBTAIN THE POSITION INDEX AS IF IA
C                                    IS SORTED IN ASCENDING ORDER
      DO 155 I = 1, 8
        IP1 = I + 1
        DO 150 J = IP1, 9
          IF (IA(IRR(I)).LE.IA(IRR(J))) GO TO 150
          ITEMP = IRR(I)
          IRR(I) = IRR(J)
          IRR(J) = ITEMP
  150   CONTINUE
  155 CONTINUE
C                                  COMPUTE LOG FACTORIAL OF THE TABLE
      IEND = 0
      X = 0.0
      FF = 0.0
      DO 170 J = 1, 9
        IRRJ = IRR(J)
        IBEG = IEND + 1
        IEND = IA(IRRJ)
        IF (IEND.LE.0) GO TO 165
        DO 160 I = IBEG, IEND
          X = X + 1.0
          FF = FF + ALOG(X)
  160   CONTINUE
  165   F(IRRJ) = FF
  170 CONTINUE
C                                  OBTAIN THE DIAGONAL ELEMENT OF THE
C                                    SMALLEST ELEMENT IN THE 2 BY 2
      Z1 = MINA
      IF (ROW.EQ.COL) GO TO 175
      Z3 = A(1,1,1)
      Z4 = A(2,2,1)
      GO TO 180
  175 Z3 = A(1,2,1)
      Z4 = A(2,1,1)
  180 Z2 = A(3,3,1) - Z1 - Z3 - Z4
C                                  COMPUTE THE PROBABILITY FOR A GIVEN
C                                    TABLE
      NUM = -F(9)
      DO 185 J = 5, 8
        NUM = NUM + F(J)
  185 CONTINUE
      DENOM = 0.0
      DO 190 J = 1, 4
        DENOM = DENOM + F(J)
  190 CONTINUE
      PF = EXP(NUM-DENOM)
C                                  COMPUTE THE PROBABILITY FOR A MORE
C                                    EXTREME TABLE
      PA = PF
      PT = 0.0
      I1 = MINA
      DO 195 I = 1, I1
        Z3 = Z3 + 1
        Z4 = Z4 + 1
        PA = PA*(Z1*Z2)/(Z3*Z4)
        Z1 = Z1 - 1
        Z2 = Z2 - 1
        PT = PT + PA
  195 CONTINUE
C                                  TOTAL EXACT PROBABILITIES
      STAT(5) = PF + PT
C                                  TOCHERS PROBABILITY
      STAT(6) = PT
 9000 CONTINUE
      IF (IER1.NE.0) CALL UERTST(IER1,'CTRBYC')
      IF (IER.NE.0) GO TO 9005
      IER = IER1
      GO TO 9010
 9005 CALL UERTST(IER,'CTRBYC')
 9010 RETURN
      END

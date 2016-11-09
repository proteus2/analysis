C   IMSL ROUTINE NAME   - RLFOR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - FIT A UNIVARIATE CURVILINEAR REGRESSION
C                           MODEL USING ORTHOGONAL POLYNOMIALS WITH
C                           OPTIONAL WEIGHTING - EASY TO USE VERSION
C
C   USAGE               - CALL RLFOR (XYW,IX,N,RSQ,MDP,ALBP,ANOVA,B,IB,
C                           PRED,IP,WK,IER)
C
C   ARGUMENTS    XYW    - INPUT/OUTPUT N BY MDP(1)+3 DATA MATRIX
C                           CONTAINING THE INDEPENDENT AND DEPENDENT
C                           (RESPONSE) VARIABLE SETTINGS IN COLUMNS 1
C                           AND 2, RESPECTIVELY. COLUMN 3 CONTAINS
C                           POSITIVE WEIGHTS FOR THE RESPONSE SETTINGS
C                           AND THE WEIGHTS SHOULD BE INVERSELY
C                           PROPORTIONAL TO THE RESPONSE VARIANCES. FOR
C                           AN UNWEIGHTED ANALYSIS, SET COLUMN 3 TO
C                           UNITY. IN EITHER CASE, COLUMN 3 MUST NOT
C                           CONTAIN ANY NON-POSITIVE ELEMENTS. THE
C                           REMAINING COLUMNS ARE WORK STORAGE.
C                         ON OUTPUT, COLUMN 1 IS SCALED TO THE INTERVAL
C                           (-2,2). COLUMN 3 CONTAINS THE SQUARE ROOTS
C                           OF THE WEIGHTS.
C                IX     - INPUT ROW DIMENSION OF THE MATRIX XYW EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                N      - INPUT NUMBER OF DATA POINTS (NUMBER OF ROWS
C                           IN XYW). N MUST EXCEED THE VALUE MDP(1)+1.
C                RSQ    - INPUT VALUE SUCH THAT RSQ IS GREATER THAN
C                           ZERO AND LESS THAN OR EQUAL TO 100 AND
C                           USED TO CONTROL THE DEGREE OF THE FITTED
C                           MODEL. ADDITION OF HIGHER ORDER TERMS IS
C                           TERMINATED WHEN THE COEFFICIENT OF
C                           DETERMINATION (PERCENT OF VARIATION
C                           EXPLAINED BY TERMS IN THE MODEL) EXCEEDS
C                           RSQ, UNLESS TERMINATED EARLIER BY PRIOR
C                           REACHING OF DEGREE MDP(1).
C                MDP    - INPUT/OUTPUT OPTION AND CONTROL VECTOR OF
C                           LENGTH 3.
C                         MDP(1) CONTAINS THE INPUT MAXIMUM DEGREE
C                           ALLOWED FOR THE FITTED MODEL.
C                         MDP(2) CONTAINS THE OUTPUT DEGREE OF THE
C                           FITTED MODEL.
C                         MDP(3) CONTAINS THE INPUT PREDICTION ANALYSIS
C                           OPTION. IF MDP(3) IS ZERO, A PREDICTION
C                           ANALYSIS IS NOT PERFORMED. IF MDP(3) IS
C                           NONZERO, THE ANALYSIS IS PERFORMED FOR THE
C                           GIVEN INDEPENDENT VARIABLE SETTINGS (COLUMN
C                           1 OF INPUT XYW), FOR THE TRUE MEAN OF THE
C                           RESPONSE. RESULTS ARE RETURNED IN PRED.
C                ALBP   - INPUT VECTOR OF LENGTH 2 CONTAINING RISK
C                           LEVELS, EACH IN THE EXCLUSIVE INTERVAL
C                           (0,1), FOR THE PARAMETER AND PREDICTION
C                           ANALYSIS INTERVAL ESTIMATES IN LOCATIONS 1
C                           AND 2, RESPECTIVELY. IF MDP(3) IS ZERO,
C                           ALBP(2) IS NOT REQUIRED. THE VALUE 0.05 IS
C                           A TYPICAL ONE FOR AN ALBP ELEMENT.
C                ANOVA  - OUTPUT VECTOR OF LENGTH 13 CONTAINING THE
C                           ANALYSIS OF VARIANCE TABLE INFORMATION AND
C                           OTHER RESULTS.
C                         ANOVA(1), ANOVA(2), AND ANOVA(3) CONTAIN
C                           THE DEGREES OF FREEDOM FOR THE
C                           REGRESSION, RESIDUAL, AND CORRECTED TOTAL
C                           SOURCES OF VARIATION, RESPECTIVELY.
C                         ANOVA(4), ANOVA(5), AND ANOVA(6) CONTAIN
C                           THE SUMS OF SQUARES CORRESPONDING
C                           TO THE DEGREES OF FREEDOM IN LOCATIONS 1,
C                           2, AND 3.
C                         ANOVA(7) AND ANOVA(8) CONTAIN
C                           THE MEAN SQUARES FOR THE REGRESSION
C                           AND RESIDUAL SOURCES OF VARIATION,
C                           RESPECTIVELY.
C                         ANOVA(9) AND ANOVA(10) CONTAIN
C                           THE COMPUTED F-VALUE AND THE
C                           CORRESPONDING TAIL AREA OF THE F
C                           DISTRIBUTION, RESPECTIVELY, FOR TESTING
C                           THE NULL HYPOTHESIS THAT ALL REGRESSION
C                           COEFFICIENT PARAMETERS ARE ZERO.
C                         ANOVA(11) CONTAINS THE PERCENTAGE OF VARIATION
C                           EXPLAINED BY THE ESTIMATED MODEL.
C                         ANOVA(12) CONTAINS THE STANDARD DEVIATION OF
C                           RESIDUALS.
C                         ANOVA(13) CONTAINS THE RESIDUAL STANDARD
C                           DEVIATION EXPRESSED AS A PERCENTAGE OF THE
C                           RESPONSE MEAN.
C                B      - OUTPUT (MDP(1)+3) BY 12 MATRIX CONTAINING IN
C                           THE FIRST MDP(2)+1 ROWS, INFERENCE RESULTS
C                           FOR THE INTERCEPT AND EACH REGRESSION
C                           COEFFICIENT IN ASCENDING DEGREE ORDER BY
C                           ROWS. INTERCEPT RESULTS APPEAR IN THE LAST
C                           ROW.
C                         B(I,1) FOR I=1,...,MDP(2)+1 CONTAIN
C                           THE VARIABLE MEANS, WITH THE RESPONSE
C                           MEAN LAST
C                         B(I,2) FOR I=1,...,MDP(2)+1 CONTAIN
C                           THE PARAMETER ESTIMATES
C                         B(I,3) FOR I=1,...,MDP(2)+1 CONTAIN
C                           THE STANDARD ERRORS OF PARAMETER
C                           ESTIMATES
C                         B(I,4) FOR I=1,...,MDP(2)+1 CONTAIN
C                           THE LOWER BOUNDS FOR INTERVAL ESTIMATES
C                         B(I,5) FOR I=1,...,MDP(2)+1 CONTAIN
C                           THE UPPER BOUNDS FOR INTERVAL ESTIMATES
C                         B(I,6) FOR I=1,...,MDP(2) CONTAIN
C                           THE ADJUSTED SUMS OF SQUARES
C                         B(I,7) FOR I=1,...,MDP(2) CONTAIN
C                           THE PARTIAL F-VALUES
C                         B(I,8) FOR I=1,...,MDP(2) CONTAIN
C                           THE F DISTRIBUTION TAIL AREAS CORRES-
C                           PONDING TO THE PARTIAL F-VALUES.
C                         THE ELEMENTS IN COLUMNS 6,7, AND 8 OF ROW
C                           MDP(2)+1 ARE WORK STORAGE AS IS THE RE-
C                           MAINDER OF B.
C                IB     - INPUT ROW DIMENSION OF THE MATRIX B EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                PRED   - OUTPUT N BY 6 MATRIX, REQUIRED AND USED ONLY
C                           IF MDP(3) IS NONZERO. COLUMN 1 CONTAINS THE
C                           PREDICTED VALUES CORRESPONDING TO THE DATA
C                           POINTS GIVEN IN COLUMN 1 OF XYW. COLUMNS 3
C                           AND 4 CONTAIN THE LOWER AND UPPER 100
C                           (1-ALBP(2)) PERCENT INTERVAL ESTIMATE
C                           LIMITS, RESPECTIVELY, FOR THE TRUE MEAN OF
C                           THE RESPONSE. COLUMNS 2,5, AND 6 ARE WORK
C                           STORAGE.
C                IP     - INPUT ROW DIMENSION OF THE MATRIX PRED EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                WK     - WORK VECTOR OF LENGTH THE MAXIMUM OF
C                           (MDP(1)+1)*(MDP(1)+3) AND 4*N.
C                           WK MUST BE TYPED DOUBLE PRECISION IN THE
C                           CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT A PERFECT FIT TO THE
C                             DATA WAS OBTAINED. IN SUCH A CASE OUTPUT
C                             ELEMENTS 7 THROUGH 13 OF ANOVA AND ALL
C                             OUTPUT COLUMNS OF B EXCEPTING 1 AND 2
C                             HAVE BEEN SET TO ZERO. A PREDICTION
C                             ANALYSIS IS NOT PERFORMED.
C                           IER=34 INDICATES A WARNING OR WARNING WITH
C                             FIX ERROR OCCURRED IN RLDCW,RLFOTW, OR
C                             RLOPDC.
C                         TERMINAL ERROR
C                           IER=131 INDICATES THAT COLUMN 3 OF XYW
C                             CONTAINED NON-POSITIVE ELEMENTS.
C                           IER=132 INDICATES N OR ALBP WAS SPECIFIED
C                             INCORRECTLY.
C                           IER=133 INDICATES A TERMINAL ERROR OCCURRED
C                             IN IMSL ROUTINES RLDCVA, RLDCW, OR RLFOTW
C                             OR THAT AN ERROR OCCURRED IN IMSL
C                             SUBROUTINE MDBETA OR MDBETI.
C                             USAGE OF RLFOR IS NOT APPROPRIATE IN THIS
C                             CASE.
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/MDBETA,MDBETI,MLGAMD=DLGAMA,
C                           RLDCVA,RLDCW,RLDOPM,RLFOTW,RLOPDC,RLPOL,
C                           RLPRDI,UERTST,UGETIO
C                       - SINGLE(H36,H48,H60)/MDBETA,MDBETI,
C                           MLGAMA=ALGAMA,RLDCVA,RLDCW,RLDOPM,RLFOTW,
C                           RLOPDC,RLPOL,RLPRDI,UERTST,UGETIO
C                       - DOUBLE/MDBETA,MDBETI,MLGAMD=DLGAMA,RLDCVA,
C                           RLDCW,RLDOPM,RLFOTW,RLOPDC,RLPOL,RLPRDI,
C                           UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IN A USER APPLICATION PROGRAM CALLING RLFOR, IT
C                MAY BE DESIRABLE, PRIOR TO CALLING RLFOR, TO UTILIZE
C                IMSL ROUTINE BDTRGI (TRANSFORMATION).
C            2.  THE USER MAY WISH TO EXERCISE MORE CONTROL OVER THE
C                THE DEGREE OF THE FITTED MODEL THAN ALLOWED BY JOINT
C                USAGE OF PARAMETERS RSQ AND MDP(1). SETTING RSQ TO
C                100.0 AND MDP(1) TO THE DESIRED DEGREE THE USER MAY
C                MAKE MULTIPLE CALLS TO RLFOR WITH MDP(1) = 1,2,...
C                AN F TEST FOR THE SIGNIFICANCE OF THE DEGREE TERM
C                JUST ADDED MAY BE PERFORMED AFTER EACH CALL USING
C                THE DIFFERENCE IN THE SUCCESSIVE ERROR SUMS OF
C                SQUARES TO OBTAIN THE CONTRIBUTION OF THE LAST TERM
C                ADDED. EQUIVALENTLY, THE PARTIAL F-TEST FOR THE
C                HIGHEST DEGREE TERM IN THE CURRENT MODEL MAY BE
C                USED. SEE OUTPUT PARAMETER B. TWO OF THE STOPPING
C                RULES THAT MAY BE CONSIDERED ARE
C                A.  TERMINATE WITH THE FIRST NON-SIGNIFICANT TERM.
C                B.  TERMINATE WHEN TWO NON-SIGNIFICANT TERMS IN A
C                    ROW ARE FOUND.
C            3.  ALTHOUGH MDP(1) MAY BE AS LARGE AS N-2, TYPICAL
C                USAGE WOULD REQUIRE N GREATER THAN MDP(1)+10.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLFOR(XYW,IX,N,RSQ,MDP,ALBP,ANOVA,B,IB,PRED,IP,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,N,MDP(3),IB,IP,IER
      REAL               XYW(IX,1),RSQ,ALBP(2),ANOVA(13),B(IB,12),
     *                   PRED(IP,6)
      DOUBLE PRECISION   WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IOPT,J,JER,K,KER,MDP1,MDP2,MDP21,MDP3
      REAL               DF,DF1,DF2,F,P,SP
      REAL               ONE,ZERO,T,TEMP,P5,HUND,TEMP2,EPS,XINF,XETA
      DOUBLE PRECISION   SUM2,SUM
      DATA               ONE/1.0/,ZERO/0./
      DATA               HUND/100./,P5/.5/
      DATA               EPS/Z3C100000/
      DATA               XINF /Z7FFFFFFF/
      DATA               XETA /Z00100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      KER = 0
      LER = 0
      IOPT = 1
      DO 5 I = 1, 2
        IF (MDP(3).EQ.0.AND.I.EQ.2) GO TO 5
        IF (ALBP(I).LE.ZERO .OR. ALBP(I).GE.ONE) GO TO 135
    5 CONTINUE
      MDP1 = MDP(1)
      MDP3 = MDP(3)
      IF (N.LE.MDP1+1) GO TO 135
      SUM = 0.D0
      SUM2 = 0.D0
      DO 10 I = 1, N
        IF (XYW(I,3).LE.ZERO) GO TO 145
C                                  SAVE ORIGINAL INDEPENDENT VARIABLES
        XYW(I,4) = XYW(I,1)
        SUM = SUM + DBLE(XYW(I,1))
        SUM2 = SUM2 + DBLE(XYW(I,2))
   10 CONTINUE
      TEMP = SUM/N
      TEMP2 = SUM2/N
      SUM = 0.0D0
      SUM2 = 0.0D0
      DO 15 I = 1, N
        SUM = SUM + DBLE(XYW(I,1)-TEMP)*DBLE(XYW(I,1)-TEMP)
        SUM2 = SUM2 + DBLE(XYW(I,2)-TEMP2)*DBLE(XYW(I,2)-TEMP2)
   15 CONTINUE
      TEMP = SUM
      TEMP2 = SUM2
      CALL RLFOTW(XYW(1,1),XYW(1,2),N,RSQ,MDP1,XYW(1,3),MDP(2),WK,B(1,
     *  9),B(1,10),B(1,11),B(1,12),LER)
      IF (LER.GT.128) GO TO 140
      IF (LER.EQ.68 .OR. LER.EQ.37) IER = 34
      IF (LER.EQ.36 .OR. LER.EQ.38) IER = 33
C                                  TEST FOR PERFECT FIT
      IF (B(MDP1+3,10).LE.EPS*HUND*TEMP2) KER = 33
      MDP2 = MDP(2)
      MDP21 = MDP2 + 1
      DF = N - MDP21
      DF1 = P5
      DF2 = DF*P5
      K = MDP1 + 3
C                                  CHECK FOR ZERO DEGREE OF POLYNOMIAL
      IF (MDP2.EQ.0) GO TO 40
      DO 25 I = 1, MDP2
        SUM = 0.D0
        DO 20 J = 1, N
          SUM = SUM + XYW(J,4)**I
   20   CONTINUE
        B(I,1) = SUM/N
   25 CONTINUE
      IF (KER.EQ.33) GO TO 40
      CALL RLDCW(B(K,10),XYW(1,1),XYW(1,3),N,MDP2,IOPT,B(1,11),B(1,12),
     *  B(1,8),XYW(1,4),IX,JER)
      IF (JER.GT.128) GO TO 140
      IF (MDP3.EQ.0 .OR. KER.EQ.33) GO TO 40
      CALL RLOPDC(XYW(1,1),N,B(1,11),B(1,12),B(1,9),MDP2,IOPT,WK,PRED(1,
     *  1),JER)
      IF (JER.EQ.33) IER = 34
      SP = ONE - ALBP(2)
      CALL MDBETI(SP,DF1,DF2,F,JER)
      IF (JER.NE.0) GO TO 140
      T = DF2*F/(DF1*(ONE-F))
      T = SQRT(T)
      DO 35 I = 1, N
        SUM = B(1,8)
        DO 30 J = 2, MDP21
          SUM = SUM + DBLE(B(J,8))*DBLE(XYW(I,J+2))*DBLE(XYW(I,J+2))
   30   CONTINUE
        PRED(I,2) = SUM
   35 CONTINUE
      CALL RLPRDI(PRED(1,1),PRED(1,2),N,T,IOPT,PRED(1,3),IP)
   40 ANOVA(1) = MDP2
      ANOVA(2) = DF
      ANOVA(3) = N - 1
      ANOVA(4) = 0.0
      IF (MDP2.EQ.0) GO TO 50
      DO 45 I = 1, MDP2
        ANOVA(4) = ANOVA(4) + B(I+2,10)
   45 CONTINUE
   50 ANOVA(5) = B(MDP1+3,10)
      ANOVA(6) = ANOVA(4) + ANOVA(5)
      IF (ANOVA(6).LT.ZERO) ANOVA(6) = ZERO
      ANOVA(4) = ANOVA(6) - ANOVA(5)
      IF (ANOVA(4).LT.ZERO) ANOVA(4) = ZERO
      IF (KER.NE.33.AND.MDP2.NE.0) GO TO 70
      DO 55 I = 7, 13
        ANOVA(I) = ZERO
   55 CONTINUE
      DO 65 I = 2, 8
        DO 60 J = 1, K
          B(J,I) = ZERO
   60   CONTINUE
   65 CONTINUE
      IF (MDP2.GT.0) GO TO 105
      B(MDP21,1) = XYW(MDP21,2)
      B(MDP21,2) = XYW(MDP21,2)
      GO TO 130
C
   70 ANOVA(7) = ZERO
      IF (ANOVA(1).NE.ZERO) ANOVA(7) = ANOVA(4)/ANOVA(1)
      ANOVA(8) = ANOVA(5)/ANOVA(2)
C
      IF (ANOVA(8).GT.ONE) GO TO 75
      IF (ANOVA(7).LT.XINF*ANOVA(8)) GO TO 80
      ANOVA(9) = XINF
      GO TO 85
   75 IF (ANOVA(7).GT.XETA*ANOVA(8)) GO TO 80
      ANOVA(9) = XETA
      GO TO 85
   80 ANOVA(9) = ANOVA(7)/ANOVA(8)
   85 CONTINUE
C
      SP = ONE - ANOVA(2)/(ANOVA(2)+ANOVA(1)*ANOVA(9))
      DF1 = P5*ANOVA(1)
      CALL MDBETA(SP,DF1,DF2,P,JER)
      IF (JER.NE.0) GO TO 140
      ANOVA(10) = ONE - P
C
      IF (ANOVA(6).GT.ONE) GO TO 90
      IF (ANOVA(4).LT.XINF*ANOVA(6)) GO TO 95
      ANOVA(11) = HUND
      GO TO 100
   90 IF (ANOVA(4).GT.XETA*ANOVA(6)) GO TO 95
      ANOVA(11) = ZERO
      GO TO 100
   95 ANOVA(11) = HUND*(ANOVA(4)/ANOVA(6))
  100 CONTINUE
C
      ANOVA(12) = SQRT(ANOVA(8))
  105 SUM = 0.D0
      DO 110 I = 1, N
        SUM = SUM + DBLE(XYW(I,2))
  110 CONTINUE
      B(MDP2+1,1) = SUM/N
      ANOVA(13) = 100.0*N*ANOVA(12)/SUM
      IF (KER.NE.0) GO TO 115
      CALL RLDCVA(B(1,8),MDP2,B(1,11),B(1,12),B(MDP2+2,9),WK,MDP21,JER)
      IF (JER.GT.128) GO TO 140
  115 CALL RLDOPM(B(1,9),MDP2,B(1,11),B(1,12),WK)
      DO 120 I = 2, MDP21
        B(I-1,2) = B(I,9)
        IF (KER.NE.33) B(I-1,3) = SQRT(B(I,8))
  120 CONTINUE
      B(MDP21,2) = B(1,9)
      IF (KER.EQ.33) GO TO 130
      B(MDP21,3) = SQRT(B(1,8))
      DF1 = P5
      SP = ONE - ALBP(1)
      CALL MDBETI(SP,DF1,DF2,F,JER)
      IF (JER.GT.128) GO TO 140
      T = DF2*F/(DF1*(ONE-F))
      T = SQRT(T)
      DO 125 I = 1, MDP2
        TEMP = T*B(I,3)
        B(I,4) = B(I,2) - TEMP
        B(I,5) = B(I,2) + TEMP
        B(I,7) = (B(I,2)/B(I,3))**2
        B(I,6) = B(I,7)*ANOVA(8)
        SP = ONE - ANOVA(2)/(ANOVA(2)+B(I,7))
        CALL MDBETA(SP,DF1,DF2,P,JER)
        IF (JER.NE.0) GO TO 140
        B(I,8) = ONE - P
  125 CONTINUE
      TEMP = T*B(MDP21,3)
      B(MDP21,4) = B(MDP21,2) - TEMP
      B(MDP21,5) = B(MDP21,2) + TEMP
      B(MDP21,6) = ZERO
      B(MDP21,7) = ZERO
      B(MDP21,8) = ZERO
  130 IER = MAX0(IER,KER)
      IF (KER.NE.0) IER = KER
      IF (IER.NE.0) GO TO 9000
      GO TO 9005
  135 IER = 132
      GO TO 9000
  140 IER = 133
      GO TO 9000
  145 IER = 131
 9000 CONTINUE
      IF (IER.NE.KER.AND.KER.NE.0) CALL UERTST(KER,6HRLFOR )
      CALL UERTST(IER,6HRLFOR )
 9005 RETURN
      END

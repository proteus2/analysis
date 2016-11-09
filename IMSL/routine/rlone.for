C   IMSL ROUTINE NAME   - RLONE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - ANALYSIS OF A SIMPLE LINEAR REGRESSION MODEL
C
C   USAGE               - CALL RLONE(XY,IX,N,IMOD,IPRED,ALBAP,DES,ANOVA,
C                           STAT,PRED,IP,NN,IER)
C
C   ARGUMENTS    XY     - INPUT N BY 2 MATRIX OF DATA CONTAINING THE
C                           INDEPENDENT AND DEPENDENT (RESPONSE) VARIA-
C                           BLE SETTINGS IN COLUMNS 1 AND 2, RESPEC-
C                           TIVELY.
C                         ON OUTPUT, XY IS DESTROYED IF IMOD EQUALS
C                           ZERO.
C                IX     - INPUT ROW DIMENSION OF MATRIX XY EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                N      - INPUT NUMBER OF DATA POINTS (NUMBER OF ROWS
C                           IN XY). N MUST EXCEED 2 OR 1 AS IMOD IS
C                           ZERO OR NONZERO, RESPECTIVELY.
C                IMOD   - INPUT MODEL OPTION.  IF IMOD IS ZERO, THE
C                           MODEL INCLUDES AN INTERCEPT PARAMETER.
C                           OTHERWISE A FIT THROUGH THE ORIGIN IS
C                           OBTAINED.
C                IPRED  - INPUT PREDICTION ANALYSIS OPTION. IF IPRED IS
C                           POSITIVE, A PREDICTION ANALYSIS IS PERFORMED
C                           FOR THE GIVEN (COLUMN 1 OF PRED) INDEPENDENT
C                           VARIABLE SETTINGS, FOR THE TRUE MEAN OF THE
C                           RESPONSE AND FOR THE AVERAGE OF IPRED
C                           FUTURE OBSERVATIONS ON THE RESPONSE. RESULTS
C                           ARE RETURNED IN PRED. FOR IPRED NON-
C                           POSITIVE, THE ANALYSIS IS NOT PERFORMED.
C                ALBAP  - INPUT VECTOR OF LENGTH 3 CONTAINING RISK
C                           LEVELS, EACH IN THE EXCLUSIVE INTERVAL
C                           (0,1), FOR THE REGRESSION COEFFICIENT,
C                           INTERCEPT, AND PREDICTION ANALYSIS INTERVAL
C                           ESTIMATES, RESPECTIVELY. IF IMOD IS
C                           NON-ZERO, ALBAP(2) IS NOT REQUIRED. IF IPRED
C                           IS NON-POSITIVE, ALBAP(3) IS NOT REQUIRED.
C                DES    - OUTPUT VECTOR OF LENGTH 5 CONTAINING BASIC
C                           DESCRIPTIVE STATISTICS FOR THE INDEPENDENT
C                           VARIABLE AND RESPONSE (CALL THEM X AND Y,
C                           RESPECTIVELY).  DES(I) CONTAINS WHEN
C                             I = 1, MEAN OF X
C                             I = 2, MEAN OF Y
C                             I = 3, STANDARD DEVIATION OF X
C                             I = 4, STANDARD DEVIATION OF Y
C                             I = 5, ESTIMATE OF CORRELATION BETWEEN X
C                                    AND Y.
C                ANOVA  - OUTPUT VECTOR OF LENGTH 14 CONTAINING STANDARD
C                           ANALYSIS OF VARIANCE TABLE ENTRIES AND
C                           OTHER RESULTS.  ANOVA(I) CONTAINS WHEN
C                             I = 1,2,3, THE DEGREES OF FREEDOM FOR RE-
C                               GRESSION, RESIDUAL, AND TOTAL SOURCES OF
C                               VARIATION, RESPECTIVELY.
C                             I = 4,5,6, THE SUMS OF SQUARES CORRES-
C                               PONDING TO THE DEGREES OF FREEDOM IN
C                               LOCATIONS 1,2, AND 3.
C                             I = 7,8, THE MEAN SQUARES FOR THE REGRES-
C                               SION AND RESIDUAL SOURCES OF VARIATION,
C                               RESPECTIVELY.
C                             I = 9,10, THE COMPUTED F-VALUE AND COR-
C                               RESPONDING TAIL AREA OF THE F-DISTRIBU-
C                               TION, RESPECTIVELY, FOR TESTING THE
C                               NULL HYPOTHESIS THAT THE REGRESSION
C                               COEFFICIENT PARAMETER IS ZERO.
C                             I = 11, PERCENTAGE OF VARIATION EXPLAINED
C                               BY THE FITTED MODEL.
C                             I = 12, STANDARD DEVIATION OF RESIDUALS.
C                             I = 13, RESIDUAL STANDARD DEVIATION EX-
C                               PRESSED AS A PERCENTAGE OF THE RESPONSE
C                               MEAN.
C                             I = 14, DURBIN-WATSON STATISTIC.
C                STAT   - OUTPUT VECTOR OF LENGTH 9 CONTAINING
C                           INFERENCES FOR THE MODEL PARAMETERS. THE
C                           LAST 5 ELEMENTS ARE DEFINED ONLY IF IMOD IS
C                           ZERO. STAT(I) CONTAINS WHEN
C                             I = 1, ESTIMATE OF THE REGRESSION COEF-
C                               FICIENT.
C                             I = 2, STANDARD ERROR OF THE ESTIMATE IN
C                               STAT(1).
C                             I = 3,4, LOWER AND UPPER 100 (1-ALBAP(1))
C                               PERCENT INTERVAL ESTIMATE LIMITS, RE-
C                               SPECTIVELY, FOR THE REGRESSION
C                               COEFFICIENT PARAMETER.
C                             I = 5, ESTIMATE OF THE INTERCEPT
C                             I = 6, STANDARD ERROR OF THE ESTIMATE IN
C                               STAT(5)
C                             I = 7,8, LOWER AND UPPER 100(1-ALBAP(2))
C                               PERCENT INTERVAL ESTIMATE LIMITS, RE-
C                               SPECTIVELY, FOR THE INTERCEPT PARAMETER
C                               AND DEFINED ONLY IF IMOD IS ZERO.
C                             I = 9, ESTIMATE OF THE COVARIANCE BETWEEN
C                               THE REGRESSION COEFFICIENT AND INTERCEPT
C                               ESTIMATORS.
C                PRED   - INPUT/OUTPUT NN BY 7 MATRIX, USED ONLY IF
C                           IPRED IS POSITIVE. ON INPUT, COLUMN 1
C                           CONTAINS THE INDEPENDENT VARIABLE SETTINGS
C                           FOR WHICH A PREDICTION ANALYSIS IS DESIRED.
C                         ON OUTPUT COLUMN 2 CONTAINS THE PREDICTED
C                           VALUES CORRESPONDING TO COLUMN 1.  COLUMNS
C                           4 AND 5 CONTAIN THE LOWER AND UPPER 100(1-
C                           ALBAP(3)) PERCENT INTERVAL ESTIMATE LIMITS,
C                           RESPECTIVELY, FOR THE TRUE MEAN OF THE RE-
C                           SPONSE.  COLUMNS 6 AND 7 CONTAIN RESULTS
C                           ANALOGOUS TO THOSE IN COLUMNS 4 AND 5, FOR
C                           THE AVERAGE OF IPRED FUTURE OBSERVA-
C                           TIONS ON THE RESPONSE.  COLUMN 3 IS WORK
C                           SPACE.
C                IP     - INPUT ROW DIMENSION OF MATRIX PRED EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                NN     - INPUT NUMBER OF POINTS FOR WHICH A PREDICTION
C                           ANALYSIS IS DESIRED (NUMBER OF ROWS IN
C                           PRED). NN MUST BE POSITIVE IF IPRED IS
C                           POSITIVE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THE FIRST COLUMN OF
C                             XY WAS CONSTANT.
C                           IER = 130 INDICATES N, ALBAP, OR NN WAS
C                             SPECIFIED INCORRECTLY.
C                           IER = 131 INDICATES THAT AN ERROR OCCURRED
C                             IN SUBROUTINE MDBETA OR MDBETI.
C                         WARNING ERROR
C                           IER = 36 INDICATES THAT A PERFECT FIT TO
C                             THE DATA WAS OBTAINED. IN SUCH A CASE,
C                             ALL OUTPUT ELEMENTS HAVE BEEN SET TO ZERO,
C                             EXCEPT ELEMENTS 1 THROUGH 6 OF ANOVA,
C                             ALL ELEMENTS OF DES, AND ELEMENTS 1 AND 5
C                             OF STAT.
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/MDBETA,MDBETI,MLGAMD=DLGAMA,
C                           RLPRDI,UERTST,UGETIO,VBLA=DSDOT,VBLA=SDOT
C                       - SINGLE(H36)/MDBETA,MDBETI,MLGAMA=ALGAMA,
C                           RLPRDI,UERTST,UGETIO,VBLA=DSDOT,VBLA=SDOT
C                       - SINGLE(H48,H60)/MDBETA,MDBETI,
C                           MLGAMA=ALGAMA,RLPRDI,UERTST,UGETIO,VBLA=SDOT
C                       - DOUBLE/MDBETA,MDBETI,MLGAMD=DLGAMA,RLPRDI,
C                           UERTST,UGETIO,VBLA=DDOT
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IN A USER APPLICATION PROGRAM CALLING RLONE, IT MAY
C                BE DESIRABLE, PRIOR TO CALLING RLONE, TO UTILIZE IMSL
C                ROUTINES RLFITI (LACK OF FIT) AND/OR BDTRGI
C                (TRANSFORMATION). AFTER THE RLONE CALL, IMSL ROUTINES
C                RLINCF (RESPONSE CONTROL), RLINPF (INVERSE PREDICTION),
C                RLRES (RESIDUAL ANALYSIS), AND USPLT (PLOTTING)
C                (USPLTD FOR DOUBLE PRECISION) MAY BE UTILIZED.
C            2.  ALTHOUGH N MAY BE AS SMALL AS 3 (OR 2 IF IMOD IS NOT
C                EQUAL TO ZERO), TYPICAL USAGE WOULD REQUIRE N GREATER
C                THAN 11.
C            3.  A SIMPLE AND TYPICAL USAGE OF RLONE WOULD HAVE
C                IMOD=IPRED=0, ALBAP(1)=ALBAP(2)=0.05, AND PRED
C                DIMENSIONED 1 BY 1 IN THE MAIN PROGRAM.
C            4.  ANOVA(11) = 100.0*ANOVA(4)/ANOVA(6)
C                FOR ALL VALUES OF IMOD.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLONE  (XY,IX,N,IMOD,IPRED,ALBAP,DES,ANOVA,STAT,PRED,
     *                   IP,NN,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,N,IMOD,IPRED,IP,NN,IER
      REAL               XY(IX,2),ALBAP(3),DES(5),ANOVA(14),STAT(9),
     *                   PRED(IP,7)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,KER
      REAL               DF1,DF2,F,P,SP
      REAL               SUM1,SUM2,SUMA,SUMB,EPS,TEN,TEPS,XBAR,YBAR,V,
     *                   P5,ZERO,ONE,HUND,W,Q
      REAL               SDOT
      DOUBLE PRECISION   DSDOT
      DOUBLE PRECISION   SUM,SUM12,DZERO
      DATA               ZERO /0.0/,ONE /1.0/,HUND /100.0/,P5 /.5/,
     *                   TEN/10.0/
      DATA               DZERO/0.0D0/
      DATA               EPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      KER = 0
      IF (N.LE.2 .AND. IMOD.EQ.0) GO TO 145
      IF (N.LE.1 .AND. IMOD.NE.0) GO TO 145
      DO 5 I=1,3
         IF (I.EQ.2 .AND. IMOD.NE.0) GO TO 5
         IF (I.EQ.3 .AND. IPRED.LE.0) GO TO 5
         IF (ALBAP(I).LE.ZERO .OR. ALBAP(I).GE.ONE) GO TO 145
    5 CONTINUE
      IF (IPRED.LE.0) GO TO 10
      IF (NN.LE.0) GO TO 145
C                                  CHECK FOR X CONSTANT
   10 TEPS = TEN*EPS*ABS(XY(1,1))
      DO 15 I=2,N
         IF (ABS(XY(I,1)-XY(1,1)).GT.TEPS) GO TO 20
   15 CONTINUE
      GO TO 140
C                                  CHECK FOR Y CONSTANT
   20 TEPS = TEN*EPS*ABS(XY(1,2))
      DO 25 I=2,N
         IF (ABS(XY(I,2)-XY(1,2)).GT.TEPS) GO TO 30
   25 CONTINUE
      KER = 36
      IF (IMOD.NE.0 .AND. XY(1,2).NE.ZERO) KER = 30
   30 K = N-1
      SUM = DZERO
      SUM12 = DZERO
      DO 35 I=1,N
         SUM12 = SUM12+XY(I,2)
         SUM = SUM+XY(I,1)
   35 CONTINUE
      XBAR = SUM/N
      YBAR = SUM12/N
      IF (IMOD.NE.0) GO TO 45
      DO 40 I=1,N
         XY(I,1) = XY(I,1)-XBAR
         XY(I,2) = XY(I,2)-YBAR
   40 CONTINUE
      K = K-1
   45 CONTINUE
      SUM1 = SDOT(N,XY(1,1),1,XY(1,1),1)
      SUM2 = SDOT(N,XY(1,2),1,XY(1,2),1)
      SUM12 = DSDOT(N,XY(1,1),1,XY(1,2),1)
      STAT(1) = SUM12/SUM1
      STAT(5) = YBAR-STAT(1)*XBAR
      DES(1) = XBAR
      DES(2) = YBAR
      DES(5) = ZERO
      IF (IMOD.EQ.0) GO TO 60
      SUMA = ZERO
      SUMB = ZERO
      SUM = DZERO
      DO 50 I=1,N
         V = XY(I,1)-XBAR
         W = XY(I,2)-YBAR
         SUMA = SUMA+V*V
         SUMB = SUMB+W*W
         SUM = SUM+DBLE(V)*DBLE(W)
   50 CONTINUE
      DES(3) = SQRT(SUMA/(N-1))
      DES(4) = SQRT(SUMB/(N-1))
      IF (KER.EQ.30) GO TO 55
      DES(5) = SUM/(DES(3)*DES(4)*(N-1))
      GO TO 65
   55 KER = 0
      GO TO 65
   60 DES(3) = SQRT(SUM1/(N-1))
      DES(4) = SQRT(SUM2/(N-1))
      IF (KER.EQ.36) GO TO 65
      DES(5) = SUM12/(DES(3)*DES(4)*(N-1))
   65 ANOVA(1) = ONE
      ANOVA(2) = K
      ANOVA(3) = K+ONE
      ANOVA(4) = STAT(1)*SUM12
      ANOVA(6) = SUM2
      ANOVA(5) = ANOVA(6)-ANOVA(4)
      IF (ANOVA(5).GT.EPS*HUND*SUM2) GO TO 90
      ANOVA(5) = ZERO
      KER = 36
      DO 70 I=7,14
         ANOVA(I) = ZERO
   70 CONTINUE
      DO 75 I=2,9
         IF (I.NE.5) STAT(I) = ZERO
   75 CONTINUE
      IF (IPRED.LE.0) GO TO 135
      DO 85 I=1,NN
         DO 80 J=1,7
            PRED(I,J) = ZERO
   80    CONTINUE
   85 CONTINUE
      GO TO 135
   90 ANOVA(7) = ANOVA(4)
      ANOVA(8) = ANOVA(5)/ANOVA(2)
      STAT(2) = SQRT(ANOVA(8)/SUM1)
      ANOVA(9) = ANOVA(7)/ANOVA(8)
      DF1 = P5
      DF2 = P5*ANOVA(2)
      SP = ONE-ANOVA(2)/(ANOVA(2)+ANOVA(9))
      CALL MDBETA(SP,DF1,DF2,P,IER)
      IF (IER.NE.0) GO TO 150
      ANOVA(10) = ONE-P
      ANOVA(11) = HUND*ANOVA(4)/ANOVA(6)
      ANOVA(12) = SQRT(ANOVA(8))
      ANOVA(13) = HUND*ANOVA(12)/DES(2)
      SP = ONE-ALBAP(1)
      CALL MDBETI(SP,DF1,DF2,F,IER)
      IF (IER.NE.0) GO TO 150
      Q = DF2*F/(DF1*(ONE-F))
      Q = STAT(2)*SQRT(Q)
      STAT(3) = STAT(1)-Q
      STAT(4) = STAT(1)+Q
      IF (IMOD.EQ.0) GO TO 100
      DO 95 I=5,9
         STAT(I) = ZERO
   95 CONTINUE
      GO TO 105
  100 STAT(6) = SQRT(ANOVA(8)*(ONE/N+XBAR*XBAR/SUM1))
      SP = ONE-ALBAP(2)
      CALL MDBETI(SP,DF1,DF2,F,IER)
      IF (IER.NE.0) GO TO 150
      Q = DF2*F/(DF1*(ONE-F))
      Q = STAT(6)*SQRT(Q)
      STAT(7) = STAT(5)-Q
      STAT(8) = STAT(5)+Q
      STAT(9) = -ANOVA(8)*XBAR/SUM1
  105 Q = YBAR-STAT(5)-STAT(1)*XBAR
      V = XY(1,2)-STAT(1)*XY(1,1)
      IF (IMOD.EQ.0) V = V+Q
      SUMA = ZERO
      SUMB = V*V
      DO 110 I=2,N
         W = XY(I,2)-STAT(1)*XY(I,1)
         IF (IMOD.EQ.0) W = W+Q
         V = W-V
         SUMA = SUMA+V*V
         SUMB = SUMB+W*W
         V = W
  110 CONTINUE
      ANOVA(14) = SUMA/SUMB
      IF (IPRED.LE.0) GO TO 135
      IF (IMOD.EQ.0) GO TO 120
      DO 115 I=1,NN
         PRED(I,2) = STAT(1)*PRED(I,1)
         PRED(I,3) = PRED(I,1)*PRED(I,1)/SUM1
  115 CONTINUE
      GO TO 130
  120 SUM2 = ONE/N
      SUM1 = ONE/SUM1
      DO 125 I=1,NN
         PRED(I,2) = STAT(1)*PRED(I,1)+STAT(5)
         V = PRED(I,1)-XBAR
         PRED(I,3) = SUM2+V*V*SUM1
  125 CONTINUE
  130 SP = ONE-ALBAP(3)
      CALL MDBETI(SP,DF1,DF2,F,IER)
      IF (IER.NE.0) GO TO 150
      Q = DF2*F/(DF1*(ONE-F))
      Q = SQRT(Q*ANOVA(8))
      CALL RLPRDI(PRED(1,2),PRED(1,3),NN,Q,IPRED,PRED(1,4),IP)
  135 IER = KER
      IF (IER.NE.0) GO TO 9000
      GO TO 9005
  140 IER = 129
      GO TO 9000
  145 IER = 130
      GO TO 9000
  150 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,6HRLONE )
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - SSRBLK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - STRATIFIED RANDOM SAMPLING WITH CONTINUOUS
C                           DATA-INFERENCES REGARDING THE POPULATION
C                           MEAN AND TOTAL USING RATIO OR REGRESSION
C                           ESTIMATION
C
C   USAGE               - CALL SSRBLK (Y,IY,IOPT,NBR,NH,IN,XBARH,ALPHA,
C                           TEMP,BR,HMUSIG,IH,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT NBR(2) BY 2 SUBMATRIX OF THE MATRIX
C                           (CALL IT YY) CONTAINING THE ENTIRE
C                           STRATIFIED RANDOM SAMPLE.  THE SUBMATRIX Y
C                           MUST BE EITHER THE SAME AS YY OR IT MUST
C                           CONTAIN ALL OR PART OF THE SAMPLE FOR A
C                           SINGLE STRATUM ONLY.  IN THE CASE WHERE Y
C                           IS A PROPER SUBSET OF YY, THE LAST SUBMATRIX
C                           FROM ANY ONE OR MORE OF THE STRATA MAY HAVE
C                           FEWER THAN NBR(2) ROWS.  IN ALL CASES THE
C                           PAIRS OF OBSERVATIONS WITHIN ANY ONE STRATUM
C                           MUST APPEAR CONTIGUOUSLY IN YY.  THE
C                           AUXILIARY VARIABLE SETTINGS ARE IN COLUMN 1
C                           WITH CORRESPONDING VARIABLE OF INTEREST
C                           SETTINGS IN COLUMN 2. SEE PROGRAMMING NOTES.
C                IY     - INPUT ROW DIMENSION OF MATRIX Y EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IOPT   - INPUT ESTIMATION OPTION.
C                           IF IOPT = 0, RATIO ESTIMATION IS TO BE USED
C                             FOR INFERENCE ABOUT THE POPULATION MEAN
C                             AND TOTAL.
C                           IF IOPT IS NEGATIVE, REGRESSION ESTIMATION,
C                             WITH REGRESSION COEFFICIENT(S) PREASSIGNED
C                             VIA BR, IS TO BE USED.
C                           IF IOPT IS POSITIVE, REGRESSION ESTIMATION,
C                             WITH REGRESSION COEFFICIENT(S) ESTIMATED
C                             FROM THE DATA, IS TO BE USED
C                NBR    - INPUT VECTOR OF LENGTH 9.  NBR(I) CONTAINS
C                           WHEN
C                             I=1, NUMBER OF STRATA INTO WHICH THE
C                               SAMPLE IS DIVIDED.
C                             I=2, NUMBER OF PAIRS OF OBSERVATIONS IN
C                               EACH SUBMATRIX Y, NOT INCLUDING THE LAST
C                               SUBMATRIX IN EACH STRATUM, WHERE THE
C                               NUMBER MAY BE LESS THAN OR EQUAL TO
C                               NBR(2).  HOWEVER NBR(2) SHOULD BE THE
C                               SAME FOR ALL CALLS.
C                             I=3, THE NUMBER OF THE SUBMATRIX STORED
C                               IN Y. SEE REMARKS.
C                             I=4, THE TEMPORARY MEAN INDICATOR.  IF
C                               NBR(4) = 0, THE USER SUPPLIES TEMPORARY
C                               MEANS IN TEMP.  THE USER MAY CHANGE TEMP
C                               ONLY WHEN SUBMATRIX Y IS THE FIRST SET
C                               OF OBSERVATIONS IN A STRATUM.  IF NBR(4)
C                               DOES NOT EQUAL ZERO, THE FIRST PAIR OF
C                               OBSERVATIONS IN EACH STRATUM IS
C                               UTILIZED.
C                             I=5, OPTION FOR A SEPARATE OR COMBINED
C                               ESTIMATION TECHNIQUE.  IF NBR(5) = 0,
C                               SEPARATE RATIO OR REGRESSION ESTIMATION
C                               IS TO BE USED.  OTHERWISE A COMBINED
C                               ESTIMATION TECHNIQUE IS UTILIZED.
C                             I=(6-9), INTEGER WORK AREA.
C                NH     - INPUT NBR(1) BY 3 MATRIX CONTAINING THE
C                           NUMBER OF PAIRS OF SAMPLE UNITS IN EACH
C                           STRATUM FOR THE SAMPLE (MATRIX YY) AND FOR
C                           THE POPULATION.  THE STRATA SIZES FOR THE
C                           SAMPLE AND THE POPULATION ARE IN COLUMNS
C                           1 AND 2 RESPECTIVELY.  THE STRATA SIZES MUST
C                           BE ORDERED IN CORRESPONDENCE WITH THE
C                           ORDERING OF STRATA IN YY.  IN THE CASE WHERE
C                           POPULATION STRATA SIZES ARE NOT KNOWN,
C                           ESTIMATES MUST BE ENTERED IN THEIR PLACE.
C                           COLUMN 3 IS WORK STORAGE.
C                IN     - INPUT ROW DIMENSION OF MATRIX NH EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                XBARH  - INPUT VECTOR OF LENGTH NBR(1) CONTAINING, FOR
C                           EACH STRATUM, THE POPULATION MEAN OF THE
C                           AUXILIARY VARIATE, PROVIDED NBR(5) = 0.
C                           OTHERWISE, ONLY XBARH(1) IS DEFINED AND IT
C                           MUST CONTAIN THE POPULATION MEAN OF THE
C                           AUXILIARY VARIATE.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVALS FOR THE MEAN AND TOTAL
C                           PARAMETERS.
C                TEMP   - INPUT VECTOR OF LENGTH 2.  IF NBR(4) = 0, TEMP
C                           MUST CONTAIN TEMPORARY MEANS FOR THE TWO
C                           COLUMNS OF YY, RESPECTIVELY.  OTHERWISE TEMP
C                           IS UNDEFINED.
C                BR     - INPUT/OUTPUT VECTOR OF LENGTH NBR(1).
C                           IF IOPT IS NEGATIVE, BR IS DEFINED ONLY ON
C                             INPUT.  WHEN NBR(5) = 0, BR CONTAINS A
C                             PREASSIGNED REGRESSION COEFFICIENT FOR
C                             EACH STRATUM.  WHEN NBR(5) IS NONZERO,
C                             ONLY BR(1) IS DEFINED AND MUST CONTAIN THE
C                             PREASSIGNED REGRESSION COEFFICIENT COMMON
C                             TO ALL STRATA.
C                           IF IOPT IS POSITIVE, BR IS DEFINED ONLY ON
C                             OUTPUT.  WHEN NBR(5) = 0, BR CONTAINS THE
C                             ESTIMATED REGRESSION COEFFICIENT FOR EACH
C                             STRATUM.  WHEN NBR(5) IS NONZERO, ONLY
C                             BR(1) IS DEFINED AND CONTAINS THE
C                             ESTIMATED REGRESSION COEFFICIENT COMMON
C                             TO ALL STRATA.
C                           IF IOPT = 0, BR IS DEFINED ONLY ON OUTPUT.
C                             WHEN NBR(5) = 0, BR CONTAINS THE ESTIMATE
C                             OF THE RATIO FOR EACH STRATUM.  WHEN
C                             NBR(5) IS NONZERO, ONLY BR(1) IS DEFINED
C                             AND CONTAINS THE COMBINED ESTIMATE OF THE
C                             RATIO.
C                HMUSIG - OUTPUT NBR(1) BY 6 MATRIX CONTAINING THE
C                           WITHIN STRATA MEAN ESTIMATES, VARIANCE
C                           ESTIMATES AND COEFFICIENT OF VARIATION OF
C                           THE MEAN ESTIMATES.  THE ESTIMATES FOR THE
C                           AUXILIARY VARIABLE ARE IN COLUMNS 1, 2, AND
C                           3 RESPECTIVELY.  THE CORRESPONDING ESTIMATES
C                           FOR THE VARIABLE OF INTEREST SETTINGS ARE
C                           IN COLUMNS 4, 5, AND 6, RESPECTIVELY.  THE
C                           ESTIMATES ARE ORDERED IN CORRESPONDENCE
C                           WITH THE ORDERING OF STRATA IN YY.
C                IH     - INPUT ROW DIMENSION OF MATRIX HMUSIG EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                STAT   - OUTPUT VECTOR OF LENGTH 20.  STAT(I) CONTAINS
C                           WHEN
C                             I=1, ESTIMATE OF THE MEAN.
C                             I=2, ESTIMATE OF THE TOTAL.
C                             I=3, VARIANCE ESTIMATE OF THE MEAN
C                               ESTIMATE.
C                             I=4, VARIANCE ESTIMATE OF THE TOTAL
C                               ESTIMATE.
C                             I=5, LOWER CONFIDENCE LIMIT FOR THE MEAN.
C                             I=6, UPPER CONFIDENCE LIMIT FOR THE MEAN.
C                             I=7, LOWER CONFIDENCE LIMIT FOR THE TOTAL.
C                             I=8, UPPER CONFIDENCE LIMIT FOR THE TOTAL.
C                             I=9, ESTIMATE OF THE COEFFICIENT OF
C                               VARIATION OF THE MEAN AND TOTAL
C                               ESTIMATES.
C                             I=10, ESTIMATE OF THE MEAN FOR THE
C                               AUXILIARY VARIATE, CONSIDERING ONLY THE
C                               DATA IN COLUMN 1 OF YY.  DEFINED ONLY
C                               WHEN NBR(5) IS NONZERO.
C                             I=11, ESTIMATE OF THE MEAN FOR THE
C                               VARIABLE OF INTEREST, CONSIDERING ONLY
C                               THE DATA IN COLUMN 2 OF YY.  DEFINED
C                               ONLY WHEN NBR(5) IS NONZERO.
C                             I=12-20, WORK AREA.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           WARNING ERROR
C                             IER = 33 INDICATES THE MAXIMUM ELEMENT IN
C                               COLUMN 3 OF HMUSIG TIMES THE SQUARE ROOT
C                               OF NBR(1) EXCEEDS 30, WHEN NBR(5) = 0.
C                               WHEN NBR(5) IS NONZERO, THIS INDICATES
C                               THE ESTIMATE OF THE COEFFICIENT OF
C                               VARIATION OF STAT(10) EXCEEDS 10.  IN
C                               EITHER CASE THE BIAS IN THE ESTIMATES
C                               STAT(1) AND STAT(2) MAY NOT BE
C                               NEGLIGIBLE.
C                           TERMINAL ERROR
C                             IER = 130 INDICATES NBR(1) IS LESS THAN 2
C                               OR THAT NBR(2) IS LESS THAN 1.
C                             IER = 131 INDICATES NBR(3) IS LESS THAN 1
C                               OR THAT NBR(2)*(NBR(3)-NBR(1)) EXCEEDS
C                               THE NUMBER OF PAIRS OF ELEMENTS IN YY
C                               (TOTAL OF COLUMN 1 OF NH).
C                             IER = 132 INDICATES THAT AT LEAST ONE
C                               ELEMENT OF COLUMN 1 OF NH IS LESS THAN 3
C                               OR THAT AT LEAST ONE ELEMENT IN COLUMN 2
C                               IS LESS THAN THE ELEMENT IN THE SAME ROW
C                               IN COLUMN 1.
C                             IER = 133 INDICATES THAT ALPHA IS NOT IN
C                               THE EXCLUSIVE INTERVAL (0,1).
C                             IER = 134 INDICATES THAT NBR(2) EXCEEDS
C                               THE TOTAL OF COLUMN 1 OF NH, THE TOTAL
C                               SAMPLE SIZE.  WHEN NBR(2) IS LESS
C                               THAN THE TOTAL SAMPLE SIZE, IER = 134
C                               INDICATES THAT NBR(2) IS GREATER
C                               THAN SOME NH(I,1),I=1,2,...,NBR(1).
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      BETWEEN THE FIRST AND LAST CALL TO SSRBLK, EXCEPT FOR
C                OPTIONALLY CHANGING TEMP BETWEEN STRATA, ONLY NBR(3)
C                MAY BE MODIFIED AND IT SHOULD FOLLOW THE PATTERN
C                1,2,... . THOUGH THIS PATTERN IS THE OBVIOUS ONE TO
C                FOLLOW, IT IS NOT NECESSARY IN ITS ENTIRETY. FOR THE
C                FIRST AND LAST SUBMATRIX WITHIN EACH STRATUM, NBR(3)
C                MUST HAVE ITS CORRECT PATTERN VALUE. FOR THE INTER-
C                MEDIATE SUBMATRICES WITHIN EACH STRATUM, NBR(3) MAY
C                TAKE ANY VALUE BETWEEN THE CORRESPONDING FIRST AND
C                LAST SUBMATRIX NBR(3) VALUES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSRBLK (Y,IY,IOPT,NBR,NH,IN,XBARH,ALPHA,TEMP,BR,HMUSIG,
     *                   IH,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IY,IOPT,NBR(1),IN,NH(IN,1),IH,IER
      REAL               Y(IY,1),XBARH(1),ALPHA,TEMP(1),BR(1),HMUSIG(IH,
     *                   1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,IJ,IL,IS,JJ,L,MM,N,N2,NBR21,NN,NS
      REAL               A,HH,HHH,HNH,HNH2,HNHM,P,PP,S2HX,S2HY,SH2,SHXY,
     *                   T(20),TH
      DOUBLE PRECISION   SUMA,SUMB,SUMC,SUMD,SUME,SUMF,SUMG,SUMH,SUMI,
     *                   SUMJ,WK,WL,WM,WN,WO,WP
      EQUIVALENCE        (T(1),SUMA),(T(3),SUMB),(T(5),SUMC),(T(7),SUMD)
     *                   ,(T(9),SUME),(T(11),SUMF),(T(13),SUMG),(T(15),
     *                   SUMH),(T(17),SUMI),(T(19),SUMJ)
      DATA               SUMA /0.D0/,SUMB /0.D0/,SUMC /0.D0/,SUMD /0.D0/
      DATA               SUME /0.D0/,SUMF /0.D0/,SUMG /0.D0/,SUMH /0.D0/
      DATA               SUMI /0.D0/,SUMJ /0.D0/
C                                  FIRST EXECUTABLE STATEMENT
      IL = 0
      L = NBR(1)
      IF (NBR(3).GT.1) GO TO 35
C                                  TEST FOR INPUT ERRORS
      IF (NBR(3).LT.1) GO TO 180
      IF (ALPHA.LE.0.0 .OR. ALPHA.GE.1.0) GO TO 190
      IF (NBR(1).LT.2 .OR. NBR(2).LT.1) GO TO 175
      N = 0
      N2 = 0
      DO 5 I=1,L
         N2 = N2+NH(I,2)
         N = N+NH(I,1)
    5 CONTINUE
      NBR(7) = N2
      NBR(9) = N
      DO 10 I=1,L
         IF (NH(I,1).LT.3) GO TO 185
         IF (NH(I,1).LT.NBR(2) .AND. NBR(2).LT.N) GO TO 195
         IF (NH(I,2).LT.NH(I,1)) GO TO 185
   10 CONTINUE
      NBR(6) = 1
      NBR(8) = 1
      DO 15 I=11,20
   15 T(I) = 0.0
      IER = 0
      MM = 1
      IF (NBR(2)-N) 20, 30, 195
   20 MM = 0
      NBR21 = NBR(2)-1
      DO 25 I=1,L
         NH(I,3) = (NH(I,1)+NBR21)/NBR(2)
         MM = MM+NH(I,3)
         IF (I.GT.1) NH(I,3) = NH(I,3)+NH(I-1,3)
   25 CONTINUE
      NBR(6) = MM
   30 IF (NBR(4)) 40, 45, 40
C                                  ALL ENTRIES EXCEPT THE FIRST
C                                  BEGIN HERE
   35 MM = NBR(6)
      IF (NBR(2)*(NBR(3)-NBR(1)).GT.NBR(9)) GO TO 180
      IF (NBR(4).EQ.0) GO TO 45
      II = NBR(8)
      IF (NBR(3).NE.(NH(II,3)+1)) GO TO 45
   40 TEMP(1) = Y(1,1)
      TEMP(2) = Y(1,2)
   45 IF (MM.NE.1) GO TO 60
   50 IF (IL.EQ.L) GO TO 120
C                                  PROCESS STRATUM IL
      IL = IL+1
      NN = NH(IL,1)
      IF (IL.EQ.1) GO TO 55
      II = JJ+1
      JJ = JJ+NN
      IF (NBR(4).EQ.0) GO TO 65
      TEMP(1) = Y(II,1)
      TEMP(2) = Y(II,2)
      GO TO 65
   55 II = 1
      JJ = NN
      GO TO 65
   60 IS = NBR(8)
      IF (NBR(3).GT.NH(IS,3)) IS = IS+1
      NBR(8) = IS
      NS = NH(IS,3)-1
      IF (IS.NE.1) NS = NS-NH(IS-1,3)
      NN = NBR(2)
      IF (NBR(3).EQ.NH(IS,3)) NN = NH(IS,1)-NS*NN
      II = 1
      JJ = NN
      IF (NBR(3).EQ.1) GO TO 65
      IF (NBR(8).EQ.1) GO TO 75
      IF (NBR(3).NE.NH(IS-1,3)+1) GO TO 75
      IF (NBR(4).NE.0) GO TO 65
      TEMP(1) = Y(1,1)
      TEMP(2) = Y(1,2)
   65 DO 70 I=1,10
   70 STAT(I) = 0.0
   75 DO 80 I=1,10
         T(I) = STAT(I)
   80 CONTINUE
C                                  COMPUTE TOTALS AND TEMPORARY SUMS
C                                  OF SQUARES AND CROSS PRODUCTS
C                                  FOR EACH STRATUM.
      DO 85 I=II,JJ
         SUMA = SUMA+DBLE(Y(I,1))
         SUMB = SUMB+DBLE(Y(I,2))
         WK = DBLE(Y(I,1))-DBLE(TEMP(1))
         WL = DBLE(Y(I,2))-DBLE(TEMP(2))
         SUMC = SUMC+WK**2
         SUMD = SUMD+WL**2
         SUME = SUME+WK*WL
   85 CONTINUE
      IF (IL.NE.0) GO TO 95
      DO 90 I=1,20
         STAT(I) = T(I)
   90 CONTINUE
      IF (NBR(3).NE.NH(IS,3)) GO TO 9005
C                                  STRATUM IS COMPLETE.  NOW ADJUST
C                                  STATISTICS FOR THIS STRATUM.
   95 IF (IL.NE.0) IS = IL
      WL = NH(IS,1)
      HMUSIG(IS,1) = SUMA/WL
      HMUSIG(IS,4) = SUMB/WL
      HNH = NH(IS,1)
      HNHM = HNH-1.
      WM = WL-1.D0
      HMUSIG(IS,2) = (SUMC-(WL*(HMUSIG(IS,1)-TEMP(1))**2))/WM
      HMUSIG(IS,5) = (SUMD-(WL*(HMUSIG(IS,4)-TEMP(2))**2))/WM
      A = HMUSIG(IS,1)
      IF (NBR(5).EQ.0) A = XBARH(IS)
      HMUSIG(IS,3) = 100.*SQRT(HMUSIG(IS,2)/HNH)/A
      HMUSIG(IS,6) = 100.*SQRT(HMUSIG(IS,5)/HNH)/HMUSIG(IS,4)
      HNH2 = NH(IS,2)
      WL = HNH2*(HNH2-HNH)/(HNH*HNHM)
      SHXY = SUME-HNH*(HMUSIG(IS,1)-TEMP(1))*(HMUSIG(IS,4)-TEMP(2))
      S2HX = HMUSIG(IS,2)*HNHM
      S2HY = HMUSIG(IS,5)*HNHM
      IF (IOPT.LE.0) GO TO 105
C                                  REGRESSION ESTIMATION WITH
C                                  COEFFICIENTS ESTIMATED
      IF (NBR(5).EQ.0) GO TO 100
C                                  ---- COMBINED
      SUMF = SUMF+WL*S2HY
      SUMG = SUMG+WL*SHXY
      SUMH = SUMH+WL*S2HX
      HHH = HNH2/NBR(7)
      WL = HHH*HHH*(1.-(HNH/HNH2))/(HNH*HNHM)
      SUMI = SUMI+SHXY*WL
      SUMJ = SUMJ+S2HX*WL
      GO TO 115
C                                  ---- SEPARATE
  100 BR(IS) = SHXY/S2HX
      WL = HNH2*(HNH2-HNH)/(HNH*(HNH-2.))
      SUMF = SUMF+WL*(S2HY-(BR(IS)*BR(IS)*S2HX))
      GO TO 115
C                                  REGRESSION ESTIMATION WITH
C                                  COEFFICIENTS PREASSIGNED
C                                  OR RATIO ESTIMATION
  105 IF (NBR(5).EQ.0) GO TO 110
C                                  ---- COMBINED
      SUMF = SUMF+WL*S2HY
      SUMG = SUMG+WL*SHXY
      SUMH = SUMH+WL*S2HX
      GO TO 115
C                                  ---- SEPARATE
  110 IF (IOPT.EQ.0) BR(IS) = HMUSIG(IS,4)/HMUSIG(IS,1)
      SUMF = SUMF+WL*(S2HY+BR(IS)*(BR(IS)*S2HX-2.0*SHXY))
  115 IF (MM.EQ.1) GO TO 50
      IF (IS.NE.L) GO TO 9005
C                                  PROCESSING OF ALL STRATA
C                                  COMPLETE
  120 SH2 = NBR(7)
      IF (NBR(5).EQ.0) GO TO 130
      WN = 0.D0
      WO = 0.D0
      DO 125 I=1,L
         WP = NH(I,2)
         WN = WN+WP*HMUSIG(I,1)
         WO = WO+WP*HMUSIG(I,4)
  125 CONTINUE
      STAT(10) = WN/SH2
      STAT(11) = WO/SH2
  130 IF (IOPT.EQ.0) GO TO 145
C                                  REGRESSION ESTIMATION
      IF (NBR(5).EQ.0) GO TO 135
C                                  ---- COMBINED
      IF (IOPT.GT.0) BR(1) = SUMI/SUMJ
      STAT(4) = SUMF-BR(1)*(SUMG+SUMG-BR(1)*SUMH)
      STAT(1) = STAT(11)+BR(1)*(XBARH(1)-STAT(10))
      GO TO 160
C                                  ---- SEPARATE
  135 STAT(4) = SUMF
      WO = 0.D0
      DO 140 I=1,L
         HNH = NH(I,2)
         WO = WO+(HNH/SH2)*(HMUSIG(I,4)+BR(I)*(XBARH(I)-HMUSIG(I,1)))
  140 CONTINUE
      STAT(1) = WO
      GO TO 160
C                                  RATIO ESTIMATION
  145 IF (NBR(5).EQ.0) GO TO 150
C                                  ---- COMBINED
      BR(1) = STAT(11)/STAT(10)
      STAT(4) = SUMF+BR(1)*(BR(1)*SUMH-SUMG-SUMG)
      STAT(1) = BR(1)*XBARH(1)
      GO TO 160
C                                  ---- SEPARATE
  150 STAT(4) = SUMF
      WO = 0.D0
      DO 155 I=1,L
         WO = WO+BR(I)*XBARH(I)*NH(I,2)/SH2
  155 CONTINUE
      STAT(1) = WO
C                                  NOW COMPUTE REMAINING STATS
  160 STAT(2) = STAT(1)*SH2
      STAT(3) = STAT(4)/(SH2*SH2)
      P = 1.-(.5*ALPHA)
      CALL MDNRIS(P,TH,IER)
      P = AMAX1(0.0,STAT(3))
      P = SQRT(P)
      PP = P*TH
      STAT(5) = STAT(1)-PP
      STAT(6) = STAT(1)+PP
      P = AMAX1(0.0,STAT(4))
      P = SQRT(P)
      PP = P*TH
      STAT(7) = STAT(2)-PP
      STAT(8) = STAT(2)+PP
      STAT(9) = (100.*P)/STAT(2)
      HNH = NBR(1)
      HNH = SQRT(HNH)
      SUMA = 0.D0
      DO 165 I=1,L
         IF (HMUSIG(I,3)*HNH.GT.30.) GO TO 170
         HHH = NH(I,2)
         SUMA = SUMA+HMUSIG(I,2)*HHH*HHH/NH(I,1)
  165 CONTINUE
      IF (NBR(5).EQ.0) GO TO 9005
      HHH = SUMA/(SH2*SH2)
      HHH = AMAX1(0.0,HHH)
      IF (((100.*SQRT(HHH))/STAT(10)).LE.10.) GO TO 9005
  170 IER = 33
      GO TO 9000
  175 IER = 130
      GO TO 9000
  180 IER = 131
      GO TO 9000
  185 IER = 132
      GO TO 9000
  190 IER = 133
      GO TO 9000
  195 IER = 134
 9000 CONTINUE
      CALL UERTST(IER,6HSSRBLK)
 9005 RETURN
      END

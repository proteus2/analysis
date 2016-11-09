C   IMSL ROUTINE NAME   - BECOVW
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - MEANS AND VARIANCE-COVARIANCE OR CORRELATION
C                           MATRIX FROM DATA POSSIBLY CONTAINING MISSING
C                           OBSERVATIONS, WITH WEIGHTING ON OPTION.
C
C   USAGE               - CALL BECOVW (X,IX,WT,NBR,XMISS,XM,VCV,INCD,WK,
C                           IER)
C
C   ARGUMENTS    X      - INPUT. X IS AN NBR(3) BY NBR(1) SUBMATRIX
C                           OF THE MATRIX (CALL IT XX) OF DATA FOR
C                           WHICH MEANS, VARIANCES AND COVARIANCES,
C                           OR CORRECTED SUMS OF SQUARES AND
C                           CROSSPRODUCTS ARE DESIRED. THE LAST
C                           SUBMATRIX IN XX MAY HAVE FEWER THAN NBR(3)
C                           ROWS. SEE EXAMPLE.
C                IX     - INPUT, ROW DIMENSION OF X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                WT     - INPUT VECTOR OF LENGTH NBR(3) IF WEIGHTING IS
C                           DESIRED.  IN THIS CASE WT CONTAINS THE
C                           NONNEGATIVE WEIGHTS CORRESPONDING TO THE
C                           ROWS OF X.  IF NO WEIGHTING (OR WEIGHTS ALL
C                           1) IS DESIRED, WT IS OF LENGTH 1 AND WT(1)
C                           MAY BE SET TO ANY NEGATIVE QUANTITY ON THE
C                           FIRST CALL TO BECOVW FOR A GIVEN PROBLEM,
C                           AND THEN OTHER ELEMENTS OF WT ARE NOT USED.
C                           IF WT(1) IS NONNEGATIVE, ALL THE ELEMENTS
C                           OF WT MUST BE NONNEGATIVE.
C                NBR    - INPUT VECTOR OF LENGTH 6. NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF VARIABLES.
C                           I=2, NUMBER OF OBSERVATIONS PER VARIABLE
C                             IN XX.
C                           I=3, NUMBER OF OBSERVATIONS PER VARIABLE IN
C                             EACH SUBMATRIX X, NOT INCLUDING THE LAST
C                             SUBMATRIX WHERE THE NUMBER MAY BE LESS
C                             THAN OR EQUAL TO NBR(3). HOWEVER, NBR(3)
C                             SHOULD BE THE SAME FOR ALL CALLS.
C                           I=4, THE NUMBER OF THE SUBMATRIX STORED IN
C                             X. SEE EXAMPLE.
C                             NBR(4) CAN ALSO BE USED AS A NO-DATA
C                             INDICATOR.  IF NBR(4) IS NEGATIVE, NO DATA
C                             FROM X ARE USED AND THE ONLY COMPUTATIONS
C                             THAT ARE PERFORMED ARE TO COMPUTE THE
C                             COVARIANCE OR CORRELATION MATRIX AS
C                             SPECIFIED BY NBR(6) FROM THE CORRECTED
C                             SUMS OF SQUARES AND CROSSPRODUCTS MATRIX
C                             WHICH IS INPUT IN VCV. IF THE CORRELATION
C                             MATRIX IS DESIRED BUT VCV CONTAINS THE
C                             VARIANCE-COVARIANCE MATRIX, THE
C                             CORRELATION MATRIX CAN BE OBTAINED BY
C                             SETTING NBR(4) TO A NEGATIVE QUANTITY
C                             AND SETTING NBR(5) TO ZERO REGARDLESS
C                             OF HOW MISSING VALUES HAD BEEN HANDLED
C                             IN THE ORIGINAL COMPUTATIONS. NBR(6) MUST
C                             ALSO BE SET APPROPRIATELY.  ALSO, THE WK
C                             VECTOR MUST BE CONTAIN THE VALUES FROM
C                             THE LAST CALL TO BECOVW.
C                           I=5, THE MISSING VALUE OPTION.
C                             IF NBR(5) = 0, ALL DATA IN X ARE ASSUMED
C                               TO BE VALID.
C                             OTHERWISE, THE VALUES IN XMISS ARE
C                               INTERPRETED AS MISSING VALUE CODES AND
C                               ANY X VALUE EQUAL TO THE CORRESPONDING
C                               VALUE IN XMISS IS EXCLUDED FROM THE
C                               COMPUTATIONS.
C                               IF NBR(5) IS NEGATIVE, THE EXCLUSION IS
C                                 LISTWISE.  (THE ENTIRE ROW OF X IS
C                                 EXCLUDED IF ANY OF THE VALUES OF THE
C                                 ROW IS EQUAL TO THE CORRESPONDING
C                                 VALUE IN XMISS.)
C                               IF NBR(5) = 1, RAW CROSS-PRODUCTS ARE
C                                 COMPUTED FROM ALL VALID PAIRS, AND
C                                 MEANS AND VARIANCES ARE COMPUTED
C                                 FROM ALL VALID DATA ON THE INDIVIDUAL
C                                 VARIABLES.  CORRECTED CROSS-PRODUCTS,
C                                 COVARIANCES AND CORRELATIONS ARE THEN
C                                 COMPUTED USING THESE QUANTITIES.
C                               IF NBR(5) = 2, RAW CROSS-PRODUCTS, MEANS
C                                 AND VARIANCES ARE COMPUTED AS IN THE
C                                 CASE OF NBR(5) = 1. HOWEVER, CORRECTED
C                                 CROSS-PRODUCTS AND COVARIANCES ARE
C                                 COMPUTED USING MEANS THAT ARE COMPUTED
C                                 ONLY FROM THE VALID PAIRS OF DATA.
C                                 CORRELATIONS ARE THEN COMPUTED
C                                 USING THESE COVARIANCES AND THE
C                                 VARIANCES FROM ALL VALID DATA ON
C                                 THE INDIVIDUAL OBSERVATIONS.
C                               IF NBR(5) = 3, RAW CROSS-PRODUCTS,
C                                 MEANS, VARIANCES AND COVARIANCES ARE
C                                 COMPUTED AS IN THE CASE OF NBR(5) = 2.
C                                 HOWEVER, CORRELATIONS ARE COMPUTED
C                                 USING THESE COVARIANCES, BUT USING
C                                 VARIANCES THAT ARE COMPUTED ONLY FROM
C                                 THE VALID PAIRS OF DATA.
C                           I=6, THE VCV OPTION.
C                             IF NBR(6) = 0, VCV CONTAINS THE VARIANCE-
C                               COVARIANCE MATRIX.
C                             IF NBR(6) IS POSITIVE, VCV CONTAINS THE
C                               CORRECTED SUMS OF SQUARES AND CROSS-
C                               PRODUCTS MATRIX.
C                             IF NBR(6) IS NEGATIVE, VCV CONTAINS THE
C                               CORRELATION MATRIX. IF NBR(6) = -1, THE
C                               DIAGONAL ELEMENTS OF VCV CONTAINS THE
C                               STANDARD DEVIATIONS, INSTEAD OF THE
C                               CORRELATIONS (WHICH ARE ALL 1.0).
C                XMISS  - INPUT VECTOR OF LENGTH NBR(1) REQUIRED IF
C                           NBR(5) .NE. 0. IN THIS CASE, ANY VALUE IN
C                           THE DATA MATRIX X EQUAL TO THE CORRESPONDING
C                           VALUE IN XMISS IS ASSUMED TO BE INVALID OR
C                           MISSING. THE INVALID OR MISSING DATA IS
C                           HANDLED BY EITHER PAIRWISE OR LISTWISE
C                           DELETION AS SPECIFIED IN NBR(5).
C                XM     - INPUT/OUTPUT VECTOR OF LENGTH NBR(1)
C                           CONTAINING THE VARIABLE MEANS. IF WT(1) IS
C                           NONNEGATIVE, THESE ARE WEIGHTED MEANS.
C                           IF NBR(4) IS NOT 1, ON INPUT XM CONTAINS
C                           MEANS FROM PREVIOUS CALLS TO BECOVW.
C                VCV    - INPUT/OUTPUT NBR(1) BY NBR(1) MATRIX STORED
C                           IN SYMMETRIC STORAGE MODE REQUIRING
C                           (NBR(1)*(NBR(1)+1))/2 STORAGE LOCATIONS.
C                           VCV CONTAINS EITHER THE CORRELATION MATRIX
C                           (POSSIBLY WITH THE STANDARD DEVIATIONS ON
C                           THE DIAGONALS), THE VARIANCE-COVARIANCE
C                           MATRIX, OR THE CORRECTED SUMS OF SQUARES AND
C                           CROSSPRODUCTS MATRIX, AS CONTROLLED BY THE
C                           VCV OPTION, NBR(6). IF WT(1) IS NONNEGATIVE,
C                           THESE QUANTITIES ARE WEIGHTED.
C                           IF NBR(4) IS NOT 1, ON INPUT VCV CONTAINS
C                           CORRECTED SUMS OF SQUARES AND CROSSPRODUCTS
C                           FROM PREVIOUS CALLS TO BECOVW.
C                INCD   - INPUT/OUTPUT INCIDENCE MATRIX USED ONLY IF
C                           NBR(5) IS NOT ZERO.
C                           IF NBR(5) IS NEGATIVE, IMPLYING LISTWISE
C                             DELETION OF MISSING VALUES, INCD IS
C                             OF LENGTH 1 AND ITS ELEMENT IS THE
C                             NUMBER OF VALID OBSERVATIONS.
C                           IF NBR(5) IS POSITIVE, INCD IS OF LENGTH
C                             (NBR(1)*(NBR(1)+1))/2+1.  THE FIRST
C                             PART OF INCD IS AN NBR(1) BY NBR(1)
C                             MATRIX IN SYMMETRIC STORAGE MODE
C                             CONTAINING THE NUMBERS OF PAIRS OF VALID
C                             OBSERVATIONS WHICH WERE USED IN
C                             CALCULATING THE VALUES IN VCV, AND THE
C                             LAST ELEMENT OF INCD IS THE TOTAL NUMBER
C                             OF OBSERVATIONS THAT CONTAINED ANY
C                             MISSING VALUES.
C                           IF NBR(5) IS ZERO, INCD
C                             MAY BE DIMENSIONED WITH LENGTH 1.
C                           IF NBR(4) IS NOT 1, ON INPUT INCD CONTAINS
C                             COUNTS FROM PREVIOUS CALLS TO BECOVW.
C                WK     - INPUT/OUTPUT WORK VECTOR (MUST NOT BE
C                           CHANGED BETWEEN CALLS TO BECOVW).
C                           IF NBR(5) IS LESS THAN 1,
C                             IF WT(1) IS NEGATIVE,
C                               WK IS OF LENGTH 3*NBR(1).
C                             IF WT(1) IS NONNEGATIVE,
C                               WK IS OF LENGTH 3*NBR(1)+1.
C                           IF NBR(5) EQUALS 1 OR 2,
C                             IF WT(1) IS NEGATIVE,
C                               WK IS OF LENGTH NBR(1)*(NBR(1)+2).
C                             IF WT(1) IS NONNEGATIVE,
C                               WK IS OF LENGTH (NBR(1)*(3*NBR(1)+5))/2.
C                           IF NBR(5) EQUALS 3,
C                             IF WT(1) IS NEGATIVE,
C                               WK IS OF LENGTH 2*NBR(1)*(NBR(1)+1).
C                             IF WT(1) IS NONNEGATIVE,
C                               WK IS OF LENGTH 5*NBR(1)*(NBR(1)+1)/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT CORRELATIONS WERE
C                             REQUESTED BUT THAT THE OBSERVATIONS ON
C                             SOME VARIABLE WERE CONSTANT. THE PERTINENT
C                             CORRELATION COEFFICIENTS ARE SET TO
C                             NEGATIVE MACHINE INFINITY.
C                           IER=34 INDICATES THAT VARIANCES AND
C                             COVARIANCES WERE REQUESTED BUT THAT
C                             FEWER THAN TWO VALID OBSERVATIONS
C                             WERE PRESENT FOR SOME VARIABLE. THE
C                             PERTINENT STATISTICS ARE SET TO NEGATIVE
C                             MACHINE INFINITY.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NBR(3)*(NBR(4)-1)
C                             EXCEEDS NBR(2) (WHEN NBR(4) IS POSITIVE).
C                           IER=130 INDICATES THAT NBR(1) IS LESS THAN
C                             1 OR NBR(2) IS LESS THAN 2 OR NBR(3)
C                             EXCEEDS NBR(2).
C                           IER=131 INDICATES THAT NBR(4)=1 AND WT(1) IS
C                             NONNEGATIVE (INPLYING WEIGHTS ARE TO BE
C                             USED) AND THAT SOME WEIGHTS ARE NEGATIVE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE WORKSPACE MAY CONTAIN STATISTICS OF INTEREST. LET
C                           M = NBR(1)
C                           M2 = 2*M
C                           M3 = 3*M     AND
C                           MM1 = M*(M+1)/2.
C                THE WORKSPACE UTILIZATION IS:
C
C                             START  LENGTH         CONTENTS
C        ALL OPTIONS            1      M    INDICATORS OF CONSTANT DATA
C                              M+1     M    FIRST NONMISSING DATA
C
C        IF NBR(5).LE.0       M2+1     M    DEVIATION FROM TEMP. MEAN
C          AND WEIGHTING      M3+1     1    SUM OF WEIGHTS
C
C        IF NBR(5).EQ.1 OR 2  M2+1    M**2  PAIRWISE MEANS
C          AND WEIGHTING    M2+M**2+1 MM1   PAIRWISE SUMS OF WEIGHTS
C
C        IF NBR(5).EQ.3       M2+1    M**2  PAIRWISE MEANS
C          AND NO WEIGHTING M2+M**2+1 M**2  PAIRWISE SUMS OF PRODUCTS
C          AND WEIGHTING    M2+M**2+1 MM1   PAIRWISE SUMS OF WEIGHTS
C                       M2+MM1+M**2+1 M**2  PAIRWISE SUMS OF PRODUCTS
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BECOVW (X,IX,WT,NBR,XMISS,XM,VCV,INCD,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,NBR(1),INCD(1),IER
      REAL               X(IX,1),WT(1),XMISS(1),XM(1),VCV(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IMISS,IMEAN,IPAIRM,IPAIRV,IPAIRW,J,JMEAN,K,
     *                   KMEAN,L,L1,L2,M1,M2,MM1,MM1P1,MSQ,NN,NNN
      REAL               DENOM, CNT,CNTM1,RAT,SUMWT,TJ,TJK,TKJ,VCVL,WI,
     *                   WTT,XINFM,XINF
      DATA               WTT,CNT/2*0.0/
      DATA               XINFM /ZFFFFFFFF/
      DATA               XINF /Z7FFFFFFF/
C
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NBR(1).GE.1 .AND. NBR(2).GE.2 .AND. NBR(3).LE.NBR(2)) GO TO 5
      IER = 130
      GO TO 9000
    5 IF (NBR(4).LT.0) GO TO 15
      IF (NBR(4).EQ.0) GO TO 10
      IF (NBR(3)*(NBR(4)-1).LE.NBR(2)) GO TO 15
C                                  TERMINAL ERROR - INVALID NBR VALUES
   10 IER = 129
      GO TO 9000
   15 NNN = (NBR(2)+NBR(3)-1)/NBR(3)
      M1 = NBR(1)
      M2 = 2*M1
      MSQ = M1*M1
      MM1 = (MSQ+M1)/2
      MM1P1 = MM1+1
      IPAIRM = M2
      IPAIRW = M2+MSQ
      IPAIRV = M2+MSQ
      IF (NBR(5).EQ.3 .AND. WT(1).GT.0.0) IPAIRV = IPAIRV+MM1
      IF (NBR(4).LT.0) GO TO 175
      IF (NBR(4).NE.1) GO TO 35
C                                  INITIALIZE XM,VCV,CNT,INCD,WK
      IF (NBR(5) .NE. 0) INCD(1) = 0
      CNT = 0.0
      WTT = WT(1)
      IF (NBR(5).LT.1 .AND. WTT.GE.0.0) WK(3*M1+1) = 0.0
      L = 0
      DO 25 I=1,M1
         WK(I) = 0.0
         WK(M1+I) = X(1,I)
         XM(I) = 0.0
         IF (NBR(5).LT.1) GO TO 25
         DO 20 J=1,M1
            L = L+1
            WK(IPAIRM+L) = 0.0
            IF (NBR(5).EQ.3) WK(IPAIRV+L) = 0.0
   20    CONTINUE
   25 CONTINUE
      DO 30 I=1,MM1
         IF (NBR(5).GT.0) INCD(I) = 0
         IF (NBR(5).GT.0 .AND. WTT.GT.0.0) WK(IPAIRW+I) = 0.0
         VCV(I) = 0.0
   30 CONTINUE
      IF (NBR(5).GT.0) INCD(MM1P1) = 0
C                                  ALL ENTRIES EXCEPT FIRST BEGIN HERE
   35 NN = NBR(3)
      IF (NBR(4).EQ.NNN) NN = NBR(2)-(NNN-1)*NBR(3)
      IF (NBR(5).EQ.0) GO TO 50
      DO 45 I=1,NN
         DO 40 J=1,M1
            IF (WK(M1+J).EQ.XMISS(J)) WK(M1+J) = X(I,J)
   40    CONTINUE
   45 CONTINUE
   50 IF (WTT.GE.0.0) GO TO 105
C                                  ACCUMULATE MEANS AND ADJUSTED SUMS
C                                    OF SQUARES AND CROSS PRODUCTS.
C                                    TEST FOR CONSTANT OBSERVATION SET.
C
C                                  NO WEIGHTING
C
      IF (NBR(5).GE.1) GO TO 80
C                                  NO MISSING VALUES OR ELSE LISTWISE
C                                    DELETION OPTION
      DO 75 I=1,NN
         IF (NBR(5).EQ.0) GO TO 60
C                                  CHECK TO SEE IF ANY MISSING VALUES
         DO 55 J=1,M1
            IF (X(I,J).EQ.XMISS(J)) GO TO 75
   55    CONTINUE
   60    CNT = CNT+1.0
         L = 0
         DO 70 J=1,M1
            IF (X(I,J).NE.WK(M1+J)) WK(J) = 1.0
            WK(M2+J) = X(I,J)-XM(J)
            TJ = WK(M2+J)
            DO 65 K=1,J
               L = L+1
               TJK = TJ*WK(M2+K)
               VCV(L) = VCV(L)+TJK-(TJK/CNT)
   65       CONTINUE
            XM(J) = XM(J)+(TJ/CNT)
   70    CONTINUE
   75 CONTINUE
      IF (NBR(5).LT.0) INCD(1) = CNT
      GO TO 170
   80 CONTINUE
C                                  PAIRWISE DELETION OPTION
      DO 100 I=1,NN
         L = 0
         IMISS = 0
         JMEAN = IPAIRM-M1
         JVAR = IPAIRV-M1
         DO 95 J=1,M1
            JMEAN = JMEAN+M1
            JVAR = JVAR+M1
            IF (X(I,J).NE.XMISS(J)) GO TO 85
            IMISS = 1
            L = L+J
            GO TO 95
   85       IF (X(I,J).NE.WK(M1+J)) WK(J) = 1.0
            KMEAN = IPAIRM-M1
            KVAR = IPAIRV-M1
            DO 90 K=1,J
               L = L+1
               KMEAN = KMEAN+M1
               KVAR = KVAR+M1
               IF (X(I,K).EQ.XMISS(K)) GO TO 90
               INCD(L) = INCD(L)+1
               TJK = X(I,J)-WK(JMEAN+K)
               TKJ = X(I,K)-WK(KMEAN+J)
               VCV(L) = VCV(L)+TKJ*TJK-(TKJ*TJK)/INCD(L)
               WK(JMEAN+K) = WK(JMEAN+K)+TJK/INCD(L)
               IF (NBR(5).EQ.3) WK(JVAR+K) =
     *                    WK(JVAR+K)+TJK*TJK-(TJK*TJK)/INCD(L)
               IF (J.EQ.K) GO TO 90
               WK(KMEAN+J) = WK(KMEAN+J)+TKJ/INCD(L)
               IF (NBR(5).EQ.3) WK(KVAR+J) =
     *                    WK(KVAR+J)+TKJ*TKJ-(TKJ*TKJ)/INCD(L)
   90       CONTINUE
            XM(J) = XM(J)+(X(I,J)-XM(J))/INCD(L)
   95    CONTINUE
         IF (IMISS.EQ.1) INCD(MM1P1) = INCD(MM1P1)+1
  100 CONTINUE
      GO TO 170
C                                  WEIGHTING
  105 CONTINUE
      IF (NBR(5).GT.0) GO TO 140
C                                  NO MISSING VALUES OR ELSE LISTWISE
C                                    DELETION OPTION
      DO 135 I=1,NN
         IF (NBR(5).EQ.0) GO TO 115
         DO 110 J=1,M1
            IF (X(I,J).EQ.XMISS(J)) GO TO 135
  110    CONTINUE
  115    WI = WT(I)
         IF (WI.GE.0.0) GO TO 120
         IER = 131
         GO TO 9000
  120    CNT = CNT+1.0
         IF (WI.EQ.0.0) GO TO 135
         WK(3*M1+1) = WK(3*M1+1)+WI
         RAT = WI/WK(3*M1+1)
         L = 0
         DO 130 J=1,M1
            IF (X(I,J).NE.WK(M1+J)) WK(J) = 1.0
            WK(M2+J) = X(I,J)-XM(J)
            TJ = WK(M2+J)
            DO 125 K=1,J
               L = L+1
               TJK = TJ*WK(M2+K)*WI
               VCV(L) = VCV(L)+TJK-(TJK*RAT)
  125       CONTINUE
            XM(J) = XM(J)+(TJ*RAT)
  130    CONTINUE
  135 CONTINUE
      IF (NBR(5).LT.0) INCD(1) = CNT
      GO TO 170
  140 CONTINUE
C                                  PAIRWISE DELETION OPTION
      DO 165 I=1,NN
         WI = WT(I)
         IF (WI.GE.0.0) GO TO 145
         IER = 131
         GO TO 9000
  145    L = 0
         IMISS = 0
         JMEAN = IPAIRM-M1
         JVAR = IPAIRV-M1
         DO 160 J=1,M1
            JMEAN = JMEAN+M1
            JVAR = JVAR+M1
            IF (X(I,J).NE.XMISS(J)) GO TO 150
            IMISS = 1
            L = L+J
            GO TO 160
  150       IF (X(I,J).NE.WK(M1+J)) WK(J) = 1.0
            KMEAN = IPAIRM-M1
            KVAR = IPAIRV-M1
            DO 155 K=1,J
               L = L+1
               KMEAN = KMEAN+M1
               KVAR = KVAR+M1
               IF (X(I,K).EQ.XMISS(K)) GO TO 155
               INCD(L) = INCD(L)+1
               IF (WI .EQ. 0.0) GO TO 155
               WK(IPAIRW+L) = WK(IPAIRW+L)+WI
               RAT = WI/WK(IPAIRW+L)
               TJK = X(I,J)-WK(JMEAN+K)
               TKJ = X(I,K)-WK(KMEAN+J)
               VCV(L) = VCV(L)+WI*(TJK*TKJ-(TJK*TKJ*RAT))
               WK(JMEAN+K) = WK(JMEAN+K)+(TJK*RAT)
               IF (NBR(5).EQ.3) WK(JVAR+K) =
     *                    WK(JVAR+K)+WI*(TJK*TJK-(TJK*TJK*RAT))
               IF (J.EQ.K) GO TO 155
               WK(KMEAN+J) = WK(KMEAN+J)+(TKJ*RAT)
               IF (NBR(5).EQ.3) WK(KVAR+J) =
     *                    WK(KVAR+J)+WI*(TKJ*TKJ-(TKJ*TKJ*RAT))
  155       CONTINUE
            IF (WI .EQ. 0.0) GO TO 160
            XM(J) = XM(J)+(X(I,J)-XM(J))*WI/WK(IPAIRW+L)
  160    CONTINUE
         IF (IMISS.EQ.1) INCD(MM1P1) = INCD(MM1P1)+1
  165 CONTINUE
  170 IF (NBR(4).NE.NNN) GO TO 9005
C                                  POSTPROCESSING BEGINS HERE
  175 IF (NBR(5).EQ.0) CNT = NBR(2)
      IF (NBR(5).LT.0) CNT = INCD(1)
      CNTM1 = CNT-1
C                                  CHECK FOR CONSTANT VALUES
      L = 0
      DO 180 I=1,M1
         L = L+I
         IF (WK(I).NE.0.0) GO TO 180
         VCV(L) = 0.0
         IF (NBR(6).LT.0) IER = 33
  180 CONTINUE
      IF (NBR(5).NE.1) GO TO 195
C                                  ADJUST SUM OF SQUARES AND
C                                    CROSS-PRODUCTS MATRIX.
      L = 0
      IMEAN = IPAIRM-M1
      DO 190 I=1,M1
         JMEAN = IPAIRM-M1
         IMEAN = IMEAN+M1
         DO 185 J=1,I
            L = L+1
            JMEAN = JMEAN+M1
            CNT = INCD(L)
            IF (WTT.GE.0.0) CNT = WK(IPAIRW+L)
            TJ = XM(I)*XM(J)+WK(JMEAN+I)*WK(IMEAN+J)-XM(I)*WK(JMEAN+I)
     *      -XM(J)*WK(IMEAN+J)
            VCV(L) = VCV(L)+CNT*TJ
  185    CONTINUE
         IF (VCV(L).LT.0.0) VCV(L) = 0.0
  190 CONTINUE
  195 IF (NBR(6).GT.0) GO TO 9005
      IF (NBR(6).LT.0) GO TO 255
C                                  COMPUTE THE VARIANCE-COVARIANCE
C                                    MATRIX
      IF (NBR(5).GT.0) GO TO 215
      IF (NBR(5).LT.0 .AND. INCD(1).LT.2) GO TO 205
      DO 200 L=1,MM1
  200 VCV(L) = VCV(L)/CNTM1
      IF (IER.EQ.0) GO TO 9005
      GO TO 9000
  205 IER = 34
      DO 210 L=1,MM1
  210 VCV(L) = XINFM
      GO TO 9000
  215 DO 225 L=1,MM1
         IF (INCD(L).GE.2) GO TO 220
         IER = 34
         VCV(L) = XINFM
         GO TO 225
  220    VCV(L) = VCV(L)/(INCD(L)-1.0)
  225 CONTINUE
      IF (NBR(5).NE.3) GO TO 250
C                                  CHANGE THE PAIRWISE SUMS OF PRODUCTS
C                                    TO COVARIANCES
      L = 0
      L2 = -M1
      DO 245 I = 1,M1
         K = I-1
         L1 = I-M1
         L2 = L2-L1+1
         IF (K .LT. 1) GO TO 236
         DO 235 J = 1,K
            L = L+1
            L1 = L1+M1
            L2 = L2+1
            IF (INCD(L) .GE. 2) GO TO 230
            WK(IPAIRV+L1) = 0.0
            WK(IPAIRV+L2) = 0.0
            GO TO 235
  230       WK(IPAIRV+L1) = WK(IPAIRV+L1)/(INCD(L)-1)
            WK(IPAIRV+L2) = WK(IPAIRV+L2)/(INCD(L)-1)
  235    CONTINUE
  236    L = L+1
         L2 = L2+1
         IF (INCD(L) .GE. 2) GO TO 240
         WK(IPAIRV+L2) = 0.0
         GO TO 245
  240    WK(IPAIRV+L2) = WK(IPAIRV+L2)/(INCD(L)-1)
  245 CONTINUE
  250 IF (IER.EQ.0) GO TO 9005
      GO TO 9000
C                                  COMPUTE THE CORRELATION MATRIX
  255 IF (M1.EQ.1) GO TO 295
      L1 = 1
      L2 = 2
      JVAR = IPAIRV
      DO 290 J=2,M1
         JVAR = JVAR+M1
         L1 = L1+J
         L = 0
         K = J-1
         IVAR = IPAIRV-M1
         DO 285 I=1,K
            IVAR = IVAR+M1
            L = L+I
            IF (VCV(L).EQ.0.0 .OR. VCV(L1).EQ.0.0) GO TO 275
            IF (NBR(5).GT.0) GO TO 260
            IF (CNTM1.LT.1) GO TO 275
            VCV(L2) = VCV(L2)/(SQRT(VCV(L1))*SQRT(VCV(L)))
            GO TO 270
  260       IF (INCD(L2).LT.2) GO TO 275
            IF (NBR(5).NE.3) VCV(L2) = (VCV(L2)*
     *         SQRT((INCD(L)-1.0)*(INCD(L1)-1.0)))/
     *         ((INCD(L2)-1.0)*SQRT(VCV(L1))*SQRT(VCV(L)))
            IF (NBR(5).NE.3) GO TO 270
            IF (WK(JVAR+I).LE.0.0 .OR. WK(IVAR+J).LE.0.0) GO TO 275
            DENOM = SQRT(WK(JVAR+I)*WK(IVAR+J))
            IF(DENOM.GT.1.0) GO TO 265
            RAT = XINF*DENOM
            IF(VCV(L2).GT.RAT) GO TO 275
  265       VCV(L2) =  VCV(L2)/DENOM
  270       IF (VCV(L2).GT.1.0) VCV(L2) = 1.0
            IF (VCV(L2).LT.-1.0) VCV(L2) = -1.0
            GO TO 280
  275       VCV(L2) = XINFM
            IF (IER.EQ.0) IER = 33
  280       L2 = L2+1
  285    CONTINUE
         L2 = L2+1
  290 CONTINUE
  295 IF (NBR(6).NE.-1) GO TO 305
C                                  COMPUTE STANDARD DEVIATIONS AND STORE
C                                    IN DIAGONALS
      L = 0
      DO 300 I=1,M1
         L = L+I
         VCVL = VCV(L)
         VCV(L) = XINFM
         IF (NBR(5).LE.0 .AND. CNTM1.GE.1.0) VCV(L) = SQRT(VCVL/CNTM1)
         IF (NBR(5).GT.0 .AND. INCD(L).GT.1) VCV(L) =
     *   SQRT(VCVL/(INCD(L)-1.0))
  300 CONTINUE
      IF (IER.EQ.0) GO TO 9005
      GO TO 9000
C                                  SET DIAGONALS TO 1.0
  305 L = 0
      DO 310 I=1,M1
         L = L+I
         VCV(L) = 1.0
  310 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'BECOVW')
 9005 RETURN
      END

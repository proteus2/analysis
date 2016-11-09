C   IMSL ROUTINE NAME   - SSSBLK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STRATIFIED RANDOM SAMPLING WITH CONTINUOUS
C                           DATA - INFERENCES REGARDING THE POPULATION
C                           MEAN AND TOTAL
C
C   USAGE               - CALL SSSBLK (Y,NBR,NH,IN,ALPHA,TEMP,HMUSIG,
C                           IH,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT SUBVECTOR OF LENGTH NBR(2) OF THE VECTOR
C                           (CALL IT YY) CONTAINING THE ENTIRE STRATI-
C                           FIED RANDOM SAMPLE. THE SUBVECTOR Y MUST BE
C                           EITHER THE SAME AS YY OR IT MUST CONTAIN
C                           ALL OR PART OF THE SAMPLE FOR A SINGLE
C                           STRATUM ONLY.  IN THE CASE WHERE Y IS A
C                           PROPER SUBSET OF YY, THE LAST SUBVECTOR
C                           FROM ANY ONE OR MORE OF THE STRATA MAY HAVE
C                           FEWER THAN NBR(2) ELEMENTS.  IN ALL CASES
C                           THE OBSERVATIONS WITHIN ANY ONE STRATUM
C                           MUST APPEAR CONTIGUOUSLY IN YY.
C                NBR    - INPUT VECTOR OF LENGTH 5.  NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF STRATA INTO WHICH THE SAMPLE
C                             IS DIVIDED.
C                           I=2, NUMBER OF OBSERVATIONS IN EACH SUBVEC-
C                             TOR Y, NOT INCLUDING THE LAST SUBVECTOR
C                             IN EACH STRATUM, WHERE THE NUMBER MAY BE
C                             LESS THAN OR EQUAL TO NBR(2).  HOWEVER,
C                             NBR(2) SHOULD BE THE SAME FOR ALL CALLS.
C                           I=3, THE NUMBER OF THE SUBVECTOR STORED
C                             IN Y. SEE REMARKS.
C                           I=4, THE TEMPORARY MEAN INDICATOR.  IF
C                             NBR(4) = 0, THE USER SUPPLIES THE TEMPOR-
C                             ARY MEAN IN TEMP.  THE USER MAY CHANGE
C                             TEMP ONLY WHEN SUBVECTOR Y IS THE FIRST
C                             SET OF OBSERVATIONS IN A STRATUM.  IF
C                             NBR(4) DOES NOT EQUAL ZERO, THE FIRST
C                             ELEMENT IN EACH STRATUM IS USED AS THE
C                             TEMPORARY MEAN FOR THE DATA IN THAT
C                             STRATUM.
C                           I=5, THE WITHIN STRATUM VARIANCE ASSUMPTION
C                             INDICATOR.
C                             IF NBR(5) = 0, THE TRUE WITHIN STRATUM
C                               VARIANCE IS ASSUMED CONSTANT, AND A
C                               POOLED ESTIMATE OF THAT VARIANCE IS RE-
C                               TURNED IN HMUSIG(1,2).
C                             IF NBR(5) IS NONZERO, SEPARATE WITHIN
C                               STRATA VARIANCE ESTIMATES ARE COMPUTED
C                               AND RETURNED IN COLUMN 2 OF HMUSIG.
C                NH     - INPUT NBR(1) BY 3 MATRIX CONTAINING THE NUM-
C                           BER OF SAMPLE UNITS IN EACH STRATUM FOR THE
C                           SAMPLE (VECTOR YY) AND FOR THE POPULATION.
C                           THE STRATA SIZES FOR THE SAMPLE AND THE
C                           POPULATION ARE IN COLUMNS 1 AND 2, RESPEC-
C                           TIVELY.  THE STRATA SIZES MUST BE ORDERED IN
C                           CORRESPONDENCE WITH THE ORDERING OF STRATA
C                           IN YY.  IN THE CASE WHERE POPULATION STRATA
C                           SIZES ARE NOT KNOWN, ESTIMATES MUST BE
C                           ENTERED IN THEIR PLACE.  COLUMN 3 IS WORK
C                           STORAGE.
C                IN     - INPUT ROW DIMENSION OF MATRIX NH EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVALS FOR THE MEAN AND TOTAL
C                           PARAMETERS.
C                TEMP   - INPUT TEMPORARY MEAN.  REQUIRED ONLY IF
C                           NBR(4) = 0.
C                HMUSIG - OUTPUT NBR(1) BY 2 MATRIX CONTAINING THE
C                           WITHIN STRATA MEAN AND VARIANCE ESTIMATES.
C                           THE MEANS AND VARIANCES ARE IN COLUMNS 1
C                           AND 2, RESPECTIVELY.  THE ESTIMATES ARE OR-
C                           DERED IN CORRESPONDENCE WITH THE ORDERING OF
C                           STRATA IN YY.  IN THE CASE NBR(5) = 0, ONLY
C                           ELEMENT HMUSIG(1,2) OF COLUMN 2 IS DEFINED,
C                           AND IT CONTAINS THE POOLED WITHIN STRATA
C                           VARIANCE ESTIMATE.
C                IH     - INPUT ROW DIMENSION OF MATRIX HMUSIG EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                STAT   - OUTPUT VECTOR OF LENGTH 11.  STAT(I) CONTAINS,
C                         WHEN
C                           I=1, ESTIMATE OF THE MEAN.
C                           I=2, ESTIMATE OF THE TOTAL.
C                           I=3, VARIANCE ESTIMATE OF THE MEAN ESTIMATE.
C                           I=4, VARIANCE ESTIMATE OF THE TOTAL ESTIMATE
C                           I=5, LOWER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=6, UPPER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=7, LOWER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=8, UPPER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=9, ESTIMATE OF THE COEFFICIENT OF VARI-
C                             ATION OF THE MEAN AND TOTAL ESTIMATES.
C                           I=10, NUMBER OF DEGREES OF FREEDOM ASSOCIA-
C                             TED WITH THE VARIANCE ESTIMATES OF THE
C                             MEAN AND TOTAL ESTIMATES.  WHEN NBR(5) IS
C                             NONZERO, STAT(10) CONTAINS AN EFFECTIVE
C                             NUMBER OF DEGREES OF FREEDOM DETERMINED
C                             ACCORDING TO THE SATTERTHWAITE APPROX-
C                             IMATION.
C                           I=11, VARIANCE ESTIMATE OF THE MEAN ESTI-
C                             MATE ASSUMING THAT SAMPLING WAS SIMPLE
C                             RANDOM INSTEAD OF STRATIFIED RANDOM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES NBR(1) IS LESS THAN 2 OR
C                             THAT NBR(2) IS LESS THAN 1.
C                           IER = 130 INDICATES NBR(3) IS LESS THAN 1 OR
C                             THAT NBR(2)*(NBR(3)-NBR(1)) EXCEEDS THE
C                             NUMBER (TOTAL OF COLUMN 1 OF NH) OF
C                             ELEMENTS IN YY.
C                           IER = 131 INDICATES THAT AT LEAST ONE
C                             ELEMENT OF COLUMN 1 OF NH IS LESS THAN 2
C                             OR THAT AT LEAST ONE ELEMENT IN COLUMN 2
C                             IS LESS THAN THE ELEMENT IN THE SAME ROW
C                             IN COLUMN 1.
C                           IER = 132 INDICATES THAT ALPHA IS NOT IN THE
C                             EXCLUSIVE INTERVAL (0,1) OR THAT A TER-
C                             MINAL ERROR OCCURRED IN IMSL ROUTINE
C                             MDSTI.
C                           IER = 133 INDICATES THAT NBR(2) IS GREATER
C                             THAN SOME NH(I,1),I=1,...,NBR(1) WHEN
C                             NBR(2) IS LESS THAN THE TOTAL SAMPLE SIZE,
C                             OR THAT NBR(2) EXCEEDS THE TOTAL SAMPLE
C                             SIZE.
C
C   REQD. IMSL ROUTINES - MDNRIS,MDSTI,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      BETWEEN THE FIRST AND LAST CALL TO SSSBLK, EXCEPT FOR
C                OPTIONALLY CHANGING TEMP BETWEEN STRATA, ONLY NBR(3)
C                MAY BE MODIFIED AND IT SHOULD FOLLOW THE PATTERN
C                1,2,... . THOUGH THIS PATTERN IS THE OBVIOUS ONE TO
C                FOLLOW, IT IS NOT NECESSARY IN ITS ENTIRETY. FOR THE
C                FIRST AND LAST SUBVECTOR WITHIN EACH STRATUM, NBR(3)
C                MUST HAVE ITS CORRECT PATTERN VALUE. FOR THE INTER-
C                MEDIATE SUBVECTORS WITHIN EACH STRATUM, NBR(3) MAY
C                TAKE ANY VALUE BETWEEN THE CORRESPONDING FIRST AND
C                LAST SUBVECTOR NBR(3) VALUES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSSBLK (Y,NBR,NH,IN,ALPHA,TEMP,HMUSIG,IH,STAT,IER)
C
      DIMENSION          Y(1),NBR(1),NH(IN,1),HMUSIG(IH,1),STAT(1),T(4)
      DOUBLE PRECISION   SUMA,SUMB,WK
      EQUIVALENCE        (T(1),SUMA),(T(3),SUMB)
      DATA               ZERO,ONE/0.0,1.0/,HUND/100.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IL = 0
      IF (NBR(1) .GE. 2 .AND. NBR(2) .GE. 1) GO TO 5
C                                  TERMINAL ERROR - INVALID NBR VALUES
      IER = 129
      GO TO 9000
    5 IF (NBR(3) .GE. 1) GO TO 15
   10 IER = 130
      GO TO 9000
   15 L = NBR(1)
      N = 0
      DO 20 I=1,L
         N = N+NH(I,1)
         IF (NH(I,1) .LT. 2) GO TO 25
         IF (NH(I,2) .LT. NH(I,1)) GO TO 25
   20 CONTINUE
      IF (NBR(2)*(NBR(3)-NBR(1)) .GT. N) GO TO 10
      GO TO 30
   25 IER = 131
      GO TO 9000
   30 IF (ALPHA .GT. ZERO .AND. ALPHA .LT. ONE ) GO TO 35
      IER = 132
      GO TO 9000
   35 IF (NBR(2) .GT. N) GO TO 50
      IF (NBR(2) .LT. N) GO TO 40
      MM = 1
      GO TO 55
   40 MM = 0
      NBR2M1 = NBR(2)-1
      DO 45 I=1,L
         IF (NH(I,1) .LT. NBR(2)) GO TO 50
C                                  COMPUTE THE NUMBER OF SUBVECTORS IN
C                                    EACH STRATUM AND COMPUTE THE
C                                    PARTIAL SUMS
         NH(I,3) = (NH(I,1)+NBR2M1)/NBR(2)
         MM = MM+NH(I,3)
         IF (I .EQ. 1) GO TO 45
         NH(I,3) = NH(I-1,3)+NH(I,3)
   45 CONTINUE
      GO TO 55
C                                  TERMINAL ERROR - NBR(2) .GT. NH(I,1)
C                                    FOR I = 1,...,NBR(1)
   50 IER = 133
      GO TO 9000
   55 IF (NBR(3) .NE. 1) GO TO 60
      STAT(11) = 1.0
      IF (NBR(4) .NE. 0) TEMP = Y(1)
   60 IF (MM .EQ. 1) GO TO 85
C                                  COMPUTE THE NUMBER OF THE STRATUM
C                                    FOR THE CURRENT SUBVECTOR AND
C                                    DETERMINE WHETHER THE SUBVECTOR
C                                    IS THE FIRST, LAST, OR AN INTER-
C                                    MEDIATE SUBVECTOR OF THAT STRATUM
      IS = STAT(11)
      IF (NBR(3) .GT. NH(IS,3)) IS = IS+1
      NS = NH(IS,3)-1
      IF (IS .NE. 1) NS = NH(IS,3)-NH(IS-1,3)-1
      NN = NBR(2)
      IF (NBR(3) .EQ. NH(IS,3)) NN = NH(IS,1)-NS*NN
      II = 1
      JJ = NN
      IF (NBR(3) .EQ. 1) GO TO 65
      IF (IS .EQ. 1) GO TO 75
      IF (NBR(3) .NE. NH(IS-1,3)+1) GO TO 75
C                                  MOVE THE FIRST ELEMENT OF Y OF
C                                    CURRENT STRATUM INTO TEMP
      IF (NBR(4) .NE. 0) TEMP = Y(1)
C                                  INITIALIZE STAT VECTOR FOR EACH
C                                    STRATUM
   65 DO 70 I=1,10
         STAT(I) = ZERO
   70 CONTINUE
   75 T(1) = STAT(1)
      T(2) = STAT(2)
      T(3) = STAT(3)
      T(4) = STAT(4)
      DO 80 I=II,JJ
         SUMA = SUMA+DBLE(Y(I))
         WK = DBLE(Y(I))-DBLE(TEMP)
         SUMB = SUMB+WK**2
   80 CONTINUE
      IF (IL .NE. 0) GO TO 95
      STAT(1) = T(1)
      STAT(2) = T(2)
      STAT(3) = T(3)
      STAT(4) = T(4)
      IF (NBR(3) .EQ. NH(IS,3)) GO TO 95
      GO TO 9005
C                                  NBR(2) IS EQUAL TO THE TOTAL SAMPLE
C                                    SIZE CASE
   85 IF (IL .EQ. L) GO TO 100
      IL = IL+1
      NN = NH(IL,1)
      IF (IL .EQ. 1) GO TO 90
      II = JJ+1
      JJ = JJ+NN
      IF (NBR(4) .NE. 0) TEMP = Y(II)
      GO TO 65
   90 II = 1
      JJ = NN
      GO TO 65
   95 IF (IL .NE. 0) IS = IL
C                                  COMPUTE MEAN AND VARIANCE ESTIMATES
      FNN = NH(IS,1)
      SUMA = SUMA/FNN
      WK = DBLE(FNN)*(SUMA-TEMP)**2
      HMUSIG(IS,1) = SUMA
      HMUSIG(IS,2) = SUMB-WK
      STAT(11) = IS
      IF (MM .EQ. 1) GO TO 85
      IF (IS .NE. L) GO TO 9005
  100 IF (NBR(5) .NE. 0) GO TO 110
      SS = ZERO
      DO 105 I=1,L
         SS = SS+HMUSIG(I,2)
  105 CONTINUE
      FNL = N-L
      STAT(10) = FNL
      SS = SS/FNL
      HMUSIG(1,2) = SS
      GO TO 120
  110 DO 115 I=1,L
         RNH = ONE/FLOAT(NH(I,1)-1)
         HMUSIG(I,2) = HMUSIG(I,2)*RNH
  115 CONTINUE
C                                  COMPUTE THE SUM OF NH(I,2),I=1,...,L
  120 NN = 0
      DO 125 I=1,L
         NN = NN+NH(I,2)
  125 CONTINUE
C                                  COMPUTE STAT(1) AND STAT(2)
      RNN = ONE/FLOAT(NN)
      STAT(1) = ZERO
      DO 130 I=1,L
         STAT(1) = STAT(1)+FLOAT(NH(I,2))*HMUSIG(I,1)*RNN
  130 CONTINUE
      STAT(2) = STAT(1)*FLOAT(NN)
      GHSH2 = ZERO
      GHSH4 = ZERO
      QA = ZERO
      QB = ZERO
      QC = ZERO
      QD = ZERO
      QE = ZERO
      QQA = ZERO
      QQB = ZERO
      QQC = ZERO
      STATA = ZERO
      STATB = ZERO
      DO 135 I=1,L
         RNH = ONE/FLOAT(NH(I,1))
         RNH1 = ONE/FLOAT(NH(I,1)-1)
         GH = FLOAT(NH(I,2)*(NH(I,2)-NH(I,1)))*RNH
         STATA = STATA+GH
         T(1) = GH*HMUSIG(I,2)
         STATB = STATB+T(1)
         GHSH2 = GHSH2+T(1)
         GHSH4 = GHSH4+T(1)**2*RNH1
         T(3) = FLOAT(NH(I,2))*RNN
         QA = QA+T(3)
         T(4) = T(3)*RNH
         QB = QB+T(4)
         QC = QC+T(3)*T(4)
         T(2) = T(3)*HMUSIG(I,1)
         QD = QD+T(2)*HMUSIG(I,1)
         QE = QE+T(2)
         TT = T(3)*HMUSIG(I,2)
         QQA = QQA+TT
         TS = TT*RNH
         QQB = QQB+TS
         QQC = QQC+TS*T(3)
  135 CONTINUE
      GHSH2 = GHSH2**2
      QE = QE**2
      STAT(4) = STATB
      IF (NBR(5) .EQ. 0) STAT(4) = STATA*SS
      STAT(3) = STAT(4)*(RNN**2)
      STAT(9) = (HUND*SQRT(STAT(4)))/STAT(2)
      STAT(10) = N-L
      IF (NBR(5) .NE. 0) STAT(10) = GHSH2/GHSH4
C                                  CALL ROUTINE MDSTI
      CALL MDSTI (ALPHA,STAT(10),X,JER)
      IF (JER .EQ. 0) GO TO 140
      IER = 132
      GO TO 9000
  140 T(1) = X*SQRT(STAT(3))
      T(2) = X*SQRT(STAT(4))
      STAT(5) = STAT(1)-T(1)
      STAT(6) = STAT(1)+T(1)
      STAT(7) = STAT(2)-T(2)
      STAT(8) = STAT(2)+T(2)
      T(3) = FLOAT(NN-N)/FLOAT(N*(NN-1))
      T(4) = ONE-ONE/FLOAT(NN)
      IF (NBR(5) .NE. 0) GO TO 145
      Q = SS*(T(4)*QA-QB+QC+QD-QE)
      STAT(11) = T(3)*Q
      GO TO 9005
  145 QQ = T(4)*QQA-QQB+QQC+QD-QE
      STAT(11) = T(3)*QQ
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HSSSBLK)
 9005 RETURN
      END

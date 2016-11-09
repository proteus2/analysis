C   IMSL ROUTINE NAME   - SSSCAN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SINGLE STAGE CLUSTER SAMPLING WITH CONTINUOUS
C                           DATA - INFERENCES REGARDING THE POPULATION
C                           MEAN AND TOTAL
C
C   USAGE               - CALL SSSCAN (Y,IOPT,NBR,MC,IM,SIZE,TSIZE,
C                                      ALPHA,TEMP,CMUSIG,IC,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT SUBVECTOR OF LENGTH NBR(2) OF THE VECTOR
C                           (CALL IT YY) CONTAINING THE ENTIRE CLUSTER
C                           SAMPLE. THE SUBVECTOR Y MUST BE EITHER THE
C                           SAME AS YY OR IT MUST CONTAIN ALL OR PART
C                           OF THE SAMPLE FOR A SINGLE CLUSTER ONLY.
C                           IN THE CASE WHERE Y IS A PROPER SUBSET OF
C                           YY, THE LAST SUBVECTOR FROM ANY ONE OR
C                           MORE OF THE CLUSTERS MAY HAVE FEWER THAN
C                           NBR(2) ELEMENTS. IN ALL CASES THE OBSERVA-
C                           TIONS WITHIN ANY ONE CLUSTER MUST APPEAR
C                           CONTIGUOUSLY IN YY.
C                IOPT   - INPUT ESTIMATION OPTION.
C                         IF IOPT IS NEGATIVE, UNBIASED ESTIMATION IS
C                           UTILIZED.
C                         IF IOPT IS ZERO, RATIO TO SIZE ESTIMATION
C                           IS UTILIZED.
C                         IF IOPT IS POSITIVE, PROBABILITY PROPORTIONAL
C                           TO SIZE ESTIMATION IS UTILIZED.
C                NBR    - INPUT VECTOR OF LENGTH 7. NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF CLUSTERS INTO WHICH THE
C                             SAMPLE IS DIVIDED.
C                           I=2, NUMBER OF OBSERVATIONS IN EACH
C                             SUBVECTOR Y, NOT INCLUDING THE LAST
C                             SUBVECTOR IN EACH CLUSTER, WHERE THE
C                             NUMBER MAY BE LESS THAN OR EQUAL TO
C                             NBR(2). HOWEVER, NBR(2) SHOULD BE THE
C                             SAME FOR ALL CALLS.
C                           I=3, THE NUMBER OF THE SUBVECTOR STORED
C                             IN Y. SEE REMARKS.
C                           I=4, THE TEMPORARY MEAN INDICATOR.
C                             IF NBR(4) = 0, THE USER SUPPLIES THE
C                             TEMPORARY MEAN IN TEMP. THE USER MAY
C                             CHANGE TEMP ONLY WHEN SUBVECTOR Y IS
C                             THE FIRST SET OF OBSERVATIONS IN A
C                             CLUSTER. IF NBR(4) DOES NOT EQUAL ZERO,
C                             THE FIRST ELEMENT IN EACH CLUSTER IS
C                             UTILIZED.
C                           I=5, NUMBER OF CLUSTERS IN THE SAMPLED
C                             POPULATION.
C                           I=6, NUMBER OF ELEMENTS IN THE POPULATION
C                             (SUM OF ALL THE CLUSTER SIZES IN THE
C                             POPULATION). NOT REQUIRED IF NBR(7)
C                             IS NONZERO.
C                           I=7, OPTION FOR PROBABILITY PROPORTIONAL TO
C                             SIZE ESTIMATION. REQUIRED ONLY WHEN IOPT
C                             IS POSITIVE.
C                           IF NBR(7) = 0, ALL CLUSTERS IN THE POPU-
C                             LATION ARE OF KNOWN SIZE.
C                           IF NBR(7) IS NONZERO, THE CLUSTER SIZES
C                             ARE KNOWN ONLY APPROXIMATELY OR A MEASURE
C                             OF CLUSTER SIZE OTHER THAN THE NUMBER OF
C                             ELEMENTS PER CLUSTER IS TO BE USED.
C                MC     - INPUT NBR(1) BY 2 MATRIX CONTAINING, IN
C                           COLUMN 1, THE NUMBER OF ELEMENTS IN EACH
C                           CLUSTER IN THE SAMPLE. THE SAMPLED CLUSTER
C                           SIZES MUST BE ORDERED IN CORRESPONDENCE
C                           WITH THE ORDERING OF CLUSTERS IN YY.
C                           COLUMN 2 IS WORK STORAGE.
C                IM     - INPUT ROW DIMENSION OF MATRIX MC EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                SIZE   - INPUT VECTOR OF LENGTH NBR(1), REQUIRED ONLY
C                           WHEN NBR(7) IS NONZERO, CONTAINING A MEASURE
C                           OF CLUSTER SIZE FOR EACH CLUSTER IN THE
C                           SAMPLE.  THE SAMPLED CLUSTER SIZE MEASURES
C                           MUST BE ORDERED IN CORRESPONDENCE WITH
C                           THE ORDERING OF CLUSTERS IN YY.
C                TSIZE  - INPUT MEASURE OF TOTAL SIZE OF ALL CLUSTERS
C                           IN THE POPULATION. REQUIRED ONLY WHEN
C                           NBR(7) IS NONZERO.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVALS FOR THE MEAN AND TOTAL
C                           PARAMETERS.
C                TEMP   - INPUT TEMPORARY MEAN.  REQUIRED ONLY IF
C                           NBR(4) = 0.
C                CMUSIG - OUTPUT NBR(1) BY 2 MATRIX CONTAINING THE
C                           WITHIN CLUSTER MEANS AND VARIANCES IN
C                           COLUMNS 1 AND 2, RESPECTIVELY. THE VALUES
C                           ARE ORDERED IN CORRESPONDENCE WITH THE
C                           ORDERING OF CLUSTERS IN YY.
C                IC     - INPUT ROW DIMENSION OF MATRIX CMUSIG EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                STAT   - OUTPUT VECTOR OF LENGTH 9. STAT(I) CONTAINS,
C                         WHEN
C                           I=1, ESTIMATE OF THE MEAN.
C                           I=2, ESTIMATE OF THE TOTAL.
C                           I=3, VARIANCE ESTIMATE OF THE MEAN ESTIMATE.
C                           I=4, VARIANCE ESTIMATE OF THE TOTAL ESTIMATE
C                           I=5, LOWER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=6, UPPER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=7, LOWER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=8, UPPER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=9, ESTIMATE (EXPRESSED AS A PERCENTAGE)
C                             OF THE COEFFICIENT OF VARIATION OF THE
C                             MEAN AND TOTAL ESTIMATES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES NBR(1) IS LESS THAN 2 OR
C                             THAT NBR(2) IS LESS THAN 1.
C                           IER = 130 INDICATES NBR(3) IS LESS THAN 1 OR
C                             THAT NBR(2)*(NBR(3)-NBR(1)) EXCEEDS THE
C                             NUMBER (TOTAL OF COLUMN 1 OF MC) OF
C                             ELEMENTS IN YY.
C                           IER = 131 INDICATES THAT AT LEAST ONE
C                             ELEMENT OF COLUMN 1 OF MC IS LESS THAN 1
C                             OR THAT ALPHA IS NOT IN THE EXCLUSIVE
C                             INTERVAL (0,1).
C                           IER = 132 INDICATES THAT NBR(2) IS GREATER
C                             THAN SOME MC(I,1),I=1,...,NBR(1) WHEN
C                             NBR(2) IS LESS THAN THE TOTAL SAMPLE
C                             SIZE, OR THAT NBR(2) EXCEEDS THE TOTAL
C                             SAMPLE SIZE.
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      BETWEEN THE FIRST AND LAST CALL TO SSSCAN, EXCEPT FOR
C                OPTIONALLY CHANGING TEMP BETWEEN CLUSTERS, ONLY
C                NBR(3) MAY BE MODIFIED AND IT SHOULD FOLLOW THE PATTERN
C                1,2,... . THOUGH THIS PATTERN IS THE OBVIOUS ONE TO
C                FOLLOW, IT IS NOT NECESSARY IN ITS ENTIRETY. FOR THE
C                FIRST AND LAST SUBVECTOR WITHIN EACH CLUSTER, NBR(3)
C                MUST HAVE ITS CORRECT PATTERN VALUE. FOR THE INTER-
C                MEDIATE SUBVECTORS WITHIN EACH STRATUM, NBR(3) MAY
C                TAKE ANY VALUE BETWEEN THE CORRESPONDING FIRST AND
C                LAST SUBVECTOR NBR(3) VALUES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSSCAN (Y,IOPT,NBR,MC,IM,SIZE,TSIZE,ALPHA,TEMP,CMUSIG,
     1                   IC,STAT,IER)
C
      DIMENSION          Y(1),NBR(1),MC(IM,1),CMUSIG(IC,1),STAT(1)
      DIMENSION          SIZE(1),T(4)
      DOUBLE PRECISION   SUMA,SUMB,WK,DZERO
      EQUIVALENCE        (T(1),SUMA),(T(3),SUMB)
      DATA               ZERO,HALF,ONE,HUND/0.0,0.5,1.0,100.0/
      DATA               DZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IN = 0
      IF (NBR(1) .GE. 2 .AND. NBR(2) .GE. 1) GO TO 5
C                                  TERMINAL ERROR - INVALID NBR VALUES
      IER = 129
      GO TO 9000
    5 IF (NBR(3) .GE. 1) GO TO 15
   10 IER = 130
      GO TO 9000
   15 N = NBR(1)
      M = 0
      DO 20 I=1,N
         M = M+MC(I,1)
         IF (MC(I,1) .LT. 1) GO TO 25
   20 CONTINUE
      IF (NBR(2)*(NBR(3)-NBR(1)) .GT. M) GO TO 10
      GO TO 30
   25 IER = 131
      GO TO 9000
   30 IF (ALPHA .LE. ZERO .OR. ALPHA .GE. ONE) GO TO 25
      IF (NBR(2) .GT. M) GO TO 45
      IF (NBR(2) .LT. M) GO TO 35
      MM = 1
      GO TO 50
   35 MM = 0
      NBR2 = NBR(2)-1
      DO 40 I=1,N
         IF (MC(I,1) .LT. NBR(2)) GO TO 45
C                                  COMPUTE THE NUMBER OF SUBVECTORS IN
C                                    EACH CLUSTER AND COMPUTE THE
C                                    PARTIAL SUMS
         MC(I,2) = (MC(I,1)+NBR2)/NBR(2)
         MM = MM+MC(I,2)
         IF (I .EQ. 1) GO TO 40
         MC(I,2) = MC(I-1,2)+MC(I,2)
   40 CONTINUE
      GO TO 50
C                                  TERMINAL ERROR - NBR(2) .GT. MC(I,1)
C                                    FOR I = 1,...,NBR(1)
   45 IER = 132
      GO TO 9000
   50 IF (NBR(3) .NE. 1) GO TO 55
      STAT(9) = 1.0
      IF (NBR(4) .NE. 0) TEMP = Y(1)
   55 IF (MM .EQ. 1) GO TO 80
C                                  COMPUTE THE NUMBER OF THE CLUSTER
C                                    FOR THE CURRENT SUBVECTOR AND
C                                    DETERMINE WHETHER THE SUBVECTOR
C                                    IS THE FIRST, LAST, OR AN INTER-
C                                    MEDIATE SUBVECTOR OF THAT CLUSTER
      IS = STAT(9)
      IF (NBR(3) .GT. MC(IS,2)) IS = IS+1
      NS = MC(IS,2)-1
      IF (IS .NE. 1) NS = NS-MC(IS-1,2)
      NN = NBR(2)
      IF (NBR(3) .EQ. MC(IS,2)) NN = MC(IS,1)-NS*NN
      II = 1
      JJ = NN
      IF (NBR(3) .EQ. 1) GO TO 60
      IF(IS.EQ.1) GO TO 70
      IF (NBR(3) .NE. MC(IS-1,2)+1) GO TO 70
C                                  MOVE THE FIRST ELEMENT OF Y OF
C                                    CURRENT CLUSTER INTO TEMP
      IF (NBR(4) .NE. 0) TEMP = Y(1)
C                                  INITIALIZE STAT VECTOR FOR EACH
C                                    CLUSTER
   60 DO 65 I=1,8
         STAT(I) = ZERO
   65 CONTINUE
   70 T(1) = STAT(1)
      T(2) = STAT(2)
      T(3) = STAT(3)
      T(4) = STAT(4)
      DO 75 I=II,JJ
         SUMA = SUMA+DBLE(Y(I))
         WK = DBLE(Y(I))-DBLE(TEMP)
         SUMB = SUMB+WK**2
   75 CONTINUE
      IF (IN .NE. 0) GO TO 95
      STAT(1) = T(1)
      STAT(2) = T(2)
      STAT(3) = T(3)
      STAT(4) = T(4)
      IF (NBR(3) .EQ. MC(IS,2)) GO TO 95
      GO TO 9005
C                                  NBR(2) IS EQUAL TO THE TOTAL SAMPLE
C                                    SIZE CASE
   80 IF (IN .EQ. N) GO TO 100
      IN = IN+1
      NN = MC(IN,1)
      IF (IN .EQ. 1) GO TO 90
      II = JJ+1
      JJ = JJ+NN
      IF (NBR(4) .NE. 0) TEMP = Y(II)
      GO TO 60
   90 II = 1
      JJ = NN
      GO TO 60
   95 IF (IN .NE. 0) IS = IN
C                                  COMPUTE MEAN AND VARIANCE ESTIMATES
      FNN = MC(IS,1)
      WK = SUMA/FNN
      WK = DBLE(FNN)*(WK-TEMP)**2
      CMUSIG(IS,1) = SUMA
      CMUSIG(IS,2) = SUMB-WK
      STAT(9) = IS
      IF (MM .EQ. 1) GO TO 80
      IF (IS .NE. N) GO TO 9005
  100 FRN1 = ONE/FLOAT(N-1)
      FRN = ONE/FLOAT(N)
      IF (IOPT .GT. 0) GO TO 130
      FRN6 = ONE/FLOAT(NBR(6))
      WK = DZERO
      DO 105 I=1,N
         WK = WK+DBLE(CMUSIG(I,1))
  105 CONTINUE
      TT = FLOAT(NBR(5)*(NBR(5)-N))*FRN
      IF (IOPT .EQ. 0) GO TO 115
C                                  IOPT IS LESS THAN 0
      STAT(2) = FLOAT(NBR(5))*FRN*WK
      STAT(1) = STAT(2)*FRN6
      WK = WK*FRN
      SUMA = DZERO
      DO 110 I=1,N
         SUMA = SUMA+(CMUSIG(I,1)-WK)**2
  110 CONTINUE
      SUMA = SUMA*FRN1
      STAT(4) = SUMA*TT
      GO TO 125
C                                  IOPT IS EQUAL TO ZERO
  115 STAT(1) = WK/FLOAT(M)
      STAT(2) = FLOAT(NBR(6))*STAT(1)
      WK = DZERO
      SUMA = DZERO
      SUMB = DZERO
      DO 120 I = 1,N
         WK = WK+DBLE(CMUSIG(I,1))**2
         SUMA = SUMA+DBLE(CMUSIG(I,1))*FLOAT(MC(I,1))
         SUMB = SUMB+FLOAT(MC(I,1))**2
  120 CONTINUE
      R = STAT(2)*FRN6
      STAT(4) = (WK+R**2*SUMB-(R+R)*SUMA)*FRN1
      STAT(4) = STAT(4)*TT
  125 STAT(3) = STAT(4)*FRN6*FRN6
      GO TO 160
C                                  IOPT IS GREATER THAN ZERO
  130 WK = DZERO
      SUMA = DZERO
      IF (NBR(7) .NE. 0) GO TO 145
C                                  IOPT .GT. 0 AND NBR(7) .EQ. 0 CASE
      DO 135 I=1,N
         CMUSIG(I,1) = CMUSIG(I,1)/FLOAT(MC(I,1))
         WK = WK+CMUSIG(I,1)
  135 CONTINUE
      STAT(1) = WK*FRN
      STAT(2) = FLOAT(NBR(6))*STAT(1)
      DO 140 I=1,N
         WK = DBLE(CMUSIG(I,1))-DBLE(STAT(1))
         SUMA = SUMA+WK**2
  140 CONTINUE
      SUMA = SUMA*FRN1
      STAT(3) = SUMA*FRN
      STAT(4) = FLOAT(NBR(6))**2*STAT(3)
      GO TO 160
C                                  NBR(7) .NE. 0 CASE
  145 DO 150 I=1,N
         WK = WK+CMUSIG(I,1)/SIZE(I)
  150 CONTINUE
      STAT(1) = WK*FRN
      STAT(2) = TSIZE*STAT(1)
      DO 155 I=1,N
         ZR = TSIZE/SIZE(I)
         SUMA = SUMA+(CMUSIG(I,1)*ZR-DBLE(STAT(2)))**2
  155 CONTINUE
      STAT(4) = SUMA*FRN*FRN1
      STAT(3) = STAT(4)/(TSIZE**2)
C                                  CALL ROUTINE MDNRIS
  160 DO 165 I=1,N
         FMC = ONE/FLOAT(MC(I,1))
         CMUSIG(I,2) = CMUSIG(I,2)*FMC
         IF (IOPT .GT. 0 .AND. NBR(7) .EQ. 0) GO TO 165
         CMUSIG(I,1) = CMUSIG(I,1)*FMC
  165 CONTINUE
      P = ONE-ALPHA*HALF
      CALL MDNRIS (P,TT,IER)
      T(1) = TT*SQRT(STAT(3))
      T(2) = TT*SQRT(STAT(4))
      STAT(5) = STAT(1)-T(1)
      STAT(6) = STAT(1)+T(1)
      STAT(7) = STAT(2)-T(2)
      STAT(8) = STAT(2)+T(2)
      STAT(9) = (HUND*SQRT(STAT(4)))/STAT(2)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HSSSCAN)
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - SSSEST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TWO-STAGE SAMPLING WITH CONTINUOUS DATA AND
C                           EQUISIZED PRIMARY UNITS - INFERENCES
C                           REGARDING THE POPULATION MEAN AND TOTAL
C
C   USAGE               - CALL SSSEST(Y,NBR,ALPHA,TEMP,SMUSIG,IS,STAT,
C                           IER)
C
C   ARGUMENTS    Y      - INPUT SUBVECTOR OF LENGTH NBR(2) OF THE VECTOR
C                           (CALL IT YY) CONTAINING THE ENTIRE TWO-STAGE
C                           SAMPLE. THE SUBVECTOR Y MUST BE EITHER THE
C                           SAME AS YY OR IT MUST CONTAIN ALL OR PART OF
C                           THE SAMPLE FOR A SINGLE PRIMARY UNIT ONLY.
C                           IN THE CASE WHERE Y IS A PROPER SUBSET OF
C                           YY, THE LAST SUBVECTOR FROM ANY ONE OR MORE
C                           OF THE PRIMARY UNITS MAY HAVE FEWER
C                           THAN NBR(2) ELEMENTS. IN ALL CASES THE
C                           OBSERVATIONS WITHIN ANY ONE PRIMARY UNIT
C                           MUST APPEAR CONTIGUOUSLY IN YY.
C                NBR    - INPUT VECTOR OF LENGTH 7. NBR(I) CONTAINS,
C                           WHEN,
C                           I=1, NUMBER OF PRIMARY UNITS INTO WHICH
C                             THE SAMPLE IS DIVIDED.
C                           I=2, NUMBER OF OBSERVATIONS IN EACH
C                             SUBVECTOR Y, NOT INCLUDING THE LAST
C                             SUBVECTOR IN EACH PRIMARY UNIT, WHERE
C                             THE NUMBER MAY BE LESS THAN OR EQUAL TO
C                             NBR(2). HOWEVER, NBR(2) SHOULD BE THE SAME
C                             FOR ALL CALLS.
C                           I=3, THE NUMBER OF THE SUBVECTOR STORED IN
C                             Y. SEE REMARKS.
C                           I=4, TEMPORARY MEAN INDICATOR. IF NBR(4)=0,
C                             THE USER SUPPLIES THE TEMPORARY MEAN IN
C                             TEMP. THE USER MAY CHANGE TEMP ONLY WHEN
C                             SUBVECTOR Y IS THE FIRST SET OF
C                             OBSERVATIONS IN A PRIMARY UNIT. IF NBR(4)
C                             IS NON-ZERO, THE FIRST ELEMENT OF EACH
C                             PRIMARY UNIT IS UTILIZED.
C                           I=5, NUMBER OF PRIMARY UNITS IN THE SAMPLED
C                             POPULATION.
C                           I=6, NUMBER OF ELEMENTS IN EACH PRIMARY
C                             UNIT IN THE POPULATION.
C                           I=7, NUMBER OF ELEMENTS IN THE SAMPLE IN
C                             EACH SAMPLED PRIMARY UNIT.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVALS FOR THE MEAN AND
C                           TOTAL PARAMETERS.
C                TEMP   - INPUT TEMPORARY MEAN. REQUIRED ONLY
C                           IF NBR(4)=0.
C                SMUSIG - OUTPUT NBR(1) BY 2 MATRIX CONTAINING THE
C                           WITHIN PRIMARY UNIT MEAN AND VARIANCE
C                           ESTIMATES IN COLUMNS 1 AND 2, RESPECTIVELY.
C                           THE ESTIMATES ARE ORDERED IN CORRESPONDENCE
C                           WITH THE ORDERING OF PRIMARY UNITS IN YY.
C                IS     - INPUT ROW DIMENSION OF MATRIX SMUSIG EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                STAT   - OUTPUT VECTOR OF LENGTH 9. STAT(I) CONTAINS,
C                           WHEN
C                           I=1, ESTIMATE OF THE MEAN.
C                           I=2, ESTIMATE OF THE TOTAL.
C                           I=3, VARIANCE ESTIMATE OF THE MEAN ESTIMATE.
C                           I=4, VARIANCE ESTIMATE OF THE TOTAL
C                             ESTIMATE.
C                           I=5, LOWER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=6, UPPER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=7, LOWER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=8, UPPER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=9, ESTIMATE (EXPRESSED AS A PERCENTAGE) OF
C                             THE COEFFICIENT OF VARIATION OF THE MEAN
C                             AND TOTAL ESTIMATES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES NBR(1) IS LESS THAN 2 OR
C                             THAT NBR(2) IS LESS THAN 1.
C                           IER = 130 INDICATES NBR(3) IS LESS THAN 1 OR
C                             THAT NBR(2)*(NBR(3)-NBR(1)) EXCEEDS THE
C                             NUMBER (NBR(1)*NBR(7)) OF ELEMENTS IN YY.
C                           IER = 131 INDICATES THAT NBR(1) EXCEEDS
C                             NBR(5) OR THAT ALPHA IS NOT IN THE
C                             EXCLUSIVE INTERVAL (0,1).
C                           IER = 132 INDICATES THAT NBR(7) IS LESS THAN
C                             2 OR THAT NBR(7) EXCEEDS NBR(6).
C                           IER = 133 INDICATES THAT NBR(2) EXCEEDS
C                             NBR(1)*NBR(7).
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      BETWEEN THE FIRST AND LAST CALL TO SSSEST, EXCEPT FOR
C                OPTIONALLY CHANGING TEMP BETWEEN PRIMARY UNITS, ONLY
C                NBR(3) MAY BE MODIFIED AND IT SHOULD FOLLOW THE PATTERN
C                1,2,... . THOUGH THIS PATTERN IS THE OBVIOUS ONE TO
C                FOLLOW, IT IS NOT NECESSARY IN ITS ENTIRETY. FOR THE
C                FIRST AND LAST SUBVECTOR WITHIN EACH PRIMARY UNIT,
C                NBR(3) MUST HAVE ITS CORRECT PATTERN VALUE. FOR THE
C                INTERMEDIATE SUBVECTORS WITHIN EACH PRIMARY UNIT,
C                NBR(3) MAY TAKE ANY VALUE BETWEEN THE CORRESPONDING
C                FIRST AND LAST SUBVECTOR NBR(3) VALUES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSSEST(Y,NBR,ALPHA,TEMP,SMUSIG,IS,STAT,IER)
C
      DIMENSION          Y(1),NBR(1),SMUSIG(IS,1),STAT(1),TA(4)
      DOUBLE PRECISION   SA,SB,YDB,QSUM,SUM,TT
      EQUIVALENCE        (TA(1),SUM),(TA(3),QSUM)
C                                  FIRST EXECUTABLE STATEMENT
      N = NBR(1)
      M = NBR(7)
      RN = N
      RM = M
      RNN = NBR(5)
      RMM = NBR(6)
      RN1 = RN-1.0
      RM1 = RM-1.0
      NM = N*M
      RNM = RNN*RMM
      IER = 0
C                                  TERMINAL ERRORS - INVALID NBR AND
C                                  ALPHA VALUES
      IF(N .GE. 2 .AND. NBR(2) .GE. 1) GO TO 5
      IER = 129
      GO TO 9000
    5 IF(NBR(3) .GE. 1) GO TO 15
   10 IER = 130
      GO TO 9000
   15 IF(NBR(2)*(NBR(3)-NBR(1)) .GT. NM) GO TO 10
      IF(NBR(1) .LE. NBR(5)) GO TO 25
   20 IER = 131
      GO TO 9000
   25 IF(ALPHA .LE. 0. .OR. ALPHA .GE. 1.) GO TO 20
      IF(NBR(7) .GE. 2) GO TO 35
   30 IER = 132
      GO TO 9000
   35 IF(NBR(7) .GT. NBR(6)) GO TO 30
      IF(NBR(2) .LE. NM)GO TO 40
      IER = 133
      GO TO 9000
C                                  COMPUTE THE NUMBER OF THE LAST
C                                  SUBVECTOR IN THE PRIMARY UNIT
   40 L = ((M+NBR(2)-1)/NBR(2))
      IF(NBR(3) .NE. 1) GO TO 45
      STAT(8) = 1.0
      STAT(9) = 1.0
C                                  DETERMINE THE NUMBER OF THE CURRENT
C                                  SUBVECTOR IN THE PRIMARY UNIT
   45 KSPOT = STAT(8)
      IF(KSPOT .NE. 1) GO TO 50
      SB = 0.D0
      STAT(1) = 0.0
      STAT(2) = 0.0
      STAT(3) = 0.0
      STAT(4) = 0.0
   50 TA(1) = STAT(1)
      TA(2) = STAT(2)
      TA(3) = STAT(3)
      TA(4) = STAT(4)
      ISPOT = STAT(9)
      IF(KSPOT .GT. 1) GO TO 55
C                                  MOVE FIRST ELEMENT OF Y INTO TEMP
      IF(NBR(4) .NE. 0) TEMP = Y(1)
      T = TEMP
   55 IF(NBR(2) .EQ. NM) GO TO 70
C                                  NBR(2) .NE. NM
      MEND = NBR(2)
      IF(KSPOT .EQ. L) MEND = M-(L-1)*NBR(2)
   60 DO 65 J=1,MEND
         TT = DBLE(Y(J))-T
         QSUM = QSUM+TT*TT
         SUM = SUM+DBLE(Y(J))
   65 CONTINUE
      STAT(1) = TA(1)
      STAT(2) = TA(2)
      STAT(3) = TA(3)
      STAT(4) = TA(4)
      GO TO 85
C                                  NBR(2) .EQ. NM
   70 JLJ = 1
      MLM = M
      DO 80 ISPOT = 1,N
         SUM = 0.D0
         QSUM = 0.D0
         IF(NBR(4) .NE. 0) TEMP = Y(JLJ)
         T = TEMP
         DO 75 J = JLJ,MLM
            TT = DBLE(Y(J))-T
            SUM = SUM+DBLE(Y(J))
            QSUM = QSUM+TT*TT
   75    CONTINUE
         SUM = SUM/RM
         SMUSIG(ISPOT,1) = SUM
         SMUSIG(ISPOT,2) = QSUM-RM*(SUM-T)*(SUM-T)
         JLJ = JLJ+M
         MLM = MLM+M
   80 CONTINUE
      GO TO 100
   85 IF(ISPOT .EQ. NBR(1) .AND. KSPOT .EQ. L) GO TO 95
      IF(KSPOT .EQ. L) GO TO 90
      STAT(8) = STAT(8)+1.0
      GO TO 9005
   90 SUM = SUM/RM
      SMUSIG(ISPOT,1) = SUM
      SMUSIG(ISPOT,2) = QSUM-RM*(SUM-T)*(SUM-T)
      STAT(8) = 1.0
      STAT(9) = STAT(9) + 1.0
      GO TO 9005
   95 SUM = SUM/RM
      SMUSIG(ISPOT,1) = SUM
      SMUSIG(ISPOT,2) = QSUM-RM*(SUM-T)*(SUM-T)
  100 FN1 = RN/RNN
      FN2 = RM/RMM
      SB = 0.D0
      SA = 0.D0
      YDB = 0.D0
      DO 105 I = 1,N
         SB = SB+DBLE(SMUSIG(I,2))
         SMUSIG(I,2) = SMUSIG(I,2)/RM1
         YDB = DBLE(SMUSIG(I,1)) + YDB
  105 CONTINUE
      SB = SB/(RN*RM1)
      YDB = YDB/RN
      DO 110 I=1,N
         SA = SA+(SMUSIG(I,1)-YDB)*(SMUSIG(I,1)-YDB)
  110 CONTINUE
      SA = SA/RN1
      P = 1.-ALPHA*.5
      CALL MDNRIS(P,T1,IER)
      STAT(1) = YDB
      STAT(2) = RNM*STAT(1)
      STAT(3) = ((1.-FN1)*SA)/RN+FN1*(1.-FN2)*SB/(RN*RM)
      STAT(4) = RNM*RNM*STAT(3)
      DELTA1 = T1*SQRT(STAT(3))
      DELTA2 = T1*SQRT(STAT(4))
      STAT(5) = STAT(1)-DELTA1
      STAT(6) = STAT(1)+DELTA1
      STAT(7) = STAT(2)-DELTA2
      STAT(8) = STAT(2) + DELTA2
      STAT(9) = 100.*SQRT(STAT(4))/STAT(2)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HSSSEST)
 9005 RETURN
      END

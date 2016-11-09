C   IMSL ROUTINE NAME   - SSPBLK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STRATIFIED RANDOM SAMPLING WITH PROPORTION
C                           DATA - INFERENCES REGARDING THE
C                           POPULATION PROPORTION AND TOTAL
C
C   USAGE               - CALL SSPBLK (NBR,NH,IN,ALPHA,PH,STAT,IER)
C
C   ARGUMENTS    NBR    - INPUT NUMBER OF STRATA INTO WHICH THE SAMPLE
C                           IS DIVIDED.
C                NH     - INPUT NBR BY 3 MATRIX CONTAINING THE OBSERVED
C                           NUMBER OF UNITS IN EACH STRATUM FOR THE
C                           CLASS OF INTEREST AND THE SAMPLE IN COLUMNS
C                           1 AND 2, RESPECTIVELY. THE NUMBER OF UNITS
C                           IN EACH STRATUM IN THE POPULATION ARE IN
C                           COLUMN 3. EACH ROW MUST CONTAIN THE INFOR-
C                           MATION FOR A SINGLE STRATUM ONLY. IN THE
C                           CASE WHERE POPULATION STRATA SIZES ARE NOT
C                           KNOWN, ESTIMATES MUST BE ENTERED IN THEIR
C                           PLACE.
C                IN     - INPUT ROW DIMENSION OF MATRIX NH EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVALS FOR THE PROPORTION AND
C                           TOTAL PARAMETERS.
C                PH     - OUTPUT VECTOR OF LENGTH NBR CONTAINING THE
C                           WITHIN STRATA PROPORTION ESTIMATES.
C                STAT   - OUTPUT VECTOR OF LENGTH 10. STAT(I) CONTAINS,
C                         WHEN
C                           I=1, ESTIMATE OF THE PROPORTION.
C                           I=2, ESTIMATE OF THE TOTAL.
C                           I=3, VARIANCE ESTIMATE OF THE PROPORTION
C                             ESTIMATE.
C                           I=4, VARIANCE ESTIMATE OF THE TOTAL
C                             ESTIMATE.
C                           I=5, LOWER CONFIDENCE LIMIT FOR THE
C                             PROPORTION.
C                           I=6, UPPER CONFIDENCE LIMIT FOR THE
C                             PROPORTION.
C                           I=7, LOWER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=8, UPPER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=9, ESTIMATE (EXPRESSED AS A PERCENTAGE)
C                             OF THE COEFFICIENT OF VARIATION OF THE
C                             TOTAL ESTIMATE.
C                           I=10, VARIANCE ESTIMATE OF THE MEAN ESTI-
C                             MATE ASSUMING THAT SAMPLING WAS SIMPLE
C                             RANDOM INSTEAD OF STRATIFIED RANDOM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES NBR IS LESS THAN 2.
C                           IER = 130 INDICATES AT LEAST ONE ROW OF NH
C                             IS SPECIFIED INCORRECTLY. THIS ERROR
C                             OCCURS, IF IN ANY ROW, THE FIRST ELEMENT
C                             EXCEEDS THE SECOND ELEMENT, OR THE SECOND
C                             ELEMENT EXCEEDS THE THIRD ELEMENT, OR THE
C                             SECOND ELEMENT IS LESS THAN 2.
C                           IER = 131 INDICATES THAT ALPHA IS NOT IN THE
C                             EXCLUSIVE INTERVAL (0,1).
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSPBLK (NBR,NH,IN,ALPHA,PH,STAT,IER)
C
      DIMENSION          NH(IN,1),PH(1),STAT(1)
      DOUBLE PRECISION   SUMC,SUMD,SUME,SUMF,SUMG,DZERO
      DATA               ZERO,HALF,ONE,TWO,HUND/0.0,0.5,1.0,2.0,100.0/
      DATA               DZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NBR .GE. 2) GO TO 5
C                                  TERMINAL ERROR - NBR .LT. 2
      IER = 129
      GO TO 9000
C                                  TERMINAL ERROR - ALPHA OUT OF RANGE
    5 IF (ALPHA .GT. ZERO .AND. ALPHA .LT. ONE) GO TO 10
      IER = 131
      GO TO 9000
   10 NN = 0
      DO 15 I=1,NBR
         NN = NN+NH(I,3)
   15 CONTINUE
      SUMB = NN
      SUMC = DZERO
      SUMD = DZERO
      SUME = DZERO
      SUMF = DZERO
      SUMG = DZERO
      SUMA = ZERO
      STAT(2) = ZERO
      STAT(4) = ZERO
      DO 30 I=1,NBR
         SAL = NH(I,1)
         SNL = NH(I,2)
         BNL = NH(I,3)
         IF (SAL .GT. SNL .OR. SNL .GT. BNL) GO TO 20
         IF (SNL .GE. TWO) GO TO 25
C                                  TERMINAL ERROR - NH SPECIFIED
C                                    INCORRECTLY
   20    IER = 130
         GO TO 9000
   25    SUMA = SUMA+SNL
         PH(I) = SAL/SNL
         PHA = PH(I)*(ONE-PH(I))
         TA =(DBLE(BNL-SNL)/(SNL-ONE))*PHA
         STAT(2) = STAT(2)+DBLE(BNL*PH(I))
         STAT(4) = STAT(4)+DBLE(TA)*BNL
         WH = BNL/SUMB
         TB = TA/SUMB
         SUMC = SUMC+TB
         TC = TB/SNL
         SUMD = SUMD+TC
         SUME = SUME+TC*WH
         TD = WH*PH(I)
         SUMF = SUMF+TD
         SUMG = SUMG+TD*PH(I)
   30 CONTINUE
      RSUMB = ONE/SUMB
      STAT(1) = STAT(2)*RSUMB
      STAT(3) = STAT(4)*RSUMB*RSUMB
      P = ONE-ALPHA*HALF
      CALL MDNRIS (P,TA,IER)
      TB = TA*SQRT(STAT(3))
      TC = SQRT(STAT(4))
      TA = TA*TC
      STAT(5) = AMAX1(ZERO,STAT(1)-TB)
      STAT(6) = AMIN1(ONE,STAT(1)+TB)
      STAT(7) = AMAX1(ZERO,STAT(2)-TA)
      STAT(8) = AMIN1(SUMB,STAT(2)+TA)
      STAT(9) = HUND*TC/STAT(2)
      TA = (SUMB-SUMA)/(SUMA*(SUMB-ONE))
      STAT(10) = TA*((ONE-RSUMB)*SUMC-SUMD+SUME+SUMG-SUMF*SUMF)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HSSPBLK)
 9005 RETURN
      END

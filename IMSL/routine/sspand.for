C   IMSL ROUTINE NAME   - SSPAND
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SIMPLE RANDOM SAMPLING WITH PROPORTION
C                           DATA-INFERENCES REGARDING THE POPULATION
C                           PROPORTION AND TOTAL
C
C   USAGE               - CALL SSPAND (IOPT,NBR,ALPHA,STAT,IDIST,IER)
C
C   ARGUMENTS    IOPT   - INPUT SUBPOPULATION AND NUMBER OF CLASSES
C                         INDICATOR. THE USUAL SITUATION FOR WHICH THIS
C                         PROGRAM IS APPLICABLE IS ONE WHERE A POPULA-
C                         TION IS SAMPLED AND EACH SAMPLE UNIT IS
C                         CATEGORIZED INTO ONE OF TWO CLASSES. EXTEN-
C                         SIONS TO MORE THAN TWO CLASSES WHERE SOME
C                         CLASSES ARE OR ARE NOT OMITTED, TO SUBPOPULA-
C                         TIONS, AND TO ANY COMBINATION OF THE PRECED-
C                         ING IS ALLOWED.
C                           IF IOPT = 0, THE SAMPLE IS FROM A POPULA-
C                             TION WHERE THERE ARE TWO OR MORE CLASSES,
C                             BUT NO CLASSES ARE OMITTED.
C                           IF IOPT IS NEGATIVE, THE SAMPLE IS FROM A
C                             POPULATION WHERE THERE ARE THREE OR MORE
C                             CLASSES AND AT LEAST ONE CLASS IS OMIT-
C                             TED, OR THE SAMPLE IS FROM A SUBPOPULA-
C                             TION WHERE THERE ARE TWO OR MORE CLASSES
C                             AND CLASSES MAY OR MAY NOT BE OMITTED.
C                             THE NUMBER OF UNITS, NBR(4), IN THE
C                             POPULATION (OR SUBPOPULATION) CORRESPOND-
C                             ING TO THE INCLUDED CLASSES IS UNKNOWN.
C                           IF IOPT IS POSITIVE, THE SITUATION IS AS
C                             FOR IOPT NEGATIVE, EXCEPT THAT NBR(4) IS
C                             KNOWN.
C                NBR    - INPUT VECTOR OF LENGTH 5. NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF SAMPLE UNITS IN THE CLASS
C                             OF INTEREST, FOR THE POPULATION (OR SUB-
C                             POPULATION) OF INTEREST.
C                           I=2, NUMBER OF SAMPLE UNITS IN THE IN-
C                             CLUDED CLASSES, FOR THE POPULATION (OR
C                             SUBPOPULATION) OF INTEREST.
C                           I=3, NUMBER OF SAMPLE UNITS IN THE ENTIRE
C                             RANDOM SAMPLE. REQUIRED ONLY WHEN IOPT IS
C                             NEGATIVE.
C                           I=4, NUMBER OF UNITS IN THE POPULATION (OR
C                             SUBPOPULATION) CORRESPONDING TO THE
C                             INCLUDED CLASSES. NOT REQUIRED WHEN IOPT
C                             IS NEGATIVE.
C                           I=5, NUMBER OF UNITS IN THE POPULATION.
C                             REQUIRED ONLY WHEN IOPT IS NEGATIVE.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVALS FOR THE PROPORTION AND
C                           TOTAL PARAMETERS FOR THE CLASS AND POPULA-
C                           TION (OR SUBPOPULATION) OF INTEREST.
C                STAT   - OUTPUT VECTOR OF LENGTH 9. STAT(I) CONTAINS,
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
C                           I=9,ESTIMATE (EXPRESSED AS A PERCENTAGE) OF
C                             THE COEFFICIENT OF VARIATION OF THE
C                             TOTAL ESTIMATE. NOT RETURNED IF NBR(1)=0.
C                IDIST  - OUTPUT INDICATOR OF THE DISTRIBUTION USED TO
C                         APPROXIMATE THE HYPERGEOMETRIC DISTRIBUTION
C                         FOR THE CONFIDENCE INTERVAL CALCULATIONS.
C                           IF IDIST = 0, THE NORMAL IS USED.
C                           IF IDIST IS NEGATIVE, THE POISSON IS USED.
C                           IF IDIST IS POSITIVE, THE BINOMIAL IS USED.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES NBR(1) IS LESS THAN 0,
C                             OR THAT NBR(2) IS LESS THAN 2, OR THAT
C                             NBR(1) EXCEEDS NBR(2).
C                           IER = 130 INDICATES, WHEN IOPT IS NEGATIVE,
C                             THAT NBR(3) IS LESS THAN NBR(2), OR THAT
C                             NBR(5) IS LESS THAN NBR(3).
C                           IER = 131 INDICATES, WHEN IOPT IS NON-NEG-
C                             ATIVE, THAT NBR(2) EXCEEDS NBR(4).
C                           IER = 132 INDICATES THAT ALPHA IS NOT IN
C                             THE EXCLUSIVE INTERVAL (0,1), OR THAT AN
C                             ERROR OCCURRED IN IMSL ROUTINE
C                             BELBIN OR BELPOS.
C
C   REQD. IMSL ROUTINES - H32/BELBIN,BELPOS,MDBETA,MDBETI,MDCH,MDCHI,
C                           MDNOR,MDNRIS,MERFI,MERRC=ERFC,MGAMAD=DGAMMA,
C                           MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/BELBIN,BELPOS,MDBETA,MDBETI,
C                           MDCH,MDCHI,MDNOR,MDNRIS,MERFI,MERRC=ERFC,
C                           MGAMA=GAMMA,MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSPAND (IOPT,NBR,ALPHA,STAT,IDIST,IER)
C
      DIMENSION          NBR(1),STAT(1)
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      RN1 = NBR(1)
      RN2 = NBR(2)
      RN3 = NBR(3)
      RN4 = NBR(4)
      RN5 = NBR(5)
      IF(RN1 .GE. 0.0 .AND. RN1 .LE. RN2) GO TO 10
C                                  TERMINAL - NBR(1) .LT. 0 OR NBR(2)
C                                  .LT. 2 OR NBR(1) .GT. NBR(2)
    5 IER = 129
      GO TO 9000
   10 IF(RN2 .LT. 2.) GO TO 5
      IF(ALPHA .GT. 0. .AND. ALPHA .LT. 1.) GO TO 15
C                                  TERMINAL ERROR - ALPHA IS OUT
C                                  OF RANGE
      IER = 132
      GO TO 9000
   15 PSTAR = AMIN1(RN1,RN2-RN1)
      PSTAR = PSTAR/RN2
      STAT(1) = RN1/RN2
      SCON = STAT(1)*(1.-STAT(1))
      IF(IOPT .GE. 0) GO TO 20
      IF(RN3 .GE. RN2 .AND. RN5 .GE. RN3) GO TO 30
C                                  TERMINAL ERROR - NBR(3) .LT. NBR(2)
C                                  OR NBR(5) .LT. NBR(3)
      IER = 130
      GO TO 9000
   20 IF(RN2 .LE. RN4) GO TO 25
C                                  TERMINAL ERROR - NBR(2) .GT. NBR(4)
      IER = 131
      GO TO 9000
   25 RN3 = RN2
      RN5 = RN4
   30 F1 = RN3/RN5
      P = RN1/RN3
      STAT(2) = RN5*P
      STAT(3) = ((1.0-F1)/(RN2-1.0))*SCON
      IF(IOPT .GE. 0) P = STAT(1)
      STAT(4) =((RN5*(RN5-RN3))/(RN3-1.0))* P*(1.0-P)
      SQS4 = SQRT(STAT(4))
      SK = (RN5*RN2)/RN3
   35 IF(PSTAR .LT. 0.05) GO TO 50
      RN2PS = RN2*PSTAR
      IF(PSTAR .GE. 0.30) GO TO 40
      IF(RN2PS .LT. (-184.0*PSTAR+79.2))GO TO 55
      GO TO 45
   40 IF(RN2PS .LT. (-40.0*PSTAR+36.0))GO TO 55
   45 IDIST = 0
      PP = 1.-ALPHA*.5
      CALL MDNRIS(PP,T,JER)
      IF(JER .NE. 0) GO TO 70
      DELTA1 = T*SQRT(STAT(3))+1./(RN2+RN2)
      DELTA2 = T*SQS4+1./(RN3+RN3)
      STAT(5) = STAT(1)-DELTA1
      STAT(6) = STAT(1)+DELTA1
      STAT(7) = STAT(2)-DELTA2
      STAT(8) = STAT(2)+DELTA2
      GO TO 65
   50 CALL BELPOS(NBR(1),1,ALPHA,RLAMHT,RLAMLR,RLAMUP,JER)
      IF(JER .NE. 0) GO TO 70
      IDIST = -1
      RLAMLR = RLAMLR/RN2
      RLAMUP = RLAMUP/RN2
      GO TO 60
   55 CALL BELBIN(NBR(2),NBR(1),ALPHA,RLAMHT,RLAMLR,RLAMUP,JER)
      IF(JER .NE. 0) GO TO 70
      IDIST = 1
   60 SQSF = SQRT(1.-F1)
      STAT(5) = STAT(1)-SQSF*(STAT(1)-RLAMLR)
      STAT(6) = STAT(1)+SQSF*(RLAMUP-STAT(1))
      STAT(7) = SK*STAT(5)
      STAT(8) = SK*STAT(6)
      STAT(9) = 0.0
      IF(RN1 .GT. 0.0 .OR. RN1 .LT. 0.0) GO TO 65
      GO TO 9005
   65 STAT(9) = 100.*SQS4/STAT(2)
      GO TO 9005
   70 IER = 132
 9000 CONTINUE
      CALL UERTST(IER,6HSSPAND)
 9005 RETURN
      END

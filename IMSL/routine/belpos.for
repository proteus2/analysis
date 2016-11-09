C   IMSL ROUTINE NAME   - BELPOS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INTERVAL ESTIMATE OF THE PARAMETER LAMBDA
C                           OF THE POISSON DISTRIBUTION
C
C   USAGE               - CALL BELPOS (NPOIS,NN,ALPHA,RLAMHT,RLAMLR,
C                           RLAMUP,IER)
C
C   ARGUMENTS    NPOIS  - INPUT VECTOR OF LENGTH NN CONTAINING A RANDOM
C                           SAMPLE OF SIZE NN FROM A POISSON
C                           DISTRIBUTION WITH PARAMETER LAMBDA.
C                NN     - INPUT VALUE CONTAINING THE LENGTH OF NPOIS.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING A 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVAL FOR PARAMETER LAMBDA.
C                RLAMHT - OUTPUT ESTIMATE OF LAMBDA.
C                RLAMLR - OUTPUT LOWER CONFIDENCE LIMIT FOR LAMBDA.
C                RLAMUP - OUTPUT UPPER CONFIDENCE LIMIT FOR LAMBDA.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES AT LEAST ONE ELEMENT OF
C                             NPOIS IS NEGATIVE.
C                           IER=130 INDICATES NN IS LESS THAN 1.
C                           IER=131 INDICATES ALPHA IS NOT IN THE
C                             EXCLUSIVE INTERVAL (0,1).
C                           IER=132 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL ROUTINE MDCHI.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDCHI,MDNOR,MDNRIS,MERFI,MERRC=ERFC,
C                           MGAMAD=DGAMMA,UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDCHI,MDNOR,MDNRIS,MERFI,
C                           MERRC=ERFC,MGAMA=GAMMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BELPOS (NPOIS,NN,ALPHA,RLAMHT,RLAMLR,RLAMUP,IER)
C
      INTEGER            NPOIS(1),NN,IER
      REAL               ALPHA,RLAMHT,RLAMLR,RLAMUP
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
      RNN = NN
      IF(NN .GE. 1) GO TO 5
C                                  TERMINAL ERROR - NN IS LESS THAN 1
      IER = 130
      GO TO 9000
    5 IF(ALPHA .GT. 0.0 .AND. ALPHA .LT. 1.) GO TO 10
C                                  TERMINAL ERROR - ALPHA IS OUT OF
C                                  RANGE
      IER = 131
      GO TO 9000
   10 SUM = 0.0
      DO 15 I=1,NN
C                                  TERMINAL ERROR - AT LEAST ONE
C                                  ELEMENT OF NPOIS IS NEGATIVE.
         IF(NPOIS(I) .LT. 0) GO TO 9000
         SUM = SUM+FLOAT(NPOIS(I))
   15 CONTINUE
      IER = 0
      DU = SUM+SUM+2.0
      PL = ALPHA*.5
      PU = 1.-PL
      DL = DU-2.0
      RLAMHT = SUM/RNN
      IF(SUM .LE. 0.0) GO TO 25
      CALL MDCHI(PU,DU,XU,JER)
      IF(JER .GT. 0) GO TO 30
      RLAMUP = XU*.5/RNN
      IF(SUM .GT. 0.9 .AND. SUM .LT. 1.1) GO TO 20
      CALL MDCHI(PL,DL,XL,JER)
      IF(JER .GT. 0) GO TO 30
      RLAMLR = XL*.5/RNN
      GO TO 9005
   20 RLAMLR = -ALOG(PU)/RNN
      GO TO 9005
   25 RLAMUP = (-ALOG(PL))/RNN
      RLAMLR = 0.0
      GO TO 9005
   30 IER = 132
 9000 CONTINUE
      CALL UERTST(IER,'BELPOS')
 9005 RETURN
      END

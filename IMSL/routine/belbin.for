C   IMSL ROUTINE NAME   - BELBIN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INTERVAL ESTIMATE OF THE PARAMETER P OF THE
C                           BINOMIAL DISTRIBUTION
C
C   USAGE               - CALL BELBIN (NTRIAL,NX,ALPHA,PHAT,PLOWER,
C                           PUPPER,IER)
C
C   ARGUMENTS    NTRIAL - INPUT NUMBER OF BERNOULLI TRIALS.
C                NX     - INPUT NUMBER OF SUCCESSES IN THE NTRIAL TRIALS
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING A 100(1-ALPHA) PERCENT
C                           CONFIDENCE INTERVAL FOR PARAMETER P, THE
C                           PROBABILITY OF A SUCCESS ON ANY ONE TRIAL.
C                PHAT   - OUTPUT ESTIMATE OF P.
C                PLOWER - OUTPUT LOWER CONFIDENCE LIMIT FOR P.
C                PUPPER - OUTPUT UPPER CONFIDENCE LIMIT FOR P.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES NTRIAL IS LESS THAN 1.
C                           IER=130 INDICATES NX IS EITHER NEGATIVE OR
C                             IT EXCEEDS NTRIAL.
C                           IER=131 INDICATES ALPHA IS NOT IN THE
C                             EXCLUSIVE INTERVAL (0,1).
C                           IER=132 INDICATES THAT A TERMINAL ERROR
C                             OCCURRED IN IMSL ROUTINE MDBETI.
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MDBETI,MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/MDBETA,MDBETI,MLGAMA=ALGAMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BELBIN (NTRIAL,NX,ALPHA,PHAT,PLOWER,PUPPER,IER)
C
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF(NTRIAL .GE. 1) GO TO 5
C                                  TERMINAL ERROR - NTRIAL IS LESS
C                                  THAN 1
      IER = 129
      GO TO 9000
    5 IF (NX .GE. 0 .AND. NX .LE. NTRIAL) GO TO 10
C                                  TERMINAL ERROR - NX IS NEGATIVE OR
C                                  NX IS GREATER THAN NTRIAL
      IER = 130
      GO TO 9000
   10 IF(ALPHA .GT. 0.0 .AND. ALPHA .LT. 1.0) GO TO 15
C                                  TERMINAL ERROR - ALPHA IS OUT OF
C                                  RANGE
      IER = 131
      GO TO 9000
   15 RN = NTRIAL
      X = NX
      RNMX = RN-X
      PHAT = X/RN
      AL = .5*ALPHA
      IF (X .LE. 0.0) GO TO 20
      CALL MDBETI(AL,X,RNMX+1.,PLOWER,JER)
      IF(JER .NE. 0) GO TO 35
      GO TO 25
   20 PLOWER = 0.0
   25 IF(X .GE.RN) GO TO 30
      CALL MDBETI(AL,RNMX,X+1.,PUPPER,JER)
      IF(JER .NE. 0) GO TO 35
      PUPPER = 1.-PUPPER
      GO TO 9005
   30 PUPPER = 1.0
      GO TO 9005
   35 IER = 132
 9000 CONTINUE
      CALL UERTST(IER,'BELBIN')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - NMKTS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - K-SAMPLE TRENDS TEST AGAINST ORDERED
C                           ALTERNATIVES
C
C   USAGE               - CALL NMKTS (X,XM,K,DSEED,XSTAT,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH XM(1)+XM(2)+...+XM(K)
C                           CONTAINING ALL SAMPLE OBSERVATIONS (SAMPLE 1
C                           FOLLOWED BY SAMPLE 2 AND SO FORTH). THE
C                           SAMPLES ARE ASSUMED TO HAVE BEEN ORDERED BY
C                           THE EXPERIMENTER TO REPRESENT THE DESIRED
C                           ALTERNATIVE HYPOTHESES.
C                XM     - INPUT VECTOR OF LENGTH K CONTAINING THE K
C                           SAMPLE SIZES. XM(I) MUST CONTAIN THE NUMBER
C                           OF OBSERVATIONS IN SAMPLE I,
C                           WHERE I = 1,2,...,K.
C                K      - INPUT NUMBER OF SAMPLES TO BE CONSIDERED.
C                           K SHOULD BE GREATER THAN OR EQUAL TO 3.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                XSTAT  - OUTPUT VECTOR OF LENGTH 17 PROVIDING THE
C                           TEST RESULTS.
C                         XSTAT(1) CONTAINS THE STATISTIC USED TO TEST
C                           THE HYPOTHESIS. (TIES WERE RANDOMIZED)
C                         XSTAT(2) CONTAINS THE CONSERVATIVE STATISTIC
C                           CONSIDERING TIES. THIS STATISTIC WOULD TEND
C                           (IF TIES WERE PRESENT) TO BE UNFAVORABLE
C                           FOR THE ALTERNATIVE HYPOTHESES IN QUESTION.
C                         XSTAT(3) CONTAINS THE PROBABILITY OF REJECTING
C                           THE NULL HYPOTHESIS IN ERROR. THIS
C                           PROBABILITY IS ASSOCIATED WITH THE STATISTIC
C                           IN XSTAT(1).
C                         XSTAT(4) CONTAINS THE PROBABILITY OF REJECTING
C                           THE NULL HYPOTHESIS IN ERROR. THIS
C                           PROBABILITY IS ASSOCIATED WITH THE STATISTIC
C                           IN XSTAT(2).
C                         XSTAT(5) CONTAINS THE SAME PROBABILITY AS
C                           XSTAT(3) BUT WITH A CONTINUITY CORRECTION
C                         XSTAT(6) CONTAINS THE SAME PROBABILITY AS
C                           XSTAT(4) BUT WITH A CONTINUITY CORRECTION
C                         XSTAT(7) CONTAINS THE EXPECTED MEAN OF THE
C                           STATISTIC
C                         XSTAT(8) CONTAINS THE EXPECTED KURTOSIS OF THE
C                           STATISTIC (SKEWNESS EXPECTATION IS ZERO)
C                         XSTAT(9) CONTAINS THE SIZE OF THE TOTAL SET OF
C                           SAMPLES (SUM OF THE VECTOR XM)
C                         XSTAT(10) CONTAINS THE COEFFICIENT OF RANK
C                           CORRELATION USING XSTAT(1)
C                         XSTAT(11) CONTAINS THE COEFFICIENT OF RANK
C                           CORRELATION USING XSTAT(2)
C                         XSTAT(12) CONTAINS THE NUMBER OF TIES (BETWEEN
C                           SAMPLES) WHICH OCCURRED
C                         XSTAT(13) CONTAINS THE T-STATISTIC ASSOCIATED
C                           WITH XSTAT(3)
C                         XSTAT(14) CONTAINS THE T-STATISTIC ASSOCIATED
C                           WITH XSTAT(4)
C                         XSTAT(15) CONTAINS THE T-STATISTIC ASSOCIATED
C                           WITH XSTAT(5)
C                         XSTAT(16) CONTAINS THE T-STATISTIC ASSOCIATED
C                           WITH XSTAT(6)
C                         XSTAT(17) CONTAINS THE DEGREES OF FREEDOM FOR
C                           THE T-STATISTIC
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT SOME SAMPLE SIZE
C                             XM(I), FOR I = 1,2,...,K, WAS SPECIFIED
C                             LESS THAN OR EQUAL TO ZERO.
C                         WARNING ERROR
C                           IER=34 INDICATES THAT THE DEGREES OF FREEDOM
C                             WERE LESS THAN 1. XSTAT(3),...,XSTAT(6)
C                             CONTAIN MEANINGLESS RESULTS.
C                           IER=35 INDICATES THAT TIES EXIST AND HAVE
C                             BEEN RANDOMIZED.
C                           IER=36 INDICATES THAT ONLY TWO SAMPLES
C                             WERE SPECIFIED, AND THAT POSSIBLY A MORE
C                             POWERFUL TWO-SAMPLE PROCEDURE MIGHT BE
C                             APPROPRIATE.
C
C   REQD. IMSL ROUTINES - H32/GGUBS,MDBETA,MDTD,MERRC=ERFC,
C                           MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/GGUBS,MDBETA,MDTD,MERRC=ERFC,
C                           MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE CLOSER XSTAT(10) AND XSTAT(11) ARE TO UNITY, THE
C                MORE ONE WOULD BE INCLINED TO REJECT THE HYPOTHESIS
C                OF RANDOMNESS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NMKTS  (X,XM,K,DSEED,XSTAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,IER
      REAL               X(1),XSTAT(17),XM(K)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBEG,IBEG2,IEND,I1,J,JER,N,KK,N1
      REAL               DARRAY(1),ONED18,DUM,RL,S1,S2,XJ,XK,XL,XX,YJ
      DATA               ONED18/.5555556E-01/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (K .EQ. 2) IER = 36
      N1 = K-1
      RL = 0.
      XK = 0.
      XL = 0.
      XJ = 0.
      S1 = 0.
      S2 = 0.
      XSTAT(12) = 0.
      DO 10 I = 1,N1
C                                  TERMINAL ERROR - SOME SAMPLE SIZE IS
C                                  LESS THAN OR EQUAL TO ZERO
         IF (XM(I) .LE. 0.) GO TO 40
C                                  FIND TOTAL NUMBER OF OBSERVATIONS
         XK = XK+XM(I)
         XX = XM(I)*XM(I)
         XL = XL+XX*(XM(I)+XM(I)+3.)
         XJ = XM(I)*XX*(6.*XX+15.*XM(I)+10.)+XJ
         I1 = I+1
         DO 5 J = I1,K
            RL = RL+XM(I)*XM(J)
    5    CONTINUE
   10 CONTINUE
      IF (XM(K) .LE. 0.) GO TO 40
      XK = XK+XM(K)
      XL = XL+XM(K)*XM(K)*(XM(K)+XM(K)+3.)
      XJ = XJ+XM(K)**3*(6.*XM(K)*XM(K)+15.*XM(K)+10.)
      XSTAT(9) = XK
      KK = XK
C                                  CALCULATE XSTAT(1) AND XSTAT(2)
      IEND = 0
      DO 30 I = 1,N1
         IBEG = IEND+1
         IEND = IEND+XM(I)
         IBEG2 = IEND+1
         DO 25 J = IBEG,IEND
            DO 20 N = IBEG2,KK
            IF (X(J) .GE. X(N)) GO TO 15
            S1 = S1+2.
            S2 = S2+2.
            GO TO 20
   15       IF (X(J) .NE. X(N)) GO TO 20
C                                  RANDOMIZE AND TALLY THE TIES
            XSTAT(12) = XSTAT(12)+1.
            CALL GGUBS (DSEED,1,DARRAY)
            DUM = DARRAY(1)
            IF (DUM .GT. .5) S1 = S1+2.
C                                  WARNING TIES EXIST
            IER = 35
   20       CONTINUE
   25    CONTINUE
   30 CONTINUE
      DUM = RL
      XSTAT(1) = S1-DUM
      XSTAT(2) = S2-DUM
      YJ = XK*XK
      S1 = YJ*(XK+XK+3.)-XL
C                                  CALCULATE EXPECTED MEAN OF XSTAT(1)
      XSTAT(7) = ONED18*S1
C                                  CALCULATE COEFFICIENTS OF RANK
C                                  CORRELATION
      XSTAT(10) = XSTAT(1)/DUM
      XSTAT(11) = XSTAT(2)/DUM
C                                  CALCULATE KURTOSIS
      XSTAT(8) = -1.44*(YJ*XK*(6.*YJ+15.*XK+10.)-XJ)/(S1*S1)
C                                  CALCULATE DEGREES OF FREEDOM
      XSTAT(17) = -3.*(2.+XSTAT(8))/XSTAT(8)
C                                  CALCULATE T-STATISTICS AND
C                                  PROBABILITIES
      DUM = (XSTAT(17)+1.)*XSTAT(7)
      XSTAT(3) = XSTAT(1)
      XSTAT(4) = XSTAT(2)
      XSTAT(5) = XSTAT(1)-1.
      XSTAT(6) = XSTAT(2)-1.
      DO 35 I = 3,6
         N = I+10
         XSTAT(N) = XSTAT(I)*SQRT(XSTAT(17)/(DUM-XSTAT(I)*XSTAT(I)))
         CALL MDTD(XSTAT(N),XSTAT(17),XSTAT(I),JER)
         IF (JER .NE. 0) IER = 34
         XSTAT(I) = .5*XSTAT(I)
   35 CONTINUE
      GO TO 50
   40 IER = 129
   50 IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HNMKTS )
 9005 RETURN
      END

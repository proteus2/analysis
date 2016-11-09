C   IMSL ROUTINE NAME   - NAFRE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - FRIEDMANS TEST FOR RANDOMIZED COMPLETE
C                           BLOCK DESIGNS
C
C   USAGE               - CALL NAFRE(Y,NB,NT,ALPHA,STAT,RJ,D,IWK,WK,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH NB*NT CONTAINING
C                           THE OBSERVATIONS.  THE FIRST NT POSITIONS
C                           OF Y CONTAIN THE OBSERVATIONS ON TREATMENTS
C                           1,2,...,NT IN THE FIRST BLOCK, THE SECOND
C                           NT POSITIONS CONTAIN THE CORRESPONDING
C                           OBSERVATIONS IN THE SECOND BLOCK, AND
C                           SO ON.
C                NB     - INPUT.  NUMBER OF BLOCKS.
C                NT     - INPUT.  NUMBER OF TREATMENTS.
C                ALPHA  - INPUT.  CRITICAL LEVEL FOR MULTIPLE
C                           COMPARISONS.  ALPHA SHOULD BE BETWEEN
C                           0.0 AND 1.0 EXCLUSIVE.  (IF THIS IS NOT
C                           THE CASE, NO ERROR IS GENERATED AND D IS
C                           NOT COMPUTED.)
C                STAT   - OUTPUT VECTOR OF LENGTH 6 CONTAINING TEST
C                           STATISTICS AND CORRESPONDING APPROXIMATE
C                           PROBABILITIES OF OBSERVING THESE VALUES
C                           OR GREATER VALUES IF THE APPROPRIATE NULL
C                           HYPOTHESES ARE TRUE.
C                         STAT(1) CONTAINS THE FRIEDMAN TEST STATISTIC
C                           FOR A TWO-SIDED ALTERNATIVE HYPOTHESIS.
C                         STAT(2) CONTAINS THE APPROXIMATE F STATISTIC
C                           FOR A TWO-SIDED ALTERNATIVE HYPOTHESIS.
C                         STAT(3) CONTAINS THE PAGE TEST STATISTIC FOR
C                           THE ORDERED ALTERNATIVE HYPOTHESIS THAT
C                           THE MEDIAN OF TREATMENT I IS LESS THAN
C                           OR EQUAL TO THE MEDIAN OF TREATMENT I+1
C                           AND THAT AT LEAST ONE MEDIAN IS STRICTLY
C                           LESS THAN SOME OTHER ONE.
C                         STAT(4) CONTAINS THE APPROXIMATE PROBABILITY
C                           CORRESPONDING TO STAT(1).  (CHI-SQUARED
C                           APPROXIMATION)
C                         STAT(5) CONTAINS THE APPROXIMATE PROBABILITY
C                           CORRESPONDING TO STAT(2).  (F APPROXI-
C                           MATION)
C                         STAT(6) CONTAINS THE APPROXIMATE PROBABILITY
C                           CORRESPONDING TO STAT(3).  (NORMAL APPROXI-
C                           MATION)
C                RJ     - OUTPUT VECTOR OF LENGTH NT CONTAINING THE
C                           SUM OF THE RANKS OF EACH TREATMENT.
C                D      - OUTPUT MINIMUM ABSOLUTE DIFFERENCE IN TWO
C                           ELEMENTS OF RJ TO INFER AT THE ALPHA LEVEL
C                           THAT THE MEDIANS OF THE CORRESPONDING
C                           TREATMENTS ARE DIFFERENT.
C                IWK    - WORK VECTOR OF LENGTH NT.
C                WK     - WORK VECTOR OF LENGTH 2*NT.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NB OR NT WAS LESS
C                             THAN 2.
C                         WARNING ERROR
C                           IER=33 INDICATES THAT THERE WAS ONE OR MORE
C                             TIES IN THE RANKS OF TREATMENTS WITHIN
C                             THE BLOCKS.
C                           IER=34 INDICATES THAT THE RANKS OF THE
C                             TREATMENTS WERE EXACTLY THE SAME IN ALL
C                             THE BLOCKS.  STAT(1) AND STAT(2) ARE SET
C                             TO MACHINE INFINITY.  STAT(4) AND STAT(5)
C                             ARE COMPUTED AS THE 1-NB POWER OF NT
C                             FACTORIAL, UNLESS THIS NUMBER IS TOO SMALL
C                             FOR SOME COMPUTATIONS, IN WHICH CASE
C                             STAT(4) AND STAT(5) ARE SET TO ZERO.
C                           IER=35 INDICATES THAT BOTH THE CONDITIONS OF
C                             IER=33 AND IER=34 OCCURRED.
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MDCH,MDFDRE,MDNOR,MDNRIS,MDSTI,
C                           MERFI,MERRC=ERFC,MGAMAD=DGAMMA,
C                           MLGAMD=DLGAMA,UERTST,UGETIO,VSRTR
C                       - H36,H48,H60/MDBETA,MDCH,MDFDRE,MDNOR,MDNRIS,
C                           MDSTI,MERFI,MERRC=ERFC,MGAMA=GAMMA,
C                           MLGAMA=ALGAMA,UERTST,UGETIO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NAFRE  (Y,NB,NT,ALPHA,STAT,RJ,D,IWK,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NB,NT,IWK(NT),IER
      REAL               Y(1),ALPHA,STAT(6),RJ(NT),D,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IS,J,KER,NTP1
      REAL               A2,B2B,BBB,FDEM,FNUM,FOUR,ONE,ONE44,RMEAN,RNB,
     *                   RNBM1,RNBNT,RNBTM1,RNT,RNTM1,RNTP1,RNTSQ,S,T,
     *                   THREE,TWLF,TWO,VAR,XETA,XINF,XJ,XMEXE,YNOR,ZERO
      INTEGER            V7I,V7J,V7JDON,V7JJ,V7JJJ,V7J2,V7K,V7K1,V7L,
     *                   V7N1,V7ISM1
      REAL               V7Y
      DATA               ZERO /0.0/,ONE /1.0/,TWO /2.0/,THREE
     *                   /3.0/,FOUR /4.0/,TWLF /12.0/,ONE44 /144.0/
      DATA               XINF/Z7FFFFFFF/
      DATA               XETA/Z00100000/
      DATA               XMEXE/-180.2182/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NB.GT.1 .AND. NT.GT.1) GO TO 5
      IER = 129
      GO TO 9000
    5 IS = 1
      RNB = NB
      RNT = NT
      RNBNT = RNB*RNT
      NTP1 = NT+1
      RNBM1 = RNB-ONE
      RNTM1 = RNT-ONE
      RNBTM1 = RNTM1*RNBM1
      RNTP1 = RNT+ONE
      RNTSQ = RNTP1*RNTP1
      DO 10 J=1,NT
         RJ(J) = ZERO
   10 CONTINUE
      A2 = ZERO
      DO 20 I=1,NB
C                                  ASSIGN RANKS TO EACH BLOCK
C                                  CODE SUBROUTINE NMRANK IN-LINE
         V7ISM1 = IS-1
         DO 1005 V7I = 1,NT
            WK(NT+V7I) = Y(V7ISM1+V7I)
 1005    IWK(V7I) = V7I
C                                  SORT ELEMENTS OF VECTOR R INTO
C                                  ASCENDING SEQUENCE SAVING
C                                  PERMUTATIONS
         CALL VSRTR (WK(NTP1),NT,IWK)
         S = 0.0
         T = 0.0
         V7N1 = NT-1
         V7L = 1
 1010    DO 1030 V7J = V7L,V7N1
            V7JJ = V7J
            V7Y = WK(NT+V7J)
            IF (ABS(V7Y-WK(NT+V7J+1)) .GT. 0.0) GO TO 1028
C                                  COUNT THE NUMBER OF TIES
            V7K = 1
            V7J2 = V7J+2
            IF (V7J2 .GT. NT) GO TO 1020
            DO 1015 V7I = V7J2,NT
               IF (ABS(V7Y-WK(NT+V7I)) .GT. 0.0) GO TO 1020
 1015       V7K = V7K+1
 1020       V7Y = V7J+.5*V7K
            V7K1 = V7K+1
            DO 1025 V7I = 1,V7K1
               V7JJ = V7J+V7I-1
               V7JJJ = IWK(V7JJ)
 1025       WK(V7JJJ) = V7Y
            V7I = V7K*(V7K+1)
            S = S+V7I
            T = T+V7I*(V7K+2)
            GO TO 1035
 1028       V7JDON = IWK(V7J)
 1030    WK(V7JDON) = V7J
 1035    V7L = V7JJ+1
         IF (V7L .LE. V7N1) GO TO 1010
         IF (V7L .NE. NT) GO TO 1040
         V7JDON = IWK(NT)
         WK(V7JDON) = NT
 1040    CONTINUE
         IF (S.GT.ZERO) IER = 33
         DO 15 J=1,NT
C                                  COMPUTE SUM OF WITHIN-BLOCKS RANKS
            RJ(J) = RJ(J)+WK(J)
            A2 = A2+WK(J)*WK(J)
   15    CONTINUE
         IS = IS+NT
   20 CONTINUE
      B2B = ZERO
      STAT(3) = ZERO
      DO 25 J=1,NT
         XJ = J
         B2B = B2B+RJ(J)*RJ(J)
C                                  COMPUTE PAGE TEST STATISTIC
         STAT(3) = STAT(3)+XJ*RJ(J)
   25 CONTINUE
      IF (ALPHA.LE.ZERO .OR. ALPHA.GE.ONE) GO TO 30
C                                  COMPUTE D
      CALL MDSTI(ALPHA,RNBTM1,T,KER)
      D = TWO*RNB*(A2-B2B/RNB)/(RNBM1*RNTM1)
      D = T*SQRT(D)
   30 BBB = RNB*RNT*RNTP1*RNTP1/FOUR
C                                  COMPUTE F APPROXIMATION
      FDEM = A2-B2B/RNB
      IF (FDEM.LE.ZERO) GO TO 45
      FNUM = RNBM1*(B2B/RNB-BBB)
      IF (FDEM.GT.ONE) GO TO 35
      IF (FNUM.GE.XINF*FDEM) GO TO 45
      GO TO 40
   35 IF (FNUM.GT.XETA*FDEM) GO TO 40
      STAT(2) = ZERO
      STAT(5) = ONE
      STAT(1) = ZERO
      STAT(4) = ONE
      GO TO 65
   40 STAT(2) = RNBM1*(B2B/RNB-BBB)/(A2-B2B/RNB)
C                                  OBTAIN APPROXIMATE PROBABILITY
C                                  CORRESPONDING TO STAT(2)
      CALL MDFDRE(STAT(2),RNTM1,RNBTM1,S,KER)
      STAT(5) = ONE-S
C                                  COMPUTE FRIEDMAN TEST STATISTIC
      STAT(1) = RNB*RNTM1*STAT(2)/(STAT(2)+RNBM1)
C                                  OBTAIN APPROXIMATE PROBABILITY
C                                  CORRESPONDING TO STAT(1)
      CALL MDCH(STAT(1),RNTM1,S,KER)
      STAT(4) = ONE-S
      GO TO 65
C                                  HANDLE CASE OF PERFECT FIT
   45 STAT(2) = XINF
      IF (RNT.GT.30.) GO TO 55
      BBB = ONE
      DO 50 I=1,NT
   50 BBB = BBB*I
      B2B = -ALOG(BBB)*RNBM1
      IF (B2B.LT.XMEXE) GO TO 55
      BBB = EXP(B2B)
      GO TO 60
   55 BBB = ZERO
   60 STAT(5) = BBB
      STAT(1) = XINF
      STAT(4) = STAT(5)
      IF (IER.EQ.33) IER = 35
      IF (IER.EQ.0) IER = 34
C                                  COMPUTE MEAN AND VARIANCE FOR PAGE
C                                  TEST STATISTIC
   65 RMEAN = RNBNT*RNTSQ/FOUR
      VAR = RNBNT*RNT*RNTSQ*RNTM1/ONE44
C                                  OBTAIN APPROXIMATE PROBABILITY
C                                  CORRESPONDING TO STAT(3)
      YNOR = (STAT(3)-RMEAN)/SQRT(VAR)
      CALL MDNOR(YNOR,S)
      STAT(6) = ONE-S
      IF (IER.EQ.0) GO TO 9005
 9000 CALL UERTST(IER,6HNAFRE )
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - GGDT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERAL DISCRETE DISTRIBUTION RANDOM DEVIATE
C                           GENERATOR USING TABLE LOOKUP METHOD
C
C   USAGE               - CALL GGDT (DSEED,NR,L1,IOPT,PF,CP,DEL,LCP,
C                                    L2,IR,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT.  A DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT.  THE NUMBER OF RANDOM DEVIATES TO BE
C                           GENERATED.
C                L1     - INPUT IF IOPT EQUALS 0.  INPUT/OUTPUT IF IOPT
C                           IS NOT EQUAL TO 0.  AS INPUT L1 IS THE LEFT
C                           HAND LIMIT OF THE DISCRETE DISTRIBUTION.
C                           IF IOPT .NE. 0, ON OUTPUT L1 WILL BE THE
C                           MINIMUM INTEGER (NOT LESS THAN THE INPUT
C                           VALUE OF L1) FOR WHICH THE CUMULATIVE
C                           DISTRIBUTION FUNCTION IS GREATER THAN DEL.
C                IOPT   - INPUT.  OPTION PARAMETER. IOPT=0 INDICATES
C                           THAT THE VECTOR CP ON INPUT CONTAINS THE
C                           CUMULATIVE PROBABILITIES, AND THE FUNCTION
C                           IN EXTERNAL SUBROUTINE PF WILL NOT BE USED.
C                           OTHERWISE GGDT CALLS THE USER-WRITTEN SUB-
C                           ROUTINE PF TO COMPUTE THE VALUES FOR CP.
C                PF     - INPUT.  A SUBROUTINE SUPPLIED BY THE USER TO
C                           COMPUTE THE PROBABILITY FUNCTION FOR THE
C                           DISCRETE DISTRIBUTION.  GGDT CALLS PF AS
C                           FOLLOWS
C                                      CALL  PF (I,P)
C                           WHERE P IS THE PROBABILITY THAT THE RANDOM
C                           VARIABLE IS EQUAL TO I.  PF MUST BE IN AN
C                           EXTERNAL STATEMENT IN THE CALLING PROGRAM.
C                           HOWEVER, PF IS NOT USED IF IOPT EQUALS ZERO.
C                CP     - VECTOR OF LENGTH LCP.  IF IOPT=0, CP IS INPUT
C                           AND CONTAINS THE CUMULATIVE PROBABILITIES
C                           (THAT IS, CP(I) EQUALS THE PROBABILITY THAT
C                           THE RANDOM VARIABLE IS LESS THAN OR EQUAL TO
C                           L1+I-1).  THE ELEMENTS OF CP SHOULD BE NON-
C                           DECREASING.  CP(LCP) SHOULD BE ONE.  IF IOPT
C                           IS NOT EQUAL TO 0,  CP IS COMPUTED FROM PF.
C                           L1 IS POSSIBLY INCREMENTED IN SUCH A FASHION
C                           THAT THE CUMULATIVE PROBABILITY ASSOCIATED
C                           WITH L IS GREATER THAN DEL. L2 IS DETERMINED
C                           SO THAT CP(L2) IS GREATER THAN (1.0 - DEL).
C                           IF (L2 - L1 + 1) IS GREATER THAN LCP, A
C                           TERMINAL ERROR IS GENERATED.
C                DEL    - INPUT.  MAXIMUM ERROR ALLOWED IN COMPUTING THE
C                           CUMULATIVE DISTRIBUTION FUNCTION FOR CP.
C                           DEL IS A (SMALL) POSITIVE NUMBER.
C                           IF IOPT = 0, DEL IS NOT USED.
C                LCP    - INPUT.  THE LENGTH OF CP.
C                L2     - OUTPUT IF IOPT IS NOT EQUAL TO ZERO.  L2 IS
C                           THE MAXIMUM INTEGER DETERMINED TO HAVE NON-
C                           ZERO PROBABILITY.  L2 IS COMPUTED CORRECTLY
C                           UNLESS IER = 130.  IF IOPT EQUALS ZERO, L2
C                           IS NOT USED.
C                IR     - OUTPUT.  VECTOR OF LENGTH NR CONTAINING
C                           THE RANDOM DEVIATES GENERATED.
C                IER    - OUTPUT.  ERROR PARAMETER.  (THESE ERRORS WILL
C                           NOT OCCUR IF IOPT = 0.)
C                         WARNING ERROR (WITH FIX).
C                           IER = 65  INDICATES THAT FOR SOME I,
C                             CP(I) WAS COMPUTED TO BE LESS THAN
C                             (1.0 - DEL), AND YET (CP(I+1) - 1.0) WAS
C                             GREATER THAN (1.0 - CP(I)).  IN THIS CASE,
C                             L2 IS CHOSEN SO THAT L2-L1+1 = I.
C                         TERMINAL ERROR.
C                           IER = 129  INDICATES THAT THE CUMULATIVE
C                             DISTRIBUTION FUNCTION VALUE IN CP(LCP)
C                             WAS NOT SUFFICIENTLY CLOSE TO ONE.  CP
C                             NEEDS TO BE ASSIGNED MORE STORAGE AND
C                             LCP SHOULD BE SET TO  L2-L1+1 .
C                           IER = 130  INDICATES THAT THE SAME ERROR
C                             OCCURED AS INDICATED BY IER=129 AND THAT
C                             L2 WAS NOT COMPUTED CORRECTLY BECAUSE THE
C                             COMPUTED CUMULATIVE PROBABILITIES DID NOT
C                             INCREASE FOR TEN SUCCESSIVE CALLS TO PF.
C
C   REQD. IMSL ROUTINES - GGUBFS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGDT (DSEED,NR,L1,IOPT,PF,CP,DEL,LCP,L2,IR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,L1,IOPT,LCP,L2,IR(NR),IER
      REAL               CP(LCP),DEL
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IT,K,J,LL,L1M1,L1M2,N,NLCP,NN(129)
      REAL               ANTDEL,CPK1,CPN,DROP,FR,RF,RLCP,ST,TOC,U,XLCP
      REAL               ZERO,ONE,TWO
      DOUBLE PRECISION   DSUM
      DATA               ZERO/0.0/,ONE/1.0/,TWO/2.0/,NN(1)/1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NLCP = LCP
      XLCP = LCP
C                                  CHECK WHETHER OR NOT CP IS PRE-SET
      IF (IOPT.EQ.0) GO TO 45
C                                  FIND L1 FOR WHICH P IS FIRST .GT. DEL
      DROP = ZERO
    5 CALL PF(L1,P)
      DROP = DROP + P
      L1 = L1 + 1
      IF (DROP.LT.DEL) GO TO 5
C                                  FILL CUMULATIVE PROBABILTY VECTOR CP
      CP(1) = P
      DSUM = DBLE(P)
      ANTDEL = ONE - (DEL+DROP-P)
      L1M2 = L1 - 2
      DO 10 K=2,LCP
         L2 = L1M2 + K
         CALL PF(L2,P)
         DSUM = DSUM + DBLE(P)
         CP(K) = DSUM
         IF (CP(K).GE.ANTDEL) GO TO 30
   10 CONTINUE
C                                  SET TERMINAL ERROR FLAG AND FIND L2
      IER = 129
      CPK = CP(LCP)
      IT = 0
   15 L2 = L2 + 1
      CALL PF(L2,P)
      CPK1 = CPK + P
      IF (CPK1.GT.CPK) GO TO 20
      IT = IT + 1
      IF (IT.LT.10) GO TO 25
      IER = 130
      GO TO 9000
   20 IT = 0
   25 CPK = CPK1
      IF (CPK.LT.ANTDEL) GO TO 15
      GO TO 9000
C                                  IF NECESSARY, SET WARNING AND FIX
   30 IF (CP(K).LT.TWO-CP(K-1)) GO TO 35
      IER = 65
      L2 = L2 - 1
C                                  IF NECESSARY, SCALE CP VECTOR
   35 NLCP = L2 - L1M2
      XLCP = NLCP
      IF (CP(NLCP).GE.ONE) GO TO 45
      RLCP = ONE/CP(NLCP)
      DO 40 I=1,NLCP
   40 CP(I) = CP(I)*RLCP
   45 CP(NLCP) = ONE
C                                  FILL SEARCH STEPS VECTOR NN
      L1M1 = L1 - 1
      J = 1
      LL = ALOG(XLCP) - ONE
      LL = MIN0(LL,7)
      FR = TWO**(-LL)
      RF = TWO**LL
      ST = FR
      DO 55 N=1,NLCP
         CPN = CP(N)
         IF(CPN.LT.ST) GO TO 55
   50    J = J + 1
         NN(J) = N
         ST = ST + FR
         IF(CPN.GE.ST) GO TO 50
   55 CONTINUE
C                                  SEARCH CP USING NN AND FILL IR
      DO 65 I=1,NR
         U = GGUBFS(DSEED)
         J = ONE + U*RF
         N = NN(J) - 1
   60    N = N + 1
         IF (U.GT.CP(N)) GO TO 60
   65 IR(I) = N + L1M1
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GGDT  ')
 9005 RETURN
      END

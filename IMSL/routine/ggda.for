C   IMSL ROUTINE NAME   - GGDA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - GENERAL DISCRETE DISTRIBUTION RANDOM DEVIATE
C                           GENERATOR USING ALIAS METHOD
C
C   USAGE               - CALL GGDA (DSEED,NR,NDMP,P,IA,WK,IR)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF RANDOM DEVIATES TO BE
C                           GENERATED.
C                NDMP   - INPUT NUMBER OF DISCRETE MASS POINTS IN THE
C                           DISTRIBUTION.
C                P      - INPUT VECTOR OF LENGTH NDMP CONTAINING THE
C                           PROBABILITIES ASSOCIATED WITH THE INDIVIDUAL
C                           MASS POINTS.  THE ELEMENTS OF P MUST BE NON-
C                           NEGATIVE AND MUST SUM TO ONE.
C                IA     - INPUT AND WORK VECTOR OF LENGTH NDMP.THE FIRST
C                           ELEMENT OF IA IS INPUT USED AS AN INITIAL
C                           ENTRY FLAG.THE FIRST TIME GGDA IS CALLED FOR
C                           A SPECIFIC PROBABILITY DISTRIBUTION IA(1)
C                           MUST BE SET TO NEGATIVE ONE. IA SHOULD NOT
C                           BE CHANGED BETWEEN SUCCESSIVE CALLS TO GGDA
C                           FOR A SPECIFIC PROBABILITY DISTRIBUTION.
C                WK     - REAL WORK VECTOR OF LENGTH NDMP. WK SHOULD NOT
C                           BE CHANGED BETWEEN SUCCESSIVE CALLS TO GGDA
C                           FOR A SPECIFIC PROBABILITY DISTRIBUTION.
C                IR     - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           RANDOM DEVIATES GENERATED.
C
C   REQD. IMSL ROUTINES - GGUBFS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGDA (DSEED,NR,NDMP,P,IA,WK,IR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,NDMP,IR(NR),IA(NDMP)
      REAL               P(NDMP),WK(NDMP)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,M
      REAL               A,ONE,RNDMP,U
      DATA               ONE /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      RNDMP = NDMP
C                                  CHECK FIRST CALL INDICATOR.
      IF (IA(1).NE.-1) GO TO 40
C                                  SET UP WORK AND ALIAS VECTORS.
      A = ONE/RNDMP
      K = 0
      M = 0
      DO 5 I=1,NDMP
    5 WK(I) = P(I)*RNDMP
   10 K = K+1
      IF (K.GT.NDMP) GO TO 30
      IF (P(K).LE.A) GO TO 10
   15 M = M+1
      IF (M.GT.NDMP) GO TO 30
      IF (P(M).GE.A) GO TO 15
      J = M
   20 IA(J) = K
      WK(K) = WK(K)-(ONE-WK(J))
      IF (WK(K).GE.ONE) GO TO 15
      J = K
   25 K = K+1
      IF (K.GT.NDMP) GO TO 30
      IF (P(K).LE.A) GO TO 25
      GO TO 20
   30 DO 35 I=1,NDMP
   35 WK(I) = WK(I)+I-1
C                                  GENERATE RANDOM DEVIATES
   40 DO 50 J=1,NR
         U = RNDMP*GGUBFS(DSEED)
         I = U+ONE
         IF (U.GT.WK(I)) GO TO 45
         IR(J) = I
         GO TO 50
   45    IR(J) = IA(I)
   50 CONTINUE
      RETURN
      END

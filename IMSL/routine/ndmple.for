C   IMSL ROUTINE NAME   - NDMPLE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NONPARAMETRIC PROBABILITY DENSITY FUNCTION
C                           (ONE DIMENSIONAL) ESTIMATION BY THE
C                           PENALIZED LIKELIHOOD METHOD
C
C   USAGE               - CALL NDMPLE (X,N,MMAX,IND,ALPHA,EPS,ITMAX,
C                           F,B,XJAC,IXJAC,ILOHI,DELC,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE RANDOM
C                           SAMPLE.  THE SAMPLE MUST BE SORTED IN
C                           ASCENDING ORDER.  ALL OF THE SAMPLE POINTS
C                           MUST FALL IN THE INTERVAL (B(1),B(2)).
C                N      - INPUT SIZE OF RANDOM SAMPLE.
C                MMAX   - INPUT NUMBER OF MESH NODES FOR THE DISCRETE
C                           PDF ESTIMATE.  MMAX MUST BE AN ODD INTEGER
C                           GREATER THAN 4.
C                IND    - INPUT OPTION PARAMETER.
C                         IND EQUAL TO ZERO IMPLIES THAT
C                           NDMPLE CALCULATES VECTOR F BY A BOOTSTRAP
C                           PROCEDURE.
C                         IND NONZERO IMPLIES INITIAL
C                           ESTIMATES CONTAINED IN VECTOR F ON INPUT.
C                ALPHA  - INPUT POSITIVE PENALTY WEIGHING FACTOR WHICH
C                           CONTROLS HOW SMOOTH THE ESTIMATE IS.
C                           SEE ALGORITHM SECTION IN THE MANUAL DOCUMENT
C                           FOR FURTHER DETAILS.
C                EPS    - INPUT POSITIVE CONVERGENCE CRITERION.
C                           NEWTONS METHOD ENDS WHEN THE EUCLIDEAN
C                           NORM OF THE CHANGE VECTOR IS LESS THAN EPS.
C                           EPS = 10**-5 IS SUGGESTED.
C                ITMAX  - INPUT MAXIMUM NUMBER OF ITERATIONS
C                           ALLOWED FOR NEWTONS METHOD. ITMAX=30 IS
C                           SUGGESTED.  ON OUTPUT ITMAX CONTAINS THE
C                           NUMBER OF ITERATIONS USED.
C                F      - INPUT/OUTPUT VECTOR OF LENGTH MMAX+2
C                           CONTAINING THE VALUES OF THE DISCRETE PDF
C                           ESTIMATE AT THE MMAX EQUALLY SPACED MESH
C                           NODES.  F IS AN INPUT VECTOR IF IND IS
C                           NONZERO, IN WHICH CASE F(1) AND F(MMAX)
C                           MUST ALWAYS BE ZERO.  ON OUTPUT F CONTAINS
C                           THE ESTIMATED VALUES AT THE NODES.
C                         THE LAST TWO ELEMENTS OF F ARE WORK STORAGE.
C                B      - INPUT/OUTPUT VECTOR OF LENGTH MMAX+2.
C                           ON INPUT, THE INTERVAL OF SUPPORT FOR THE
C                           ESTIMATE IS CONTAINED IN (B(1),B(2)).
C                           B IS DESTROYED ON OUTPUT AND REPLACED BY
C                           MESH INFORMATION. SEE REMARKS.
C                         THE REMAINDER OF B IS USED AS WORK STORAGE.
C                XJAC   - WORK MATRIX OF DIMENSION (MMAX-2) BY 5,
C                           CONTAINING THE JACOBIAN IN BAND MATRIX FORM.
C                IXJAC  - INPUT ROW DIMENSION OF THE MATRIX XJAC EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                ILOHI  - WORK VECTOR OF LENGTH 2*(MMAX+2).
C                DELC   - OUTPUT VECTOR OF LENGTH 6*(MMAX-2).
C                           ON OUTPUT, DELC(1) AND DELC(2) CONTAIN THE
C                           LOG LIKELIHOOD AND LOG PENALTY TERMS,
C                           RESPECTIVELY.  DELC(3) AND DELC(4) CONTAIN
C                           THE EXACT MEAN AND VARIANCE OF THE
C                           ESTIMATE, RESPECTIVELY.
C                         THE REMAINDER OF DELC IS USED AS WORK STORAGE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THE JACOBIAN WAS
C                             ALGORITHMICALLY SINGULAR.  TRY DIFFERENT
C                             INITIAL F ESTIMATES OR TRY IND=0.
C                           IER=130 INDICATES THAT ONE OF THE INPUT
C                             PARAMETERS MMAX,ALPHA,EPS,B(1),B(2), OR X
C                             WAS INCORRECTLY SPECIFIED.
C                           IER=131 INDICATES ITMAX WAS EXCEEDED.
C
C
C   REQD. IMSL ROUTINES - IQHSCU,LEQT1B,NDEST,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF IND IS NONZERO ON INPUT, THEN NONNEGATIVE INITIAL
C                GUESSES FOR F(I), I=1,...,MMAX ARE ALSO INPUT. F(1)
C                AND F(MMAX) MUST BE ZERO. NO OTHER VALUES OF F(I)
C                SHOULD BE ZERO. NEWTONS METHOD NEVER CHANGES VALUES
C                THAT ARE INITIALLY ZERO (THUS PROVIDING A FALSE
C                CONVERGENCE).
C            2.  IF IND IS ZERO ON INPUT, THEN A BOOTSTRAP ALGORITHM
C                IS USED TO PROVIDE INITIAL GUESSES FOR F(I),
C                I=1,...,MMAX. THE SMALLEST PROBLEM MMAX=5 IS SOLVED.
C                THIS ESTIMATE PROVIDES VALUES FOR MMAX=9 WHICH
C                PROVIDES VALUES FOR MMAX=17 AND SO ON.
C            3.  ON OUTPUT, INSPECT THE LOG LIKELIHOOD OF THE
C                ESTIMATE IN DELC(1). IF A FALSE CONVERGENCE OCCURRED
C                AS MENTIONED IN REMARK 1, THEN THE ESTIMATED DENSITY
C                WILL BE ZERO AT A SAMPLE DATA POINT, SO THAT THE LOG
C                LIKELIHOOD WILL BE A LARGE NEGATIVE NUMBER (CLOSE TO
C                NEGATIVE MACHINE INFINITY). OFTEN IMSL ROUTINE
C                LEQT1B RETURNS IER=129 AS XJAC BECOMES SINGULAR IN
C                THIS SITUATION.
C            4.  IN CASES OF INCORRECT OPERATION OF NDMPLE, CONSIDER
C                THE FOLLOWING POINTS
C                (A)  X(I) MUST BE LESS THAN OR EQUAL TO X(I+1) FOR
C                     I=1,...,N-1.
C                (B)  MMAX MUST BE GREATER THAN OR EQUAL TO FIVE AND
C                     ODD.
C                (C)  B(1) MUST BE LESS THAN X(1) AND
C                     X(N) MUST BE LESS THAN B(2)
C                (D)  ALPHA MUST BE GREATER THAN ZERO.
C            5.  UNDERFLOWS MAY OCCUR IN IMSL ROUTINE LEQT1B AND
C                SHOULD BE IGNORED.
C            6.  USING IND NONZERO WITH POOR INITIAL GUESSES IN F
C                MAY RESULT IN DIVISIONS BY ZERO AND FAILURE OF
C                NDMPLE. GENERALLY, IND EQUAL TO ZERO IS RECOMMENDED.
C                HOWEVER, WHEN THE EXPERIMENTER IS TRYING SEVERAL
C                ALPHA VALUES, HE WILL HAVE GOOD INITIAL GUESSES FOR
C                F FROM PREVIOUS RUNS AND MAY WISH TO USE IND NONZERO
C                FOR EFFICIENCY REASONS.
C            7.  THE DENSITY ESTIMATES AT THE MESH NODES ARE RETURNED
C                IN F. THE CORRESPONDING MESH NODE VALUES ARE B(1),
C                B(1)+D, B(1)+2.0*D, ..., B(2) WHERE D IS
C                (B(2)-B(1))/(MMAX-1).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NDMPLE (X,N,MMAX,IND,ALPHA,EPS,ITMAX,F,B,
     1                   XJAC,IXJAC,ILOHI,DELC,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MMAX,IND,ITMAX,IXJAC,ILOHI(1),IER
      REAL               X(N),ALPHA,EPS,XJAC(IXJAC,5),DELC(1),B(1),F(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            MINCR(9),I1,I2,I3,I4,I5,I,IH,IL,IMPTR,IPTR,
     1                   ITER,K,KM1,KM2,M,MM1,MM2,MM4,MM4T2,MM4T31,
     2                   MM7,MOLD,IOPT
      REAL               BBIG,BDEL,BK,BKM1,BLEFT,BRIGHT,BSMALL,
     1                   CK,CKMCM1,CKM1,CKM2,CKP1,CKP2,CONS,DELTA1,
     2                   FACTOR,FK,FKM1,FKM2,FKP1,H,HINV,H2,H3,SUM,
     3                   XNORM,TEMP,WK(1),BEND(2)
      REAL               FOUR,HALF,P10M4,ONE,SIX,THREE,TWELVE,TWO,ZERO
      DOUBLE PRECISION   SUM1,SUM2,SUM3
      DATA               FOUR,HALF,P10M4,ONE,SIX,THREE,TWELVE,TWO,ZERO/
     1                   4.,.5,1.E-4,1.,6.,3.,12.,2.,0./
      DATA               MINCR/3,5,9,17,33,65,129,253,100001/
      DATA               I1/1/,I2/2/,I3/3/,I4/4/,I5/5/,IOPT/1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 130
C                                  ERROR CHECKS
      IF (MMAX.LE.4.OR.(MMAX/2)*2.EQ.MMAX.OR.ALPHA.LE.ZERO.OR.EPS.LE.ZER
     1O.OR.N.LE.0.OR.B(I1).GE.X(1).OR.X(N).GE.B(I2)) GO TO 9000
      IF (N.EQ.1) GO TO 10
      DO 5 I=2,N
         IF (X(I).LT.X(I-1)) GO TO 9000
    5 CONTINUE
   10 BLEFT = B(I1)
      BRIGHT = B(I2)
      BEND(1) = BLEFT
      BEND(2) = BRIGHT
      IMPTR = 1
      BDEL = BRIGHT-BLEFT
      M = MINCR(IMPTR)+2
      IF (IND.NE.0) GO TO 15
      F(I1) = ZERO
      F(I2) = ZERO
      F(I3) = TWO/BDEL
      F(I4) = ZERO
      F(I5) = ZERO
C                                  NOTICE EXTRA ZERO VALUES AROUND
C                                  F VECTOR -- (M = MMAX+2)
C                                  INTERNALLY F = SQUARE ROOT OF F
C                                  REFINE MESH -- USE OLD ESTIMATES
   15 IF (MMAX.LE.MINCR(IMPTR)) GO TO 135
      MOLD = M-2
      IMPTR = IMPTR+1
      M = MIN0(MMAX,MINCR(IMPTR))+2
      IF (IND.NE.0) M = MMAX+2
      MM1 = M-1
      MM2 = M-2
      MM4 = M-4
      MM7 = M-7
      MM4T2 = MM4+MM4
      MM4T31 = MM4T2+MM4+1
C                                  CALCULATE MESH INTERVAL WIDTH
      H = BDEL/(M-3)
      HINV = ONE/H
      H2 = H*H
      H3 = H2*H
      DO 20 I=1,M
         B(I) = BLEFT+(I-2)*H
         ILOHI(I) = 0
         ILOHI(I+M) = -1
   20 CONTINUE
C                                  CALCULATE ARRAYS FOR THE MAPPING
C                                  * THAT ARE I-I WHERE *(I) IS
C                                  SUCH THAT B(*(I)).LE.X(I).LT.
C                                  B(*(I)+1),I=1,N*(I) = K FOR I =
C                                  ILO(K) TO IHI(K) UNLESS ILO(K).GT.
C                                  ILI(K)
      IPTR = 1
      DO 30 K=2,M
         KM1 = K-1
         I = KM1+M
         ILOHI(KM1) = IPTR
         ILOHI(I) = IPTR-1
   25    IF (X(IPTR).GE.B(K)) GO TO 30
         ILOHI(I) = ILOHI(I)+1
         IPTR = IPTR+1
         IF (IPTR.GT.N) GO TO 35
         GO TO 25
   30 CONTINUE
   35 FACTOR = TWO*ALPHA/H3
C                                  INITIALIZE MESH NODE VALUES
      IF (IND.NE.0) GO TO 45
      CALL NDEST (B(I3),MM4,F(I2),MOLD,IOPT,BEND,WK,DELC(I3),JER)
      TEMP = ONE/(M*M)
      DO 40 I=3,MM2
         F(I) = AMAX1(TEMP,SQRT(DELC(I)))
   40 CONTINUE
      GO TO 55
   45 DO 50 I=2,MM2
         F(M-I) = SQRT(F(MM1-I))
   50 CONTINUE
   55 F(MM1) = ZERO
      F(M) = ZERO
C                                  FORM NEWTON S EQUATION FOR ZEROES
      ITER = 0
   60 ITER = ITER+1
C
C                                  EVALUATE DELC.XJAC - LAGRANGIAN
      XJAC(1,1) = ZERO
      XJAC(1,2) = ZERO
      XJAC(2,1) = ZERO
      BSMALL = ZERO
      SUM = ZERO
C                                  CK.. ARE TRUE ESTIMATES = FK**2
      DO 95 K=3,MM2
         KM1 = K-1
         KM2 = K-2
         FK = F(K)
         FKM1 = F(KM1)
         FKM2 = F(KM2)
         CKM2 = FKM2**2
         CKM1 = FKM1**2
         CK = FK**2
         CKP1 = F(K+1)**2
         CKP2 = F(K+2)**2
         BK = B(K)
         BKM1 = B(KM1)
         SUM = SUM+CK
         IF (KM2.LT.3) GO TO 65
         TEMP = FOUR*FK*FKM2*FACTOR
         XJAC(KM2,1) = TEMP
   65    SUM1 = ZERO
         SUM2 = ZERO
         SUM3 = ZERO
         IL = ILOHI(K)
         IH = ILOHI(K+M)
         IF (IL.GT.IH) GO TO 75
         DO 70 I=IL,IH
            TEMP = (X(I)-BK)*HINV
            CONS = (ONE-TEMP)/(CK+(CKP1-CK)*TEMP)
            SUM1 = SUM1-CONS
            SUM2 = SUM2+CONS*CONS
   70    CONTINUE
   75    IL = ILOHI(KM1)
         IH = ILOHI(KM1+M)
         IF (IL.GT.IH) GO TO 85
         CKMCM1 = CK-CKM1
         DO 80 I=IL,IH
            CONS = (X(I)-BKM1)*HINV
            TEMP = CKM1+CKMCM1*CONS
            SUM1 = SUM1-CONS/TEMP
            TEMP = TEMP*TEMP
            SUM2 = SUM2+(CONS*CONS)/TEMP
            SUM3 = SUM3+CONS*(ONE-CONS)/TEMP
   80    CONTINUE
   85    TEMP = FACTOR*(CKM2+CKP2-FOUR*(CKM1+CKP1)+SIX*CK)+SUM1
         TEMP = TEMP+TEMP
         BSMALL = BSMALL+TWO*CK*TEMP
         XJAC(KM2,3) = TEMP+FOUR*CK*(SIX*FACTOR+SUM2)
         IF (KM2.EQ.1) GO TO 90
         XJAC(KM2,2) = FOUR*FK*FKM1*(-FOUR*FACTOR+SUM3)
   90    DELC(KM2) = FK*TEMP
         DELC(KM2+MM4) = -TWO*FK
   95 CONTINUE
      BSMALL = HINV-SUM+BSMALL
      BBIG = FOUR*SUM
C                                  SAVE PORTION OF DELC
      TEMP = -TWO*BSMALL/BBIG
      DO 100 K=1,MM4
         DELC(K+MM4T2) = DELC(K+MM4)
         XJAC(K,3) = XJAC(K,3)+TEMP
  100 CONTINUE
C                                  FILL OUT SYMMETRIC BAND STRUCTURE
      DO 105 I=1,MM4
         K = 1+MOD(I+MM7,MM4)
         XJAC(K,5) = XJAC(I,1)
         K = 1+MOD(K,MM4)
         XJAC(K,4) = XJAC(I,2)
  105 CONTINUE
C                                  SOLVE A SYMMETRIC BAND LINEAR SYSTEM
      CALL LEQT1B (XJAC,MM4,2,2,IXJAC,DELC,2,MM4,0,DELC(MM4T31),IER)
      IF (IER.NE.0) GO TO 9000
      CONS = ZERO
      TEMP = ZERO
      DO 110 I=1,MM4
         BK = DELC(I+MM4T2)
         TEMP = TEMP+BK*DELC(I)
         CONS = CONS+BK*DELC(MM4+I)
  110 CONTINUE
      DELTA1 = (HINV-SUM-TEMP)/CONS
      XNORM = ZERO
      DO 115 I=1,MM4
         TEMP = DELC(I)+DELTA1*DELC(MM4+I)
         F(I+2) = F(I+2)-TEMP
         XNORM = XNORM+TEMP*TEMP
  115 CONTINUE
      XNORM = SQRT(XNORM)
      IF (XNORM.LT.EPS) GO TO 125
      IER = 131
      IF (ITER.GE.ITMAX) GO TO 9000
C                                  AD HOC PROJECTION TO PLUS QUADRANT
      DO 120 I=3,MM2
         IF (F(I).GT.-P10M4) GO TO 120
         F(I) = P10M4
  120 CONTINUE
      GO TO 60
C                                  REPLACE F(*) WITH TRUE ESTIMATES
  125 IPTR = 2
      DO 130 I=3,MM2
         F(I) = F(I)**2
  130 CONTINUE
      GO TO 15
C                                  EVALUATE LOG LIKELIHOOD AND PENALTY
  135 IER = 0
      SUM1 = ZERO
      DO 140 K=2,MM1
         SUM1 = SUM1+(F(K-1)-TWO*F(K)+F(K+1))**2
  140 CONTINUE
      DELC(I2) = -HALF*FACTOR*SUM1
      SUM2 = ZERO
      DO 145 I=1,N
         CALL NDEST (X(I),I1,F(I2),MMAX,IOPT,BEND,WK,DELC(I5),JER)
         SUM2 = SUM2+ALOG(DELC(I5))
  145 CONTINUE
      DELC(I1) = SUM2
C                                  EVALUATE M.L.P.E. MEAN AND VARIANCE
      SUM1 = ZERO
      SUM2 = ZERO
      DO 150 K=2,MM2
         FK = F(K)
         FKP1 = F(K+1)
         BK = B(K)
         CONS = FK+FKP1
         TEMP = CONS+FKP1
         SUM1 = SUM1+H2*TEMP/SIX+HALF*H*BK*CONS
         SUM2 = SUM2+H3*(TEMP+FKP1)/TWELVE+H2*BK*TEMP/THREE+HALF*H*BK
     1   *BK*CONS
  150 CONTINUE
      DELC(3) = SUM1
      DELC(4) = SUM2-SUM1*SUM1
      DO 155 I=1,MMAX
         F(I) = F(I+1)
  155 CONTINUE
      ITMAX = ITER
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HNDMPLE)
 9005 RETURN
      END

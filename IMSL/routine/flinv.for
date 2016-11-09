C   IMSL ROUTINE NAME   - FLINV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - INVERSE LAPLACE TRANSFORM OF A USER
C                           SUPPLIED COMPLEX FUNCTION
C
C   USAGE               - CALL FLINV (F,N,T,ALPHA,NSIG,KMAX,FINV,IER)
C
C   ARGUMENTS    F      - A USER SUPPLIED COMPLEX FUNCTION, F(S),
C                           SPECIFYING THE LAPLACE TRANSFORM WHOSE IN-
C                           VERSE IS TO BE CALCULATED (INPUT).
C                           THE CALLING SEQUENCE OF THIS FUNCTION MUST
C                           BE OF THE FORM F(S) WHERE S IS A COMPLEX
C                           VARIABLE. F SHOULD BE TYPE EXTERNAL IN THE
C                           MAIN PROGRAM.
C                N      - THE NUMBER OF POINTS AT WHICH THE INVERSE
C                           LAPLACE TRANSFORM IS TO BE CALCULATED
C                           (INPUT).
C                T      - A VECTOR OF LENGTH N CONTAINING THE POINTS
C                           AT WHICH THE CALCULATION OF THE INVERSE
C                           LAPLACE TRANSFORM IS DESIRED (INPUT).
C                           T(I) MUST BE GREATER THAN ZERO FOR ALL I.
C                ALPHA  - THE MAXIMUM OF THE REAL PARTS OF THE
C                           SINGULARITIES OF F(S), OR AN ESTIMATED
C                           VALUE GREATER THAN THIS MAXIMUM (INPUT).
C                NSIG   - INTEGER VALUE SPECIFYING THE NUMBER OF
C                           SIGNIFICANT DIGITS DESIRED IN THE OUPUT
C                           VECTOR FINV (INPUT).
C                KMAX   - THE ALGORITHM IS PERMITTED TO USE 3*KMAX
C                           FUNCTION EVALUATIONS FOR EACH T(I) (INPUT).
C                FINV   - OUTPUT VECTOR OF LENGTH N. FINV(I) CONTAINS
C                           THE VALUE OF THE INVERSE LAPLACE TRANSFORM
C                           OF THE USER SUPPLIED FUNCTION AT T(I).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           TERMINAL ERROR
C                             IER = 129 INDICATES THAT THE ALGORITHM
C                               WAS NOT ABLE TO ACHIEVE THE ACCURACY
C                               REQUESTED WITHIN KMAX FUNCTION
C                               EVALUATIONS FOR SOME T(I).
C                             IER = 130 INDICATES THAT OVERFLOW WOULD
C                               HAVE OCCURRED FOR A PARTICULAR VALUE
C                               OF T.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FLINV (F,N,T,ALPHA,NSIG,KMAX,FINV,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSIG,KMAX,IER,MMMM
      REAL               T(N),ALPHA,FINV(N)
      COMPLEX            F
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II,ILOOP,ITER,ITMAX,I,IH,JMAX,JODD,J,KK,K,LL,
     *                   NCNT,NEWTRY,NT
      REAL               AL,AX,A,DEFINT,DELA,EREST,ERJ,ER,EXPAT,
     *                   HALF,ONE,PIK,PI,QEPALJ,QEPAL,RATIO,REPALI,
     *                   REPALJ,REPAL,RERR,RPEXE,SERR,T2,TCAP,TEN,ZERO
      COMPLEX            DIFF,FK,FSUM,S(99),SI,SVI,W,WK,ZTEMP,ZK,Z
      COMPLEX            TEMP1,TEMP2,EPSUM(99)
      DATA               PI/3.141593E0/
      DATA               ZERO/0.0/,ONE/1.0/,HALF/.5/
      DATA               RPEXE/174.0/
      DATA               TEN/10.0/,AX/0.1/
      DATA               DEFINT/2.92E0/
      DATA               JMAX/99/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      RERR = TEN**(-NSIG)
      KK = 0
      ITMAX = 8
C                                  BEGIN MAIN LOOP
      DO 90 I = 1, N
        NEWTRY = 0
        AL = ALPHA
        T2 = DEFINT*T(I)
        TCAP = T2*.5E0
        ITER = 1
    5   A = AL - ALOG(RERR)/T2
        RATIO = AMIN1(ONE,T(I)/TCAP)
        TEMP1 = CEXP(CMPLX(ZERO,PI*T(I)/TCAP))
        TEMP2 = CEXP(CMPLX(ZERO,-PI*T(I)/TCAP))
        Z = (TEMP1+TEMP2)/CMPLX(2.E0,ZERO)
        W = (TEMP1-TEMP2)/CMPLX(ZERO,2.E0)
   10   ZK = Z
        WK = W
        K = 1
        FSUM = HALF*F(CMPLX(A,ZERO))
        EXPAT = EXP(A*T(I))/TCAP
        JODD = 1
        J = 1
        S(1) = FSUM
        QEPAL = EXPAT*REAL(FSUM)
        QEPALJ = QEPAL
        PIK = PI
   15   FK = F(CMPLX(A,PIK/TCAP))
C                                  GENERATE PARTIAL SUM
        FSUM = FSUM + REAL(FK)*ZK - AIMAG(FK)*WK
        ZTEMP = ZK*Z - WK*W
        WK = WK*Z + ZK*W
        ZK = ZTEMP
        PIK = PIK + PI
        K = K + 1
        IF (J.EQ.1) REPALJ = EXPAT*REAL(FSUM)
C                                  BEGIN EPSILON ALGORITHM
        JODD = -JODD
        J = J + 1
        S(J) = FSUM
        SVI = CMPLX(ZERO,ZERO)
        II = J
        DO 20 ILOOP = 2, J
          II = II - 1
          DIFF = S(II+1) - S(II)
          SI = SVI
          IF (CABS(DIFF).NE.ZERO) SI = SI + CMPLX(ONE,ZERO)/DIFF
          SVI = S(II)
          S(II) = SI
   20   CONTINUE
        IF (JODD.LT.0) GO TO 30
        REPALJ = EXPAT*REAL(S(1))
        IF (J+1.GE.JMAX) GO TO 25
C                                  CONVERGENCE CHECK WITHIN ONE
C                                  SAWTOOTH OF EPSILON ALGORITHM
        ERJ = ABS(REPALJ-QEPALJ)/AMAX1(AX,ABS(REPALJ))
        QEPALJ = REPALJ
        IF (ERJ.GT.RERR*RATIO) GO TO 30
   25   J = 1
        S(1) = FSUM
        QEPALJ = EXPAT*REAL(FSUM)
        REPAL = REPALJ
C                                  CONVERGENCE CHECK BETWEEN TWO
C                                  SAWTEETH
        ER = ABS(REPAL-QEPAL)/AMAX1(AX,ABS(REPAL))
        IF (ER.LE.RATIO*RERR) GO TO 35
        QEPAL = REPAL
   30   IF (K.LT.KMAX) GO TO 15
        IER = 129
        GO TO 9000
   35   CONTINUE
        KK = KK + K
        IF (MOD(ITER,2).EQ.0) GO TO 45
        ITER = ITER + 1
        A = A + AMAX1(1.15/TCAP,0.1*ABS(A))
        IF (A*T(I).LE.RPEXE) GO TO 40
        IER = 130
        GO TO 9000
   40   REPALI = REPAL
        GO TO 10
   45   CONTINUE
        ITER = ITER + 1
C                                  CONVERGENCE CHECK BETWEEN RESULTS
C                                  FROM THE USE OF TWO SUCCESSIVE A
        EREST = ABS(REPAL-REPALI)/AMAX1(AX,ABS(REPAL),ABS(REPALI))
        IF (EREST.LE.RERR) GO TO 85
        DELA = 0.1*ABS(RERR/EREST)
        IF (ITER.GT.3) RATIO = RATIO*DELA
        IF (ITER.LE.ITMAX) GO TO 80
C                                  GET KMAX PARTIAL SUMS
        SERR = RERR*AX
        AL = ALPHA
        T2 = DEFINT*T(I)
        TCAP = T2*.5E0
        A = AL - ALOG(SERR)/T2
        TEMP1 = CEXP(CMPLX(ZERO,PI*T(I)/TCAP))
        TEMP2 = CEXP(CMPLX(ZERO,-PI*T(I)/TCAP))
        Z = (TEMP1+TEMP2)/CMPLX(2.E0,ZERO)
        W = (TEMP1-TEMP2)/CMPLX(ZERO,2.E0)
   50   ZK = Z
        WK = W
        FSUM = HALF*F(CMPLX(A,ZERO))
        EXPAT = EXP(A*T(I))/TCAP
        PIK = PI
        JMAX = MAX0(99,KMAX)
        IF (MOD(JMAX,2).EQ.0) JMAX = JMAX - 1
        IBEG = JMAX - 98
        II = 1
        LL = JMAX/99
        DO 55 NT = 1, JMAX
          FK = F(CMPLX(A,PIK/TCAP))
          FSUM = FSUM + REAL(FK)*ZK - AIMAG(FK)*WK
          ZTEMP = ZK*Z - WK*W
          WK = WK*Z + ZK*W
          ZK = ZTEMP
          PIK = PIK + PI
          NCNT = NT/LL
          IF (NCNT.GT.0) S(NCNT) = FSUM
   55   CONTINUE
C                                     INITIALIZE EPSUM TO ZERO AND SAVE
C                                       PSUMS
        DO 60 IH = 1, 99
          EPSUM(IH) = ZERO
   60   CONTINUE
C                                     GET ACCELLARATED SUM
        MMMM = -1
        DO 70 IEND = 2, 99
          DO 65 IH = 99, IEND, MMMM
            DIFF = S(IH) - S(IH-1)
            S(IH) = S(IH-1)
            IF (CABS(DIFF).NE.ZERO) S(IH) = ONE/DIFF + EPSUM(IH)
            EPSUM(IH) = S(IH-1)
   65     CONTINUE
   70   CONTINUE
        EPSUM(1) = S(99)
        IF (NEWTRY.EQ.0) REPAL = EXPAT*REAL(EPSUM(1))
        IF (NEWTRY.EQ.1) REPALJ = EXPAT*REAL(EPSUM(1))
        IF (NEWTRY.EQ.1) GO TO 75
        NEWTRY = NEWTRY + 1
        A = A + AMAX1(1.15/TCAP,0.1E0*ABS(A))
        GO TO 50
   75   EREST = ABS(REPAL-REPALJ)/AMAX1(AX,ABS(REPAL),ABS(REPALJ))
        IF (EREST.LE.RERR) GO TO 85
        IER = 129
        GO TO 9000
   80   AL = AL - ALOG(DELA)/T2
        GO TO 5
   85   FINV(I) = REPAL
   90 CONTINUE
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,'FLINV ')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - ZXSSQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - MINIMUM OF THE SUM OF SQUARES OF M FUNCTIONS
C                           IN N VARIABLES USING A FINITE DIFFERENCE
C                           LEVENBERG-MARQUARDT ALGORITHM
C
C   USAGE               - CALL ZXSSQ(FUNC,M,N,NSIG,EPS,DELTA,MAXFN,IOPT,
C                           PARM,X,SSQ,F,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
C
C   ARGUMENTS    FUNC   - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
C                           THE RESIDUAL VECTOR F(1),F(2),...,F(M) FOR
C                           GIVEN PARAMETER VALUES X(1),X(2),...,X(N).
C                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
C                             CALL FUNC(X,M,N,F)
C                             WHERE X IS A VECTOR OF LENGTH N AND F IS
C                               A VECTOR OF LENGTH M.
C                             FUNC MUST APPEAR IN AN EXTERNAL STATEMENT
C                               IN THE CALLING PROGRAM.
C                             FUNC MUST NOT ALTER THE VALUES OF
C                               X(I),I=1,...,N, M, OR N.
C                M      - THE NUMBER OF RESIDUALS OR OBSERVATIONS
C                           (INPUT)
C                N      - THE NUMBER OF UNKNOWN PARAMETERS (INPUT).
C                NSIG   - FIRST CONVERGENCE CRITERION. (INPUT)
C                           CONVERGENCE CONDITION SATISFIED IF ON TWO
C                           SUCCESSIVE ITERATIONS, THE PARAMETER
C                           ESTIMATES AGREE, COMPONENT BY COMPONENT,
C                           TO NSIG DIGITS.
C                EPS    - SECOND CONVERGENCE CRITERION. (INPUT)
C                           CONVERGENCE CONDITION SATISFIED IF, ON TWO
C                           SUCCESSIVE ITERATIONS THE RESIDUAL SUM
C                           OF SQUARES ESTIMATES HAVE RELATIVE
C                           DIFFERENCE LESS THAN OR EQUAL TO EPS. EPS
C                           MAY BE SET TO ZERO.
C                DELTA  - THIRD CONVERGENCE CRITERION. (INPUT)
C                           CONVERGENCE CONDITION SATISFIED IF THE
C                           (EUCLIDEAN) NORM OF THE APPROXIMATE
C                           GRADIENT IS LESS THAN OR EQUAL TO DELTA.
C                           DELTA MAY BE SET TO ZERO.
C                             NOTE, THE ITERATION IS TERMINATED, AND
C                             CONVERGENCE IS CONSIDERED ACHIEVED, IF
C                             ANY ONE OF THE THREE CONDITIONS IS
C                             SATISFIED.
C                MAXFN  - INPUT MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C                           (I.E., CALLS TO SUBROUTINE FUNC) ALLOWED.
C                           THE ACTUAL NUMBER OF CALLS TO FUNC MAY
C                           EXCEED MAXFN SLIGHTLY.
C                IOPT   - INPUT OPTIONS PARAMETER.
C                         IOPT=0 IMPLIES BROWN'S ALGORITHM WITHOUT
C                           STRICT DESCENT IS DESIRED.
C                         IOPT=1 IMPLIES STRICT DESCENT AND DEFAULT
C                           VALUES FOR INPUT VECTOR PARM ARE DESIRED.
C                         IOPT=2 IMPLIES STRICT DESCENT IS DESIRED WITH
C                           USER PARAMETER CHOICES IN INPUT VECTOR PARM.
C                PARM   - INPUT VECTOR OF LENGTH 4 USED ONLY FOR
C                         IOPT EQUAL TWO.  PARM(I) CONTAINS, WHEN
C                           I=1, THE INITIAL VALUE OF THE MARQUARDT
C                             PARAMETER USED TO SCALE THE DIAGONAL OF
C                             THE APPROXIMATE HESSIAN MATRIX, XJTJ,
C                             BY THE FACTOR (1.0 + PARM(1)).  A SMALL
C                             VALUE GIVES A NEWTON STEP, WHILE A LARGE
C                             VALUE GIVES A STEEPEST DESCENT STEP.
C                             THE DEFAULT VALUE FOR PARM(1) IS 0.01.
C                           I=2, THE SCALING FACTOR USED TO MODIFY THE
C                             MARQUARDT PARAMETER, WHICH IS DECREASED
C                             BY PARM(2) AFTER AN IMMEDIATELY SUCCESSFUL
C                             DESCENT DIRECTION, AND INCREASED BY THE
C                             SQUARE OF PARM(2) IF NOT.  PARM(2) MUST
C                             BE GREATER THAN ONE, AND TWO IS DEFAULT.
C                           I=3, AN UPPER BOUND FOR INCREASING THE
C                             MARQUARDT PARAMETER.  THE SEARCH FOR A
C                             DESCENT POINT IS ABANDONED IF PARM(3) IS
C                             EXCEEDED.  PARM(3) GREATER THAN 100.0 IS
C                             RECOMMENDED.  DEFAULT IS 120.0.
C                           I=4, VALUE FOR INDICATING WHEN CENTRAL
C                             RATHER THAN FORWARD DIFFERENCING IS TO BE
C                             USED FOR CALCULATING THE JACOBIAN.  THE
C                             SWITCH IS MADE WHEN THE NORM OF THE
C                             GRADIENT OF THE SUM OF SQUARES FUNCTION
C                             BECOMES SMALLER THAN PARM(4).  CENTRAL
C                             DIFFERENCING IS GOOD IN THE VICINITY
C                             OF THE SOLUTION, SO PARM(4) SHOULD BE
C                             SMALL.  THE DEFAULT VALUE IS 0.10.
C                X      - VECTOR OF LENGTH N CONTAINING PARAMETER
C                           VALUES.
C                         ON INPUT, X SHOULD CONTAIN THE INITIAL
C                           ESTIMATE OF THE LOCATION OF THE MINIMUM.
C                         ON OUTPUT, X CONTAINS THE FINAL ESTIMATE
C                           OF THE LOCATION OF THE MINIMUM.
C                SSQ    - OUTPUT SCALAR WHICH IS SET TO THE RESIDUAL
C                           SUMS OF SQUARES, F(1)**2+...+F(M)**2, FOR
C                           THE FINAL PARAMETER ESTIMATES.
C                F      - OUTPUT VECTOR OF LENGTH M CONTAINING THE
C                           RESIDUALS FOR THE FINAL PARAMETER ESTIMATES.
C                XJAC   - OUTPUT M BY N MATRIX CONTAINING THE
C                           APPROXIMATE JACOBIAN AT THE OUTPUT VECTOR X.
C                IXJAC  - INPUT ROW DIMENSION OF MATRIX XJAC EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                XJTJ   - OUTPUT VECTOR OF LENGTH (N+1)*N/2 CONTAINING
C                           THE N BY N MATRIX (XJAC-TRANSPOSED) * (XJAC)
C                           IN SYMMETRIC STORAGE MODE.
C                WORK   - WORK VECTOR OF LENGTH 5*N + 2*M + (N+1)*N/2.
C                         ON OUTPUT, WORK(I) CONTAINS FOR
C                           I=1, THE NORM OF THE GRADIENT DESCRIBED
C                             UNDER INPUT PARAMETERS DELTA AND PARM(4).
C                           I=2, THE NUMBER OF FUNCTION EVALUATIONS
C                             REQUIRED DURING THE WORK(5) ITERATIONS.
C                           I=3, THE ESTIMATED NUMBER OF SIGNIFICANT
C                             DIGITS IN OUTPUT VECTOR X.
C                           I=4, THE FINAL VALUE OF THE MARQUARDT
C                             SCALING PARAMETER DESCRIBED UNDER PARM(1).
C                           I=5, THE NUMBER OF ITERATIONS (I.E., CHANGES
C                             TO THE X VECTOR) PERFORMED.
C                           SEE PROGRAMMING NOTES FOR DESCRIPTION OF
C                             THE LATTER ELEMENTS OF WORK.
C                INFER  - AN INTEGER THAT IS SET, ON OUTPUT, TO
C                           INDICATE WHICH CONVERGENCE CRITERION WAS
C                           SATISFIED.
C                         INFER = 0 INDICATES THAT CONVERGENCE FAILED.
C                           IER GIVES FURTHER EXPLANATION.
C                         INFER = 1 INDICATES THAT THE FIRST CRITERION
C                           WAS SATISFIED.
C                         INFER = 2 INDICATES THAT THE SECOND CRITERION
C                           WAS SATISFIED.
C                         INFER = 4 INDICATES THAT THE THIRD CRITERION
C                           WAS SATISFIED.
C                         IF MORE THAN ONE OF THE CONVERGENCE CRITERIA
C                           WERE SATISFIED ON THE FINAL ITERATION,
C                           INFER CONTAINS THE CORRESPONDING SUM.
C                           (E.G., INFER = 3 IMPLIES FIRST AND SECOND
C                           CRITERIA SATISFIED SIMULTANEOUSLY).
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 IMPLIES A SINGULARITY WAS DETECTED
C                             IN THE JACOBIAN AND RECOVERY FAILED.
C                           IER=130 IMPLIES AT LEAST ONE OF M, N, IOPT,
C                             PARM(1), OR PARM(2) WAS SPECIFIED
C                             INCORRECTLY.
C                           IER=132 IMPLIES THAT AFTER A SUCCESSFUL
C                             RECOVERY FROM A SINGULAR JACOBIAN, THE
C                             VECTOR X HAS CYCLED BACK TO THE
C                             FIRST SINGULARITY.
C                           IER=133 IMPLIES THAT MAXFN WAS EXCEEDED.
C                         WARNING ERROR
C                           IER=38 IMPLIES THAT THE JACOBIAN IS ZERO.
C                             THE SOLUTION X IS A STATIONARY POINT.
C                           IER=39 IMPLIES THAT THE MARQUARDT
C                             PARAMETER EXCEEDED PARM(3).  THIS
C                             USUALLY MEANS THAT THE REQUESTED
C                             ACCURACY WAS NOT ACHIEVED.
C
C   REQD. IMSL ROUTINES - LEQT1P,LUDECP,LUELMP,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXSSQ  (FUNC,M,N,NSIG,EPS,DELTA,MAXFN,IOPT,PARM,
     *                   X,SSQ,F,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,NSIG,MAXFN,IOPT,IXJAC,INFER,IER
      REAL               EPS,DELTA,PARM(1),X(N),SSQ,F(M),XJAC(1),
     *                   XJTJ(1),WORK(1)
C                                  XJAC USED INTERNALLY IN PACKED FORM
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IMJC,IGRAD1,IGRADL,IGRADU,IDELX1,IDELXL,
     *                   IDELXU,ISCAL1,ISCALL,ISCALU,IXNEW1,IXNEWL,
     *                   IXBAD1,IFPL1,IFPL,IFPU,IFML1,IFML,IEVAL,
     *                   IBAD,ISW,ITER,J,IJAC,I,K,L,IS,JS,LI,LJ,ICOUNT,
     *                   IZERO,LEVEL,LEVOLD
      REAL               AL,CONS2,DNORM,DSQ,
     *                   ERL2,ERL2X,F0,F0SQ,F0SQS4,G,HALF,
     *                   HH,ONE,ONEP10,ONEP5,ONESF0,AX,
     *                   PREC,REL,RHH,SIG,SQDIF,SSQOLD,SUM,TEN,
     *                   TENTH,XDIF,XHOLD,UP,ZERO,
     *                   XDABS,RELCON,P01,TWO,HUNTW,DELTA2
      DATA               SIG/6.3/
      DATA               AX/0.1/
      DATA               P01,TENTH,HALF,ZERO,ONE,ONEP5,TWO,
     *                   TEN,HUNTW,ONEP10/0.01,0.1,0.5,0.0,
     *                   1.,1.5,2.,10.0,1.2E2,1.E10/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      LEVEL = 0
      CALL UERSET (LEVEL,LEVOLD)
      IF (M.LE.0.OR.M.GT.IXJAC.OR.N.LE.0.OR.IOPT.LT.0.OR.IOPT.GT.2)
     *   GO TO 305
      IMJC = IXJAC-M
      IF (IOPT.NE.2) GO TO 5
      IF (PARM(2).LE.ONE.OR.PARM(1).LE.ZERO) GO TO 305
C                                  MACHINE DEPENDENT CONSTANTS
    5 PREC = TEN**(-SIG-ONE)
      REL = TEN**(-SIG*HALF)
      RELCON = TEN**(-NSIG)
C                                  WORK VECTOR IS CONCATENATION OF
C                                  SCALED HESSIAN,GRADIENT,DELX,SCALE,
C                                  XNEW,XBAD,F(X+DEL),F(X-DEL)
      IGRAD1 = ((N+1)*N)/2
      IGRADL = IGRAD1+1
      IGRADU = IGRAD1+N
      IDELX1 = IGRADU
      IDELXL = IDELX1+1
      IDELXU = IDELX1+N
      ISCAL1 = IDELXU
      ISCALL = ISCAL1+1
      ISCALU = ISCAL1+N
      IXNEW1 = ISCALU
      IXNEWL = IXNEW1+1
      IXBAD1 = IXNEW1+N
      IFPL1 = IXBAD1+N
      IFPL = IFPL1+1
      IFPU = IFPL1+M
      IFML1 = IFPU
      IFML = IFML1+1
      IMJC = IXJAC - M
C                                  INITIALIZE VARIABLES
      AL = ONE
      CONS2 = TENTH
      IF (IOPT.EQ.0) GO TO 20
      IF (IOPT.EQ.1) GO TO 10
      AL = PARM(1)
      F0 = PARM(2)
      UP = PARM(3)
      CONS2 = PARM(4)
      GO TO 15
   10 AL = P01
      F0 = TWO
      UP = HUNTW
   15 ONESF0 = ONE/F0
      F0SQ = F0*F0
      F0SQS4 = F0SQ**4
   20 IEVAL = 0
      DELTA2 = DELTA*HALF
      ERL2 = ONEP10
      IBAD = -99
      ISW = 1
      ITER = -1
      INFER = 0
      IER = 0
      DO 25 J=IDELXL,IDELXU
         WORK(J) = ZERO
   25 CONTINUE
      GO TO 165
C                                  MAIN LOOP
   30 SSQOLD = SSQ
C                                  CALCULATE JACOBIAN
      IF (INFER.GT.0.OR.IJAC.GE.N.OR.IOPT.EQ.0.OR.ICOUNT.GT.0) GO TO 55
C                                  RANK ONE UPDATE TO JACOBIAN
      IJAC = IJAC+1
      DSQ = ZERO
      DO 35 J=IDELXL,IDELXU
         DSQ = DSQ+WORK(J)*WORK(J)
   35 CONTINUE
      IF (DSQ.LE.ZERO) GO TO 55
      DO 50 I=1,M
         G = F(I)-WORK(IFML1+I)
         K = I
         DO 40 J=IDELXL,IDELXU
            G = G+XJAC(K)*WORK(J)
            K = K+IXJAC
   40    CONTINUE
         G = G/DSQ
         K = I
         DO 45 J=IDELXL,IDELXU
            XJAC(K) = XJAC(K)-G*WORK(J)
            K = K+IXJAC
   45    CONTINUE
   50 CONTINUE
      GO TO 80
C                                  JACOBIAN BY INCREMENTING X
   55 IJAC = 0
      K = -IMJC
      DO 75 J=1,N
         K = K+IMJC
         XDABS = ABS(X(J))
         HH = REL*(AMAX1(XDABS,AX))
         XHOLD = X(J)
         X(J) = X(J)+HH
         CALL FUNC (X,M,N,WORK(IFPL))
         IEVAL = IEVAL+1
         X(J) = XHOLD
         IF (ISW.EQ.1) GO TO 65
C                                  CENTRAL DIFFERENCES
         X(J) = XHOLD-HH
         CALL FUNC (X,M,N,WORK(IFML))
         IEVAL = IEVAL+1
         X(J) = XHOLD
         RHH = HALF/HH
         DO 60 I=IFPL,IFPU
            K = K+1
            XJAC(K) = (WORK(I)-WORK(I+M))*RHH
   60    CONTINUE
         GO TO 75
C                                  FORWARD DIFFERENCES
   65    RHH = ONE/HH
         DO 70 I=1,M
            K = K+1
            XJAC(K) = (WORK(IFPL1+I)-F(I))*RHH
   70    CONTINUE
   75 CONTINUE
C                                  CALCULATE GRADIENT
   80 ERL2X = ERL2
      ERL2 = ZERO
      K = -IMJC
      DO 90 J=IGRADL,IGRADU
         K = K+IMJC
         SUM = ZERO
         DO 85 I=1,M
            K = K+1
            SUM = SUM+XJAC(K)*F(I)
   85    CONTINUE
         WORK(J) = SUM
         ERL2 = ERL2+SUM*SUM
   90 CONTINUE
      ERL2 = SQRT(ERL2)
C                                  CONVERGENCE TEST FOR NORM OF GRADIENT
      IF (IJAC.GT.0) GO TO 95
      IF (ERL2.LE.DELTA2) INFER = INFER+4
      IF (ERL2.LE.CONS2) ISW = 2
C                                  CALCULATE THE LOWER SUPER TRIANGE OF
C                                  JACOBIAN (TRANSPOSED) * JACOBIAN
   95 L = 0
      IS = -IXJAC
      DO 110 I=1,N
         IS = IS+IXJAC
         JS = -IXJAC
         DO 105 J=1,I
            JS = JS+IXJAC
            L = L+1
            SUM = ZERO
            DO 100 K=1,M
               LI = IS+K
               LJ = JS+K
               SUM = SUM+XJAC(LI)*XJAC(LJ)
  100       CONTINUE
            XJTJ(L) = SUM
  105    CONTINUE
  110 CONTINUE
C                                  CONVERGENCE CHECKS
      IF (INFER.GT.0) GO TO 315
      IF (IEVAL.GE.MAXFN) GO TO 290
C                                  COMPUTE SCALING VECTOR
      IF (IOPT.EQ.0) GO TO 120
      K = 0
      DO 115 J=1,N
         K = K+J
         WORK(ISCAL1+J) = XJTJ(K)
  115 CONTINUE
      GO TO 135
C                                  COMPUTE SCALING VECTOR AND NORM
  120 DNORM = ZERO
      K = 0
      DO 125 J=1,N
         K = K+J
         WORK(ISCAL1+J) = SQRT(XJTJ(K))
         DNORM = DNORM+XJTJ(K)*XJTJ(K)
  125 CONTINUE
      DNORM = ONE/SQRT(DNORM)
C                                  NORMALIZE SCALING VECTOR
      DO 130 J=ISCALL,ISCALU
         WORK(J) = WORK(J)*DNORM*ERL2
  130 CONTINUE
C                                  ADD L-M FACTOR TO DIAGONAL
  135 ICOUNT = 0
  140 K = 0
      DO 150 I=1,N
         DO 145 J=1,I
            K = K+1
            WORK(K) = XJTJ(K)
  145    CONTINUE
         WORK(K) = WORK(K)+WORK(ISCAL1+I)*AL
         WORK(IDELX1+I) = WORK(IGRAD1+I)
  150 CONTINUE
C                                  CHOLESKY DECOMPOSITION
  155 CALL LEQT1P (WORK,1,N,WORK(IDELXL),N,0,G,XHOLD,IER)
      IF (IER.EQ.0) GO TO 160
      IER = 0
      IF (IJAC.GT.0) GO TO 55
      IF (IBAD.LE.0) GO TO 240
      IF (IBAD.GE.2) GO TO 310
      GO TO 190
  160 IF (IBAD.NE.-99) IBAD = 0
C                                  CALCULATE SUM OF SQUARES
  165 DO 170 J=1,N
         WORK(IXNEW1+J) = X(J)-WORK(IDELX1+J)
  170 CONTINUE
      CALL FUNC (WORK(IXNEWL),M,N,WORK(IFPL))
      IEVAL = IEVAL+1
      SSQ = ZERO
      DO 175 I=IFPL,IFPU
         SSQ = SSQ+WORK(I)*WORK(I)
  175 CONTINUE
      IF (ITER.GE.0) GO TO 185
C                                  SSQ FOR INITIAL ESTIMATES OF X
      ITER = 0
      SSQOLD = SSQ
      DO 180 I=1,M
         F(I) = WORK(IFPL1+I)
  180 CONTINUE
      GO TO 55
  185 IF (IOPT.EQ.0) GO TO 215
C                                  CHECK DESCENT PROPERTY
      IF (SSQ.LE.SSQOLD) GO TO 205
C                                  INCREASE PARAMETER AND TRY AGAIN
  190 ICOUNT = ICOUNT+1
      AL = AL*F0SQ
      IF (IJAC.EQ.0) GO TO 195
      IF (ICOUNT.GE.4.OR.AL.GT.UP) GO TO 200
  195 IF (AL.LE.UP) GO TO 140
      IF (IBAD.EQ.1) GO TO 310
      IER = 39
      GO TO 315
  200 AL = AL/F0SQS4
      GO TO 55
C                                  ADJUST MARQUARDT PARAMETER
  205 IF (ICOUNT.EQ.0) AL = AL/F0
      IF (ERL2X.LE.ZERO) GO TO 210
      G = ERL2/ERL2X
      IF (ERL2.LT.ERL2X) AL = AL*AMAX1(ONESF0,G)
      IF (ERL2.GT.ERL2X) AL = AL*AMIN1(F0,G)
  210 AL = AMAX1(AL,PREC)
C                                  ONE ITERATION CYCLE COMPLETED
  215 ITER = ITER+1
      DO 220 J=1,N
         X(J) = WORK(IXNEW1+J)
  220 CONTINUE
      DO 225 I=1,M
         WORK(IFML1+I) = F(I)
         F(I) = WORK(IFPL1+I)
  225 CONTINUE
C                                  RELATIVE CONVERGENCE TEST FOR X
      IF (AL.GT.5.0) GO TO 30
      DO 230 J=1,N
         XDIF = ABS(WORK(IDELX1+J))/AMAX1(ABS(X(J)),AX)
         IF (XDIF.GT.RELCON) GO TO 235
  230 CONTINUE
      INFER = 1
C                                  RELATIVE CONVERGENCE TEST FOR SSQ
  235 SQDIF = ABS(SSQ-SSQOLD)/AMAX1(SSQOLD,AX)
      IF (SQDIF.LE.EPS) INFER = INFER+2
      GO TO 30
C                                  SINGULAR DECOMPOSITION
  240 IF (IBAD) 255,245,265
C                                  CHECK TO SEE IF CURRENT
C                                  ITERATE HAS CYCLED BACK TO
C                                  THE LAST SINGULAR POINT
  245 DO 250 J=1,N
         XHOLD = WORK(IXBAD1+J)
         IF (ABS(X(J)-XHOLD).GT.RELCON*AMAX1(AX,ABS(XHOLD))) GO TO 255
  250 CONTINUE
      GO TO 295
C                                  UPDATE THE BAD X VALUES
  255 DO 260 J=1,N
         WORK(IXBAD1+J) = X(J)
  260 CONTINUE
      IBAD = 1
C                                  INCREASE DIAGONAL OF HESSIAN
  265 IF (IOPT.NE.0) GO TO 280
      K = 0
      DO 275 I=1,N
         DO 270 J=1,I
            K = K+1
            WORK(K) = XJTJ(K)
  270    CONTINUE
         WORK(K) = ONEP5*(XJTJ(K)+AL*ERL2*WORK(ISCAL1+I))+REL
  275 CONTINUE
      IBAD = 2
      GO TO 155
C                                  REPLACE ZEROES ON HESSIAN DIAGONAL
  280 IZERO = 0
      DO 285 J=ISCALL,ISCALU
         IF (WORK(J).GT.ZERO) GO TO 285
         IZERO = IZERO+1
         WORK(J) = ONE
  285 CONTINUE
      IF (IZERO.LT.N) GO TO 140
      IER = 38
      GO TO 315
C                                  TERMINAL ERROR
  290 IER = IER+1
  295 IER = IER+1
      IER = IER+1
  305 IER = IER+1
  310 IER = IER+129
      IF (IER.EQ.130) GO TO 335
C                                  OUTPUT ERL2,IEVAL,NSIG,AL, AND ITER
  315 G = SIG
      DO 320 J=1,N
         XHOLD = ABS(WORK(IDELX1+J))
         IF (XHOLD.LE.ZERO) GO TO 320
         G = AMIN1(G,-ALOG10(XHOLD)+ALOG10(AMAX1(AX,ABS(X(J)))))
  320 CONTINUE
      IF(N.GT.2) GO TO 330
      DO 325 J = 1,N
  325 WORK(J+5) = WORK(J+IGRAD1)
  330 WORK(1) = ERL2+ERL2
      WORK(2) = IEVAL
      SSQ = SSQOLD
      WORK(3) = G
      WORK(4) = AL
      WORK(5) = ITER
  335 CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HZXSSQ )
 9005 RETURN
      END

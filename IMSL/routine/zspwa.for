C   IMSL ROUTINE NAME   - ZSPWA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2,ZSPWB,ZSPWC,ZSPWD,ZSPWE,
C                           ZSPWF,ZSPWG
C                       - DOUBLE/VBLA=DNRM2,ZSPWB,ZSPWC,ZSPWD,ZSPWE,
C                           ZSPWF,ZSPWG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSPWA (FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,DIAG,MODE,
     *                   FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,R,LR,QTF,
     *                   WA1,WA2,WA3,WA4,PAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR
      REAL               X(N),FVEC(N),XTOL,EPSFCN,DIAG(N),FACTOR,
     *                   FJAC(LDFJAC,N),R(LR),QTF(N),WA1(N),WA2(N),
     *                   WA3(N),WA4(N),PAR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IFLAG,ITER,IWA(1),I,JM1,J,L,MSUM,NCFAIL,NCSUC,
     *                   NSLOW1,NSLOW2
      REAL               ACTRED,DELTA,EPSMCH,FNORM1,FNORM,ONE,P0001,
     *                   P001,P1,P5,PNORM,PRERED,RATIO,SPMPAR,SUM,TEMP,
     *                   XNORM,ZERO
      REAL               SNRM2
      LOGICAL            JEVAL,SING
      EXTERNAL           FCN
      DATA               SPMPAR /Z3C100000/
      DATA               ONE,P1,P5,P001,P0001,ZERO /1.0E0,1.0E-1,5.0E-1,
     *                   1.0E-3,1.0E-4,0.0E0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
      INFO = 0
      IFLAG = 0
      NFEV = 0
C                                  CHECK THE INPUT PARAMETERS FOR
C                                  ERRORS.
      IF (N.LE.0 .OR. XTOL.LT.ZERO .OR. MAXFEV.LE.0 .OR. ML.LT.0 .OR.
     *MU.LT.0 .OR. FACTOR.LE.ZERO .OR. LDFJAC.LT.N .OR.
     *LR.LT.(N*(N+1))/2) GO TO 150
      IF (MODE.NE.2) GO TO 10
      DO 5 J=1,N
         IF (DIAG(J).LE.ZERO) GO TO 150
    5 CONTINUE
   10 CONTINUE
C                                  EVALUATE THE FUNCTION AT THE STARTING
C                                  POINT AND CALCULATE ITS NORM.
      IFLAG = 1
      CALL FCN(X,FVEC,N,PAR)
      NFEV = 1
      IF (IFLAG.LT.0) GO TO 150
      FNORM = SNRM2(N,FVEC,1)
C                                  DETERMINE THE NUMBER OF CALLS TO FCN
C                                  NEEDED TO COMPUTE THE JACOBIAN
C                                  MATRIX.
C
      MSUM = MIN0(ML+MU+1,N)
C
C                                  INITIALIZE ITERATION COUNTER AND
C                                  MONITORS.
      ITER = 1
      NCSUC = 0
      NCFAIL = 0
      NSLOW1 = 0
      NSLOW2 = 0
C                                  BEGINNING OF THE OUTER LOOP.
   15 CONTINUE
      JEVAL = .TRUE.
C                                  CALCULATE THE JACOBIAN MATRIX.
      IFLAG = 2
      CALL ZSPWB(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,WA2,
     *PAR)
      NFEV = NFEV+MSUM
      IF (IFLAG.LT.0) GO TO 150
C                                  COMPUTE THE QR FACTORIZATION OF THE
C                                  JACOBIAN.
      CALL ZSPWG(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)
C                                  ON THE FIRST ITERATION AND IF MODE IS
C                                  1, SCALE ACCORDING TO THE NORMS OF
C                                  THE COLUMNS OF THE INITIAL JACOBIAN.
      IF (ITER.NE.1) GO TO 35
      IF (MODE.EQ.2) GO TO 25
      DO 20 J=1,N
         DIAG(J) = WA2(J)
         IF (WA2(J).EQ.ZERO) DIAG(J) = ONE
   20 CONTINUE
   25 CONTINUE
C                                  ON THE FIRST ITERATION, CALCULATE THE
C                                  NORM OF THE SCALED X AND INITIALIZE
C                                  THE STEP BOUND DELTA.
      DO 30 J=1,N
         WA3(J) = DIAG(J)*X(J)
   30 CONTINUE
      XNORM = SNRM2(N,WA3,1)
      DELTA = FACTOR*XNORM
      IF (DELTA.EQ.ZERO) DELTA = FACTOR
   35 CONTINUE
C                                  FORM (Q TRANSPOSE)*FVEC AND STORE IN
C                                  QTF.
      DO 40 I=1,N
         QTF(I) = FVEC(I)
   40 CONTINUE
      DO 60 J=1,N
         IF (FJAC(J,J).EQ.ZERO) GO TO 55
         SUM = ZERO
         DO 45 I=J,N
            SUM = SUM+FJAC(I,J)*QTF(I)
   45    CONTINUE
         TEMP = -SUM/FJAC(J,J)
         DO 50 I=J,N
            QTF(I) = QTF(I)+FJAC(I,J)*TEMP
   50    CONTINUE
   55    CONTINUE
   60 CONTINUE
C                                  COPY THE TRIANGULAR FACTOR OF THE QR
C                                  FACTORIZATION INTO R.
      SING = .FALSE.
      DO 75 J=1,N
         L = J
         JM1 = J-1
         IF (JM1.LT.1) GO TO 70
         DO 65 I=1,JM1
            R(L) = FJAC(I,J)
            L = L+N-I
   65    CONTINUE
   70    CONTINUE
         R(L) = WA1(J)
         IF (WA1(J).EQ.ZERO) SING = .TRUE.
   75 CONTINUE
C                                  ACCUMULATE THE ORTHOGONAL FACTOR IN
C                                  FJAC.
      CALL ZSPWF(N,N,FJAC,LDFJAC,WA1)
C                                  RESCALE IF NECESSARY.
      IF (MODE.EQ.2) GO TO 85
      DO 80 J=1,N
         DIAG(J) = AMAX1(DIAG(J),WA2(J))
   80 CONTINUE
   85 CONTINUE
C                                  BEGINNING OF THE INNER LOOP.
   90 CONTINUE
C                                  IF REQUESTED, CALL FCN TO ENABLE
C                                  PRINTING OF ITERATES.
      IF (NPRINT.LE.0) GO TO 95
      IFLAG = 0
      IF (IFLAG.LT.0) GO TO 150
   95 CONTINUE
C                                  DETERMINE THE DIRECTION P.
      CALL ZSPWC(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C                                  STORE THE DIRECTION P AND X + P.
C                                  CALCULATE THE NORM OF P.
      DO 100 J=1,N
         WA1(J) = -WA1(J)
         WA2(J) = X(J)+WA1(J)
         WA3(J) = DIAG(J)*WA1(J)
  100 CONTINUE
      PNORM = SNRM2(N,WA3,1)
C                                  ON THE FIRST ITERATION, ADJUST THE
C                                  INITIAL STEP BOUND.
      IF (ITER.EQ.1) DELTA = AMIN1(DELTA,PNORM)
C                                  EVALUATE THE FUNCTION AT X + P AND
C                                  CALCULATE ITS NORM.
      IFLAG = 1
      CALL FCN(WA2,WA4,N,PAR)
      NFEV = NFEV+1
      IF (IFLAG.LT.0) GO TO 150
      FNORM1 = SNRM2(N,WA4,1)
C                                  COMPUTE THE SCALED ACTUAL REDUCTION.
      ACTRED = -ONE
      IF (FNORM1.LT.FNORM) ACTRED = ONE-(FNORM1/FNORM)**2
C                                  COMPUTE THE SCALED PREDICTED
C                                  REDUCTION.
      L = 1
      DO 110 I=1,N
         SUM = ZERO
         DO 105 J=I,N
            SUM = SUM+R(L)*WA1(J)
            L = L+1
  105    CONTINUE
         WA3(I) = QTF(I)+SUM
  110 CONTINUE
      TEMP = SNRM2(N,WA3,1)
      PRERED = ONE
      IF (TEMP.LT.FNORM) PRERED = ONE-(TEMP/FNORM)**2
C                                  COMPUTE THE RATIO OF THE ACTUAL TO
C                                  THE PREDICTED REDUCTION.
      RATIO = ZERO
      IF (PRERED.GT.ZERO) RATIO = ACTRED/PRERED
C                                  UPDATE THE STEP BOUND.
      IF (RATIO.GE.P1) GO TO 115
      NCSUC = 0
      NCFAIL = NCFAIL+1
      DELTA = P5*DELTA
      GO TO 120
  115 CONTINUE
      NCFAIL = 0
      NCSUC = NCSUC+1
      IF (RATIO.GE.P5 .OR. NCSUC.GT.1) DELTA = AMAX1(DELTA,PNORM/P5)
      IF (ABS(RATIO-ONE).LE.P1) DELTA = PNORM/P5
  120 CONTINUE
C                                  TEST FOR SUCCESSFUL ITERATION.
      IF (RATIO.LT.P0001) GO TO 130
C                                  SUCCESSFUL ITERATION. UPDATE X, FVEC,
C                                  AND THEIR NORMS.
      DO 125 J=1,N
         X(J) = WA2(J)
         WA2(J) = DIAG(J)*X(J)
         FVEC(J) = WA4(J)
  125 CONTINUE
      XNORM = SNRM2(N,WA2,1)
      FNORM = FNORM1
      ITER = ITER+1
  130 CONTINUE
C                                  DETERMINE THE PROGRESS OF THE
C                                  ITERATION.
      NSLOW1 = NSLOW1+1
      IF (ACTRED.GE.P001) NSLOW1 = 0
      IF (JEVAL) NSLOW2 = NSLOW2+1
      IF (ACTRED.GE.P1) NSLOW2 = 0
C                                  TEST FOR CONVERGENCE.
      IF (DELTA.LE.XTOL*XNORM .OR. FNORM.EQ.ZERO) INFO = 1
      IF (INFO.NE.0) GO TO 150
C                                  TESTS FOR TERMINATION AND STRINGENT
C                                  TOLERANCES.
      IF (NFEV.GE.MAXFEV) INFO = 2
      IF (P1*AMAX1(P1*DELTA,PNORM).LE.EPSMCH*XNORM) INFO = 3
      IF (NSLOW2.EQ.5) INFO = 4
      IF (NSLOW1.EQ.10) INFO = 5
      IF (INFO.NE.0) GO TO 150
C                                  CRITERION FOR RECALCULATING JACOBIAN
C                                  APPROXIMATION BY FORWARD DIFFERENCES.
      IF (NCFAIL.EQ.2) GO TO 145
C                                  CALCULATE THE RANK ONE MODIFICATION
C                                  TO THE JACOBIAN AND UPDATE QTF IF
C                                  NECESSARY.
      DO 140 J=1,N
         SUM = ZERO
         DO 135 I=1,N
            SUM = SUM+FJAC(I,J)*WA4(I)
  135    CONTINUE
         WA2(J) = (SUM-WA3(J))/PNORM
         WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
         IF (RATIO.GE.P0001) QTF(J) = SUM
  140 CONTINUE
C                                  COMPUTE THE QR FACTORIZATION OF THE
C                                  UPDATED JACOBIAN.
      CALL ZSPWE(N,N,R,LR,WA1,WA2,WA3,SING)
      CALL ZSPWD(N,N,FJAC,LDFJAC,WA2,WA3)
      CALL ZSPWD(1,N,QTF,1,WA2,WA3)
C                                  END OF THE INNER LOOP.
      JEVAL = .FALSE.
      GO TO 90
  145 CONTINUE
C                                  END OF THE OUTER LOOP.
      GO TO 15
  150 CONTINUE
C                                  TERMINATION, EITHER NORMAL OR USER
C                                  IMPOSED.
      IF (IFLAG.LT.0) INFO = IFLAG
      IFLAG = 0
      RETURN
      END

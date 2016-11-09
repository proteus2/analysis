C   IMSL ROUTINE NAME   - ZX4LR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZX4LP
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO,VBLA=ISAMAX,VBLA=SASUM,
C                           VBLA=SAXPY,VBLA=SCOPY,VBLA=SDOT,VBLA=SROTM,
C                           VBLA=SROTMG,VBLA=SSCAL,VBLA=SSWAP,ZX4LQ
C                       - DOUBLE/UERTST,UGETIO,VBLA=IDAMAX,VBLA=DASUM,
C                           VBLA=DAXPY,VBLA=DCOPY,VBLA=DDOT,VBLA=DROTM,
C                           VBLA=DROTMG,VBLA=DSCAL,VBLA=DSWAP,ZX4LQ
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZX4LR  (W,MDW,M,N,PRGOPT,WS,Y,MODE,D,BND,X,B,MDB,IPTR,
     *                   IND,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MDW,M,N,MODE,MDB,IPTR(1),IND(1),IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               W(MDW,1),B(MDB,MDB),D(1),X(1),Y(1),PRGOPT(1)
      REAL               WS(1),BND(1),DUMMY(2),DON
      REAL               AD,AN,ANORM,BIG,BIISQR,BLJSQR,BMIN,BNORM,BOUND,
     *                   BSCAL,CNORM,COST,COSTMN,COSTSC,EPS,FAC,ONE,
     *                   PNORM,QUOT,REPS,RERSQ,RMAX,RWNMSQ,SIZE,T,TEMPL,
     *                   THETA,TOL,XMAX,XMIN,XNORM,YMAX,YNORM,ZERO
      LOGICAL            ARTIF,FINITE,FOUND,FIRST,ONCE
      LOGICAL            COLSC,RHSSC,CSTSC
      LOGICAL            DROP,CONVRG,INDEP
      DATA               REPS/Z3C100000/
      DATA               ZERO,ONE /0.,1.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  THE FOLLOWING CODE WAS GENERATED
C                                  BY A PREPROCESSOR FROM A HIGHER
C                                  LEVEL LANGUAGE
      EPS = REPS*0.01
      GO TO 35
    5 GO TO 165
   10 GO TO 455
   15 IF (.NOT.(ARTIF)) GO TO 25
      GO TO 660
   20 CONTINUE
   25 GO TO 760
   30 GO TO 1285
   35 CONTINUE
      MODE = 1
      IF (M.LE.0 .OR. N.LE.0) GO TO 9005
      MP1 = M+1
      MP2 = M+2
      NP1 = N+1
      NP2 = N+2
      NREDC = 0
      NITER = 0
      DO 50 I=1,N
         IND(I) = IPTR(I)
         IF (.NOT.(IND(I).NE.0)) GO TO 45
         IND(I) = 1
         IF (.NOT.(WS(I).LT.ZERO)) GO TO 40
         IER = 131
         GO TO 9000
   40    BND(I) = WS(I)
   45    IPTR(I) = I
   50 CONTINUE
      IND(NP1) = 0
      IPTR(NP1) = 0
      GO TO 125
   55 CONTINUE
      N3 = 3*N
      IF (.NOT.(MDW.LT.MP2)) GO TO 60
      IER = 129
      GO TO 9000
   60 IF (.NOT.(M.GT.N)) GO TO 65
      IER = 130
      GO TO 9000
   65 IF (COLSC) GO TO 85
      DO 80 J=1,N
         IMAX = ISAMAX(M,W(1,J),1)
         T = ABS(W(IMAX,J))
         IF (.NOT.(T.NE.ZERO)) GO TO 70
         WS(J) = ONE/T
         GO TO 75
   70    WS(J) = ONE
   75    CONTINUE
   80 CONTINUE
   85 W(MP1,NP1) = ZERO
      W(MP2,1) = ZERO
      CALL SCOPY (NP1,W (MP2,1),0,W (MP2,1),MDW)
      B(1,MP1) = ZERO
      CALL SCOPY (MP2,B (1,MP1),0,B (1,MP1),1)
      D(MP1) = ONE
      IF (CSTSC) GO TO 105
      CNORM = ZERO
      DO 90 J=1,N
         CNORM = AMAX1(CNORM,ABS(WS(J)*W(MP1,J)))
   90 CONTINUE
      IF (.NOT.(CNORM.NE.ZERO)) GO TO 95
      COSTSC = ONE/CNORM
      GO TO 100
   95 COSTSC = ONE
  100 CONTINUE
  105 CALL SSCAL (N,COSTSC,W (MP1,1),MDW)
      CALL SCOPY (MP1,W (1,NP1),1,W (1,NP2),1)
      IF (RHSSC) GO TO 120
      IMAX = ISAMAX(M,W(1,NP2),1)
      T = ABS(W(IMAX,NP2))
      IF (.NOT.(T.NE.ZERO)) GO TO 110
      BSCAL = ONE/T
      GO TO 115
  110 BSCAL = ONE
  115 CONTINUE
  120 CALL SSCAL (M,BSCAL,W (1,NP2),1)
      FAC = M+M
      DON = M*M
      TOL = EPS*DON
      RERSQ = TOL
      TOL = TOL+TOL
      BNORM = ZERO
      GO TO 5
  125 CONTINUE
      KPRINT = 0
      IDIGIT = -4
      COLSC = .FALSE.
      CSTSC = .FALSE.
      COSTSC = ONE
      RHSSC = .FALSE.
      BSCAL = ONE
      NLINK = 10000
      NOPT = 1000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (.NOT.(LINK.LE.0 .OR. LINK.GT.NLINK)) GO TO 130
      IER = 138
      GO TO 9000
  130 IF (.NOT.(LINK.GT.1)) GO TO 160
      NTIMES = NTIMES+1
      IF (.NOT.(NTIMES.GT.NOPT)) GO TO 135
      IER = 138
      GO TO 9000
  135 KEY = PRGOPT(LAST+1)
      IF (KEY.EQ.13) KPRINT = PRGOPT(LAST+2)
      IF (KEY.EQ.14) IDIGIT = PRGOPT(LAST+2)
      IF (.NOT.(KEY.EQ.15)) GO TO 140
      COLSC = .TRUE.
      CALL SCOPY (N,PRGOPT (LAST+2),1,WS,1)
  140 IF (.NOT.(KEY.EQ.16)) GO TO 145
      RHSSC = .TRUE.
      T = PRGOPT(LAST+2)
      IF (T.NE.ZERO) BSCAL = T
  145 IF (.NOT.(KEY.EQ.17)) GO TO 150
      CSTSC = .TRUE.
      T = PRGOPT(LAST+2)
      IF (T.NE.ZERO) COSTSC = T
  150 NEXT = PRGOPT(LINK)
      IF (.NOT.(NEXT.LE.0 .OR. NEXT.GT.NLINK)) GO TO 155
      IER = 138
      GO TO 9000
  155 LAST = LINK
      LINK = NEXT
      GO TO 130
  160 GO TO 55
  165 CONTINUE
      J = 1
      JM1 = 0
      IBAD = N
      ITRY = N
      DROP = .FALSE.
      GO TO 175
  170 IF (.NOT.(J.LE.M .AND. IBAD.GT.M)) GO TO 325
  175 K = IPTR(J)
      JBIG = J
      BIG = ZERO
      I = J
  180 IF (.NOT.(I.LE.IBAD)) GO TO 195
      K = IPTR(I)
      IF (.NOT.(K.GT.0)) GO TO 190
      T = ABS(W(J,K)*WS(K))
      IF (.NOT.(T.GT.BIG)) GO TO 185
      BIG = T
      JBIG = I
  185 CONTINUE
  190 I = I+1
      GO TO 180
  195 ITEMP = IPTR(J)
      IPTR(J) = IPTR(JBIG)
      IPTR(JBIG) = ITEMP
      K = IPTR(J)
      D(MP2) = WS(K)**2
      CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
      CALL ZX4LQ (B,M,JM1,MDB,D,IFLAG)
      B(1,J) = ZERO
      CALL SCOPY (JM1,B (1,J),0,B (1,J),1)
      CALL SCOPY (MP2-JM1,B (J,MP2),1,B (J,J),1)
      D(J) = D(MP2)
      DO 210 L=J,M
         RWNMSQ = ZERO
         IF (.NOT.(JM1.GT.0)) GO TO 205
         DO 200 I=1,JM1
            RWNMSQ = RWNMSQ+D(I)*B(L,I)**2
  200    CONTINUE
  205    BLJSQR = D(J)*B(L,J)**2
         IF (BLJSQR+RWNMSQ.NE.RWNMSQ) GO TO 215
  210 CONTINUE
      L = 0
  215 CONTINUE
      K = IABS(IPTR(J))
      IF (.NOT.(L.GT.J)) GO TO 290
      LMAX = L
      RMAX = RWNMSQ/BLJSQR
  220 IF (.NOT.(L.LE.M)) GO TO 240
      RWNMSQ = ZERO
      IF (.NOT.(JM1.GT.0)) GO TO 230
      DO 225 I=1,JM1
         RWNMSQ = RWNMSQ+D(I)*B(L,I)**2
  225 CONTINUE
  230 BLJSQR = D(J)*B(L,J)**2
      IF (.NOT.(BLJSQR*RMAX.GT.RWNMSQ)) GO TO 235
      RMAX = RWNMSQ/BLJSQR
      LMAX = L
  235 L = L+1
      GO TO 220
  240 L = LMAX
      IF (.NOT.(IPTR(L).GT.0)) GO TO 245
      ITEMP = IPTR(J)
      IPTR(J) = IPTR(L)
      IPTR(L) = ITEMP
      IPTR(L) = -IABS(IPTR(L))
      GO TO 275
  245 ITRY = MIN0(ITRY,IBAD)
      ITRYST = MAX0(ITRY,J)
  250 IF (.NOT.(IPTR(ITRY).LT.0 .AND. ITRY.GT.J)) GO TO 255
      ITRY = ITRY-1
      GO TO 250
  255 IF (.NOT.(ITRY.LE.J)) GO TO 270
      ITRY = IBAD
  260 IF (.NOT.(IPTR(ITRY).LT.0 .AND. ITRY.GT.ITRYST)) GO TO 265
      ITRY = ITRY-1
      GO TO 260
  265 IF (.NOT.(ITRY.EQ.ITRYST)) GO TO 270
      GO TO 400
  270 CONTINUE
      ITEMP = IPTR(J)
      IPTR(J) = IPTR(ITRY)
      IPTR(ITRY) = ITEMP
      IPTR(ITRY) = -IABS(IPTR(ITRY))
  275 CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
      D(MP2) = -WS(K)**2
      DROP = .TRUE.
      CALL ZX4LQ (B,M,JM1,MDB,D,IFLAG)
      IF (.NOT.(IFLAG.GT.0)) GO TO 285
      I99918 = 1
      GO TO 410
  280 CONTINUE
  285 GO TO 320
  290 IF (.NOT.(L.EQ.J)) GO TO 305
      JM1 = J
      J = J+1
      I = IBAD
  295 IF (.NOT.(I.GT.J)) GO TO 300
      IPTR(I) = IABS(IPTR(I))
      I = I-1
      GO TO 295
  300 BNORM = AMAX1(BNORM,ABS(WS(K)*SASUM(M,W(1,K),1)))
      GO TO 320
  305 CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
      D(MP2) = -WS(K)**2
      DROP = .TRUE.
      ITEMP = IPTR(J)
      IPTR(J) = IPTR(IBAD)
      IPTR(IBAD) = ITEMP
      CALL ZX4LQ (B,M,JM1,MDB,D,IFLAG)
      IBAD = IBAD-1
      IF (.NOT.(IFLAG.GT.0)) GO TO 315
      I99918 = 2
      GO TO 410
  310 CONTINUE
  315 CONTINUE
  320 GO TO 170
  325 DO 330 I=1,NP1
         IPTR(I) = IABS(IPTR(I))
  330 CONTINUE
      IF (.NOT.(IBAD.LT.M)) GO TO 335
      GO TO 400
  335 IF (.NOT.(J.LE.M)) GO TO 365
      JM1 = M
      J = MP1
      DROP = .TRUE.
      I99918 = 3
      GO TO 410
  340 INDEP = .TRUE.
      DO 355 I=1,M
         RWNMSQ = ZERO
         IM1 = I-1
         IF (.NOT.(IM1.GT.0)) GO TO 350
         DO 345 K=1,IM1
            RWNMSQ = RWNMSQ+D(K)*B(I,K)**2
  345    CONTINUE
  350    BIISQR = D(I)*B(I,I)**2
         INDEP = INDEP .AND. (BIISQR+RWNMSQ.NE.RWNMSQ)
  355 CONTINUE
      IF (INDEP) GO TO 360
      GO TO 400
  360 CONTINUE
  365 IF (.NOT.(DROP)) GO TO 375
      I99918 = 4
      GO TO 410
  370 GO TO 395
  375 RERSQ = RERSQ/BNORM
      AN = D(1)*B(1,1)**2
      AD = AN
      I = 2
  380 IF (.NOT.(I.LE.M)) GO TO 385
      T = D(I)*B(I,I)**2
      AN = AMAX1(AN,T)
      AD = AMIN1(AD,T)
      I = I+1
      GO TO 380
  385 AN = SQRT(AN)
      AD = SQRT(AD)
      I99895 = 1
      GO TO 405
  390 CONTINUE
  395 GO TO 10
  400 CONTINUE
      IER = 132
      GO TO 9000
  405 CONTINUE
      GO TO 1340
  410 CONTINUE
      IF (.NOT.(DROP)) GO TO 450
      BNORM = ZERO
      IF (.NOT.(JM1.GT.0)) GO TO 430
      DO 425 I=1,JM1
         IM1 = I-1
         K = IABS(IPTR(I))
         IF (.NOT.(K.EQ.NP1)) GO TO 415
         T = ONE
         GO TO 420
  415    T = WS(K)
  420    BNORM = AMAX1(BNORM,ABS(T*SASUM(M,W(1,K),1)))
         CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
         D(MP2) = T**2
         CALL ZX4LQ (B,M,IM1,MDB,D,IFLAG)
         B(1,I) = ZERO
         CALL SCOPY (IM1,B (1,I),0,B (1,I),1)
         CALL SCOPY (MP2-IM1,B (I,MP2),1,B (I,I),1)
         D(I) = D(MP2)
  425 CONTINUE
  430 NREDC = NREDC+1
      DROP = .FALSE.
      AN = D(1)*B(1,1)**2
      AD = AN
      I = 2
  435 IF (.NOT.(I.LE.JM1)) GO TO 440
      T = D(I)*B(I,I)**2
      AN = AMAX1(AN,T)
      AD = AMIN1(AD,T)
      I = I+1
      GO TO 435
  440 DON = M*M
      TOL = EPS*DON
      TOL = TOL+TOL
      RERSQ = TOL/BNORM
      AN = SQRT(AN)
      AD = SQRT(AD)
      I99895 = 2
      GO TO 405
  445 CONTINUE
  450 GO TO 1335
  455 CONTINUE
      I99883 = 1
      GO TO 615
  460 IMIN = 1
      XMIN = X(1)
      DO 470 I=1,M
         IF (.NOT.(X(I).LT.XMIN)) GO TO 465
         IMIN = I
         XMIN = X(I)
  465    CONTINUE
  470 CONTINUE
      IMAX = ISAMAX(M,X,1)
      XNORM = ABS(X(IMAX))
      ARTIF = AD*ABS(XMIN).GT.XNORM*TOL*AN .AND. XMIN.LT.ZERO
      DO 480 I=1,M
         L = IPTR(I)
         IF (.NOT.(IND(L).NE.0)) GO TO 475
         BOUND = BND(L)*BSCAL/WS(L)
         XMAX = BOUND-X(I)
         ARTIF = ARTIF .OR. (AD*ABS(XMAX).GT.XNORM*TOL*AN .AND.
     *   XMAX.LT.ZERO)
  475    CONTINUE
  480 CONTINUE
      IF (.NOT.(ARTIF)) GO TO 605
      BMIN = ONE
      FIRST = .TRUE.
      DO 500 I=1,M
         L = IPTR(I)
         IF (.NOT.(IND(L).NE.0)) GO TO 495
         BOUND = BND(L)*BSCAL/WS(L)
         IF (.NOT.(FIRST)) GO TO 485
         FIRST = .FALSE.
         BMIN = BOUND
         GO TO 490
  485    BMIN = AMIN1(BMIN,BOUND)
  490    CONTINUE
  495    CONTINUE
  500 CONTINUE
      L = 1
  505 IF (.NOT.(L.LE.M)) GO TO 585
      CALL SCOPY (MP1,W (1,NP2),1,W (1,NP1),1)
      DO 510 I=1,M
         K = IPTR(I)
         IF (I.NE.L) CALL SAXPY (MP1,-BMIN*WS (K),W (1,K),1,W (1,NP1),1)
  510 CONTINUE
      W(MP1,NP1) = ZERO
      W(MP2,NP1) = ONE
      CALL SCOPY (MP2,W (1,NP1),1,B (1,MP2),1)
      D(MP2) = ONE
      CALL ZX4LQ (B,M,M,MDB,D,IFLAG)
      I99868 = 1
      GO TO 1280
  515 I = NP1
  520 IF (IPTR(I).EQ.0) GO TO 525
      I = I-1
      GO TO 520
  525 IPTR(I) = NP1
      K = IPTR(L)
      CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
      D(MP2) = -WS(K)**2
      DROP = .TRUE.
      ITEMP = IPTR(L)
      IPTR(L) = IPTR(I)
      IPTR(I) = ITEMP
      CALL ZX4LQ (B,M,M,MDB,D,IFLAG)
      I99918 = 5
      GO TO 410
  530 INDEP = .TRUE.
      DO 545 I=1,M
         RWNMSQ = ZERO
         IM1 = I-1
         IF (.NOT.(IM1.GT.0)) GO TO 540
         DO 535 K=1,IM1
            RWNMSQ = RWNMSQ+D(K)*B(I,K)**2
  535    CONTINUE
  540    BIISQR = D(I)*B(I,I)**2
         INDEP = INDEP .AND. (BIISQR+RWNMSQ.NE.RWNMSQ)
  545 CONTINUE
      IF (INDEP) GO TO 590
      I = NP1
  550 IF (IPTR(I).EQ.ITEMP) GO TO 555
      I = I-1
      GO TO 550
  555 IPTR(L) = IPTR(I)
      IPTR(I) = 0
      K = IPTR(L)
      CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
      D(MP2) = WS(K)**2
      CALL ZX4LQ (B,M,M,MDB,D,IFLAG)
      I99868 = 2
      GO TO 1280
  560 CALL SCOPY (MP2,W (1,NP1),1,B (1,MP2),1)
      D(MP2) = -ONE
      DROP = .TRUE.
      CALL ZX4LQ (B,M,M,MDB,D,IFLAG)
      IF (.NOT.(0.LT.IFLAG)) GO TO 570
      I99918 = 6
      GO TO 410
  565 GO TO 580
  570 I99868 = 3
      GO TO 1280
  575 CONTINUE
  580 L = L+1
      GO TO 505
  585 IER = 133
      GO TO 9000
  590 CONTINUE
      I99883 = 2
      GO TO 615
  595 GO TO 610
  600 CONTINUE
  605 GO TO 15
  610 CONTINUE
      GO TO 600
  615 CONTINUE
      CALL SCOPY (M,W (1,NP2),1,X,1)
      DO 620 I=1,M
         X(I) = X(I)/B(I,I)
         IF (I.EQ.M) GO TO 620
         CALL SAXPY (M-I,-X (I),B (I+1,I),1,X (I+1),1)
  620 CONTINUE
      DO 625 I=1,M
         X(I) = X(I)/D(I)
  625 CONTINUE
      DO 630 II=1,M
         I = MP1-II
         X(I) = X(I)/B(I,I)
         CALL SAXPY (I-1,-X (I),B (I,1),MDB,X,1)
  630 CONTINUE
      CALL SCOPY (M,X,1,B (1,MP2),1)
      DO 645 I=1,M
         K = IPTR(I)
         IF (.NOT.(K.EQ.NP1)) GO TO 635
         T = ONE
         GO TO 640
  635    T = WS(K)
  640    X(I) = SDOT(M,W(1,K),1,B(1,MP2),1)*T
  645 CONTINUE
      IMAX = ISAMAX(M,X,1)
      XNORM = ABS(X(IMAX))
      GO TO 655
  650 GO TO 1315
  655 CONTINUE
      GO TO 650
  660 CONTINUE
      CALL SSWAP (NP1,W (MP1,1),MDW,W (MP2,1),MDW)
      CALL SSWAP (MP1,B (MP1,1),MDB,B (MP2,1),MDB)
      CONVRG = .FALSE.
      FOUND = .TRUE.
      ONCE = .FALSE.
      FINITE = .TRUE.
  665 IF (.NOT.(ARTIF .AND. FINITE .AND. FOUND)) GO TO 725
      I99834 = 1
      GO TO 880
  670 IF (.NOT.(.NOT.FOUND)) GO TO 700
      IF (ONCE) GO TO 695
      I99918 = 7
      GO TO 410
  675 I99883 = 3
      GO TO 615
  680 DO 685 I=1,NP1
         IPTR(I) = IABS(IPTR(I))
  685 CONTINUE
      I99834 = 2
      GO TO 880
  690 ONCE = .TRUE.
  695 CONTINUE
  700 IF (.NOT.(FOUND)) GO TO 720
      I99824 = 1
      GO TO 1055
  705 IF (.NOT.(FINITE)) GO TO 715
      ARTIF = .NOT.IPTR(ILEAVE).EQ.NP1
      ONCE = .FALSE.
      I99820 = 1
      GO TO 1135
  710 CONTINUE
  715 CONTINUE
  720 GO TO 665
  725 IF (FINITE) GO TO 730
      IER = 134
      GO TO 9000
  730 IF (.NOT.(ARTIF)) GO TO 735
      GO TO 755
  735 DO 740 I=1,NP1
         IPTR(I) = IABS(IPTR(I))
  740 CONTINUE
      CALL SSWAP (NP1,W (MP1,1),MDW,W (MP2,1),MDW)
      CALL SSWAP (MP1,B (MP1,1),MDB,B (MP2,1),MDB)
      GO TO 750
  745 GO TO 20
  750 CONTINUE
      GO TO 745
  755 CONTINUE
      IER = 135
      GO TO 9000
  760 CONTINUE
      CONVRG = .FALSE.
      FOUND = .TRUE.
      FINITE = .TRUE.
      ONCE = .FALSE.
  765 IF (.NOT.(FINITE .AND. FOUND)) GO TO 870
      I99834 = 3
      GO TO 880
  770 IF (.NOT.(.NOT.FOUND)) GO TO 800
      IF (ONCE) GO TO 795
      I99918 = 8
      GO TO 410
  775 I99883 = 4
      GO TO 615
  780 DO 785 I=1,NP1
         IPTR(I) = IABS(IPTR(I))
  785 CONTINUE
      I99834 = 4
      GO TO 880
  790 ONCE = .TRUE.
  795 CONTINUE
  800 IF (.NOT.(FOUND)) GO TO 865
      I99824 = 2
      GO TO 1055
  805 ONCE = .FALSE.
      IF (.NOT.(FINITE)) GO TO 815
      I99820 = 2
      GO TO 1135
  810 GO TO 860
  815 I99918 = 9
      GO TO 410
  820 I99883 = 5
      GO TO 615
  825 DO 830 I=1,NP1
         IPTR(I) = IABS(IPTR(I))
  830 CONTINUE
      I99834 = 5
      GO TO 880
  835 FINITE = .TRUE.
      IF (.NOT.(FOUND)) GO TO 855
      I99824 = 3
      GO TO 1055
  840 IF (.NOT.(FINITE)) GO TO 850
      I99820 = 3
      GO TO 1135
  845 CONTINUE
  850 CONTINUE
  855 CONTINUE
  860 CONTINUE
  865 GO TO 765
  870 IF (.NOT.(.NOT.FINITE)) GO TO 875
      IER = 134
  875 IF (.NOT.FINITE) GO TO 9000
      GO TO 30
  880 CONTINUE
      CALL SCOPY (M,B (MP1,1),MDB,Y,1)
      DO 885 II=1,M
         I = MP1-II
         Y(I) = Y(I)/B(I,I)
         CALL SAXPY (I-1,-Y (I),B (I,1),MDB,Y,1)
  885 CONTINUE
      IMAX = ISAMAX(M,Y,1)
      PNORM = ABS(Y(IMAX))
      FIRST = .TRUE.
      I = NP1
      COSTMN = ZERO
      T = ZERO
      IENTER = 0
  890 IF (.NOT.(I.GT.M)) GO TO 915
      K = IPTR(I)
      IF (.NOT.(K.GT.0)) GO TO 910
      COST = W(MP1,K)-SDOT(M,W(1,K),1,Y,1)
      COST = COST*WS(K)
      IF (IND(K).NE.0 .AND. MOD(IND(K),2).EQ.0) COST = -COST
      IF (.NOT.(FIRST)) GO TO 895
      COSTMN = COST
      IENTER = I
      FIRST = .FALSE.
      GO TO 905
  895 IF (.NOT.(COST.LT.COSTMN)) GO TO 900
      COSTMN = COST
      IENTER = I
  900 CONTINUE
  905 CONTINUE
  910 I = I-1
      GO TO 890
  915 IF (.NOT.(FIRST)) GO TO 920
      CONVRG = .TRUE.
      FOUND = .FALSE.
      GO TO 925
  920 T = AN*RERSQ
      T = T*FAC
      CONVRG = T.GE.-COSTMN*AD
      K = IPTR(IENTER)
      IMAX = ISAMAX(M,W(1,K),1)
      ANORM = ABS(W(IMAX,K)*WS(K))
      SIZE = ANORM*PNORM+ABS(W(MP1,K)*WS(K))
      T = SIZE*TOL*AN
      T = T*FAC
      FOUND = .NOT.(T.GE.-COSTMN*AD)
      DON = M
      IF (.NOT.FOUND) FOUND = .NOT.(T/DON.GE.-COSTMN*AD)
  925 GO TO 935
  930 GO TO 1325
  935 CONTINUE
      IF (.NOT.(KPRINT.GT.1)) GO TO 950
      IF (.NOT.(FOUND)) GO TO 940
      TEMPL = ONE
      GO TO 945
  940 TEMPL = ZERO
  945 CONTINUE
      IF (AD.GT.ZERO) T = T/AD
      DUMMY(1) = COSTMN
      DUMMY(2) = T
  950 GO TO 930
  955 CONTINUE
      K = IPTR(IENTER)
      CALL SCOPY (M,W (1,K),1,Y,1)
      CALL SSCAL (M,WS (K),Y,1)
      DO 960 I=1,M
         Y(I) = Y(I)/B(I,I)
         IF (I.EQ.M) GO TO 960
         CALL SAXPY (M-I,-Y (I),B (I+1,I),1,Y (I+1),1)
  960 CONTINUE
      DO 965 I=1,M
         Y(I) = Y(I)/D(I)
  965 CONTINUE
      DO 970 II=1,M
         I = MP1-II
         Y(I) = Y(I)/B(I,I)
         CALL SAXPY (I-1,-Y (I),B (I,1),MDB,Y,1)
  970 CONTINUE
      CALL SCOPY (M,Y,1,B (1,MP2),1)
      DO 985 I=1,M
         L = IPTR(I)
         IF (.NOT.(L.EQ.NP1)) GO TO 975
         T = ONE
         GO TO 980
  975    T = WS(L)
  980    Y(I) = SDOT(M,W(1,L),1,B(1,MP2),1)*T
  985 CONTINUE
      IF (.NOT.(IND(K).NE.0 .AND. MOD(IND(K),2).EQ.0 .AND. FOUND)) GO
     *TO 995
      DO 990 I=1,M
         Y(I) = -Y(I)
  990 CONTINUE
  995 IMAX = ISAMAX(M,Y,1)
      YNORM = ABS(Y(IMAX))
      T = YNORM*TOL*AN
      DON = M
      T = T*DON
      FIRST = .TRUE.
      ILEAVE = 1
      L = IPTR(IENTER)
      IF (.NOT.(IND(L).NE.0)) GO TO 1000
      BOUND = BND(L)*BSCAL/WS(L)
      FIRST = .FALSE.
      THETA = BOUND
      ILEAVE = IENTER
 1000 DO 1050 I=1,M
         IF (.NOT.(AD*Y(I).GT.T)) GO TO 1020
         QUOT = X(I)/Y(I)
         IF (.NOT.(FIRST)) GO TO 1005
         THETA = QUOT
         ILEAVE = I
         FIRST = .FALSE.
         GO TO 1015
 1005    IF (.NOT.(QUOT.LT.THETA)) GO TO 1010
         THETA = QUOT
         ILEAVE = I
 1010    CONTINUE
 1015    GO TO 1045
 1020    L = IPTR(I)
         IF (.NOT.(IND(L).NE.0 .AND. AD*ABS(Y(I)).GT.T)) GO TO 1040
         BOUND = BND(L)*BSCAL/WS(L)
         QUOT = (X(I)-BOUND)/Y(I)
         IF (.NOT.(FIRST)) GO TO 1025
         THETA = QUOT
         ILEAVE = -I
         FIRST = .FALSE.
         GO TO 1035
 1025    IF (.NOT.(QUOT.LT.THETA)) GO TO 1030
         THETA = QUOT
         ILEAVE = -I
 1030    CONTINUE
 1035    CONTINUE
 1040    CONTINUE
 1045    CONTINUE
 1050 CONTINUE
      FINITE = .NOT.FIRST
      GO TO 1320
 1055 CONTINUE
      I99768 = 1
      GO TO 955
 1060 I = IABS(ILEAVE)
      IF (.NOT.(ARTIF .AND. IPTR(I).NE.NP1 .AND. FINITE)) GO TO 1110
      IARTIF = 1
 1065 IF (IPTR(IARTIF).EQ.NP1) GO TO 1070
      IARTIF = IARTIF+1
      GO TO 1065
 1070 T = X(IARTIF)-THETA*Y(IARTIF)
      SIZE = (AN*XNORM*FAC)**2
      IF (.NOT.((AD*T)**2+SIZE.EQ.SIZE)) GO TO 1105
      I99918 = 10
      GO TO 410
 1075 I99883 = 6
      GO TO 615
 1080 I99768 = 2
      GO TO 955
 1085 CALL SCOPY (M,W (1,NP2),1,B (1,MP2),1)
      K = IPTR(IENTER)
      CALL SAXPY (M,-THETA*WS (K),W (1,K),1,B (1,MP2),1)
      DO 1095 I=1,M
         K = IPTR(I)
         IF (.NOT.(I.NE.IARTIF)) GO TO 1090
         T = WS(K)*AMAX1(ZERO,X(I)-THETA*Y(I))
         CALL SAXPY (M,-T,W (1,K),1,B (1,MP2),1)
 1090    CONTINUE
 1095 CONTINUE
      T = SASUM(M,B(1,MP2),1)
      SIZE = BNORM*(XNORM+THETA*YNORM)+SASUM(M,W(1,NP2),1)
      IF (.NOT.(T**2+SIZE**2.EQ.SIZE**2)) GO TO 1100
      ILEAVE = IARTIF
      ARTIF = .FALSE.
 1100 CONTINUE
 1105 CONTINUE
 1110 GO TO 1310
 1115 CONTINUE
      IF (.NOT.(KPRINT.GT.1)) GO TO 1130
      IF (.NOT.(FINITE)) GO TO 1120
      TEMPL = ONE
      GO TO 1125
 1120 TEMPL = ZERO
 1125 CONTINUE
 1130 GO TO 1185
 1135 CONTINUE
      T = AD*THETA*YNORM
      IF (.NOT.(T.LE.RERSQ*AN*XNORM)) GO TO 1165
      I = IABS(ILEAVE)
      IF (.NOT.(IABS(IPTR(I)).NE.NP1)) GO TO 1160
      FIRST = .TRUE.
      DO 1155 I=1,M
         T = YNORM*AN*TOL*FAC
         IF (.NOT.(AD*Y(I).GE.T .AND. AD*ABS(X(I)).LE.XNORM*AN*TOL*FAC))
     *   GO TO 1150
         IF (.NOT.(FIRST)) GO TO 1140
         YMAX = Y(I)
         ILEAVE = I
         FIRST = .FALSE.
 1140    IF (.NOT.(YMAX.LT.Y(I))) GO TO 1145
         YMAX = Y(I)
         ILEAVE = I
 1145    CONTINUE
 1150    CONTINUE
 1155 CONTINUE
 1160 I = IABS(ILEAVE)
      IPTR(I) = -IABS(IPTR(I))
      GO TO 1180
 1165 I = NP1
 1170 IF (.NOT.(I.GT.M)) GO TO 1175
      IPTR(I) = IABS(IPTR(I))
      I = I-1
      GO TO 1170
 1175 CONTINUE
 1180 GO TO 1115
 1185 CALL SAXPY (M,-THETA,Y,1,X,1)
      IF (.NOT.(ILEAVE.LT.0)) GO TO 1190
      ILEAVE = -ILEAVE
      L = IABS(IPTR(ILEAVE))
      BOUND = BND(L)*BSCAL
      CALL SAXPY (M,-BOUND,W (1,L),1,W (1,NP2),1)
      IND(L) = IND(L)+1
 1190 IF (.NOT.(IENTER.NE.ILEAVE)) GO TO 1245
      X(ILEAVE) = THETA
      K = IPTR(IENTER)
      IF (.NOT.(IND(K).GT.0 .AND. MOD(IND(K),2).EQ.0)) GO TO 1195
      BOUND = BND(K)*BSCAL
      CALL SAXPY (M,BOUND,W (1,K),1,W (1,NP2),1)
      X(ILEAVE) = BOUND/WS(K)-THETA
      IND(K) = IND(K)+1
 1195 D(MP2) = WS(K)**2
      CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
      BNORM = AMAX1(BNORM,SASUM(M,W(1,K),1))
      CALL ZX4LQ (B,M,M,MDB,D,IFLAG)
      I99868 = 4
      GO TO 1280
 1200 K = IABS(IPTR(ILEAVE))
      T = WS(K)
      IF (.NOT.(K.EQ.NP1)) GO TO 1205
      T = ONE
      IPTR(ILEAVE) = 0
      ARTIF = .FALSE.
 1205 ITEMP = IPTR(ILEAVE)
      IPTR(ILEAVE) = IPTR(IENTER)
      IPTR(IENTER) = ITEMP
      D(MP2) = -T**2
      DROP = .TRUE.
      CALL SCOPY (MP2,W (1,K),1,B (1,MP2),1)
      CALL ZX4LQ (B,M,M,MDB,D,IFLAG)
      IF (.NOT.(0.LT.IFLAG)) GO TO 1220
      I99918 = 11
      GO TO 410
 1210 I99883 = 7
      GO TO 615
 1215 GO TO 1240
 1220 I99868 = 5
      GO TO 1280
 1225 IMAX = ISAMAX(M,X,1)
      XNORM = AMAX1(XNORM,ABS(X(IMAX)))
      AN = D(1)*B(1,1)**2
      AD = AN
      I = 2
 1230 IF (.NOT.(I.LE.M)) GO TO 1235
      T = D(I)*B(I,I)**2
      AN = AMAX1(AN,T)
      AD = AMIN1(AD,T)
      I = I+1
      GO TO 1230
 1235 AN = SQRT(AN)
      AD = SQRT(AD)
      T = D(MP1)*B(MP1,MP1)**2
      DON = M*M
      RERSQ = AMAX1(T,(EPS*DON)**2)/BNORM**2
      RERSQ = SQRT(RERSQ)
      DON = M+M
      TOL = TOL+EPS*DON
 1240 GO TO 1255
 1245 L = IABS(IPTR(IENTER))
      IF (.NOT.(IND(L).GT.0)) GO TO 1250
      T = BND(L)*BSCAL
      IF (MOD(IND(L),2).EQ.1) T = -T
      CALL SAXPY (M,T,W (1,L),1,W (1,NP2),1)
      IND(L) = IND(L)+1
 1250 CONTINUE
 1255 NITER = NITER+1
      IF (.NOT.(MOD(NITER,2*N/3).EQ.0)) GO TO 1270
      DON = M
      FAC = (FAC+FAC)*DON
      I99918 = 12
      GO TO 410
 1260 I99883 = 8
      GO TO 615
 1265 CONTINUE
 1270 IF (.NOT.(NITER.GT.N3)) GO TO 1275
      IER = 131
      GO TO 9000
 1275 GO TO 1330
 1280 CONTINUE
      T = SQRT(D(MP1))
      B(MP1,MP1) = T*B(MP1,MP1)
      B(MP2,MP1) = T*B(MP2,MP1)
      T = B(MP1,MP1)**2+D(MP2)*B(MP1,MP2)**2
      T = AMAX1(ZERO,T)
      B(MP1,MP1) = SQRT(T)
      T = B(MP2,MP1)**2+D(MP2)*B(MP2,MP2)**2
      T = AMAX1(ZERO,T)
      B(MP2,MP1) = SQRT(T)
      GO TO 1345
 1285 CONTINUE
      CALL SSCAL (M,ONE/COSTSC,Y,1)
      DO 1290 I=1,M
         T = X(I)
         IF (AD*T.LE.TOL*AN*XNORM) T = ZERO
         L = IPTR(I)
         X(I) = T*WS(L)/BSCAL
 1290 CONTINUE
      WS(1) = ZERO
      CALL SCOPY (N,WS,0,WS,1)
      DO 1295 I=1,M
         L = IPTR(I)
         WS(L) = X(I)
 1295 CONTINUE
      DO 1305 I=1,N
         IF (.NOT.(IND(I).NE.0 .AND. MOD(IND(I),2).EQ.0)) GO TO 1300
         WS(I) = AMAX1(ZERO,BND(I)-WS(I))
 1300    CONTINUE
 1305 CONTINUE
      IPTR(NP1) = NITER
      IPTR(NP2) = NREDC
      GO TO 9005
 1310 GO TO (705, 805, 840), I99824
 1315 GO TO (460, 595, 680, 780, 825, 1080, 1215, 1265), I99883
 1320 GO TO (1060, 1085), I99768
 1325 GO TO (670, 690, 770, 790, 835), I99834
 1330 GO TO (710, 810, 845), I99820
 1335 GO TO (280, 310, 340, 370, 530, 565, 675, 775, 820, 1075, 1210,
     *1260), I99918
 1340 GO TO (390, 445), I99895
 1345 GO TO (515, 560, 575, 1200, 1225), I99868
 9000 IF (IER.EQ.0) GO TO 9005
      CALL UERTST (IER,6HZX4LP )
 9005 RETURN
      END

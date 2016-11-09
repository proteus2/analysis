C   IMSL ROUTINE NAME   - ZXMWE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY ZXMWD
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZXMJN
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXMWE (FUNCT,N,NSIG,MAXFN,IOPT,X,H,G,F,W,A,B,XX,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSIG,MAXFN,IOPT,IER
      REAL               X(N),H(1),G(N),F,W(1),A(1),B(1),XX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDIFF,IFN,IGG,IG,II,IJ,IM1,IR,IS,ITN,I,JB,JJ,
     *                   JNT,JP1,J,KJ,K,LINK,L,NJ,NM1,NP1
      REAL               AEPS,ALPHA,AX,DF,DGS,DIFF,EPS,F11,F12,F1,F21,
     *                   F22,F2,FF,FIVE,GHH,GNRM,GS0,GYS,H2,HALF,HHH,HH,
     *                   HJJ,HMAX,HMIN,ONE,P1,RELX,REPS,SEVEN,SIG,TEN,
     *                   TOT,TWELVE,V,ZERO,ZZ,Z
      DATA               REPS /Z3C100000/,AX /0.1/
      DATA               ZERO /0.0/,ONE /1.0/,HALF /0.5/,SEVEN
     *                   /7.0/,FIVE /5.0/,TWELVE /12.0/,
     *                   TEN /10.0/,P1 /0.1/
C                                  INITIALIZATION
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      HH = SQRT(REPS)
      H2 = SQRT(HH)
      EPS = TEN**(-NSIG)
      IG = N
      IGG = N+N
      IS = IGG
      IDIFF = 1
      IR = N
      W(1) = -ONE
      W(2) = ZERO
      W(3) = ZERO
C                                  EVALUATE FUNCTION AT STARTING POINT
      DO 5 I=1,N
         G(I) = X(I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
    5 CONTINUE
      CALL FUNCT(N,XX,F)
      IFN = 1
      IF (IOPT.EQ.1) GO TO 50
C                                  SET OFF-DIAGONAL ELEMENTS OF H TO 0.0
      IF (N.EQ.1) GO TO 20
      IJ = 2
      DO 15 I=2,N
         DO 10 J=2,I
            H(IJ) = ZERO
            IJ = IJ+1
   10    CONTINUE
         IJ = IJ+1
   15 CONTINUE
   20 IF (IOPT.NE.0) GO TO 30
C                                  SET DIAGONAL ELEMENTS OF H TO ONE
      IJ = 0
      DO 25 I=1,N
         IJ = IJ+I
         H(IJ) = ONE
   25 CONTINUE
      GO TO 95
C                                  GET DIAGONAL ELEMENTS OF HESSIAN
   30 IM1 = 1
      NM1 = 1
      NP1 = N+1
      DO 35 I=2,NP1
         HHH = H2*AMAX1(ABS(X(IM1)),AX)
         G(IM1) = X(IM1)+HHH
         XX(IM1) = A(IM1)+(B(IM1)-A(IM1))*SIN(G(IM1))**2
         CALL FUNCT(N,XX,F2)
         G(IM1) = X(IM1)-HHH
         XX(IM1) = A(IM1)+(B(IM1)-A(IM1))*SIN(G(IM1))**2
         CALL FUNCT(N,XX,FF)
         H(NM1) = (FF-F+F2-F)/(HHH*HHH)
         G(IM1) = X(IM1)
         XX(IM1) = A(IM1)+(B(IM1)-A(IM1))*SIN(G(IM1))**2
         IM1 = I
         NM1 = I+NM1
   35 CONTINUE
      IFN = IFN+N+N
      IF (IOPT.NE.3 .OR. N.EQ.1) GO TO 50
C                                  GET THE REST OF THE HESSIAN
      JJ = 1
      II = 2
      DO 45 I=2,N
         GHH = H2*AMAX1(ABS(X(I)),AX)
         DO 40 J=1,JJ
            HHH = H2*AMAX1(ABS(X(J)),AX)
            G(I) = X(I)+GHH
            XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
            G(J) = X(J)+HHH
            XX(J) = A(J)+(B(J)-A(J))*SIN(G(J))**2
            CALL FUNCT(N,XX,F22)
            G(I) = X(I)-GHH
            XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
            CALL FUNCT(N,XX,F12)
            G(J) = X(J)-HHH
            XX(J) = A(J)+(B(J)-A(J))*SIN(G(J))**2
            CALL FUNCT(N,XX,F11)
            G(I) = X(I)+GHH
            XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
            CALL FUNCT(N,XX,F21)
            H(II) = (F22-F21-F12+F11)/(4.*HHH*GHH)
            G(J) = X(J)
            XX(J) = A(J)+(B(J)-A(J))*SIN(G(J))**2
            II = II+1
   40    CONTINUE
         G(I) = X(I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
         JJ = JJ+1
         II = II+1
   45 CONTINUE
      IFN = IFN+((N*N-N)*2)
C                                  ADD MULTIPLE OF IDENTITY TO
C                                  MAKE DIAGONAL ELEMENTS POSITIVE
   50 HMIN = H(1)
      HMAX = H(1)
      NM1 = 1
      DO 55 I=1,N
         HMIN = AMIN1(HMIN,H(NM1))
         HMAX = AMAX1(HMAX,H(NM1))
         NM1 = NM1+I+1
   55 CONTINUE
      HMIN = AMAX1(0.01*(ABS(HMAX)+ABS(HMIN))-HMIN,0.0)
      NM1 = 1
      DO 60 I=1,N
         H(NM1) = H(NM1)+HMIN
         NM1 = NM1+I+1
   60 CONTINUE
C                                  FACTOR H TO L*D*L-TRANSPOSE
      IR = N
      IF (N.GT.1) GO TO 65
      IF (H(1).GT.ZERO) GO TO 95
      H(1) = ZERO
      IR = 0
      GO TO 90
   65 NM1 = N-1
      JJ = 0
      DO 85 J=1,N
         JP1 = J+1
         JJ = JJ+J
         HJJ = H(JJ)
         IF (HJJ.GT.ZERO) GO TO 70
         H(JJ) = ZERO
         IR = IR-1
         GO TO 85
   70    IF (J.EQ.N) GO TO 85
         IJ = JJ
         L = 0
         DO 80 I=JP1,N
            L = L+1
            IJ = IJ+I-1
            V = H(IJ)/HJJ
            KJ = IJ
            DO 75 K=I,N
               H(KJ+L) = H(KJ+L)-H(KJ)*V
               KJ = KJ+K
   75       CONTINUE
            H(IJ) = V
   80    CONTINUE
   85 CONTINUE
   90 IF (IR.EQ.N) GO TO 95
      IER = 129
      GO TO 9000
   95 ITN = 0
      DF = -ONE
C                                  EVALUATE GRADIENT W(IG+I),I=1,...,N
  100 LINK = 1
      GO TO 275
  105 CONTINUE
C                                  BEGIN ITERATION LOOP
      IF (IFN.GE.MAXFN) GO TO 240
      ITN = ITN+1
      DO 110 I=1,N
         W(I) = -W(IG+I)
  110 CONTINUE
C                                  DETERMINE SEARCH DIRECTION W
C                                    BY SOLVING H*W = -G WHERE
C                                    H = L*D*L-TRANSPOSE
      IF (IR.LT.N) GO TO 140
C                                  N .EQ. 1
      G(1) = W(1)
      IF (N.GT.1) GO TO 115
      W(1) = W(1)/H(1)
      GO TO 140
C                                  N .GT. 1
  115 II = 1
C                                  SOLVE L*W = -G
      DO 125 I=2,N
         IJ = II
         II = II+I
         V = W(I)
         IM1 = I-1
         DO 120 J=1,IM1
            IJ = IJ+1
            V = V-H(IJ)*W(J)
  120    CONTINUE
         G(I) = V
         W(I) = V
  125 CONTINUE
C                                  SOLVE (D*LT)*Z = W WHERE
C                                  LT = L-TRANSPOSE
      W(N) = W(N)/H(II)
      JJ = II
      NM1 = N-1
      DO 135 NJ=1,NM1
C                                  J = N-1,N-2,...,1
         J = N-NJ
         JP1 = J+1
         JJ = JJ-JP1
         V = W(J)/H(JJ)
         IJ = JJ
         DO 130 I=JP1,N
            IJ = IJ+I-1
            V = V-H(IJ)*W(I)
  130    CONTINUE
         W(J) = V
  135 CONTINUE
C                                  DETERMINE STEP LENGTH ALPHA
  140 RELX = ZERO
      GS0 = ZERO
      DO 145 I=1,N
         W(IS+I) = W(I)
         DIFF = ABS(W(I))/AMAX1(ABS(X(I)),AX)
         RELX = AMAX1(RELX,DIFF)
         GS0 = GS0+W(IG+I)*W(I)
  145 CONTINUE
      IF (RELX.EQ.ZERO) GO TO 245
      AEPS = EPS/RELX
      IER = 130
      IF (GS0.GE.ZERO) GO TO 245
      IF (DF.EQ.ZERO) GO TO 245
      IER = 0
      ALPHA = (-DF-DF)/GS0
      IF (ALPHA.LE.ZERO) ALPHA = ONE
      ALPHA = AMIN1(ALPHA,ONE)
      IF (IDIFF.EQ.2) ALPHA = AMAX1(P1,ALPHA)
      FF = F
      TOT = ZERO
      JNT = 0
C                                  SEARCH ALONG  X+ALPHA*W
  150 IF (IFN.GE.MAXFN) GO TO 240
      DO 155 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(W(I))**2
  155 CONTINUE
      CALL FUNCT(N,XX,F1)
      IFN = IFN+1
      IF (F1.GE.F) GO TO 180
      F2 = F
      TOT = TOT+ALPHA
  160 IER = 0
      F = F1
      DO 165 I=1,N
         X(I) = W(I)
  165 CONTINUE
      IF (JNT-1) 170, 200, 205
  170 IF (IFN.GE.MAXFN) GO TO 240
      DO 175 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(W(I))**2
  175 CONTINUE
      CALL FUNCT(N,XX,F1)
      IFN = IFN+1
      IF (F1.GE.F) GO TO 205
      IF (F1+F2.GE.F+F .AND. SEVEN*F1+FIVE*F2.GT.TWELVE*F) JNT = 2
      TOT = TOT+ALPHA
      ALPHA = ALPHA+ALPHA
      GO TO 160
  180 CONTINUE
      IF (F.EQ.FF .AND. IDIFF.EQ.2 .AND. RELX.GT.EPS) IER = 130
      IF (ALPHA.LT.AEPS) GO TO 245
      IF (IFN.GE.MAXFN) GO TO 240
      ALPHA = HALF*ALPHA
      DO 185 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(W(I))**2
  185 CONTINUE
      CALL FUNCT(N,XX,F2)
      IFN = IFN+1
      IF (F2.GE.F) GO TO 195
      TOT = TOT+ALPHA
      IER = 0
      F = F2
      DO 190 I=1,N
         X(I) = W(I)
  190 CONTINUE
      GO TO 200
  195 Z = P1
      IF (F1+F.GT.F2+F2) Z = ONE+HALF*(F-F1)/(F+F1-F2-F2)
      Z = AMAX1(P1,Z)
      ALPHA = Z*ALPHA
      JNT = 1
      GO TO 150
  200 IF (TOT.LT.AEPS) GO TO 245
  205 ALPHA = TOT
C                                  SAVE OLD GRADIENT
      DO 210 I=1,N
         W(I) = W(IG+I)
  210 CONTINUE
C                                  EVALUATE GRADIENT W(IG+I), I=1,...,N
      LINK = 2
      GO TO 275
  215 IF (IFN.GE.MAXFN) GO TO 240
      GYS = ZERO
      DO 220 I=1,N
         GYS = GYS+W(IG+I)*W(IS+I)
         W(IGG+I) = W(I)
  220 CONTINUE
      DF = FF-F
      DGS = GYS-GS0
      IF (DGS.LE.ZERO) GO TO 105
      IF (DGS+ALPHA*GS0.GT.ZERO) GO TO 230
C                                  UPDATE HESSIAN H USING
C                                    COMPLEMENTARY DFP FORMULA
      SIG = ONE/GS0
      IR = -IR
      CALL ZXMJN(H,N,W,SIG,G,IR,0,ZERO)
      DO 225 I=1,N
         G(I) = W(IG+I)-W(IGG+I)
  225 CONTINUE
      SIG = ONE/(ALPHA*DGS)
      IR = -IR
      CALL ZXMJN(H,N,G,SIG,W,IR,0,ZERO)
      GO TO 105
C                                  UPDATE HESSIAN USING
C                                    DFP FORMULA
  230 ZZ = ALPHA/(DGS-ALPHA*GS0)
      SIG = -ZZ
      CALL ZXMJN(H,N,W,SIG,G,IR,0,REPS)
      Z = DGS*ZZ-ONE
      DO 235 I=1,N
         G(I) = W(IG+I)+Z*W(IGG+I)
  235 CONTINUE
      SIG = ONE/(ZZ*DGS*DGS)
      CALL ZXMJN(H,N,G,SIG,W,IR,0,ZERO)
      GO TO 105
  240 IER = 131
C                                  MAXFN FUNCTION EVALUATIONS
      GO TO 250
  245 IF (IDIFF.EQ.2) GO TO 250
C                                  CHANGE TO CENTRAL DIFFERENCES
      IDIFF = 2
      GO TO 100
  250 IF (RELX.GT.EPS .AND. IER.EQ.0) GO TO 100
C                                  MOVE GRADIENT TO G AND RETURN
      GNRM = ZERO
      DO 255 I=1,N
         G(I) = W(IG+I)
         GNRM = GNRM+G(I)*G(I)
  255 CONTINUE
      GNRM = SQRT(GNRM)
      W(1) = GNRM
      W(2) = IFN
      W(3) = -ALOG10(AMAX1(REPS,RELX))
C                                  COMPUTE H = L*D*L-TRANSPOSE
      IF (N.EQ.1) GO TO 9000
      NP1 = N+1
      NM1 = N-1
      JJ = (N*(NP1))/2
      DO 270 JB=1,NM1
         JP1 = NP1-JB
         JJ = JJ-JP1
         HJJ = H(JJ)
         IJ = JJ
         L = 0
         DO 265 I=JP1,N
            L = L+1
            IJ = IJ+I-1
            V = H(IJ)*HJJ
            KJ = IJ
            DO 260 K=I,N
               H(KJ+L) = H(KJ+L)+H(KJ)*V
               KJ = KJ+K
  260       CONTINUE
            H(IJ) = V
  265    CONTINUE
         HJJ = H(JJ)
  270 CONTINUE
      GO TO 9000
C                                  EVALUATE GRADIENT
  275 IF (IDIFF.EQ.2) GO TO 290
C                                  FORWARD DIFFERENCES
C                                    GRADIENT = W(IG+I), I=1,...,N
      DO 280 I=1,N
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  280 CONTINUE
      DO 285 I=1,N
         Z = HH*AMAX1(ABS(X(I)),AX)
         ZZ = X(I)
         X(I) = ZZ+Z
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
         CALL FUNCT(N,XX,F1)
         W(IG+I) = (F1-F)/Z
         X(I) = ZZ
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  285 CONTINUE
      IFN = IFN+N
      GO TO (105, 215), LINK
C                                  CENTRAL DIFFERENCES
C                                    GRADIENT = W(IG+I), I=1,...,N
  290 DO 295 I=1,N
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  295 CONTINUE
      DO 300 I=1,N
         Z = HH*AMAX1(ABS(X(I)),AX)
         ZZ = X(I)
         X(I) = ZZ+Z
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
         CALL FUNCT(N,XX,F1)
         X(I) = ZZ-Z
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
         CALL FUNCT(N,XX,F2)
         W(IG+I) = (F1-F2)/(Z+Z)
         X(I) = ZZ
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  300 CONTINUE
      IFN = IFN+N+N
      GO TO (105, 215), LINK
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HZXMWD )
 9005 RETURN
      END

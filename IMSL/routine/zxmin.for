C   IMSL ROUTINE NAME   - ZXMIN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MINIMUM OF A FUNCTION OF N VARIABLES USING
C                           A QUASI-NEWTON METHOD
C
C   USAGE               - CALL ZXMIN (FUNCT,N,NSIG,MAXFN,IOPT,X,H,G,F,
C                           W,IER)
C
C   ARGUMENTS    FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
C                           THE FUNCTION F FOR GIVEN PARAMETER VALUES
C                           X(1),X(2),...,X(N).
C                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
C                           CALL FUNCT(N,X,F)
C                           WHERE X IS A VECTOR OF LENGTH N.
C                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT
C                           IN THE CALLING PROGRAM. FUNCT MUST NOT
C                           ALTER THE VALUES OF X(I),I=1,...,N OR N.
C                N      - THE NUMBER OF PARAMETERS (I.E., THE LENGTH
C                           OF X) (INPUT)
C                NSIG   - CONVERGENCE CRITERION. (INPUT). THE NUMBER
C                           OF DIGITS OF ACCURACY REQUIRED IN THE
C                           PARAMETER ESTIMATES.
C                           THIS CONVERGENCE CONDITION IS SATISIFIED IF
C                           ON TWO SUCCESSIVE ITERATIONS, THE PARAMETER
C                           ESTIMATES (I.E.,X(I), I=1,...,N) AGREE,
C                           COMPONENT BY COMPONENT, TO NSIG DIGITS.
C                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E.,
C                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT)
C                IOPT   - OPTIONS SELECTOR. (INPUT)
C                         IOPT = 0 CAUSES ZXMIN TO INITIALIZE THE
C                           HESSIAN MATRIX H TO THE IDENTITY MATRIX.
C                         IOPT = 1 INDICATES THAT H HAS BEEN INITIALIZED
C                           BY THE USER TO A POSITIVE DEFINITE MATRIX.
C                         IOPT = 2 CAUSES ZXMIN TO COMPUTE THE DIAGONAL
C                           VALUES OF THE HESSIAN MATRIX AND SET H TO
C                           A DIAGONAL MATRIX CONTAINING THESE VALUES.
C                         IOPT = 3 CAUSES ZXMIN TO COMPUTE AN ESTIMATE
C                           OF THE HESSIAN IN H.
C                X      - VECTOR OF LENGTH N CONTAINING PARAMETER
C                           VALUES.
C                         ON INPUT, X MUST CONTAIN THE INITIAL
C                           PARAMETER ESTIMATES.
C                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER
C                           ESTIMATES AS DETERMINED BY ZXMIN.
C                H      - VECTOR OF LENGTH N*(N+1)/2 CONTAINING AN
C                           ESTIMATE OF THE HESSIAN MATRIX
C                           D**2F/(DX(I)DX(J)), I,J=1,...,N.
C                           H IS STORED IN SYMMETRIC STORAGE MODE.
C                         ON INPUT, IF IOPT = 0, 2, OR 3 ZXMIN INITIA-
C                           LIZES H. AN INITIAL SETTING OF H BY THE
C                           USER IS INDICATED BY IOPT=1.
C                           H MUST BE POSITIVE DEFINITE. IF IT IS NOT,
C                           A TERMINAL ERROR OCCURS.
C                         ON OUTPUT, H CONTAINS AN ESTIMATE OF THE
C                           HESSIAN AT THE FINAL PARAMETER ESTIMATES
C                           (I.E., AT X(1),X(2),...,X(N))
C                G      - A VECTOR OF LENGTH N CONTAINING AN ESTIMATE
C                           OF THE GRADIENT DF/DX(I),I=1,...,N AT THE
C                           FINAL PARAMETER ESTIMATES. (OUTPUT)
C                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION
C                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)
C                W      - A VECTOR OF LENGTH 3*N USED AS WORKING SPACE.
C                         ON OUTPUT, WORK(I), CONTAINS FOR
C                           I = 1, THE NORM OF THE GRADIENT (I.E.,
C                             SQRT(G(1)**2+G(2)**2+...+G(N)**2))
C                           I = 2, THE NUMBER OF FUNCTION EVALUATIONS
C                             PERFORMED.
C                           I = 3, AN ESTIMATE OF THE NUMBER OF
C                             SIGNIFICANT DIGITS IN THE FINAL
C                             PARAMETER ESTIMATES.
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 IMPLIES THAT THE INITIAL HESSIAN
C                             USED BY ZXMIN IS NOT POSITIVE DEFINITE,
C                             EVEN AFTER ADDING A MULTIPLE OF THE
C                             IDENTITY TO MAKE ALL DIAGONAL ELEMENTS
C                             POSITIVE.
C                           IER = 130 IMPLIES THAT THE ITERATION WAS
C                             TERMINATED DUE TO ROUNDING ERRORS
C                             BECOMING DOMINANT. THE PARAMETER
C                             ESTIMATES HAVE NOT BEEN DETERMINED TO
C                             NSIG DIGITS.
C                           IER = 131 IMPLIES THAT THE ITERATION WAS
C                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZXMJN
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXMIN (FUNCT,N,NSIG,MAXFN,IOPT,X,H,G,F,W,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSIG,MAXFN,IOPT,IER
      REAL               X(N),G(N),H(1),F,W(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IG,IGG,IS,IDIFF,IR,IJ,I,J,NM1,JJ,JP1,L,KJ,K,
     *                   IFN,LINK,ITN,II,IM1,JNT,NP1,JB,NJ
      REAL               REPS,AX,ZERO,ONE,HALF,SEVEN,FIVE,TWELVE,TEN,HH,
     *                   EPS,HJJ,V,DF,RELX,GS0,DIFF,AEPS,ALPHA,FF,TOT,
     *                   F1,F2,Z,GYS,DGS,SIG,ZZ,GNRM,P1,HHH,GHH,H2,F11,
     *                   F12,F21,F22,HMAX,HMIN
      DATA               REPS/Z3C100000/,AX/0.1/
      DATA               ZERO/0.0/,ONE/1.0/,HALF/0.5/,
     *                   SEVEN/7.0/,FIVE/5.0/,TWELVE/12.0/,
     *                   TEN/10.0/,P1/0.1/
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
    5 CONTINUE
      CALL FUNCT(N,G,F)
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
         CALL FUNCT(N,G,F2)
         G(IM1) = X(IM1)-HHH
         CALL FUNCT(N,G,FF)
         H(NM1) = (FF-F+F2-F)/(HHH*HHH)
         G(IM1) = X(IM1)
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
            G(J) = X(J)+HHH
            CALL FUNCT(N,G,F22)
            G(I) = X(I)-GHH
            CALL FUNCT(N,G,F12)
            G(J) = X(J)-HHH
            CALL FUNCT(N,G,F11)
            G(I) = X(I)+GHH
            CALL FUNCT(N,G,F21)
            H(II) = (F22-F21-F12+F11)/(4.*HHH*GHH)
            G(J) = X(J)
            II = II+1
   40    CONTINUE
         G(I) = X(I)
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
      GO TO 280
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
  155 CONTINUE
      CALL FUNCT(N,W,F1)
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
  175 CONTINUE
      CALL FUNCT(N,W,F1)
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
  185 CONTINUE
      CALL FUNCT(N,W,F2)
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
      GO TO 280
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
  250 IF (IER.NE.0) GO TO 255
      IF (RELX.LE.EPS) GO TO 255
      GO TO 100
C                                  MOVE GRADIENT TO G AND RETURN
  255 GNRM = ZERO
      DO 260 I=1,N
         G(I) = W(IG+I)
         GNRM = GNRM+G(I)*G(I)
  260 CONTINUE
      GNRM = SQRT(GNRM)
      W(1) = GNRM
      W(2) = IFN
      W(3) = -ALOG10(AMAX1(REPS,RELX))
C                                  COMPUTE H = L*D*L-TRANSPOSE
      IF (N.EQ.1) GO TO 9000
      NP1 = N+1
      NM1 = N-1
      JJ = (N*(NP1))/2
      DO 275 JB=1,NM1
         JP1 = NP1-JB
         JJ = JJ-JP1
         HJJ = H(JJ)
         IJ = JJ
         L = 0
         DO 270 I=JP1,N
            L = L+1
            IJ = IJ+I-1
            V = H(IJ)*HJJ
            KJ = IJ
            DO 265 K=I,N
               H(KJ+L) = H(KJ+L)+H(KJ)*V
               KJ = KJ+K
  265       CONTINUE
            H(IJ) = V
  270    CONTINUE
         HJJ = H(JJ)
  275 CONTINUE
      GO TO 9000
C                                  EVALUATE GRADIENT
  280 IF (IDIFF.EQ.2) GO TO 290
C                                  FORWARD DIFFERENCES
C                                    GRADIENT = W(IG+I), I=1,...,N
      DO 285 I=1,N
         Z = HH*AMAX1(ABS(X(I)),AX)
         ZZ = X(I)
         X(I) = ZZ+Z
         CALL FUNCT(N,X,F1)
         W(IG+I) = (F1-F)/Z
         X(I) = ZZ
  285 CONTINUE
      IFN = IFN+N
      GO TO (105, 215), LINK
C                                  CENTRAL DIFFERENCES
C                                    GRADIENT = W(IG+I), I=1,...,N
  290 DO 295 I=1,N
         Z = HH*AMAX1(ABS(X(I)),AX)
         ZZ = X(I)
         X(I) = ZZ+Z
         CALL FUNCT(N,X,F1)
         X(I) = ZZ-Z
         CALL FUNCT(N,X,F2)
         W(IG+I) = (F1-F2)/(Z+Z)
         X(I) = ZZ
  295 CONTINUE
      IFN = IFN+N+N
      GO TO (105, 215), LINK
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HZXMIN )
 9005 RETURN
      END

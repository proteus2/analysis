C   IMSL ROUTINE NAME   - DVCPS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DVCPR
C
C   REQD. IMSL ROUTINES - DVCPT,DVCPU,DVCPV,DVCPW,DVCPX,DVCPY,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPS (MMAX,M,NMAX,N,P,R,MTNMAX,MMAX2,A,B,PAR,TOL,X,Y,
     *                   ABT,ALPHA,A1,B1,EJ,A2,C2,DEL,UU,RES,F,HX,SK,
     *                   GRADF,AUX,XAU,IC,IR,IQJ,ICA,FCNI,FCNJ,FCNB,
     *                   JERROR)
C                                  SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER            MMAX,M,NMAX,N,P,R,MTNMAX,MMAX2,JERROR,IC(NMAX,
     *                   MMAX),IR(NMAX,MMAX),IQJ(NMAX),ICA(MMAX)
      REAL               A,B,PAR(5),TOL,X(NMAX),Y(MTNMAX),ABT(MMAX),
     *                   ALPHA(MMAX),A1(MMAX,MMAX),B1(MMAX,MMAX),
     *                   EJ(NMAX),A2(MTNMAX,MMAX),C2(MTNMAX,MMAX),
     *                   DEL(MMAX,MTNMAX),UU(MTNMAX),RES(MTNMAX),
     *                   F(MTNMAX),HX(NMAX),SK(MTNMAX),GRADF(MTNMAX),
     *                   AUX(MMAX2,MMAX),XAU(MTNMAX)
C
C                                  SPECIFICATIONS FOR COMMON /NEWT /
      COMMON /NEWT/      INWT,NU,CASI
      INTEGER            INWT,NU
      LOGICAL            CASI
C                                  SPECIFICATIONS FOR COMMON /C1 /
      COMMON /C1/        EPSNU,CONT
      REAL               EPSNU
      LOGICAL            CONT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            I1,I7,I8,ICON,IERROR,IFIN,IFI,IIM,IIP,II,IPM,
     *                   IPRINT,IP,IQ,ITEMP,I,JI,J,K11,KII1,KII,KIJ,KI,
     *                   KK,KMAX,K,L,MPNM,MPN,N1,N2,NMA,NOLD,NPU,NP,
     *                   NTEM,NTOP,NU2,NIN,NOUT
      REAL               AA(50),ALG,AUXI,BB(50),BMA,C3,C(50),DELEPS,
     *                   EPBAR,EPS1,EPSMAC,EPS,ERRNEW,ERROLD,E,HCUA,HI,
     *                   H,RABS,REPS,SIG02,TEM,TE,U1,UUN,VAIN,XN,XTE,
     *                   XXN,YKI,Z1,Z2,AM16,U2
      LOGICAL            LIN,SING
      EXTERNAL           FCNI,FCNJ,FCNB
      DATA               REPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      CASI = .FALSE.
      SING = .FALSE.
      IF (M.GT.1 .AND. M.LE.MMAX .AND. N.GT.3 .AND. N.LE.NMAX .AND.
     *P.GE.0 .AND. P.LT.M .AND. R.GE.0 .AND. P+R.LE.M) GO TO 5
      JERROR = 1
      RETURN
C                                  INITIALIZATION
    5 MPN = M*N
C
      EPSMAC = 10.*REPS
C
      EPSMAC = AMAX1(EPSMAC,1.E-4*TOL)
      IF (PAR(1).GT.0.) GO TO 10
      DELEPS = 0.
      IPRINT = 0
      VAIN = 0.
      LIN = .FALSE.
      GO TO 15
   10 DELEPS = PAR(2)
      IPRINT = PAR(3)+0.5
      VAIN = PAR(4)
      LIN = PAR(5).NE.0.
C
   15 N1 = N-1
      EPSNU = 0.
      IF (DELEPS.EQ.0.) EPSNU = 1.
      RABS = 1.
      JERROR = 0
      EPBAR = 1.E10
      CONT = .FALSE.
      ICON = 3
      NOLD = 1
      NMA = NMAX
      AM16 = NMAX-16
      NTOP = AMAX1(.9*NMAX,AM16)
      BMA = 0.07
      BB(1) = 1.
      DO 20 I=2,50
   20 BB(I) = 0.
      DO 25 NU=1,20
         NU2 = 2*NU+1
         AA(NU2) = -NU/(2.**(2*NU-1)*NU2)
         AA(NU2+1) = 0.
   25 CONTINUE
      IF (VAIN.GT.0.) GO TO 40
C                                  FIRST APPROXIMATION FOR Y AND X
      X(1) = A
      X(N) = B
      H = (B-A)/N1
      DO 30 I=2,N1
   30 X(I) = A+(I-1)*H
      DO 35 I=1,MPN
   35 Y(I) = 0.
C
   40 H = 0.
      DO 45 I=1,N1
         HX(I) = X(I+1)-X(I)
         IF (HX(I).GT.H) H = HX(I)
   45 CONTINUE
      HCUA = H**2
      IF (IPRINT.NE.0) WRITE (NOUT,50) TOL, H
   50 FORMAT (11H TOLERANCE=, E12.2, 8H MAX. H=, E12.2)
C                                  MAIN BODY OF PASVA3
      NU = 0
      EPS = AMAX1(EPSMAC,H**2*.1)
C                                  ENTER AFTER MESH CHANGE
   55 ERROLD = 1.0E20
      KMAX = MAX0(15,(N-2)/2)
      MPN = M*N
      N1 = N-1
      MPNM = M*N1
      IPM = MPNM+P+1
      C3 = 0.8
      DO 60 I=1,MPN
   60 SK(I) = 0.
      IF (NU.EQ.0) GO TO 80
C                                  AFTER MESH CHANGE WE HAVE TO
C                                    INITIALIZE SK IF NU .GT. 0
      DO 65 I7=1,N
         I8 = (I7-1)*M+1
         CALL FCNI(M,X(I7),Y(I8),F(I8))
   65 CONTINUE
      CALL FCNB(M,Y(1),Y(I8),ALPHA)
      CALL DVCPT(NU,2,2,N1,M,AA,X,F,RES,IERROR)
      DO 75 I=1,N1
         KI = (I-1)*M
         DO 70 J=1,M
            KIJ = KI+J
            SK(KIJ) = HX(I)*RES(KIJ)
   70    CONTINUE
   75 CONTINUE
      IF (NU.LT.KMAX) GO TO 80
      NU = NU-1
      GO TO 230
C                                  NU IS TOO LARGE GO TO REFINE THE
C                                    MESH. NEWTON ITERATION
   80 IF (EPSNU.GE.1.) GO TO 90
      EPS1 = EPS
   85 EPS = AMAX1(RABS,EPS)
   90 CALL DVCPV(M,N,P,R,ALPHA,A1,B1,X,Y,LIN,A2,C2,DEL,JERROR,IPRINT,
     *EPS,IR,IC,UU,RES,MMAX,MTNMAX,NMAX,MMAX2,F,HX,SK,GRADF,AUX,ICA,XAU,
     *FCNI,FCNJ,FCNB)
      IF (JERROR.NE.3) GO TO 95
      RETURN
   95 IF (EPSNU.GE.1.) GO TO 150
C                                  CONTINUATION-- CHOOSE STEP AND NEW
C                                    INITIAL PROFILE
      CONT = .TRUE.
      DO 100 I7=1,N
         I8 = (I7-1)*M+1
         CALL FCNI(M,X(I7),Y(I8),F(I8))
  100 CONTINUE
      CALL FCNB(M,Y(1),Y(I8),ALPHA)
      CONT = .FALSE.
      DO 105 I=1,P
  105 UU(I) = -ALPHA(I)
      DO 115 I=1,N1
         II = (I-1)*M
         IIP = II+P
         IIM = II+M
         DO 110 J=1,M
  110    UU(IIP+J) = .5*HX(I)*(F(IIM+J)+F(II+J))
  115 CONTINUE
      IP = P+1
      DO 120 I=IP,M
  120 UU(N1*M+I) = -ALPHA(I)
      CALL DVCPY(A2,C2,DEL,UU,M,N,P,R,IR,IC,UU,MTNMAX,MMAX,NMAX,XAU)
      IF (INWT-1) 125, 125, 130
  125 DELEPS = 2.*DELEPS
      GO TO 135
  130 DELEPS = AMAX1(.01,DELEPS/(INWT-1))
  135 DO 140 I=1,MPN
  140 Y(I) = Y(I)+DELEPS*UU(I)
      EPSNU = AMIN1(EPSNU+DELEPS,1.)
      IF (IPRINT.NE.0) WRITE (NOUT,145) DELEPS, EPSNU, EPS
  145 FORMAT (11H DEL-EPSNU=, E15.7, 7H EPSNU=, E15.7, 5H EPS=, E15.7)
      IF (EPSNU.LT.1.) GO TO 85
      EPS = EPS1
      GO TO 90
C                                  CORRECTION AND ERROR CONTROL STARTS
  150 CALL DVCPT(NU+1,2,2,N1,M,AA,X,F,RES,IERROR)
      IF (IERROR.EQ.1) GO TO 230
      DO 160 I=1,N1
         II = (I-1)*M
         DO 155 J=1,M
            KI = II+J
            AUXI = RES(KI)*HX(I)
            RES(KI) = SK(KI)-AUXI
            SK(KI) = AUXI
  155    CONTINUE
  160 CONTINUE
      IF (ICON.LE.12) GO TO 275
  165 IF (P.EQ.0) GO TO 175
      DO 170 I=1,P
  170 UU(I) = 0.
  175 DO 180 I=IPM,MPN
  180 UU(I) = 0.
      DO 185 I=1,MPNM
  185 UU(I+P) = RES(I)
      CALL DVCPW(M,N,P,R,X,Y,A1,B1,A2,C2,DEL,.TRUE.,SING,IR,IC,UU,UU,
     *LIN,MMAX,MTNMAX,NMAX,MMAX2,HX,GRADF,AUX,ICA,XAU,FCNJ,FCNB)
C
C                                  ESTIMATE FOR MAX. ABSOLUTE ERROR (BY
C                                    COMPONENTS)
      ICON = 15
      ERRNEW = 0.
      DO 190 J=1,M
  190 ABT(J) = 0.
      DO 200 I=1,N
         KK = (I-1)*M
         DO 195 J=1,M
            U1 = ABS(UU(KK+J))
            U2 = U1/AMAX1(1.0,ABS(Y(KK+J)))
            IF (U1.GT.ABT(J)) ABT(J) = U1
            IF (U2.GT.ERRNEW) ERRNEW = U2
  195    CONTINUE
  200 CONTINUE
      K = NU+1
      IF (IPRINT.EQ.0) GO TO 215
      WRITE (NOUT,205) (ABT(J),J=1,M)
      WRITE (NOUT,210) ERRNEW, NU
  205 FORMAT (30H ESTIMATED ERROR BY COMPONENTS/(10E12.3))
  210 FORMAT (17H ESTIMATED ERROR=, E12.3, 14H IN CORRECTION, I4)
  215 IF (ERRNEW.GT.TOL) GO TO 220
      JERROR = 0
      RETURN
C                                  PRECISION ACHIEVED IF NOT ENOUGH
C                                    POINTS , REFINE THE MESH
  220 IF (NU+1.GE.KMAX) GO TO 275
      IF (ERRNEW.LE.0.1*ERROLD) GO TO 225
      IF (ERRNEW.GT.C3*ERROLD) GO TO 230
      C3 = 0.5*C3
C                                  EITHER KEEP CORRECTING .
  225 ERROLD = ERRNEW
      EPS = AMAX1(EPSMAC,1.E-3*ERROLD)
      NU = NU+1
      GO TO 90
C                                  OR REFINE THE MESH, UNLESS JERROR = 4
  230 IF (JERROR.NE.4) GO TO 235
      RETURN
  235 IF (N.LE.NTOP) GO TO 240
C                                  TOO MANY GRID POINTS
      JERROR = 2
      RETURN
  240 EPS = AMAX1(EPSMAC,1.E-3*ERROLD)
      IF (ERROLD.LE.1.) GO TO 245
      EPBAR = AMIN1(1.0,1.E-2*ERROLD)
      GO TO 250
  245 EPBAR = .01*ERROLD
  250 ICON = 0
      NOLD = N
      BMA = 1.
      IF (NU.LT.1) GO TO 275
      DO 260 I=1,N1
         II = (I-1)*M
         DO 255 J=1,M
            KI = II+J
            SK(KI) = RES(KI)+SK(KI)
  255    CONTINUE
  260 CONTINUE
      NU = NU-1
      IF (NU.EQ.0) GO TO 275
      CALL DVCPT(NU,2,2,N1,M,AA,X,F,RES,IERROR)
      DO 270 I=1,N1
         II = (I-1)*M
         DO 265 J=1,M
            KI = II+J
            RES(KI) = RES(KI)*HX(I)-SK(KI)
  265    CONTINUE
  270 CONTINUE
C                                  MESH VARIATION EQUIDISTRIBUTION OF
C                                    THE L2 NORM OF THE ERROR FOR THE
C                                    O(H**(2*NU+2)) METHOD.
  275 ICON = ICON+1
      IF (IPRINT.NE.0) WRITE (NOUT,280)
  280 FORMAT (18H ---STEP CHANGE---)
      ALG = 1.5
      SIG02 = 1./(2*NU+2)
      TEM = 0.
      UUN = 0.
      DO 290 I=1,N1
         TE = 0.
         KI = (I-1)*M
         DO 285 J=1,M
            K11 = KI+J
            Z1 = ABS(Y(K11))
            Z2 = ABS(RES(K11))
            IF (Z1.GT.TEM) TEM = Z1
            IF (Z2.GT.TE) TE = Z2
  285    CONTINUE
         EJ(I) = (TE/HX(I))**SIG02
         UUN = UUN+EJ(I)
  290 CONTINUE
C
      IF (ICON.GT.1 .AND. NOLD.GT.1) GO TO 295
      EPBAR = AMAX1(AMIN1(EPBAR,TEM*BMA),TOL)
      E = EPBAR**SIG02
  295 IQ = 0
      N2 = N-2
      I = 0
  300 I = I+1
      IQJ(I) = EJ(I)/E-0.33
      IQ = IQ+IQJ(I)
      IF (I.LE.N2) GO TO 300
      IFIN = .04*N
      NMA = MIN0(NMAX-N,70)
      IF (IPRINT.NE.0) WRITE (NOUT,305) IQ, NMA
  305 FORMAT (12H NEW POINTS=, I4, 6H NMAX=, I4)
      IF (IQ.LE.IFIN .OR. NMA.LE.0) GO TO 385
      IF (IQ.LE.NMA) GO TO 315
C                                  WE ATTEMPT TO DIMINISH THE NUMBER OF
C                                    POINTS TO BE INTRODUCED
      IF (ALG.LT.1.09) GO TO 385
      ALG = ALG-.1
      XN = N+IQ
      XTE = N+NMA
      E = E*XN/AMIN1(ALG*N,XTE)
      IF (IPRINT.EQ.0) GO TO 295
      WRITE (NOUT,310) E, ALG
  310 FORMAT (7H LEVEL=, E15.5, 5H ALG=, E15.5)
      GO TO 295
C                                  CONSTRUCT NEW MESH
  315 J = 2
      DO 320 I=1,MPN
  320 UU(I) = Y(I)
      SK(1) = A
      NPU = 2*NU+3
      DO 370 I=1,N1
         KII1 = I*M
         IFI = J+IQJ(I)
         HI = HX(I)/(IQJ(I)+1)
         DO 365 L=J,IFI
            SK(L) = X(I)+HI*(L-J+1)
            KII = (L-1)*M
            IF (IQJ(I).GT.0) GO TO 330
            DO 325 K=1,M
  325       Y(KII+K) = UU(KII1+K)
            GO TO 365
  330       IF (I.LE.NU+2) GO TO 340
            IF (I.GT.N-(NU+3)) GO TO 345
            NP = NU+2
  335       ITEMP = I
            CALL DVCPU(ITEMP,NPU,NP,C,BB,X,SK(L))
            GO TO 350
  340       NP = I
            GO TO 335
  345       NP = I-N+NPU
            GO TO 335
  350       DO 360 K=1,M
               YKI = 0.
               DO 355 JI=1,NPU
  355          YKI = YKI+C(JI)*UU((I-NP+JI-1)*M+K)
               KI = KII+K
               Y(KI) = YKI
  360       CONTINUE
  365    CONTINUE
         J = IFI+1
  370 CONTINUE
C
      CASI = .FALSE.
      N = N+IQ
      N1 = N-1
      KI = M*N1
      DO 375 I=1,M
  375 Y(KI+I) = UU(MPNM+I)
      H = 0.
      DO 380 I=2,N1
         I1 = I-1
         HX(I1) = SK(I)-SK(I1)
         IF (HX(I1).GT.H) H = HX(I1)
         KI = I1*M+1
         X(I) = SK(I)
  380 CONTINUE
C
      X(N) = B
      HX(N1) = X(N)-X(N1)
      IF (HX(N1).GT.H) H = HX(N1)
      HCUA = H**2
      IF (ABS(EPBAR-TEM*BMA).LT.1.E-5) EPBAR = 1.E10
      IF (ICON.GE.5) ICON = 15
      GO TO 55
C                                  END OF MESH SELECTION
C
  385 IF (NOLD.EQ.1 .AND. ALG.LT.1.5 .AND. 2*N-1.LE.NMAX) GO TO 405
      IF (N.GT.NOLD) GO TO 165
      IF (ALG.LT.1.45) GO TO 390
      ALG = 1.4
      XN = N+IQ
      XXN = NMAX
      E = E*XN/AMIN1(ALG*N,XXN)
      IF (IPRINT.EQ.0) GO TO 295
      WRITE (NOUT,310) E, ALG
      GO TO 295
  390 IF (NU.GT.0) GO TO 415
C                                  BISECTION
      IF (IPRINT.EQ.0) GO TO 400
      WRITE (NOUT,395)
  395 FORMAT (23H BISECTION IS PERFORMED)
  400 NTEM = 2*N-1
      IF (NTEM.LE.NMAX) GO TO 405
C                                  TOO MANY GRID POINTS
      JERROR = 2
      RETURN
  405 DO 410 I=1,N1
  410 IQJ(I) = 1
      IQ = N1
      ICON = 0
      GO TO 315
  415 NU = NU-1
      ICON = 0
C                                  SINCE WE ARE GOING BACK WE MUST
C                                    RECOMPUTE THE ESTIMATE OF THE
C                                    ERROR.
      MPNM = M*N1
      DO 420 I=1,MPNM
  420 SK(I) = RES(I)+SK(I)
      IF (NU.GT.0) GO TO 435
      DO 430 I=1,N1
         KI = (I-1)*M
         DO 425 J=1,M
            KIJ = KI+J
            RES(KIJ) = -SK(KIJ)
  425    CONTINUE
  430 CONTINUE
      GO TO 275
  435 CALL DVCPT(NU,2,2,N1,M,AA,X,F,RES,IERROR)
      DO 445 I=1,N1
         HI = HX(I)
         KI = (I-1)*M
         DO 440 J=1,M
            KIJ = KI+J
            RES(KIJ) = HI*RES(KIJ)-SK(KIJ)
  440    CONTINUE
  445 CONTINUE
      GO TO 275
C                                  END OF MESH VARIATION
      END

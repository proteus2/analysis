C   IMSL ROUTINE NAME   - DRVTE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - CALCULATE FIRST, SECOND OR THIRD DERIVATIVE
C                           OF A USER SUPPLIED FUNCTION
C
C   USAGE               - FUNCTION DRVTE (F,N,X,H,ERR,IER)
C
C   ARGUMENTS    DRVTE  - ESTIMATE OF THE FIRST (IF N=1), SECOND
C                           (N=2) OR THIRD (N=3) DERIVATIVE OF F(X).
C                           (OUTPUT)
C                F      - A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM
C                           SUPPLIED BY THE USER (INPUT).
C                           F MUST BE DECLARED EXTERNAL IN THE
C                           CALLING PROGRAM.
C                N      - ORDER OF THE DERIVATIVE DESIRED (1,2 OR 3).
C                           (INPUT)
C                X      - POINT AT WHICH THE DERIVATIVE IS TO BE
C                           EVALUATED (INPUT).
C                H      - INITIAL STEP SIZE USED IN FINITE
C                           DIFFERENCE SCHEME (INPUT). H WILL BE
C                           DECREASED UNTIL THE FINITE DIFFERENCE
C                           ESTIMATES STABILIZE.  IF H=0.0,
C                           A VALUE OF 1.0 WILL BE USED.
C                ERR    - RELATIVE ERROR DESIRED IN DERIVATIVE ESTIMATE
C                           (INPUT).  CONVERGENCE IS ASSUMED WHEN
C                         ABS(D1-D2)/AMAX1(1.0,SQRT(D1**2+D2**2)).LT.ERR
C                           WHERE D1 AND D2 ARE SUCCESSIVE ESTIMATES
C                           OF THE DERIVATIVE.
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 N IS NOT EQUAL TO 1,2 OR 3
C                           IER = 130 UNABLE TO ACHIEVE REQUIRED
C                               ACCURACY IN DERIVATIVE ESTIMATION.
C                               INCREASE ERR OR SUPPLY A MORE
C                               REASONABLE VALUE OF H.  IF IER=130
C                               STILL, F MAY NOT HAVE A DERIVATIVE
C                               AT X.
C                           IER = 131 ROUNDOFF ERROR BECAME
C                               DOMINANT BEFORE ESTIMATES CONVERGED.
C                               INCREASE H.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IT IS IMPORTANT THAT THE HIGHEST AVAILABLE
C                PRECISION BE USED IN DERIVATIVE CALCULATIONS.
C                THUS ONLY A DOUBLE PRECISION VERSION OF DRVTE
C                IS AVAILABLE FOR CERTAIN COMPUTERS.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DRVTE (F,N,X,H,ERR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      DOUBLE PRECISION   F,X,H,ERR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IBIAS,I,J,KCOLSS,KCOL,K,LIM,L,MOD1,MOD2,MOD,
     *                   MSSGE,M,NCOLC,NCOMT,NDERV1,NDERV,NGOBCN,NGOBC,
     *                   NGOBRN,NGOBR,NOCN11,NOCNV1,NOCNVC,NOCNVR,
     *                   NOCONC,NOCONR,NRHVER,NROWC,NROWR,NUM1EX,NUM1,
     *                   NUM2EX,NUM2,NUMDRT,NUMDRV,NUMT,NUM
      DOUBLE PRECISION   CD(3,33),C(4),D12(33,3),DENOM,DIF1,DIF2,DIFMAX,
     *                   DIF,FIXT,HHH,HII,HI,HSTOP,PRE(3),QI,Q,QL,RELID,
     *                   RHDERV,RHH,SAV(3),SDIF1,SDIF2,SDIF,STBDX,S,XXX
      DATA               D12 /99*0.0D0/
C                                  INSULATION AND CHECK OF PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 131
      HHH = X + H
      IF (HHH.EQ.X) GO TO 9000
      IER = 0
      NDERV = N
      XXX = X
      IF (NDERV.EQ.1 .OR. NDERV.EQ.2 .OR. NDERV.EQ.3) GO TO 5
      IER = 129
      GO TO 9000
    5 HHH = DABS(H)
      NDERV1 = NDERV+1
C                                  INITIALIZATION
      QI = .5D0
      IF (HHH.EQ.0.D0) HHH = 1.0D0
      HII = HHH
      Q = QI
      L = 2
      FIXT = F(XXX)
C                                  THE FOLLOWING 10 FLAGS DESCRIBE THE
C                                  BEHAVIOR OF THE CD AND RH SCHEME
      NGOBC = 0
      NGOBR = 0
      NGOBCN = 0
      NGOBRN = 0
      NOCNVC = 0
      NOCNVR = 0
      NOCONC = 0
      NOCONR = 0
      NRHVER = 0
      MSSGE = 0
   10 MOD = 0
   15 DO 20 J=1,NDERV1
   20 C(J) = FIXT
      IF (MOD.EQ.2) GO TO 30
      IF (L.EQ.1) Q = QI/2.D0
      NOCN11 = 0
      HI = HHH
      NUMT = 0
      NUM = 0
      NUM1 = 0
      NUM2 = 0
      IBIAS = 9
      LIM = 33
      KCOL = 0
      KCOLSS = 0
      MOD = 0
      MOD1 = 0
      MOD2 = 0
      NOCNV1 = 0
C                                  START SCHEMES
      I = IBIAS
   25 IF (LIM.EQ.12 .AND. I.GT.8) GO TO 50
   30 K = 1
      M = NDERV1
      IF ((NDERV+L.NE.5) .OR. (I.EQ.1 .OR. I.EQ.9)) GO TO 35
      C(1) = C(2)
      C(4) = C(3)
      K = 2
      M = 3
C                                  STORE REQUIRED FUNCTION VALUES IN
C                                  C(J)
   35 J = K
   40 S = NDERV1-J
      IF (L.EQ.1 .OR. NDERV+J.EQ.2) GO TO 45
      S = S-1.D0
      IF (NDERV+J.GT.5) S = S-1.D0
   45 IF (S.NE.0.D0) C(J) = F(XXX+S*HHH)
      J = J+1
      IF (J.LE.M) GO TO 40
C                                  FORMULAE
      QL = L
      IF (NDERV.EQ.1) D12(I,1) = (C(1)-C(2))/QL
      IF (NDERV.EQ.2) D12(I,1) = C(1)-2.D0*C(2)+C(3)
      IF (NDERV.EQ.3 .AND. L.EQ.2) D12(I,1) = (C(1)-C(4))/2.D0-C(2)+C(3)
      IF (NDERV.EQ.3 .AND. L.EQ.1) D12(I,1) = C(1)-3.D0*(C(2)-C(3))-C(4)
      D12(I,1) = D12(I,1)/HHH**NDERV
C                                  IF SENT HERE FOR VERIFICATION, THEN
C                                  RETURN
      IF (MOD.EQ.2) GO TO 145
      NUM = NUM+1
C                                  MONITORING OF CENTRAL DIFFERENCES OR
C                                  RIGHT HAND SCHEME
      IF (I.EQ.IBIAS) GO TO 105
   50 SDIF = DABS(PRE(1)-D12(I,1))
      IF (I.EQ.IBIAS+1) GO TO 60
      IF (SDIF.EQ.0.D0 .AND. DIF.EQ.0.D0) KCOL = 1
      IF (SDIF.GE.DIF .AND. MOD.EQ.1) KCOL = 1
      IF (KCOL.NE.0) GO TO 115
C                                  A ZERO DIFFERENCE FOLLOWED BY A
C                                  NON-ZERO ONE ASSUMED TO MEAN H IS
C                                  TOO SMALL.
      IF (SDIF.NE.0.D0 .AND. DIF.EQ.0.D0) GO TO 110
      IF (D12(I,1).NE.0.D0) GO TO 55
      IF (SDIF.EQ.0.D0 .AND. MOD.EQ.-1) GO TO 110
   55 MOD = -1
      IF (DIF.GT.SDIF) MOD = 1
   60 DIF = SDIF
      IF (LIM.EQ.12 .AND. I.GT.9) GO TO 65
C                                  CALCULATION AND MONITORING OF 1ST
C                                  EXTRAPOLATIONS
      D12(I,2) = (D12(I,1)-Q**L*PRE(1))/(1.D0-Q**L)
      NUM1 = NUM1+1
      IF (I.EQ.IBIAS+1) GO TO 100
   65 SDIF1 = DABS(PRE(2)-D12(I,2))
      IF (I.EQ.IBIAS+2) GO TO 75
      IF (SDIF1.EQ.0.D0 .AND. DIF1.EQ.0.D0) KCOL = 2
      IF (SDIF1.GE.DIF1 .AND. MOD1.EQ.1) KCOL = 2
      IF (KCOL.NE.0) GO TO 115
   70 IF (SDIF1.LT.DIF1) MOD1 = 1
   75 DIF1 = SDIF1
      IF (LIM.EQ.12 .AND. I.GT.10) GO TO 80
C                                  CALCULATION AND MONITORING OF 2ND
C                                  EXTRAPOLATIONS
      D12(I,3) = (D12(I,2)-Q**(2*L)*PRE(2))/(1.D0-Q**(2*L))
      NUM2 = NUM2+1
      IF (I.EQ.IBIAS+2) GO TO 95
   80 SDIF2 = DABS(PRE(3)-D12(I,3))
      IF (I.EQ.IBIAS+3) GO TO 90
      IF (SDIF2.EQ.0.D0 .AND. DIF2.EQ.0.D0) KCOL = 3
      IF (SDIF2.GE.DIF2 .AND. MOD2.EQ.1) KCOL = 3
      IF (KCOL.NE.0) GO TO 115
   85 IF (SDIF2.LT.DIF2) MOD2 = 1
   90 DIF2 = SDIF2
      SAV(3) = PRE(3)
   95 PRE(3) = D12(I,3)
      SAV(2) = PRE(2)
  100 PRE(2) = D12(I,2)
      SAV(1) = PRE(1)
  105 PRE(1) = D12(I,1)
      HHH = Q*HHH
      I = I+1
      IF (I.LE.LIM) GO TO 25
      NOCNV1 = 1
  110 IF (LIM.EQ.12) GO TO 160
C                                  ROUNDOFF ERROR DOMINANT
      IER = 131
      GO TO 9000
C                                  AT THIS POINT CONVERGENCE HAS
C                                  OCCURRED FOR WHICHEVER SCHEME WAS RUN
  115 IF (L.EQ.2) GO TO 120
      RHDERV = (SAV(KCOL)+PRE(KCOL))/2.D0
      NROWR = I-11
      NOCONR = NOCN11
      GO TO 150
C                                  CHECK FOR RELATIVE AGREEMENT TO
C                                  PREVENT ACCIDENTAL CONVERGENCE
  120 STBDX = DABS(SAV(KCOL)-PRE(KCOL))
      STBDX = STBDX/DMAX1(1.0D0,DSQRT(SAV(KCOL)**2+PRE(KCOL)**2))
      IF (STBDX.LE.0.001D0) GO TO 125
      KCOLSS = KCOL
      MOD = 0
      MOD1 = 0
      MOD2 = 0
      KCOL = 0
      GO TO (55, 70, 85), KCOLSS
C                                  STORE CENTRAL DIFFERENCE INFORMATION
  125 DRVTE = (SAV(KCOL)+PRE(KCOL))/2.D0
  130 NUMDRV = NUM
      NUMDRT = NUMT
      NUM1EX = NUM1
      NUM2EX = NUM2
      NCOLC = KCOL
      NROWC = I-11
      NOCONC = NOCN11
      DO 140 J=1,33
         DO 135 K=1,3
  135    CD(K,J) = D12(J,K)
  140 CONTINUE
      HSTOP = HHH
C                                  CALCULATE A SUITABLE H TO USE IN A
C                                  RIGHT HAND SCHEME.
      HHH = DMIN1(DSQRT(HHH),HHH**2)
      HHH = DMAX1(HHH,.0001D0*HII)
      IF (NOCNVC.EQ.1) GO TO 155
      RHH = HHH
      I = 1
      L = 1
      MOD = 2
      GO TO 15
  145 RHDERV = D12(1,1)
      L = 2
C                                  VERIFICATION PROCESS
  150 DENOM = DSQRT(DRVTE**2+RHDERV**2)
      RELID = DABS(DRVTE-RHDERV)/DMAX1(1.0D0,DENOM)
      IF (MOD.EQ.3) GO TO 175
      IF (RELID.LE.0.05D0) GO TO 215
      IF (MOD.LT.2) GO TO 220
  155 HHH = HHH*(2.D0/Q)**3
      RHH = HHH
      L = 1
      GO TO 10
  160 IF (L.EQ.2) GO TO 185
C                                  LOOK FOR AGREEMENT OF ANY 2 ELEMENTS
C                                  OF RIGHT HAND SCHEME WITH DERIVATIVE
C                                  SELECTED BY CENTRAL DIFFERENCING
      NOCNVR = 1
      MSSGE = MSSGE+2
      IF (NOCNV1.EQ.1) GO TO 165
      NGOBR = 1
      NGOBRN = I-9
  165 MOD = 3
      M = 0
      J = NUMT+8
      K = 9
      IF (NUM-NUMT.GT.7) K = 1
  170 I = K
      RHDERV = D12(I,1)
      GO TO 150
  175 IF (RELID.LE.0.05D0) M = M+1
      IF (M.EQ.2) GO TO 180
      I = I+1
      IF (I.LE.J) GO TO 170
      IF (K.EQ.1) GO TO 220
      J = NUM-NUMT
      K = 1
      GO TO 170
  180 RHH = RHH*Q**(I-12)
      GO TO 215
C
C                                  FIND CLOSEST THREE IN A ROW IN THE CD
C                                  SCHEME AND CHOOSE AVE AS DERV
C
  185 NOCNVC = 1
      MSSGE = 4
      IF (NOCNV1.EQ.1) GO TO 190
      NGOBC = 1
      NGOBCN = I-9
  190 J = NUMT+8
      K = 11
      IF (NUM-NUMT.EQ.8) K = 3
      STBDX = 2
  195 SDIF = DABS(D12(K-2,1)-D12(K-1,1))
      DO 205 M=K,J
         SDIF1 = DABS(D12(M-1,1)-D12(M,1))
         SDIF2 = DABS(D12(M,1)-D12(M-2,1))
         DIFMAX = DMAX1(SDIF,SDIF1,SDIF2)
C                                  DIFMAX CAN NOT BE ZERO OTHERWISE
C                                  SCHEME WOULD HAVE CONVERGED
         DIF = D12(M-2,1)
         DIF1 = D12(M-1,1)
         IF (DIFMAX.EQ.SDIF1) DIF = D12(M,1)
         IF (DIFMAX.EQ.SDIF2) DIF1 = D12(M,1)
         DIFMAX = DIFMAX/DSQRT(DIF*DIF+DIF1*DIF1)
         IF (DIFMAX.GE.STBDX) GO TO 200
         STBDX = DIFMAX
         I = M
  200    SDIF = SDIF1
  205 CONTINUE
      IF (K.EQ.3) GO TO 210
      J = NUM-NUMT
      K = 3
      GO TO 195
  210 DRVTE = (D12(I,1)+D12(I-1,1)+D12(I-2,1))/3.D0
      KCOL = 1
      HHH = HI*Q**(I-9)
      GO TO 130
  215 NRHVER = 3
      MSSGE = MSSGE-1
  220 STBDX = 2.D0*STBDX
      NCOMT = 1
      IF (MSSGE.NE.-1) GO TO 225
      IF (STBDX.GE.0.8D0*ERR) NCOMT = 2
      IF (STBDX.GT.2.D0*ERR) NCOMT = 0
  225 IF (MSSGE.GE.0 .OR. NCOMT.EQ.0) IER = 130
      IF (MSSGE.EQ.0 .OR. MSSGE.EQ.3 .OR. MSSGE.EQ.4) DRVTE = RHDERV
 9000 IF (IER.GT.0) CALL UERTST(IER,6HDRVTE )
      RETURN
      END

C   IMSL ROUTINE NAME   - ZPOLR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZEROS OF A POLYNOMIAL WITH REAL
C                           COEFFICIENTS (LAGUERRE)
C
C   USAGE               - CALL ZPOLR (A,NDEG,Z,IER)
C
C   ARGUMENTS    A      - REAL VECTOR OF LENGTH NDEG+1 CONTAINING THE
C                           COEFFICIENTS IN ORDER OF DECREASING
C                           POWERS OF THE VARIABLE. (INPUT)
C                NDEG   - INTEGER DEGREE OF THE POLYNOMIAL. (INPUT)
C                           NDEG MUST BE GREATER THAN 0 AND LESS
C                           THAN 101.
C                Z      - COMPLEX VECTOR OF LENGTH NDEG CONTAINING
C                           THE COMPUTED ROOTS OF THE POLYNOMIAL.
C                           (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, INDICATES THAT THE DEGREE OF THE
C                             POLYNOMIAL IS GREATER THAN 100 OR LESS
C                             THAN 1.
C                           IER = 130, INDICATES THAT THE LEADING
C                             COEFFICIENT IS ZERO. THIS RESULTS IN AT
C                             LEAST ONE ROOT, Z(NDEG), BEING SET TO
C                             POSITIVE MACHINE INFINITY.
C                           IER = 131, INDICATES THAT ZPOLR FOUND
C                             FEWER THAN NDEG ZEROS. IF ONLY M ZEROS
C                             ARE FOUND Z(J),J=M+1,...,NDEG ARE SET TO
C                             POSITIVE MACHINE INFINITY.
C                             ZPOLR WILL TERMINATE WITH THIS ERROR
C                             IF IT CANNOT FIND ANY ONE ROOT WITHIN
C                             200*NDEG ITERATIONS OR IF IT DETERMINES,
C                             WITHIN THOSE 200*NDEG ITERATIONS, THAT
C                             IT CANNOT FIND THE ROOT IN QUESTION.
C
C   REQD. IMSL ROUTINES - SINGLE/UERSET,UERTST,UGETIO,ZQADC,ZQADR
C                       - DOUBLE/UERSET,UERTST,UGETIO,VXADD,VXMUL,VXSTO,
C                           ZQADC,ZQADR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZPOLR  (A,NDEG,Z,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NDEG,IER
      REAL               A(1)
      COMPLEX            Z(NDEG)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LEVEL,LEVOLD,I,NP1,N,J,K,NLM,N1,NTOM,LSW,ISC,
     1                   ITER,IIER,IKER,IJER
      REAL               ACF1(2),ACF2(2),ACF(2),ACDIR(2),AC(2),ACL(2),
     1                   RADIX,SINF,SDEPS,TWOD3,RNLGRX,F0,GAMA,THETA,
     2                   PHI,ZERO,ONE,TWO,TENM3,THREE,RNINE,HALF,OPTFM,
     3                   SINFSQ,SINFT,SETASQ,SC,XN,XN1,XN2,X2N,X2N1,
     4                   XN2N,RTN,G,R,ABDIRO,T,S,GO,RO,FEJER,ABDIR,
     5                   FN,V,ABSCL,S1,T1,X2,E,X,ABX,W,F
      REAL               FINITY
      DOUBLE PRECISION   DA(101),DZ(100),DZNR,DZNI,DZ0R,DZ0I,DXT,
     1                   DX,DR,DSC,DY,DX2,DV,DT,DT1,DZERO,DTWO
      COMPLEX            CF1,CF2,CF,CDIRO,CSPIR,CDIR,C,CL
      LOGICAL            STARTD,SPIRAL
      EQUIVALENCE        (CF1,ACF1(1)),(CF2,ACF2(1)),(CF,ACF(1)),
     1                   (CDIR,ACDIR(1)),(C,AC(1)),(CL,ACL(1))
      DATA               RADIX/16.0/
      DATA               SINF/Z7FFFFFFF/
      DATA               SDEPS/Z34100000/
      DATA               TWOD3/.6666667/
      DATA               RNLGRX/2.772589/
      DATA               FINITY/Z7FFFFFFF/
      DATA               F0/0.0/,GAMA/0.5/,THETA/1.0/,PHI/0.2/,
     1                   ZERO/0.0/,ONE/1.0/,TWO/2.0/,TENM3/1.0E-3/,
     2                   DZERO/0.0D0/,THREE/3.0/,RNINE/9.0/,DTWO/2.0D0/,
     3                   HALF/0.5/,OPTFM/-1.25/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CALL UERSET TO SILENCE ZQADR AND
C                                    ZQADC
      LEVEL = 0
      CALL UERSET (LEVEL,LEVOLD)
C                                  INITIALIZE CONSTANTS
      SINFSQ = SQRT(SINF)
      SINFT = SINFSQ
      SETASQ = ONE/SINFSQ
      IER = 0
      N = NDEG
      NP1 = N+1
      IF (N.GT.0) GO TO 5
      IER = 129
      GO TO 9000
    5 IF (N.LE.100) GO TO 10
      IER = 129
      GO TO 9000
   10 CONTINUE
C                                  MOVE THE COEFFICIENTS A(I) TO
C                                    DA(I).
      DO 15 I=1,NP1
         DA(I) = A(I)
   15 CONTINUE
C                                  SCALING (WHEN N .GT. 2)
      ASSIGN 85 TO LSW
      SC = ZERO
      DO 20 I=1,NP1
         SC = AMAX1(SC,ABS(SNGL(DA(I))))
   20 CONTINUE
      IF (SC.EQ.ZERO) GO TO 30
      ISC = ALOG(SC)/RNLGRX
      SC = RADIX**ISC
      DSC = SC
      IF (SC.GE.SINFT) GO TO 30
      SC = SINFT/SC
      DSC = SC
C                                  SCALE BY SC TO HAVE MAX(DA(I),I=1,
C                                    NP1) APPROACH SINFT
      DO 25 I=1,NP1
         DA(I) = DSC*DA(I)
   25 CONTINUE
C                                  FIND NUMBER I OF CONSECUTIVE LEADING
C                                    COEFFICIENTS EQUAL TO ZERO.
   30 DO 35 I=1,N
         IF (SNGL(DA(I)).NE.ZERO) GO TO 40
C                                  EACH VANISHED LEADING COEFFICIENT
C                                    YIELDS AN INFINITE ZERO.
         J = NP1-I
         Z(J) = CMPLX(FINITY,ZERO)
   35 CONTINUE
      IER = 130
      GO TO 9000
   40 IF (I.EQ.1) GO TO 65
C                                  SLIDE COEFFICIENTS BACK
      IER = 130
      DO 45 K=I,NP1
         J = K-I
         DA(J+1) = DA(K)
   45 CONTINUE
      N = N-I+1
      GO TO 60
C                                  RE-ENTRIES FOR CURRENT (REDUCED)
C                                    POLYNOMIAL.
   50 N = N1
   55 N = N-1
   60 NP1 = N+1
      ASSIGN 85 TO LSW
   65 CONTINUE
      ITER = 0
      IF (N-2) 70,75,80
   70 Z(1) = CMPLX(SNGL(-DA(2)/DA(1)),ZERO)
      GO TO 9000
   75 CALL ZQADR (SNGL(DA(1)),SNGL(DA(2)),SNGL(DA(3)),Z(2),Z(1),IIER)
      IF (IIER.NE.0) IER = 130
      GO TO 9000
C                                  CHECK FOR ZEROS = (0.,0.)
   80 IF (SNGL(DA(NP1)).NE.ZERO) GO TO LSW, (85,100)
      Z(N) = CMPLX(ZERO,ZERO)
      GO TO 55
C                                  HENCEFORTH N .GT. 2, DA(1) .NE. 0.0
C                                    AND DA(NP1) .NE. 0.0. INITIALIZE
C                                    SOME USEFUL CONSTANTS.
   85 CONTINUE
      XN = N
      XN1 = XN-ONE
      XN2 = XN1-ONE
      X2N = TWO/XN
      X2N1 = X2N/XN1
      XN2N = XN2/XN
      N1 = N-1
      RTN = SQRT(XN)
C                                  CALCULATE G , AN UPPER BOUND FOR THE
C                                    NEAREST ZERO. INITIALLY G =
C                                    CABS(GEOMETRIC MEAN OF THE ZEROS).
C
      G = EXP((ALOG(ABS(SNGL(DA(NP1))))-ALOG(ABS(SNGL(DA(1)))))/
     1        XN+TENM3)
C                                  CALCULATE LAGUERRE-STEP CDIR AND
C                                    FEJER-BOUND FOR G.
      R = SNGL(DA(N))/SNGL(DA(NP1))
      CALL ZQADR (X2N1*SNGL(DA(N-1)),X2N*SNGL(DA(N)),SNGL(DA(NP1)),C,
     1            CF1,IKER)
      R = XN2N*R
      IF (IKER .EQ. 65) CDIRO = CMPLX(AC(1)/XN1,ZERO)
      IF (IKER .NE. 65) CDIRO = C/CMPLX(R*AC(1)+XN1,R*AC(2))
      ABDIRO = ABS(REAL(CDIRO))+ABS(AIMAG(CDIRO))
      G = AMIN1(G,1.0001*AMIN1(ABS(AC(1))+ABS(AC(2)),RTN*ABDIRO))
C                                  CALCULATE THE CAUCHY-LOWER BOUND R
C                                    FOR THE SMALLEST ZERO BY SOLVING
C                                    ABS(DA(NP1)) = SUM(ABS(DA(I))
C                                    *R**(NP1-I),I=1,N)
C                                    USING NEWTON*S METHOD.
      R = G
   90 T = ABS(SNGL(DA(1)))
      S = ZERO
      DO 95 I=2,N
         S = R*S+T
         T = R*T+ABS(SNGL(DA(I)))
   95 CONTINUE
      S = R*S+T
      T = (R*T-ABS(SNGL(DA(NP1))))/S
      S = R
      R = SNGL(DBLE(R)-DBLE(T))
      IF (R.LT.S) GO TO 90
C                                  R/(2**(1/N) - 1 ) .LT. 1.445*N*R IS
C                                    ANOTHER UPPER BOUND.
      GO = AMIN1(1.445*XN*R,G)
      RO = 0.99999*S
      ASSIGN 100 TO LSW
C                                  INITIALIZE THE ITERATION TO BEGIN AT
C                                    THE ORIGIN.
  100 CONTINUE
      FEJER = GO
      G = GO
      CDIR = CDIRO
      DZNR = DZERO
      ABDIR = ABDIRO
      DZNI = DZERO
      FN = ABS(SNGL(DA(NP1)))
      SPIRAL = .FALSE.
      STARTD = .FALSE.
C                                  RE-ENTRY POINT TO ACCEPT, MODIFY,
C                                    OR REJECT THE LAGUERRE STEP.
  105 CONTINUE
C                                  ACCEPT CDIR IF CABS(CDIR) .LE.
C                                    GAMA*G.
      IF (ABDIR .LE. G*GAMA) GO TO 110
C                                  REJECT CDIR IF CABS(CDIR) .GT.
C                                    THETA*G.
      IF (ABDIR .GT. G*THETA) GO TO 215
C                                  MODIFY CDIR SO THAT CABS(CDIR) =
C                                    GAMA*G.
      IF (.NOT.(STARTD.OR.SPIRAL).AND.RO.GT.GAMA*G) GO TO 110
      V = GAMA*(G/ABDIR)
      CDIR = CMPLX(V*ACDIR(1),V*ACDIR(2))
      ABDIR = ABDIR*V
C                                  ACCEPT PREVIOUS ITERATE. SAVE DATA
C                                    ASSOCIATED WITH CURRENT ITERATE.
  110 CONTINUE
      G = FEJER
      CL = CDIR
      ABSCL = ABDIR
      F0 = FN
      DZ0R = DZNR
      DZ0I = DZNI
C                                  CDIR AT THE ORIGIN IS IN THE
C                                    DIRECTION OF DECREASING FUNCTION
C                                    VALUE
      STARTD = .TRUE.
C                                  NEXT ITERATE IS ZN=CMPLX(DZNR,DZNI).
  115 DZNR = DZ0R+ACL(1)
      DZNI = DZ0I+ACL(2)
C                                  IS ZN CLOSE TO THE REAL AXIS
C                                    RELATIVE TO STEP SIZE.
  120 CONTINUE
      ITER = ITER+1
      IF (ITER.GT.200*NDEG) GO TO 220
      IF (ABS(SNGL(DZNI)).LE.PHI*ABSCL) GO TO 175
C                                  ZN IS COMPLEX .
C                                  FACTORIZATION OF THE POLYNOMIAL BY
C                                    QUADRATIC FACTOR (Z**2-X2*Z+R)
C                                    SUM(DA(I)*Z**(N-I)) =
C                                    (Z**2-X2*Z+R)*SUM(Z(I)*Z**(N-I-2))
C                                    + Z(N-1)*(Z-X) + Z(N) FOR ALL Z ,
C                                    THE VALUE OF THE POLYNOMIAL AT
C                                    (X,Y) IS CF , FIRST DERIVATIVE OF
C                                    THE POLYNOMIAL AT (X,Y) IS CF1,
C                                    AND THE SECOND DERIVATIVE OF THE
C                                    POLYNOMIAL AT (X,Y) IS 2.*CF2,
C                                    WHERE (X,Y) IS A ZERO OF
C                                    Z**2-X2*Z+R.
C                                    E IS ERROR BOUND FOR THE VALUE OF
C                                    THE POLYNOMIAL AND DZ(I) ARE THE
C                                    COEFFICIENTS OF THE QUOTIENT
C                                    POLYNOMIAL.
C                                  INITIALIZATIONS FOR EVALUATION LOOPS
      S = ZERO
      S1 = ZERO
      DT1 = DZERO
      T1 = ZERO
      DT = DA(1)
C                                  INDEX J IS USED TO CHANGE DX ON THE
C                                    LAST ITERATION.
      J = 3
C                                  SET Z(X,Y) TO ZN(ZNR,ZNI).
      DX = DZNR
      DY = DZNI
      SC = CABS(CMPLX(SNGL(DX),SNGL(DY)))
      DSC = SC
      IF (SC.GE.SINFSQ.OR.SC.LE.SETASQ) GO TO 140
      DX2 = DX+DX
      X2 = DX2
      DR = DX*DX+DY*DY
      R = DR
      DZ(1) = DA(2)+DX2*DA(1)
      DZ(2) = DA(3)+(DX2*DZ(1)-DR*DA(1))
      IF (J.LT.N) GO TO 130
  125 DX2 = DX
      X2 = DX2
      J = N
  130 NLM = MAX0(N1,J)
      DO 135 I=J,NLM
         V = S1*R
         S1 = S
         S = T1+(X2*S-V)
         DV = DT1*DR
         DT1 = DT
         T1 = DT1
         DT = (DX2*DT-DV)+DZ(I-2)
         DZ(I) = DA(I+1)+(DX2*DZ(I-1)-DR*DZ(I-2))
  135 CONTINUE
      IF (J.LT.N) GO TO 125
      GO TO 160
  140 DX = DX/DSC
      DY = DY/DSC
      DR = (DX*DX+DY*DY)*DSC
      R = DR
      DX2 = DX+DX
      X2 = DX2
      DZ(1) = DA(2)+(DX2*DA(1))*DSC
      DZ(2) = DA(3)+(DX2*DZ(1)-DR*DA(1))*DSC
      IF (J.LT.N) GO TO 150
  145 DX2 = DX
      X2 = DX2
      J = N
  150 NLM = MAX0(N1,J)
      DO 155 I=J,NLM
         V = S1*R
         S1 = S
         S = T1+(X2*S-V)*DSC
         DV = DT1*DR
         DT1 = DT
         T1 = DT1
         DV = DX2*DT-DV
         DT = DZ(I-2)+DV*DSC
         DZ(I) = DA(I+1)+(DX2*DZ(I-1)-DR*DZ(I-2))*DSC
  155 CONTINUE
      IF (J.LT.N) GO TO 145
  160 CF = CMPLX(SNGL(DZ(NLM)),SNGL(DZNI*DZ(NLM-1)))
      FN = CABS(CF)
      E = ABS(SNGL(DA(1)))
      DO 165 I=2,N1
         E = ABS(SNGL(DZ(I-1)))+SC*E
  165 CONTINUE
      E = SDEPS*((RNINE*E*SC+THREE*ABS(SNGL(DZ(N-1))))*SC+
     1           ABS(SNGL(DZ(N))))
C                                  HAS AN ACCEPTABLE ZERO BEEN FOUND
      IF (FN.LE.E) GO TO 195
      IF (FN.GE.F0.AND.STARTD) GO TO 215
      DV = DTWO*DZNI
      V = DV
      T = DT
      CF1 = CMPLX(SNGL(DZ(NLM-1)-(DV*(DT1*DZNI))),SNGL(DV*DT))
      CF2 = CMPLX(T-V*(V*S),SNGL(DZNI)*(3.*T1-V*(V*S1)))
C                                  FIND THE LAGUERRE STEP AT ZN.
C                                  X = AMAX1(ABS(ACF(1)),ABS(ACF(2)))
C                                  C = CMPLX(ACF1(1)/X,ACF1(2)/X) /
C                                    CMPLX(ACF(1)/X,ACF(2)/X)
      IF (CABS(CF).GE.ONE) GO TO 170
C                                  IF CABS(CF1/CF) .GT. SINF, THERE IS
C                                    A ZERO WITHIN N*(1.0/SINF) OF ZN.
      IF (CABS(CF1).GT.CABS(CF)*SINF) GO TO 195
  170 C = CF1/CF
C                                  COMPUTE THE LAGUERRE STEP CDIR AND
C                                    THE BOUND FEJER AT ZN.
      CALL ZQADC (CMPLX(X2N1*ACF2(1),X2N1*ACF2(2)),CMPLX(X2N*ACF1(1),
     1            X2N*ACF1(2)),CF,CDIR,CF1,IJER)
      FEJER = ABS(ACDIR(1))+ABS(ACDIR(2))
      C = CMPLX(XN2N*AC(1),XN2N*AC(2))
      C = C*CDIR
      C = CMPLX(AC(1)+XN1,AC(2))
      CDIR = CDIR/C
      ABDIR = ABS(ACDIR(1))+ABS(ACDIR(2))
      FEJER = AMIN1(RTN*ABDIR,FEJER)
      DX = DABS(DZNR)+DABS(DZNI)
      DXT = DX+ABDIR
C                                  IS THE STEP SIZE NEGLIGIBLE
      IF (SNGL(DXT-DX).EQ.ZERO) GO TO 195
      GO TO 105
  175 CONTINUE
C                                  ZN IS REAL
C                                  FACTORIZATION OF POLYNOMIAL BY
C                                    LINEAR FACTOR (Z-X) AS FOLLOWS
C                                    SUM(DU(I)*Z**(N-I)) =
C                                    (Z-X)*SUM(Z(I)*Z**(N-I-1)) + Z(N)
C                                    FOR ALL Z ,
C                                    SO Z(N) IS VALUE OF POLYNOMIAL AT
C                                    Z=X , FIRST DERIVATIVE OF
C                                    POLYNOMIAL AT Z=X IS V , AND
C                                    SECOND DERIVATIVE OF POLYNOMIAL AT
C                                    Z=X IS 2*W . E IS ERROR BOUND FOR
C                                    THE VALUE OF POLYNOMIAL AND DZ(I)
C                                    ARE THE COEFFICIENTS OF QUOTIENT
C                                    POLYNOMIAL.
      DX = DZNR
      X = DX
      DZNI = DZERO
      ABX = ABS(X)
      DV = DA(1)
      V = DV
      W = ZERO
      DZ(1) = DA(2)+DX*DA(1)
      DO 180 I=2,N
         W = V+X*W
         DV = DZ(I-1)+DX*DV
         V = DV
         DZ(I) = DA(I+1)+DX*DZ(I-1)
  180 CONTINUE
      F = SNGL(DZ(N))
      FN = ABS(F)
      E = ABS(SNGL(DA(1)))*TWOD3
      DO 185 I=1,N1
         E = ABS(SNGL(DZ(I)))+ABX*E
  185 CONTINUE
      E = SDEPS*(THREE*ABX*E+ABS(SNGL(DZ(N))))
C                                  HAS AN ACCEPTABLE ZERO BEEN FOUND.
      IF (FN.LE.E) GO TO 205
C                                  HAS THE FUNCTION VALUE DECREASED.
      IF (FN.GE.F0.AND.STARTD) GO TO 215
C                                  FIND THE LAGUERRE STEP AT DZNR.
      IF (ABS(F).GE.ONE) GO TO 190
C                                  IF ABS(V/F) .GT. SINF, THERE IS
C                                    A ZERO WITHIN N*(1.0/SINF) OF ZN.
      IF (ABS(V).GT.ABS(F)*SINF) GO TO 205
  190 R = V/F
      CALL ZQADR (X2N1*W,X2N*V,F,C,CF1,IKER)
C                                  CALCULATE THE FEJER BOUND FOR
C                                    SMALLEST ZERO.
      FEJER = ABS(AC(1))+ABS(AC(2))
      R = XN2N*R
      CDIR = C/CMPLX(R*AC(1)+XN1,R*AC(2))
      ABDIR = ABS(ACDIR(1))+ABS(ACDIR(2))
      FEJER = AMIN1(RTN*ABDIR,FEJER)
      DX = DABS(DZNR)
C                                  IS THE STEP SIZE NEGLIGIBLE.
      DXT = DX+ABDIR
      IF (SNGL(DXT-DX).EQ.ZERO) GO TO 205
      GO TO 105
C                                  CZN IS A COMPLEX ZERO.
C                                  STORE COEFFICIENTS OF QUOTIENT
C                                    POLYNOMIAL IN DA ARRAY.
C                                  DA(1) IS UNCHANGED FOR THE DEFLATED
C                                    POLYNOMIAL.
  195 DO 200 I=3,N
         DA(I-1) = DZ(I-2)
  200 CONTINUE
      Z(N) = CMPLX(SNGL(DZNR),SNGL(DZNI))
      Z(N-1) = CONJG(Z(N))
      GO TO 50
C                                  ZN IS A REAL ZERO.
C                                  STORE COEFFICIENTS OF QUOTIENT
C                                    POLYNOMIAL IN DA ARRAY.
C                                  DA(1) IS UNCHANGED FOR THE DEFLATED
C                                    POLYNOMIAL.
  205 DO 210 I=2,N
         DA(I) = DZ(I-1)
  210 CONTINUE
      Z(N) = CMPLX(SNGL(DZNR),ZERO)
      GO TO 55
C                                  CURRENT LAGUERRE STEP IS
C                                    UNACCEPTABLE.
  215 CONTINUE
      IF (.NOT.STARTD) GO TO 245
C                                  REDUCE PREVIOUS LAGUERRE STEP BY
C                                    HALF.
      ABSCL = HALF*ABSCL
      CL = CMPLX(HALF*ACL(1),HALF*ACL(2))
      DX = DABS(DZNR)+DABS(DZNI)
      DXT = DX+ABSCL
      IF (SNGL(DXT-DX).NE.ZERO) GO TO 115
      IF (FN.LT.E*XN*XN) GO TO 240
  220 CONTINUE
      IF (N.EQ.NDEG) GO TO 230
      DO 225 I=NP1,NDEG
         Z(I-N) = Z(I)
  225 CONTINUE
  230 NTOM = NDEG-N+1
      DO 235 I=NTOM,NDEG
         Z(I) = CMPLX(FINITY,ZERO)
  235 CONTINUE
      IER = 131
      GO TO 9000
  240 IF (DZNI) 195,205,195
  245 CONTINUE
C                                  IF .NOT. STARTD, HAS CZN BEEN ON THE
C                                    INNER CAUCHY RADIUS.
      IF (SPIRAL) GO TO 250
C                                  SET SPIRAL TO .TRUE.. PUT ZN ON THE
C                                    INNER CIRCLE OF THE ANNULUS
C                                    CONTAINING THE SMALLEST ZERO IN
C                                    THE DIRECTION OF THE LAGUERRE STEP.
      SPIRAL = .TRUE.
      CSPIR = CMPLX(OPTFM/XN,ONE)
      ABSCL = RO/(XN*XN)
      C = CMPLX((ACDIR(1)/ABDIR)*RO,(ACDIR(2)/ABDIR)*RO)
      GO TO 255
C                                  SET ZN TO ANOTHER POINT ON THE
C                                    SPIRAL.
  250 C = CSPIR*CMPLX(SNGL(DZNR),SNGL(DZNI))
  255 DZNR = AC(1)
      DZNI = AC(2)
      GO TO 120
 9000 CONTINUE
C                                  RESET ERROR MESSAGE LEVEL TO
C                                    PREVIOUS LEVEL
      CALL UERSET (LEVOLD,LEVEL)
      IF (IER.GT.0) CALL UERTST (IER,6HZPOLR )
 9005 RETURN
      END

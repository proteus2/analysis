C   IMSL ROUTINE NAME   - IRATCU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - RATIONAL WEIGHTED CHEBYSHEV APPROXIMATION
C                           OF A CONTINUOUS FUNCTION.
C
C   USAGE               - CALL IRATCU (F,PHI,G,A,B,L,M,P,Q,WK,IER)
C
C   ARGUMENTS    F      - USER SUPPLIED EXTERNAL FUNCTION
C                           SUBPROGRAM F(X). (INPUT)
C                           F IS TO BE APPROXIMATED OVER THE
C                           INTERVAL (A,B).
C                PHI    - VARIABLE TRANSFORMATION WHICH MUST BE
C                           CONTINUOUS AND MONOTONIC. (INPUT)
C                           IT IS SPECIFIED BY A USER SUPPLIED
C                           EXTERNAL FUNCTION SUBPROGRAM.
C                G      - WEIGHT FUNCTION WHICH MUST BE CONTINUOUS
C                           AND NON VANISHING IN (A,B). (INPUT)
C                           IT IS SPECIFIED BY A USER SUPPLIED EXTERNAL
C                           FUNCTION SUBPROGRAM.
C                A      - LOWER END POINT OF INTERVAL (A,B). (INPUT)
C                B      - UPPER END POINT OF INTERVAL (A,B). (INPUT)
C                L      - DEGREE OF NUMERATOR OF RATIONAL APPROXIMATION
C                           (INPUT)
C                M      - DEGREE OF DENOMINATOR OF RATIONAL
C                           APPROXIMATION. (INPUT)
C                P      - OUTPUT VECTOR OF LENGTH L+1 CONTAINING
C                           THE NUMERATOR COEFFICIENTS.
C                Q      - OUTPUT VECTOR OF LENGTH M+1 CONTAINING
C                           THE DENOMINATOR COEFFICIENTS.
C                         NOTE - IRATCU ALWAYS RETURNS COEFFICIENTS IN
C                           P AND Q FOR THE RATIONAL APPROXIMATION.
C                           HOWEVER, WHEN CONVERGENCE FAILS (IER IS
C                           RETURNED GREATER THAN ZERO) THE
C                           COEFFICIENTS MAY NOT PROVIDE A GOOD
C                           APPROXIMATION.
C                WK     - WORK VECTOR OF LENGTH (NP1+6)*NP1
C                           WHERE NP1 = L+M+2.
C                         ON OUTPUT, THE MIN MAX ERROR OF THE
C                            APPROXIMATION, XLB, IS EQUAL TO THE
C                            ABSOLUTE VALUE OF WK(1).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IRATCU FAILED TO
C                             CONVERGE. A SEQUENCE OF CRITICAL POINTS
C                             WHICH IS NOT MONOTONIC WAS GENERATED.
C                             THIS INDICATES THE POSSIBILITY OF A
C                             DEGENERATE APPROXIMATION (SEE REMARKS).
C                           IER = 130 INDICATES THAT IRATCU FAILED TO
C                             CONVERGE. THE VALUE OF THE ERROR CURVE AT
C                             SOME CRITICAL POINT IS TOO LARGE. THIS
C                             INDICATES THE POSSIBILITY OF POLES IN THE
C                             RATIONAL FUNCTION (SEE REMARKS).
C                           IER = 131 INDICATES THAT THE ALGORITHM DID
C                             NOT CONVERGE WITHIN 20 ITERATIONS.
C                           IER = 132 INDICATES THAT IRATCU FAILED TO
C                             CONVERGE. THE LINEAR SYSTEM THAT DEFINES
C                             THE COEFFICIENTS P(I) AND Q(I) WAS
C                             FOUND TO BE ALGORITHMICALLY SINGULAR BY
C                             LUDATF. THIS INDICATES THE POSSIBILITY OF
C                             A DEGENERATE APPROXIMATION.
C                         WARNING ERROR
C                           IER = 33 INDICATES THAT IRATCU TERMINATED
C                             THE ITERATION BECAUSE THE ERROR WAS
C                             REDUCED AS FAR AS NUMERICALLY POSSIBLE.
C                             A GOOD APPROXIMATION IS RETURNED IN
C                             ARRAYS P AND Q (BUT THIS DOES NOT
C                             NECESSARILY GIVE THE CHEBYSHEV
C                             APPROXIMATION).
C
C   REQD. IMSL ROUTINES - LEQT1F,LUDATN,LUELMN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IRATCU ATTEMPTS TO DETERMINE THE MIN MAX ERROR XLB
C                CORRECT TO 3 DIGITS. IF THIS ERROR IS ZERO OR VERY
C                SMALL (WITH RESPECT TO WORKING PRECISION) IRATCU MAY
C                TERMINATE WITH IER=33. THIS IS AN INDICATION THAT A
C                GOOD APPROXIMATION HAS BEEN GENERATED BUT XLB
C                (MIN MAX ERROR) HAS NOT BEEN ACHIEVED CORRECT TO
C                3 DIGITS.
C            2.  IRATCU MAY FAIL TO CONVERGE IF THE ERROR CURVE
C                IS NON-STANDARD (SEE DOCUMENT REFERENCE FOR
C                FURTHER DISCUSSION). IER=129 USUALLY INDICATES THAT
C                THE APPROXIMATION IS DEGENERATE IN THE SENSE THAT
C                THE NUMERATOR AND DENOMINATOR HAVE COMMON FACTORS.
C                IN THIS CASE, AN APPROXIMATION WITH L AND M EACH
C                REDUCED BY 1 SHOULD BE ATTEMPTED. IER=130 USUALLY
C                INDICATES POLES IN THE RATIONAL APPROXIMATION. THIS
C                CAN USUALLY BE HANDLED BY MODIFYING THE INTERVAL
C                (A,B) OR BY USING SMALLER VALUES OF L AND M. IT IS
C                QUITE OFTEN THE CASE THAT EVEN WHEN IER=129 OR 130,
C                A GOOD APPROXIMATION HAS BEEN GENERATED BUT XLB
C                (MIN MAX ERROR) HAS NOT BEEN ACHIEVED CORRECT TO
C                3 DIGITS.
C            3.  THE TRANSFORMATION PHI AND WEIGHT FUNCTION G ARE
C                CHOSEN BY THE USER AS APPROPRIATE FOR THE PROBLEM
C                TO BE SOLVED. IN MANY CASES, PHI(X)=X IS USED.
C                AS A WEIGHT FUNCTION, G(X)=1 WILL CAUSE THE
C                ABSOLUTE ERROR TO BE MINIMIZED. IF IT IS DESIRED
C                TO MINIMIZE THE RELATIVE ERROR, THEN G(X) SHOULD
C                APPROXIMATE THE MAGNITUDE OF F(X) BUT MUST NOT BE
C                ZERO IN THE INTERVAL (A,B).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IRATCU (F,PHI,G,A,B,L,M,P,Q,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L,M,IER
      REAL               F,PHI,G,A,B,P(1),Q(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IB1,IB2,IB3,IB3P2,IB4,IB5,IB5P2,IB6,IC1,IC2,
     1                   IDON,II,IJK,IJK1,IJK2,IJK3,IJK4,IJK5,IJKN,
     2                   ITNO,IV,J,KD,KK,KRET,LP1,LPM,MP1,N,NP1
      REAL               AEPS,ATEST,AXLB,BB1,DEL,DN,DNP1,EPS,
     1                   HALF,H1,ONE,PI,PZ,P01,P015,P0625,SGN,
     2                   SUM,TEMPL,TEN,TEST,T1,T2,T3,T4,U,XLB,XXI,
     3                   XXM1,Y,Y2,Y3,Z,ZERO,ZZ,Z1,Z2,Z3
      DATA               EPS/1.0E-2/,ONE/1.0/,
     1                   ZERO/0.0/,HALF/0.5/,P01/0.01/,
     2                   P015/0.015/,TEN/10.0/,P0625/0.0625/
      DATA               PI/3.141593/
      DATA               AEPS/1.0E-4/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  INITIALIZATION
      XLB = ZERO
      AXLB = ONE
      ITNO = 1
      ATEST = ZERO
      LP1 = L+1
      LPM = L+M
      N = LPM+1
      NP1 = N+1
      MP1 = M+1
      DN = N
      DNP1 = DN+ONE
C                                  SET UP INDEXING FOR WORK AREA
C                                    WK(I) = WORK SPACE FOR LEQT1F
C                                    WK(IB1+I) = B(I)
C                                    WK(IB2+I) = ERR(I)
C                                    WK(IB3+I) = H(I)
C                                    WK(IB4+I) = X(I)
C                                    WK(IB5+I) = XVAL(I)
C                                    WK(IB6+I+(J-1)*NP1) = A(I,J)
C                                    FOR I=1,...,NP1 AND J=1,...,NP1.
      IB1 = NP1
      IB2 = IB1+NP1
      IB3 = IB2+NP1
      IB4 = IB3+NP1
      IB5 = IB4+NP1
      IB6 = IB5+NP1
      IB3P2 = IB3+2
      IB5P2 = IB5+2
      IDON = IB1+LP1
      DO 5 I=1,LP1
         P(I) = ZERO
    5 CONTINUE
      DO 10 I=1,MP1
         Q(I) = ZERO
   10 CONTINUE
      Q(1) = ONE
C                                  COMPUTE INITIAL GUESSES
C                                  OF CRITICAL POINTS
      XXI = ZERO
      KK = N/2
      BB1 = (B-A)*HALF
      XXM1 = (B+A)*HALF
      IF (KK.LT.1) GO TO 20
      KD = IB4+1
      DO 15 I=1,KK
         XXI = XXI+ONE
         ZZ = COS(PI*(ONE-XXI/DN))*BB1
         IJK4 = KD+I
         WK(IJK4) = ZZ+XXM1
         IJK4 = IB5-I
         WK(IJK4) = XXM1-ZZ
   15 CONTINUE
   20 WK(IB4+1) = A
      WK(IB5) = B
   25 DO 30 I=1,NP1
         IJK4 = IB4+I
         IJK5 = IB5+I
         WK(IJK5) = PHI(WK(IJK4))
   30 CONTINUE
C                                  SET UP AND SOLVE LINEAR EQUATIONS
      KRET = 1
      GO TO 50
   35 XLB = (WK(IB2)+XLB*DN)/DNP1
      IF (M.EQ.0) GO TO 90
      I = 1
   40 CONTINUE
      KRET = 2
      GO TO 50
   45 TEST = ABS(XLB-WK(IB2))
      XLB = (WK(IB2)+XLB*DN)/DNP1
      IF (TEST.LE.P01*ABS(XLB).OR.I.GE.3) GO TO 80
      I = I+1
      GO TO 40
   50 KK = L+2
C                                  SET UP LINEAR EQUATIONS
      SGN = ONE
      DO 75 IV=1,NP1
         SGN = -SGN
         IJK = IB6+IV
         WK(IJK) = ONE
         IJK5 = IB5+IV
         TEST = WK(IJK5)
         IF (L.LE.0) GO TO 60
         IC1 = IB6+IV
         IC2 = (LP1-2)*NP1+IC1
         DO 55 J=IC1,IC2,NP1
            WK(J+NP1) = WK(J)*TEST
   55    CONTINUE
   60    IJK1 = IB1+IV
         IJK4 = IB4+IV
         WK(IJK1) = F(WK(IJK4))
         IF (M.LE.0) GO TO 70
         TEMPL = XLB*G(WK(IJK4))*SGN-WK(IJK1)
         IC1 = IB6+(KK-1)*NP1+IV
         WK(IC1) = TEST*TEMPL
         IF (LPM .LT. KK) GO TO 70
         IC2 = (LPM-KK)*NP1+IC1
         DO 65 J=IC1,IC2,NP1
            WK(J+NP1) = WK(J)*TEST
   65    CONTINUE
   70    IC1 = IB6+NP1*(NP1-1)+IV
         WK(IC1) = SGN*G(WK(IJK4))
   75 CONTINUE
C                                  SOLVE LINEAR EQUATIONS
      CALL LEQT1F (WK(IB6+1),1,NP1,NP1,WK(IB1+1),0,WK(1),IER)
      IF (IER.NE.0) GO TO 275
      GO TO (35,45), KRET
   80 DO 85 I=1,M
         IJK1 = IDON+I
         Q(I+1) = WK(IJK1)
   85 CONTINUE
   90 CONTINUE
      DO 95 I=1,LP1
         IJK1 = IB1+I
         P(I) = WK(IJK1)
   95 CONTINUE
      U = ONE
      IF (XLB.NE.ZERO) U = SIGN(ONE,XLB)
      Z1 = ZERO
C                                  SEARCH FOR NEW CRITICAL POINTS
      IF (N.GT.1) GO TO 100
      WK(IB3+1) = P015*(WK(IB4+2)-WK(IB4+1))
      WK(IB3+2) = -WK(IB3+1)
      GO TO 110
  100 DO 105 I=2,N
         IJK3 = IB3+I
         IJK4 = IB4+I+1
         IJKN = IJK4-2
         WK(IJK3) = P015*(WK(IJK4)-WK(IJKN))
  105 CONTINUE
      WK(IB3+1) = WK(IB3+2)*HALF
      IJK3 = IB3+N
      WK(IB4) = -WK(IJK3)*HALF
  110 I = 1
      AXLB = ZERO
  115 IJK4 = IB4+I
      Y2 = WK(IJK4)
      IJK3 = IB3+I
      H1 = WK(IJK3)
      Y3 = Y2+H1
      PZ = Y2
      KRET = 1
C                                  COMPUTE APPROXIMATION ERROR AT PZ
  120 T1 = F(PZ)
      T2 = PHI(PZ)
      T3 = WK(IDON)
      IF (L.LT.1) GO TO 130
      DO 125 II=1,L
         IJK1 = IDON-II
         T3 = T3*T2+WK(IJK1)
  125 CONTINUE
  130 CONTINUE
      T4 = ZERO
      IF (M.LT.1) GO TO 140
      DO 135 II=1,M
         IJK1 = IB2-II
         T4 = (T4+WK(IJK1))*T2
  135 CONTINUE
  140 CONTINUE
      T4 = T4+ONE
      T2 = T3/T4
C                                  ESTIMATE MAXIMUM VALUE OF F(X)/G(X)
      T3 = G(PZ)
      T4 = ABS(T1/T3)
      ATEST = AMAX1(ATEST,T4)
      DEL = (T2-T1)/T3
C                                  DEL IS THE APPROXIMATION ERROR AT PZ
C                                    RETURN TO THE CALLING SEGEMENT
      GO TO (145,150,170,190,215,245), KRET
C                                  DIRECT SEARCH FOR ADJUSTED CRITICAL
C                                    POINTS
  145 Z2 = DEL*U
      PZ = Y3
      KRET = 2
      GO TO 120
  150 Z3 = DEL*U
      IF (Z3.GT.Z2) GO TO 155
      H1 = -H1
      Z = Z3
      Z3 = Z2
      Z2 = Z
      Y = Y3
      Y3 = Y2
      Y2 = Y
  155 Y = Y3+H1
      IF (Y.GE.A) GO TO 160
      Y = A
      GO TO 185
  160 IF (Y.LE.B) GO TO 165
      Y = B
      GO TO 185
  165 PZ = Y
      KRET = 3
      GO TO 120
  170 Z = DEL*U
      IF (Z.LE.Z3) GO TO 175
      Y2 = Y3
      Y3 = Y
      Z2 = Z3
      Z3 = Z
      GO TO 155
  175 Y = -Z3-Z3+Z+Z2
      IF (Y.EQ.ZERO) GO TO 180
      Y = (Y2+Y3)*HALF+H1*(Z2-Z3)/Y
      GO TO 185
  180 Y = Y3
  185 WK(IJK4) = Y
      PZ = Y
      KRET = 4
      GO TO 120
  190 IJK2 = IB2+I
      WK(IJK2) = DEL
      U = -U
      Z = ABS(DEL)
      AXLB = AMAX1(AXLB,Z)
C                                  CHECK FOR PROPER ORDERING OF
C                                    ADJUSTED CRITICAL POINTS
      IF (I.EQ.1) GO TO 195
      IF (I.EQ.NP1) GO TO 195
      IF (WK(IJK4).GT.WK(IJK4-1)) GO TO 195
      IF (ABS(AXLB).LE.AEPS*ATEST) GO TO 285
      IER = 129
      GO TO 290
C                                  CHECK FOR POLES IN TRIAL
C                                    RATIONAL APPROXIMANT
  195 CONTINUE
      IF (Z.LT.TEN*ATEST) GO TO 200
      IER = 130
      GO TO 290
  200 CONTINUE
      Y = ABS(XLB)
      Z2 = ONE
      IF (Y.NE.ZERO) Z2 = ABS(Z-Y)/Y
      IF (Z1.LT.Z2) Z1 = Z2
      IF (I.GE.NP1) GO TO 205
      I = I+1
      GO TO 115
  205 CONTINUE
C                                  SEARCH FOR ONE EXTRA EXTREMAL POINT
C                                    BETWEEN THE ENDPOINTS OF THE
C                                    INTERVAL AND THE PRESENT CRITICAL
C                                    POINTS.
      IJK4 = IB4+1
      IF (WK(IJK4).LE.A) GO TO 230
      H1 = (WK(IJK4)-A)*P0625
      U = -ONE
      IF (XLB.NE.ZERO) U = -SIGN(ONE,XLB)
      Z3 = ZERO
      Y = A
      I = 1
  210 PZ = Y
      KRET = 5
      GO TO 120
  215 Z = DEL*U
      IF (Z.LE.Z3) GO TO 220
      Z3 = Z
      Z2 = Y
  220 Y = Y+H1
      I = I+1
      IF (I.LE.16) GO TO 210
      Z = ABS(XLB)
      IF (Z3.LE.Z) GO TO 230
      DO 225 I=2,NP1
         IJK2 = IB3P2-I
         WK(IJK2) = WK(IJK2-1)
         IJK4 = IB5P2-I
         WK(IJK4) = WK(IJK4-1)
  225 CONTINUE
      IJK4 = IB4+1
      WK(IJK4) = Z2
      IJK2 = IB2+1
      WK(IJK2) = Z3*U
      AXLB = AMAX1(AXLB,Z3)
      GO TO 260
  230 IF (WK(IB5).GT.WK(IB5-1)) GO TO 235
      IF (ABS(AXLB).LE.AEPS*ATEST) GO TO 285
      IER = 129
      GO TO 290
  235 CONTINUE
      IF (B.LE.WK(IB5)) GO TO 265
      H1 = P0625*(B-WK(IB5))
      Z3 = ZERO
      Y = B
      U = -ONE
      IF (WK(IB3).NE.ZERO) U = -SIGN(ONE,WK(IB3))
      I = 1
  240 PZ = Y
      KRET = 6
      GO TO 120
  245 Z = DEL*U
      IF (Z3.GE.Z) GO TO 250
      Z3 = Z
      Z2 = Y
  250 Y = Y-H1
      I = I+1
      IF (I.LE.16) GO TO 240
      Z = ABS(XLB)
      IF (Z.GE.Z3) GO TO 265
      DO 255 I=1,N
         IJK2 = IB2+I
         WK(IJK2) = WK(IJK2+1)
         IJK4 = IB4+I
         WK(IJK4) = WK(IJK4+1)
  255 CONTINUE
      WK(IB5) = Z2
      WK(IB3) = Z3*U
      AXLB = AMAX1(AXLB,Z3)
  260 XLB = -XLB
      Y = Z
      Z = ONE
      IF (Y.NE.ZERO) Z = ABS(Z3-Y)/Y
      IF (Z1.LT.Z) Z1 = Z
  265 CONTINUE
C                                  CHECK FOR RELATIVE CONVERGENCE OF
C                                    XLB
      IF (Z1.LE.EPS) GO TO 290
C                                  CHECK FOR ABSOLUTE CONVERGENCE OF
C                                    XLB
      IF (ABS(AXLB).LE.AEPS*ATEST) GO TO 285
C                                  CHECK FOR NO CONVERGENCE
      IF (ITNO.GE.20) GO TO 280
      SUM = ZERO
      SGN = ONE
      DO 270 I=1,NP1
         IJK2 = IB2+I
         SUM = SUM+SGN*WK(IJK2)
         SGN = -SGN
  270 CONTINUE
      XLB = SUM/DNP1
      ITNO = ITNO+1
      GO TO 25
  275 CONTINUE
      IER = 132
      GO TO 290
  280 CONTINUE
      IER = 131
      GO TO 290
  285 CONTINUE
      IER = 33
  290 CONTINUE
      WK(1) = SIGN(AXLB,XLB)
      WK(2) = AEPS
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HIRATCU)
 9005 RETURN
      END

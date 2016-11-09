C   IMSL ROUTINE NAME   - MDGCI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - INVERSE OF A GENERAL CUMULATIVE PROBABILITY
C                           DISTRIBUTION FUNCTION, GIVEN ORDINATES OF
C                           THE DENSITY
C
C   USAGE               - CALL MDGCI (P,F,M,IOPT,B,C,X,IER)
C
C   ARGUMENTS    P      - INPUT PROBABILITY IN THE EXCLUSIVE RANGE
C                           (0,1).
C                F      - INPUT ARRAY OF LENGTH M CONTAINING THE
C                           PROBABILITY DENSITY FUNCTION ORDINATES.
C                M      - INPUT NUMBER OF ORDINATES SUPPLIED. M MUST BE
C                           GREATER THAN 1 IF LINEAR INTERPOLATION IS
C                           DESIRED AND GREATER THAN 3 IF A CURVE
C                           FITTED THROUGH THE ORDINATES IS DESIRED.
C                IOPT   - INPUT OPTION PARAMETER.
C                           IOPT = 1 IMPLIES LINEAR INTERPOLATION IS
C                             DESIRED AND THE DATA ARE EQUALLY SPACED.
C                           IOPT = 2 IMPLIES LINEAR INTERPOLATION IS
C                             DESIRED AND THE DATA ARE NOT NECESSARILY
C                             EQUALLY SPACED.
C                           IOPT = 3 IMPLIES A CURVE FITTED THROUGH THE
C                             DATA IS DESIRED AND THE DATA ARE EQUALLY
C                             SPACED.
C                           IOPT = 4 IMPLIES A CURVE FITTED THROUGH THE
C                             DATA IS DESIRED AND THE DATA ARE NOT
C                             NECESSARILY EQUALLY SPACED.
C                B      - INPUT ARRAY OF LENGTH M IF IOPT = 2, 3, OR 4.
C                           OTHERWISE, B IS OF LENGTH 2. IF IOPT = 2 OR
C                           4, B(I) IS THE ABSCISSA CORRESPONDING TO
C                           F(I). THE ABSCISSAE MUST BE SPECIFIED IN
C                           INCREASING ORDER. IF IOPT = 1 OR 3, B(1)
C                           IS THE LOWER ENDPOINT OF THE SUPPORT OF THE
C                           DISTRIBUTION (CORRESPONDING TO F(1)) AND
C                           B(2) IS THE UPPER ENDPOINT (CORRESPONDING
C                           TO F(M)).
C                C      - WORK AREA. IF IOPT = 3 OR 4, THEN C MUST BE
C                           OF LENGTH AT LEAST 3*M. OTHERWISE, ONLY
C                           C(1) IS USED.
C                X      - OUTPUT VALUE SUCH THAT P IS THE PROBABILITY
C                           THAT THE RANDOM VARIABLE IS LESS THAN OR
C                           EQUAL TO X.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES AN ERROR OCCURRED IN
C                             IMSL SUBROUTINE IQHSCU.
C                           IER = 130 INDICATES AN ELEMENT OF F WAS
C                             NEGATIVE.
C                           IER = 131 INDICATES THE ELEMENTS OF B ARE
C                             NOT IN INCREASING ORDER.
C                           IER = 132 INDICATES M IS OUT OF RANGE.
C                           IER = 133 INDICATES IOPT IS OUT OF RANGE.
C                           IER = 134 INDICATES P IS NOT IN THE
C                             EXCLUSIVE RANGE (0,1).
C
C   REQD. IMSL ROUTINES - DCSQDU,IQHSCU,MDGD,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDGCI  (P,F,M,IOPT,B,C,X,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,IOPT,M
      REAL               B(M),C(1),F(M),P,X
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,IM1,J,MM1
      REAL               A1,A2,A,BB,BP3,BP,CC,CP,DELTA,DP,D,EPS,E,
     1                   FOURTH,FOUR,H26,H2,HALF,H,OLDA,PRED,QRED,RK,
     2                   RL,SIXTH,THIRD,TWNT7,XI,NINTH,RT3D2,CI2,CIM12
      COMPLEX            AL,BT,BZ,CZ,Z,ZZ
      DATA               HALF/.5E0/,FOURTH/.25E0/,FOUR/4.0E0/
      DATA               NINTH/.1111111E0/
      DATA               RT3D2/.8660254E0/
      DATA               SIXTH/.1666667E0/
      DATA               THIRD/.3333333E0/
      DATA               TWNT7/.0370370E0/
      DATA               EPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (P.LT.1.0.AND.P.GT.0.0) GO TO 5
      IER = 134
      GO TO 9000
    5 IF (IOPT.GT.0.AND.IOPT.LT.5) GO TO 10
      IER = 133
      GO TO 9000
   10 DO 15 I=1,M
         IF (F(I).GE.0.0) GO TO 15
         IER = 130
         GO TO 9000
   15 CONTINUE
      GO TO (45,20,120,95), IOPT
C                                  UNEQUALLY SPACED - LINEAR
   20 IF (M.GT.1) GO TO 25
      IER = 132
      GO TO 9000
   25 C(1) = B(M)
      CALL MDGD (F,M,IOPT,B,C,A1,IER)
      IF (IER.NE.0) GO TO 9000
      A1 = P*A1
      A2 = 0.0
      I = 2
      IM1 = 1
      OLDA = 0.0
      IF (M.LT.3) GO TO 40
      DO 30 II=2,M
         I = II
         A2 = A2+HALF*(B(I)-B(IM1))*(F(IM1)+F(I))
         IF (A2.GE.A1) GO TO 35
         IM1 = I
         OLDA = A2
   30 CONTINUE
      X = B(I)
      GO TO 9005
   35 IF (A2.GT.A1) GO TO 40
      X = B(I)
      GO TO 9005
   40 H = B(I)-B(IM1)
      X = B(IM1)
      GO TO 75
C                                  EQUALLY SPACED - LINEAR
   45 IF (M.GT.1) GO TO 50
      IER = 132
      GO TO 9000
   50 C(1) = B(2)
      IF (B(2).GE.B(1)) GO TO 55
      IER = 131
      GO TO 9000
   55 CALL MDGD (F,M,IOPT,B,C,A1,IER)
      IF (IER.NE.0) GO TO 9000
      A1 = P*A1
      A2 = 0.0
      OLDA = 0.0
      I = 2
      IM1 = 1
      MM1 = M-1
      H = (B(2)-B(1))/MM1
      H2 = HALF*H
      IF (M.LT.3) GO TO 70
      DO 60 II=2,M
         I = II
         A2 = A2+H2*(F(IM1)+F(I))
         IF (A2.GE.A1) GO TO 65
         OLDA = A2
         IM1 = I
   60 CONTINUE
      X = B(I)
      GO TO 9005
   65 IF (A2.GT.A1) GO TO 70
      X = B(I)
      GO TO 9005
   70 X = B(IM1)
   75 E = A1-OLDA
      D = F(I)-F(IM1)
      IF (F(IM1).GE.EPS) GO TO 80
      DELTA = SQRT((H+H)*E/F(I))
      GO TO 90
   80 IF (D.GE.EPS) GO TO 85
      DELTA = E/F(IM1)
      GO TO 90
   85 DELTA = -H*F(IM1)/D*(1.0-SQRT(1.0+(E+E)*D/(H*F(IM1)*F(IM1))))
   90 X = DELTA+X
      GO TO 9005
C                                  UNEQUALLY SPACED - FITTED CURVE
   95 IF (M.GT.3) GO TO 100
      IER = 132
      GO TO 9000
  100 C(1) = B(M)
      CALL MDGD (F,M,IOPT,B,C,A1,IER)
      IF (IER.NE.0) GO TO 9000
      A1 = P*A1
      A2 = 0.0
      OLDA = 0.0
      IM1 = 1
      MM1 = M-1
      CIM12 = C(M+1)
      DO 105 II=2,MM1
         I = II
         CI2 = C(M+I)
         H = B(I)-B(IM1)
         A2 = A2+HALF*H*(F(I)+F(IM1)-(CI2+CIM12)*H*H*SIXTH)
         IF (A2.GE.A1) GO TO 110
         OLDA = A2
         IM1 = I
         CIM12 = CI2
  105 CONTINUE
      XI = B(I)
      GO TO 145
  110 IF (A2.GT.A1) GO TO 115
      X = B(I)
      GO TO 9005
  115 XI = B(IM1)
      GO TO 145
C                                  EQUALLY SPACED - CURVE
  120 IF (M.GT.3) GO TO 125
      IER = 132
      GO TO 9000
  125 C(1) = B(2)
      CALL MDGD (F,M,IOPT,B,C,A1,IER)
      IF (IER.NE.0) GO TO 9000
      A1 = P*A1
      A2 = 0.0
      OLDA = 0.0
      MM1 = M-1
      B(2) = B(M)
      H = (B(2)-B(1))/MM1
      IM1 = 1
      MM1 = M-1
      H26 = (H*H)*SIXTH
      CIM12 = C(M+1)
      DO 130 II=2,MM1
         I = II
         CI2 = C(M+I)
         A2 = A2+HALF*H*(F(I)+F(IM1)-(CI2+CIM12)*H26)
         IF (A2.GE.A1) GO TO 135
         OLDA = A2
         IM1 = I
         CIM12 = CI2
  130 CONTINUE
      XI = B(I)
      GO TO 145
  135 IF (A2.GT.A1) GO TO 140
      X = B(I)
      GO TO 9005
  140 XI = B(IM1)
  145 E = A1-OLDA
      IF (C(M+M+IM1).NE.0.0) GO TO 165
C                                  NOT A QUARTIC
      IF (C(M+IM1).NE.0) GO TO 160
C                                  NOT A CUBIC
      IF (C(IM1).NE.0) GO TO 150
C                                  MUST BE A STRAIGHT LINE
      X = XI
      IF (F(IM1).NE.0.0) X = E/F(IM1)+XI
      GO TO 9005
C                                  SOLVE A QUADRATIC
  150 A = C(IM1)*HALF
      BB = F(IM1)
      CC = -E
      D = SQRT(ABS(BB*BB-FOUR*A*CC))
      RK = 1.0E0
      OLDA = H
      A1 = H
      DO 155 I=1,2
         X = ((RK*D-BB)/C(IM1))+XI
         IF (X.GE.XI.AND.X.LE.(XI+H)) GO TO 9005
         OLDA = AMIN1(ABS(X-XI),OLDA)
         A1 = AMIN1(ABS(X-(XI+H)),A1)
         RK = -RK
  155 CONTINUE
      X = XI+H-A1
      IF (OLDA.LT.A1) X = XI+OLDA
      GO TO 9005
C                                  SOLVE A CUBIC
  160 A = 1.0/C(M+IM1)*THIRD
      BB = C(IM1)*HALF*A
      CC = F(IM1)*A
      D = -E*A
      QRED = THIRD*CC-NINTH*BB*BB
      PRED = SIXTH*(CC*BB-3.0*D)-TWNT7*BB*BB*BB
      Z = CSQRT(CMPLX(QRED*QRED*QRED+PRED+PRED,0.0))
      AL = CEXP(CMPLX(THIRD,0.0)*CLOG(PRED+Z))
      BT = CEXP(CMPLX(THIRD,0.0)*CLOG(PRED-Z))
      OLDA = H
      A1 = H
      X = REAL(AL+BT) - BB*THIRD + XI
      IF (X.GE.XI.AND.X.LE.(XI+H)) GO TO 9005
      OLDA = AMIN1(ABS(X-XI),OLDA)
      A1 = AMIN1(ABS(X-(XI+H)),A1)
      X = XI-HALF*REAL(AL+BT)-BB-THIRD+REAL(CMPLX(0.0,1.0)*
     *  CMPLX(RT3D2,0.0)*(AL-BT))
      IF (X.GE.XI.AND.X.LE.(XI+H)) GO TO 9005
      OLDA = AMIN1(ABS(X-XI),OLDA)
      A1 = AMIN1(ABS(X-(XI+H)),A1)
      X = XI-HALF*REAL(AL+BT)-BB-THIRD-REAL(CMPLX(0.0,1.0)*
     *  CMPLX(RT3D2,0.0)*(AL-BT))
      IF (X.GE.XI.AND.X.LE.(XI+H)) GO TO 9005
      OLDA = AMIN1(ABS(X-XI),OLDA)
      A1 = AMIN1(ABS(X-(XI+H)),A1)
      X = XI+H-A1
      IF (OLDA.LT.A1) X = XI+OLDA
      GO TO 9005
C                                  SOLVE QUARTIC
  165 A = 1.0/C(M+M+IM1)*FOURTH
      BB = C(M+IM1)*THIRD*A
      CC = C(IM1)*HALF*A
      D = F(IM1)*A
      E = -E*A
      BP = -HALF*CC
      CP = FOURTH*(BB*D-FOUR*E)
      DP = .125E0*(FOUR*CC*E-BB*BB*E-D*D)
      PRED = CP-(BP*BP)*THIRD
      BP3 = BP*BP*BP
      QRED = DP-(BP*CP*THIRD)+TWNT7*(BP3+BP3)
      Z = CMPLX(-QRED*HALF,0.0E0)
      Z = Z+CSQRT(CMPLX(PRED*PRED*PRED*TWNT7+QRED*QRED*FOURTH,0.0E0))
      Z = CEXP(CMPLX(THIRD,0.0)*CLOG(Z))
      Z = Z-CMPLX(PRED,0.0E0)/(3.0E0*Z)+CMPLX(CC*SIXTH,0.0E0)
      AL = CSQRT(Z+Z+CMPLX(FOURTH*BB*BB-CC,0.0E0))
      BT = (BB*Z-CMPLX(D,0.0E0))/(AL+AL)
      RK = 1.0
      RL = 1.0
      OLDA = H
      A1 = H
      DO 175 I=1,2
         BZ = CMPLX(HALF*BB,0.0E0)+RK*AL
         CZ = Z+RK*BT
         DO 170 J=1,2
            ZZ = HALF*(-BZ+RL*CSQRT(BZ*BZ-FOUR*CZ))+CMPLX(XI,0.0E0)
            OLDA = AMIN1(CABS(ZZ-CMPLX(XI,0.0)),OLDA)
            A1 = AMIN1(CABS(ZZ-CMPLX(XI+H,0.0)),A1)
            X = AIMAG(ZZ)
            IF (ABS(X).GT.EPS) GO TO 170
            X = REAL(ZZ)
            IF (X.GE.XI.AND.X.LE.(XI+H)) GO TO 9005
            RL = -RL
  170    CONTINUE
         RK = -RK
  175 CONTINUE
      X = XI+H-A1
      IF (OLDA.LT.A1) X = XI+OLDA
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMDGCI )
 9005 RETURN
      END

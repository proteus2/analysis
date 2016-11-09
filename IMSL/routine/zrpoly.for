C   IMSL ROUTINE NAME   - ZRPOLY
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZEROS OF A POLYNOMIAL WITH REAL
C                           COEFFICIENTS (JENKINS-TRAUB)
C
C   USAGE               - CALL ZRPOLY (A,NDEG,Z,IER)
C
C   ARGUMENTS    A      - INPUT REAL VECTOR OF LENGTH NDEG+1
C                           CONTAINING THE COEFFICIENTS IN ORDER OF
C                           DECREASING POWERS OF THE VARIABLE.
C                NDEG   - INPUT INTEGER DEGREE OF POLYNOMIAL.
C                           NDEG MUST BE GREATER THAN 0 AND LESS
C                           THAN 101.
C                Z      - OUTPUT COMPLEX VECTOR OF LENGTH NDEG
C                           CONTAINING THE COMPUTED ROOTS OF THE
C                           POLYNOMIAL.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*NDEG. AN APPROPRIATE
C                           EQUIVALENCE STATEMENT MAY BE REQUIRED.
C                           SEE DOCUMENT EXAMPLE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129, INDICATES THAT THE DEGREE OF THE
C                             POLYNOMIAL IS GREATER THAN 100 OR LESS
C                             THAN 1.
C                           IER=130, INDICATES THAT THE LEADING
C                             COEFFICIENT IS ZERO.
C                           IER=131, INDICATES THAT ZRPOLY FOUND FEWER
C                             THAN NDEG ZEROS. IF ONLY M ZEROS ARE
C                             FOUND, Z(J),J=M+1,...,NDEG ARE SET TO
C                             POSITIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZRPQLB,ZRPQLC,ZRPQLD,ZRPQLE,
C                           ZRPQLF,ZRPQLG,ZRPQLH,ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPOLY (A,NDEG,Z,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NDEG,IER
      REAL               A(1),Z(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,JJ,I,NM1,ICNT,N2,L,NZ,NPI
      REAL               ETA,RMRE,RINFP,REPSP,RADIX,RLO,XX,YY,SINR,
     1                   COSR,RMAX,RMIN,X,SC,XM,FF,DX,DF,BND,XXX,ARE
      REAL               PT(101)
      REAL               TEMP(101),P(101),QP(101),RK(101),QK(101),
     1                   SVK(101)
      REAL               SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   T,AA,BB,CC,FACTOR,REPSR1,ZERO,ONE,FN
      LOGICAL            ZEROK
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
C                                  THE FOLLOWING STATEMENTS SET MACHINE
C                                    CONSTANTS USED IN VARIOUS PARTS OF
C                                    THE PROGRAM. THE MEANING OF THE
C                                    FOUR CONSTANTS ARE - REPSR1 THE
C                                    MAXIMUM RELATIVE REPRESENTATION
C                                    ERROR WHICH CAN BE DESCRIBED AS
C                                    THE SMALLEST POSITIVE FLOATING
C                                    POINT NUMBER SUCH THAT 1.+REPSR1 IS
C                                    GREATER THAN 1
C                                  RINFP THE LARGEST FLOATING-POINT
C                                    NUMBER
C                                  REPSP THE SMALLEST POSITIVE
C                                    FLOATING-POINT NUMBER IF THE
C                                    EXPONENT RANGE DIFFERS IN SINGLE
C                                    AND DOUBLE PRECISION THEN REPSP
C                                    AND RINFP SHOULD INDICATE THE
C                                    SMALLER RANGE
C                                  RADIX THE BASE OF THE FLOATING-POINT
C                                    NUMBER SYSTEM USED
      DATA               RINFP/Z7FFFFFFF/
      DATA               REPSP/Z00100000/
      DATA               RADIX/16.0/
      DATA               REPSR1/Z3C100000/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  ZRPOLY USES SINGLE PRECISION
C                                    CALCULATIONS FOR SCALING, BOUNDS
C                                    AND ERROR CALCULATIONS.
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NDEG .GT. 100 .OR. NDEG .LT. 1) GO TO 165
      ETA = REPSR1
      ARE = ETA
      RMRE = ETA
      RLO = REPSP/ETA
C                                  INITIALIZATION OF CONSTANTS FOR
C                                    SHIFT ROTATION
      XX = .7071068
      YY = -XX
      SINR = .9975641
      COSR = -.06975647
      N = NDEG
      NN = N+1
C                                  ALGORITHM FAILS IF THE LEADING
C                                    COEFFICIENT IS ZERO.
      IF (A(1).NE.ZERO) GO TO 5
      IER = 130
      GO TO 9000
C                                  REMOVE THE ZEROS AT THE ORIGIN IF
C                                    ANY
    5 IF (A(NN).NE.ZERO) GO TO 10
      J = NDEG-N+1
      JJ = J+NDEG
      Z(J) = ZERO
      Z(JJ) = ZERO
      NN = NN-1
      N = N-1
      IF (NN.EQ.1) GO TO 9005
      GO TO 5
C                                  MAKE A COPY OF THE COEFFICIENTS
   10 DO 15 I=1,NN
         P(I) = A(I)
   15 CONTINUE
C                                  START THE ALGORITHM FOR ONE ZERO
   20 IF (N.GT.2) GO TO 30
      IF (N.LT.1) GO TO 9005
C                                  CALCULATE THE FINAL ZERO OR PAIR OF
C                                    ZEROS
      IF (N.EQ.2) GO TO 25
      Z(NDEG) = -P(2)/P(1)
      Z(NDEG+NDEG) = ZERO
      GO TO 145
   25 CALL ZRPQLI (P(1),P(2),P(3),Z(NDEG-1),Z(NDEG+NDEG-1),Z(NDEG),
     1   Z(NDEG+NDEG))
      GO TO 145
C                                  FIND LARGEST AND SMALLEST MODULI OF
C                                    COEFFICIENTS.
   30 RMAX = 0.
      RMIN = RINFP
      DO 35 I=1,NN
         X = ABS(P(I))
         IF (X.GT.RMAX) RMAX = X
         IF (X.NE.0..AND.X.LT.RMIN) RMIN = X
   35 CONTINUE
C                                  SCALE IF THERE ARE LARGE OR VERY
C                                    SMALL COEFFICIENTS COMPUTES A
C                                    SCALE FACTOR TO MULTIPLY THE
C                                    COEFFICIENTS OF THE POLYNOMIAL.
C                                    THE SCALING IS DONE TO AVOID
C                                    OVERFLOW AND TO AVOID UNDETECTED
C                                    UNDERFLOW INTERFERING WITH THE
C                                    CONVERGENCE CRITERION.
C                                  THE FACTOR IS A POWER OF THE BASE
      SC = RLO/RMIN
      IF (SC.GT.1.0) GO TO 40
      IF (RMAX.LT.10.) GO TO 55
      IF (SC.EQ.0.) SC = REPSP*RADIX*RADIX
      GO TO 45
   40 IF (RINFP/SC.LT.RMAX) GO TO 55
   45 L = ALOG(SC)/ALOG(RADIX)+.5
      IF (L .EQ. 0) GO TO 55
      FACTOR = RADIX**L
      DO 50 I=1,NN
   50 P(I) = FACTOR*P(I)
C                                  COMPUTE LOWER BOUND ON MODULI OF
C                                    ZEROS.
   55 DO 60 I=1,NN
   60 PT(I) = ABS(P(I))
      PT(NN) = -PT(NN)
C                                  COMPUTE UPPER ESTIMATE OF BOUND
      X = EXP((ALOG(-PT(NN))-ALOG(PT(1)))/N)
      IF (PT(N).EQ.0.) GO TO 65
C                                  IF NEWTON STEP AT THE ORIGIN IS
C                                    BETTER, USE IT.
      XM = -PT(NN)/PT(N)
      IF (XM.LT.X) X = XM
C                                  CHOP THE INTERVAL (0,X) UNTIL FF.LE.0
   65 XM = X*.1
      FF = PT(1)
      DO 70 I=2,NN
   70 FF = FF*XM+PT(I)
      IF (FF.LE.0.) GO TO 75
      X = XM
      GO TO 65
   75 DX = X
C                                  DO NEWTON ITERATION UNTIL X
C                                    CONVERGES TO TWO DECIMAL PLACES
   80 IF (ABS(DX/X).LE..005) GO TO 90
      FF = PT(1)
      DF = FF
      DO 85 I=2,N
         FF = FF*X+PT(I)
         DF = DF*X+FF
   85 CONTINUE
      FF = FF*X+PT(NN)
      DX = FF/DF
      X = X-DX
      GO TO 80
   90 BND = X
C                                  COMPUTE THE DERIVATIVE AS THE INTIAL
C                                    K POLYNOMIAL AND DO 5 STEPS WITH
C                                    NO SHIFT
      NM1 = N-1
      FN = ONE/N
      DO 95 I=2,N
   95 RK(I) = (NN-I)*P(I)*FN
      RK(1) = P(1)
      AA = P(NN)
      BB = P(N)
      ZEROK = RK(N).EQ.ZERO
      DO 115 JJ=1,5
         CC = RK(N)
         IF (ZEROK) GO TO 105
C                                  USE SCALED FORM OF RECURRENCE IF
C                                    VALUE OF K AT 0 IS NONZERO
         T = -AA/CC
         DO 100 I=1,NM1
            J = NN-I
            RK(J) = T*RK(J-1)+P(J)
  100    CONTINUE
         RK(1) = P(1)
         ZEROK = ABS(RK(N)).LE.ABS(BB)*ETA*10.
         GO TO 115
C                                  USE UNSCALED FORM OF RECURRENCE
  105    DO 110 I=1,NM1
            J = NN-I
            RK(J) = RK(J-1)
  110    CONTINUE
         RK(1) = ZERO
         ZEROK = RK(N).EQ.ZERO
  115 CONTINUE
C                                  SAVE K FOR RESTARTS WITH NEW SHIFTS
      DO 120 I=1,N
  120 TEMP(I) = RK(I)
C                                  LOOP TO SELECT THE QUADRATIC
C                                    CORRESPONDING TO EACH NEW SHIFT
      DO 140 ICNT=1,20
C                                  QUADRATIC CORRESPONDS TO A DOUBLE
C                                    SHIFT TO A NON-REAL POINT AND ITS
C                                    COMPLEX CONJUGATE. THE POINT HAS
C                                    MODULUS BND AND AMPLITUDE ROTATED
C                                    BY 94 DEGREES FROM THE PREVIOUS
C                                    SHIFT
         XXX = COSR*XX-SINR*YY
         YY = SINR*XX+COSR*YY
         XX = XXX
         SR = BND*XX
         SI = BND*YY
         U = -SR-SR
         V = BND*BND
C                                  SECOND STAGE CALCULATION, FIXED
C                                    QUADRATIC
         CALL ZRPQLB (20*ICNT,NZ)
         IF (NZ.EQ.0) GO TO 130
C                                  THE SECOND STAGE JUMPS DIRECTLY TO
C                                    ONE OF THE THIRD STAGE ITERATIONS
C                                    AND RETURNS HERE IF SUCCESSFUL.
C                                  DEFLATE THE POLYNOMIAL, STORE THE
C                                    ZERO OR ZEROS AND RETURN TO THE
C                                    MAIN ALGORITHM.
         J = NDEG-N+1
         JJ = J+NDEG
         Z(J) = SZR
         Z(JJ) = SZI
         NN = NN-NZ
         N = NN-1
         DO 125 I=1,NN
  125    P(I) = QP(I)
         IF (NZ.EQ.1) GO TO 20
         Z(J+1) = RLZR
         Z(JJ+1) = RLZI
         GO TO 20
C                                  IF THE ITERATION IS UNSUCCESSFUL
C                                    ANOTHER QUADRATIC IS CHOSEN AFTER
C                                    RESTORING K
  130    DO 135 I=1,N
  135    RK(I) = TEMP(I)
  140 CONTINUE
C                                  RETURN WITH FAILURE IF NO
C                                    CONVERGENCE WITH 20 SHIFTS
      IER = 131
C                                  CONVERT ZEROS (Z) IN COMPLEX FORM
  145 DO 150 I=1,NDEG
         NPI= NDEG+I
         P(I) = Z(NPI)
  150 CONTINUE
      N2 = NDEG+NDEG
      J = NDEG
      DO 155 I=1,NDEG
         Z(N2-1) = Z(J)
         Z(N2) = P(J)
         N2 = N2-2
         J = J-1
  155 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
C                                  SET UNFOUND ROOTS TO MACHINE INFINITY
      N2 = 2*(NDEG-NN)+3
      DO 160 I=1,N
         Z(N2) = RINFP
         Z(N2+1) = RINFP
         N2 = N2+2
  160 CONTINUE
      GO TO 9000
  165 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HZRPOLY)
 9005 RETURN
      END

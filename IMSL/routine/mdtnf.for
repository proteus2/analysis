C   IMSL ROUTINE NAME   - MDTNF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - INTEGRAL RELATED TO CALCULATION OF NON-
C                           CENTRAL T AND BIVARIATE NORMAL PROBABILITY
C                           DISTRIBUTION FUNCTIONS
C
C   USAGE               - CALL MDTNF (Y,Z,EPS,T)
C
C   ARGUMENTS    Y      - INPUT PARAMETER.  SEE REMARKS.
C                Z      - INPUT.  INTEGRATION IS FROM 0 TO Z.
C                EPS    - INPUT.  ACCURACY SHOULD NOT BE LESS THAN EPS.
C                           IF EPS=0.0 IS ENTERED, EPS=.000001 IS USED.
C                T      - OUTPUT VALUE OF THE INTEGRAL.
C
C   REQD. IMSL ROUTINES - MDNOR,MERRC=ERFC
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      MDTNF COMPUTES THE FUNCTION T(Y,Z) WHERE
C                T(Y,Z) = THE INTEGRAL, FROM 0 TO Z, OF
C                (EXP((-Y**2/2)*(1+X**2))/(2*PI*(1+X**2)))DX
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDTNF  (Y,Z,EPS,T)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               Y,Z,EPS,T
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               C,EXPOV,EP1,B,A,TA,HSQB,BEXP,ASQ,A4,B4,A4B4,
     *                   AHSQB,AB4,F,SUM,G,G1,BER,TER,D1,D2,D,AEPS
      DATA               C/.1591549/,EXPOV/174.673/
C                                  FIRST EXECUTABLE STATEMENT
      EP1 = EPS
      IF(EPS .EQ. 0.) EP1 = .000001
      T = 0.0
      B = ABS(Y)
      A = ABS(Z)
      IF(A .EQ. 0.) GO TO 35
    5 TA = ATAN(A)
      IF (A*B .LE. 4.0) GO TO 10
C                                  APPROXIMATION FOR SMALL Y*Z
      CALL MDNOR(B,T)
      T = C*(TA+ATAN(1.0/A)) - .5*(T-.5)
      GO TO 30
   10 HSQB = .5*B*B
      IF (HSQB .GT. EXPOV) GO TO 35
      BEXP = EXP(-HSQB)
      ASQ = A*A
      A4 = ASQ*ASQ
      B4 = HSQB * HSQB
      A4B4 = A4 * B4
      AHSQB = A * HSQB
      AB4 = A*B4*.5
      F = 1.0
      SUM = 0.0
      G = 3.0
C                                  BEGIN SERIES EXPANSION
   15 G1 = G
      BER = 0.0
      TER = AB4
   20 BER = BER+TER
      IF(TER .LE. BER*EP1) GO TO 25
C                                  DEVELOP COEFFICIENT SERIES
      TER = TER*(HSQB/G1)
      G1 = G1+1.0
      GO TO 20
   25 D1 = (BER+AHSQB)/F
      D2 = BER*ASQ/(F+2.0)
      D = D1-D2
      SUM = SUM+D
      T = TA-SUM*BEXP
      AEPS = EP1*T
      AHSQB = AHSQB*A4B4/((G-1.0)*G)
      AB4 = AB4*A4B4/((G +1.0)*G)
      F = F+4.0
      G = G+2.0
C                                  SHOULD SERIES EXPANSION BE TERMINATED
      IF (D2*BEXP .GE. AEPS) GO TO 15
      T = T * C
   30 IF (Z .LT. 0.0) T = -T
   35 RETURN
      END

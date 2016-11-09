C   IMSL ROUTINE NAME   - MMBSJ0
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - BESSEL FUNCTION OF THE FIRST KIND OF ORDER
C                           ZERO
C
C   USAGE               - FUNCTION MMBSJ0 (ARG,IER)
C
C   ARGUMENTS    MMBSJ0 - OUTPUT VALUE OF THE FUNCTION AT ARG. MMBSJ0
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.)
C                ARG    - INPUT ARGUMENT. THE ABSOLUTE VALUE OF ARG MUST
C                           BE LESS THAN OR EQUAL TO XMAX, WHICH IS OF
C                           THE ORDER OF 10**8. THE EXACT VALUE OF XMAX
C                           MAY ALLOW LARGER RANGES FOR ARG ON SOME
C                           COMPUTERS. SEE THE PROGRAMMING NOTES IN THE
C                           MANUAL FOR THE EXACT VALUES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ABSOLUTE VALUE
C                             OF ARG IS GREATER THAN XMAX. MMBSJ0 IS
C                             SET TO ZERO.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMBSJ0 (ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   AX,DOWN,FUDGE,FUDGEX,PI2,PROD,
     1                   PROD1,PRP,PRQ,P0,P1,QRP,QRQ,Q0,Q1,R0,R1,TWOPI1,
     2                   TWOPI2,U,UP,W,X,XDEN,XMAX,XMAX1,XNUM,XSMALL,
     3                   X01,X01P,X02,X11,X12,Z,ZSQ
      DOUBLE PRECISION   DABS,DCOS,DSIN,DSQRT
      DIMENSION          PRP(7),PRQ(8),P0(6),P1(6),QRP(5),QRQ(7),Q0(5),
     1                   Q1(5)
C                                  MACHINE DEPENDENT CONSTANTS
C                                     FUDGE = 7 * 16**(-15)
C                                     FUDGEX = 2 * 16**(-15)
C                                     PI2 = 2 / PI
C                                     TWOPI1 + TWOPI2 = 2*PI TO EXTRA
C                                     PRECISION
C                                     XMAX = 16**9, LARGEST ACCEPTABLE
C                                     ARGUMENT
C                                     XMAX1 = SMALLEST FLOATING-POINT
C                                     CONSTANT WITH ENTIRELY INTEGER
C                                     REPRESENTATION
C                                     XSMALL = 16**(-10), ARGUMENT
C                                     BELOW WHICH J0 MAY BE REPRESENTED
C                                     BY ONE TERM IN THE ASCENDING
C                                     SERIES
C                                     X01 + X02 = FIRST ZERO OF J-SUB-0
C                                     TO EXTRA PRECISION
C                                     X11 + X12 = SECOND ZERO OF
C                                     J-SUB-0 TO EXTRA PRECISION
C
      DATA PI2/Z40A2F9836E4E4415/,X01P/Z414CF454BA5C6CFF/,
     1     XMAX/Z4A10000000000000/,XSMALL/Z3710000000000000/,
     2     X01/Z41267A2A5D2E367F/,X02/Z33785631412ED80C/,
     3     X11/Z4158523D6CB0B914/,X12/Z335D415335829085/,
     4     FUDGE/Z3F00000000000007/,FUDGEX/Z3F00000000000002/,
     5     XMAX1/Z4E10000000000000/,TWOPI1/Z416487ED00000000/,
     6     TWOPI2/Z3B5110B4611A6263/
C
C                                  COEFFICIENTS FOR RATIONAL
C                                    APPROXIMATION OF J-0(X) / (X**2 -
C                                    X0**2), XSMALL .LT. ABS(X) .LE.
C                                    4.0
C
      DATA PRP/Z46CA573794BA3CCD, ZC84A13D638788682, Z49CB45410CD363E2,
     1         ZCAC04FD96CA3D13C, ZC03E0A0AA9070ECC, Z42CEE26BBA894F91,
     2         ZC511E2BA11A7964F/
      DATA QRP/Z461BD6FBD95604E5, Z482157D376E3FA41, Z4A188521B2737048,
     1         Z4B8B059E13844EA9, Z43EA08F96DC34160/
C
C                                  COEFFICIENTS FOR RATIONAL
C                                    APPROXIMATION OF J-0(X) / (X**2 -
C                                    X1**2), 4.0 .LT. ABS(X) .LE. 8.0
C
      DATA PRQ/Z44228357665284F0, Z445B9A17B1FAF674, Z4450CBD23F9D18F8,
     1         ZC438EFF0C424693A, ZC45FBC2804C4AD24, ZC3E4FE12D107D8FE,
     2         Z42612EF3BE05308B, Z435CE6C86E9C94C4/
      DATA QRQ/Z4329A256E60E8861, ZC41703C0D8DD844D, Z4491F1FAE2602AF6,
     1         ZC5290AE200BF3249, Z45781CE0B959DFC0, ZC5AEB958F7E70373,
     2         ZC232842291AC3C0F/
C                                  COEFFICIENTS FOR HART
C                                    APPROXIMATION, ABS(X) .GT. 8.0
C
      DATA P0/Z43D98A60D8DF24C3, Z4452B285FC49D24A, Z44A18162FACFC2BA,
     1        Z4458FB17172BA9CE, Z40E3BDD722D5C557, Z4299C313AEDED625/
      DATA Q0/Z43DAEDF9E9A18CC3, Z4452DF59BE6C6252, Z44A19A69994E46FC,
     1        Z4458FB17172BA9CE, Z429D1C91B97EDABA/
      DATA P1/ZC2164CDDF2D4A1E8, ZC26FD594A1F7CFF9, ZC2B9EB66BD880301,
     1        ZC2593A027883B414, ZBF240EF5FC15AF42, ZC113E7D834C33D98/
      DATA Q1/Z435D0B91E9A98D9F, Z441C60472C1DEFDE, Z442EAF21ACD49DAB,
     1        Z44164E809E20ED05, Z425A980148BFDACA/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      AX = DABS(ARG)
      IF (AX.GT.XMAX) GO TO 35
      IF (AX.GT.XSMALL) GO TO 5
      MMBSJ0 = 1.0D0
      GO TO 9005
    5 IF (AX.GT.8.0D0) GO TO 25
      IF (AX.GT.4.0D0) GO TO 15
C                                  XSMALL .LT. ABS(X) .LE. 4.0,
C                                    CALCULATION SCALED TO AVOID LOSS
C                                    OF ACCURACY ASSOCIATED WITH
C                                    HEXADECIMAL ARITHMETIC
      ZSQ = AX*AX
      XNUM = (PRP(5)*ZSQ+PRP(6))*ZSQ+PRP(7)
      XDEN = 4.0D0*ZSQ+QRP(5)
      DO 10 I=1,4
         XNUM = XNUM*ZSQ+PRP(I)
         XDEN = XDEN*ZSQ+QRP(I)
   10 CONTINUE
C                                  CALCULATION TO PRESERVE ACCURACY
C                                    NEAR THE FIRST ZERO OF J-0
      PROD = (AX-X01)-X02
      PROD1 = AX+AX+X01P
      R0 = (PROD1*XNUM)/XDEN
      R0 = R0-FUDGE*PROD1
      MMBSJ0 = PROD*R0
      GO TO 9005
C                                  4.0 .LT. ABS(X) .LE. 8.0,
C                                    CALCULATION SCALED TO AVOID LOSS
C                                    OF ACCURACY ASSOCIATED WITH
C                                    HEXADECIMAL ARITHMETIC
   15 ZSQ = 1.0D0-(AX*AX)/64.0D0
      XNUM = PRQ(7)*ZSQ+PRQ(8)
      XDEN = 2.0D0*ZSQ+QRQ(7)
      DO 20 I=1,6
         XNUM = XNUM*ZSQ+PRQ(I)
         XDEN = XDEN*ZSQ+QRQ(I)
   20 CONTINUE
      R0 = XNUM/XDEN
      R0 = R0+FUDGEX
C                                  CALCULATION TO PRESERVE ACCURACY
C                                    NEAR THE SECOND ZERO OF J-0
      PROD = (AX+X11)
      PROD = (AX-X11)*PROD-(PROD*X12)
      MMBSJ0 = R0*PROD
      GO TO 9005
C                                  ABS(X) .GT. 8.0
   25 Z = 8.0D0/AX
      W = AX/TWOPI1
      W = ((W+XMAX1)-XMAX1)+0.125D0
      U = (AX-W*TWOPI1)-W*TWOPI2
      ZSQ = Z*Z
      XNUM = P0(5)*ZSQ+P0(6)
      XDEN = ZSQ+Q0(5)
      UP = P1(5)*ZSQ+P1(6)
      DOWN = ZSQ+Q1(5)
      DO 30 I=1,4
         XNUM = XNUM*ZSQ+P0(I)
         XDEN = XDEN*ZSQ+Q0(I)
         UP = UP*ZSQ+P1(I)
         DOWN = DOWN*ZSQ+Q1(I)
   30 CONTINUE
      R0 = XNUM/XDEN
      R1 = UP/DOWN
      MMBSJ0 = DSQRT(PI2/AX)*(R0*DCOS(U)-Z*R1*DSIN(U))
      GO TO 9005
C                                  ERROR RETURN FOR ABS(X) .GT. XMAX
   35 MMBSJ0 = 0.0D0
      IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMMBSJ0)
 9005 RETURN
      END

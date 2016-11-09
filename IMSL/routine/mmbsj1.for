C   IMSL ROUTINE NAME   - MMBSJ1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - BESSEL FUNCTION OF THE FIRST KIND OF ORDER
C                           ONE
C
C   USAGE               - FUNCTION MMBSJ1 (ARG,IER)
C
C   ARGUMENTS    MMBSJ1 - OUTPUT VALUE OF THE FUNCTION AT ARG. MMBSJ1
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.)
C                ARG    - INPUT ARGUMENT. THE ABSOLUTE VALUE OF ARG
C                           MUST BE LESS THAN OR EQUAL TO XMAX, WHICH
C                           IS OF THE ORDER OF 10**8. THE EXACT VALUE OF
C                           XMAX MAY ALLOW LARGER RANGES FOR ARG ON SOME
C                           COMPUTERS. SEE THE PROGRAMMING NOTES IN THE
C                           MANUAL FOR THE EXACT VALUES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ABSOLUTE VALUE
C                             OF ARG IS GREATER THAN XMAX. MMBSJ1 IS
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
      DOUBLE PRECISION FUNCTION MMBSJ1 (ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   AP1,AP2,AQ1,AQ2,AX,B,D,DEN,DEN2,
     1                   P,Q,RTPI2,TWOPI1,TWOPI2,V,X,XC,XMAX,XMAX1,
     2                   XMIN,XNUM,XNUM2,XSMALL,X01,X02,X11,X12,Y,Z,ZSQ
      DOUBLE PRECISION   DABS,DCOS,DSIN,DSQRT
      DIMENSION          P(6),Q(5),B(8),D(6),AP1(6),AQ1(5),AP2(6),
     1                   AQ2(5)
C                                  MACHINE DEPENDENT CONSTANTS
C                                     RTPI2 = SQRT(2/PI)
C                                     TWOPI1 + TWOPI2 = 2*PI TO EXTRA
C                                     PRECISION
C                                     XMAX = 16**9, LARGEST ACCEPTABLE
C                                     ARGUMENT
C                                     XMAX1 = SMALLEST FLOATING-POINT
C                                     CONSTANT WITH ENTIRELY INTEGER
C                                     REPRESENTATION
C                                     XMIN = 2*16**(-65), ARGUMENT
C                                     BELOW WHICH J1 MAY BE REPRESENTED
C                                     BY ONE TERM IN THE ASCENDING
C                                     SERIES
C                                     X01 + X02 = FIRST ZERO OF J-SUB-1
C                                     TO EXTRA PRECISION
C                                     X11 + X12 = SECOND ZERO OF
C                                     J-SUB-1 TO EXTRA PRECISION
C
      DATA RTPI2/Z40CC42299EA1B284/,XMAX1/Z4E10000000000000/,
     1     TWOPI1/Z416487ED00000000/,TWOPI2/Z3B5110B4611A6263/,
     2     XMAX/Z4A10000000000000/,XMIN/Z0020000000000000/,
     3     X01/Z413D4EAAEB5EDE14/,X02/ZB42B00AAD4E8D905/,
     4     X11/Z41703FD7CED1C942/,X12/ZB426C89B67490DC3/,
     5     XSMALL/Z3710000000000000/
C
C                                  COEFFICIENTS FOR RATIONAL
C                                     APPROXIMATION OF J-1(ARG) / (ARG*
C                                     (ARG**2 - X0**2)), XSMALL .LT.
C                                     ABS(ARG) .LE. 4.0
C
      DATA P/Z4536DD559B0B7C0C,ZC71D043A883CE853,Z4871DE158C64455A,
     1       ZC99F464358B426A7,ZC331E92085C3C546,Z41121DAB242E9189/
      DATA Q/Z45521EA694F69585,Z4756D6B7BF375C36,Z49391785AE89FEA2,
     1       Z4B1244ED1117566A,Z4331AB2176B94DAB/
C                                  COEFFICIENTS FOR RATIONAL
C                                     APPROXIMATION OF J-1(ARG) / (ARG*
C                                     (ARG**2 - X1**2)), 4.0 .LT.
C                                     ABS(ARG) .LE. 8.0 NUMERATOR IN
C                                     MINI-NEWTON FORM
C
      DATA B/Z3F16D4006B0D22BD,ZC180DE0853CAFA1B,Z4412D36A5DC52FD4,
     1       ZC61695D15E8962F9,Z47DED253A2219B17,ZC93AB436A6806296,
     2       Z4A188AAE4A7E0460,Z4ADF87033A5B9DA8/
      DATA D/Z43278F21F064AB34,Z456E52DAAD0A7BF8,Z47D0D51111DA72EA,
     1       Z4A1115950FE9E4F1,Z4BE5172B5C19336D,Z4D5E5B25E8D16D8C/
C
C                                  COEFFICIENTS FOR HART
C                                     APPROXIMATION, ABS(ARG) .GT. 8.0
C
      DATA AP1/Z42D32725A71ABCBC,Z4413797BB3646E7F,Z447A79F68E6FDE75,
     1         Z44F526D8621E32BD,Z448998AA37C33226,Z41141D6010A865E8/
      DATA AQ1/Z42CB13D8478E3359,Z44134265806168F6,Z447A086808A4FBD0,
     1         Z44F4E658D2540F48,Z448998AA37C33225/
      DATA AP2/Z4149173B22770030,Z4253309D024E9158,Z431A9FCBFC4E8E2C,
     1        Z432D10A07BB8D1D0,Z4315F2CD95877EA5,Z3F90B4834E29AFD8/
      DATA AQ2/Z4267D19A2C564CE7,Z437132FCB9B7E7E8,Z4423C03B50CA65EA,
     1         Z443C362D661B27E3,Z441D43BCC75F5387/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      AX = DABS(ARG)
      IF (AX.GT.XMAX) GO TO 40
      IF (AX.GE.XMIN) GO TO 5
      MMBSJ1 = 0.0D0
      GO TO 9005
    5 IF (AX.GT.XSMALL) GO TO 10
      MMBSJ1 = ARG/2.0D0
      GO TO 9005
   10 IF (AX.GT.8.0D0) GO TO 30
      Y = AX*AX
      IF (AX.GT.4.0D0) GO TO 20
C                                  XSMALL .LT. ABS(ARG) .LE. 4.0
      XNUM = P(6)*Y+P(5)
      DEN = Y+Q(5)
      DO 15 I=1,4
         XNUM = XNUM*Y+P(I)
         DEN = DEN*Y+Q(I)
   15 CONTINUE
      Z = (AX-X01)-X02
      MMBSJ1 = (XNUM / DEN) * ARG * Z * (AX + X01)
      GO TO 9005
C                                  4.0 .LT. ABS(ARG) .LE. 8.0
   20 XNUM = 0.0D0
      DEN = 0.5D0
      DO 25 I=1,6
         XNUM = XNUM*Y+B(I)
         DEN = DEN*Y+D(I)
   25 CONTINUE
      XNUM = XNUM*(AX-8.0D0)*(AX+8.0D0)+B(7)
      XNUM = XNUM*(AX-4.0D0)*(AX+4.0D0)+B(8)
      Z = (AX-X11)-X12
      MMBSJ1 = (XNUM/DEN)*ARG*Z*(AX+X11)
      GO TO 9005
C                                  ABS(ARG) .GT. 8.0
   30 XC = RTPI2/DSQRT(AX)
      IF (ARG.LT.0.0D0) XC = -XC
      Z = 8.0D0/AX
      ZSQ = Z*Z
      V = ((AX/TWOPI1+XMAX1)-XMAX1)+0.375D0
      V = (AX-V*TWOPI1)-V*TWOPI2
      XNUM = AP1(6)
      DEN = 1.0D0
      XNUM2 = AP2(6)
      DEN2 = 1.0D0
      DO 35 I=1,5
         XNUM = XNUM*ZSQ+AP1(I)
         DEN = DEN*ZSQ+AQ1(I)
         XNUM2 = XNUM2*ZSQ+AP2(I)
         DEN2 = DEN2*ZSQ+AQ2(I)
   35 CONTINUE
      MMBSJ1 = XC*((XNUM/DEN)*DCOS(V)-Z*(XNUM2/DEN2)*DSIN(V))
      GO TO 9005
   40 MMBSJ1 = 0.0D0
      IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMMBSJ1)
 9005 RETURN
      END

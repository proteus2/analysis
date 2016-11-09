C   IMSL ROUTINE NAME   - MMBSI1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MODIFIED BESSEL FUNCTION OF THE FIRST KIND OF
C                           ORDER ONE
C
C   USAGE               - FUNCTION MMBSI1 (IOPT,ARG,IER)
C
C   ARGUMENTS    MMBSI1 - OUTPUT VALUE OF THE FUNCTION AT ARG. MMBSI1
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.)
C                IOPT   - INPUT OPTION PARAMETER.
C                         IF IOPT = 1, THE MODIFIED BESSEL FUNCTION OF
C                           THE FIRST KIND OF ORDER ONE FOR ARGUMENT
C                           ARG IS EVALUATED. THE ABSOLUTE VALUE OF ARG
C                           MUST BE LESS THAN OR EQUAL TO XMAX, WHICH
C                           IS AT LEAST 91. THE EXACT VALUE OF XMAX
C                           MAY ALLOW LARGER RANGES FOR ARG ON SOME
C                           COMPUTERS. SEE THE PROGRAMMING NOTES
C                           IN THE MANUAL FOR THE EXACT VALUES.
C                         IF IOPT = 2, EXP(-ABS(ARG))*THE MODIFIED
C                           BESSEL FUNCTION OF THE FIRST KIND OF ORDER
C                           ONE FOR ARGUMENT ARG IS EVALUATED.
C                ARG    - INPUT ARGUMENT. SEE DESCRIPTION OF IOPT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IOPT WAS NOT 1 OR
C                             2. MMBSI1 IS SET TO MACHINE INFINITY.
C                           IER = 130 INDICATES THAT IOPT = 1 AND THE
C                             ABSOLUTE VALUE OF ARG IS OUT OF RANGE.
C                             MMBSI1 IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMBSI1 (IOPT,ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   A,B,CONST,DEXP40,P(15),PBAR,PP(8),XCON,
     1                   Q(5),QQ(6),SUMP,SUMQ,X,XINF,XMAX,XMIN,
     1                   XSMALL,XX
      DOUBLE PRECISION   DABS,DEXP,DSQRT
C
C                                  MACHINE DEPENDENT CONSTANTS
C                                     CONST = 1/15
C                                     DEXP40 = DEXP(40)
C                                     XINF = LARGEST POSITIVE MACHINE
C                                     NUMBER
C                                     XMAX = 178.186, LARGEST ARGUMENT
C                                     ACCEPTABLE TO MMBSI1, IOPT=1
C                                     XMIN = 2 * 16**(-65), TWICE THE
C                                     SMALLEST POSITIVE MACHINE NUMBER
C                                     XSMALL = 16**(-14), ARGUMENT
C                                     BELOW WHICH MMBSI1 MAY BE
C                                     REPRESENTED BY ONE TERM IN THE
C                                     ASCENDING SERIES.
C
      DATA XCON/170.0D0/
      DATA XMAX/Z42B22F8049182478/,DEXP40/Z4F34441A72F2E5D5/
      DATA XSMALL/Z3310000000000000/
      DATA XINF/Z7FFFFFFFFFFFFFFF/,CONST/Z4011111111111111/
      DATA XMIN/Z0020000000000000/
C                                  COEFFICIENTS FOR XSMALL .LE.
C                                    ABS(ARG) .LT. 15.0
C
      DATA P/ZB13A28E5C4538F90,ZB42F03AC5C068141,ZB714FC3FB99F5D81,
     1       ZB965EC90277FCFF6,ZBC1697F59935C261,ZBE3C207FBC053D52,
     2       ZC078D9A3986D0D81,ZC2B6426C62701866,ZC4CAB6178C270CBF,
     3       ZC6A19186B98A7BE3,ZC858622408A946CD,ZCA1F19A73D76C143,
     4       ZCB65AF2101A01025,ZCCA1459E2A9BEBF0,ZCD52DC96D3B05CA6/
      DATA Q/ZC3FA7AFBC5DC099B,Z467226E209215FEA,ZC91DD31219BDCAA8,
     1       Z4B46A450B7F235C0,ZCD4B22F76A767700/
C                                  COEFFICIENTS FOR 15.0 .LE. ABS(ARG)
C
      DATA PP/ZBFF78CF459BF88B1,Z40751D7D652C037D,ZC06DAE1750DF14F2,
     1        Z4018EC52A4B7A96E,ZBED4B7064C1B5C14,ZBE17DA1CB77839E0,
     2        Z3D110C66BC67BFDC,ZBB6191DD9D969001/
      DATA QQ/ZC13E172D8C3D8EB7,Z41342662B7E5B310,ZC0D9A50DA9FACE19,
     1        Z4012FF8EEF69BCA1,ZBE95A7D2095F5300,Z3D27552155E95C8E/
      DATA PBAR/Z4066000000000000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IOPT.NE.1.AND.IOPT.NE.2) GO TO 40
      X = DABS(ARG)
      IF (X.LT.XSMALL) GO TO 30
      IF ((IOPT.EQ.1).AND.(X.GT.XMAX)) GO TO 35
      IF (X.GE.15.0D0) GO TO 10
C                                  XSMALL .LE. ABS(ARG) .LT. 15.0
      XX = X*X
      SUMP = P(1)
      DO 5 J=2,15
         SUMP = SUMP*XX+P(J)
    5 CONTINUE
      XX = XX-225.0D0
      SUMQ = ((((XX+Q(1))*XX+Q(2))*XX+Q(3))*XX+Q(4))*XX+Q(5)
      MMBSI1 = (SUMP/SUMQ)*X
      IF (IOPT.EQ.2) MMBSI1 = MMBSI1*DEXP(-X)
      GO TO 45
C                                  15.0 .LE. ABS(ARG)
   10 XX = 1.0D0/X-CONST
      SUMP = ((((((PP(1)*XX+PP(2))*XX+PP(3))*XX+PP(4))*XX+PP(5))*XX+
     *PP(6))*XX+PP(7))*XX+PP(8)
      SUMQ = (((((XX+QQ(1))*XX+QQ(2))*XX+QQ(3))*XX+QQ(4))*XX+QQ(5))*XX+
     *QQ(6)
      MMBSI1 = SUMP/SUMQ
      IF (IOPT.EQ.1) GO TO 15
      MMBSI1 = (MMBSI1+PBAR)/DSQRT(X)
      GO TO 45
C                                  CALCULATION REFORMULATED TO PRESERVE
C                                    ACCURACY ON IBM EQUIPMENT AND TO
C                                    AVOID PREMATURE OVERFLOW
   15 IF (X.GT.XCON) GO TO 20
      A = DEXP(X)
      B = 1.0D0
      GO TO 25
   20 A = DEXP(X-40.0D0)
      B = DEXP40
   25 MMBSI1 = ((MMBSI1*A+PBAR*A)/DSQRT(X))*B
      GO TO 45
C                                  RETURN FOR ABS(ARG) .LT. XSMALL
   30 IF (X.LT.XMIN) X = 0.0D0
      MMBSI1 = 0.5D0*X
      GO TO 45
C                                  ERROR RETURN FOR ABS(ARG) .GT. XMAX
   35 MMBSI1 = XINF
      IER = 130
      GO TO 45
   40 IER = 129
      MMBSI1 = XINF
   45 IF (ARG.LT.0.0D0) MMBSI1 = -MMBSI1
 9000 CONTINUE
      IF (IER.GT.0) CALL UERTST (IER,6HMMBSI1)
 9005 RETURN
      END

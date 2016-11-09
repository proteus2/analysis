C   IMSL ROUTINE NAME   - MMBSI0
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MODIFIED BESSEL FUNCTION OF THE FIRST KIND OF
C                           ORDER ZERO
C
C   USAGE               - FUNCTION MMBSI0 (IOPT,ARG,IER)
C
C   ARGUMENTS    MMBSI0 - OUTPUT VALUE OF THE FUNCTION AT ARG. MMBSI0
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.)
C                IOPT   - INPUT OPTION PARAMETER.
C                         IF IOPT = 1, THE MODIFIED BESSEL FUNCTION OF
C                           THE FIRST KIND OF ORDER ZERO FOR ARGUMENT
C                           ARG IS EVALUATED. THE ABSOLUTE VALUE OF ARG
C                           MUST BE LESS THAN OR EQUAL TO XMAX, WHICH
C                           IS AT LEAST 91. THE EXACT VALUE OF XMAX
C                           MAY ALLOW LARGER RANGES FOR ARG ON SOME
C                           COMPUTERS. SEE THE PROGRAMMING NOTES
C                           IN THE MANUAL FOR THE EXACT VALUES.
C                         IF IOPT = 2, EXP(-ABS(ARG))*THE MODIFIED
C                           BESSEL FUNCTION OF THE FIRST KIND OF ORDER
C                           ZERO FOR ARGUMENT ARG IS EVALUATED.
C                ARG    - INPUT ARGUMENT. SEE THE DESCRIPTION OF IOPT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IOPT WAS NOT 1 OR
C                             2. MMBSI0 IS SET TO MACHINE INFINITY.
C                           IER = 130 INDICATES THAT IOPT = 1 AND THE
C                             ABSOLUTE VALUE OF ARG IS OUT OF RANGE.
C                             MMBSI0 IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMBSI0 (IOPT,ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   A,B,CONST,DEXP40,P(15),PP(8),Q(5),
     1                   QQ(7),SUMP,SUMQ,X,XINF,XMAX,XSMALL,XX,XCON
      DOUBLE PRECISION   DABS,DEXP,DSQRT
C
C                                  MACHINE DEPENDENT CONSTANTS.
C                                     CONST = 1/15
C                                     DEXP40 = DEXP(40)
C                                     XINF = LARGEST POSITIVE MACHINE
C                                     NUMBER.
C                                     XMAX = 178.183, LARGEST ARGUMENT
C                                     ACCEPTABLE TO MMBSI0 WITH IOPT = 1
C                                     XSMALL = 16**(-14), ARGUMENT
C                                     BELOW WHICH MMBSI0 MAY BE
C                                     REPRESENTED BY ONE TERM IN
C                                     THE ASCENDING SERIES.
C
      DATA XMAX/Z42B22EC75954DE83/,DEXP40/Z4F34441A72F2E5D5/
      DATA XSMALL/Z3310000000000000/
      DATA XINF/Z7FFFFFFFFFFFFFFF/,CONST/Z4011111111111111/
C
C                                  COEFFICIENTS FOR XSMALL .LE.
C                                    ABS(ARG) .LT. 15.0
C
      DATA  P/ZB260D2AD564A3548,ZB547FA3E12BC74C2,ZB81D83C15EDB4C04,
     1        ZBA831207CE0A6915,ZBD1A653A43F35559,ZBF3F4CA9FEFE99CA,
     2        ZC1717F1B9BAB101B,ZC396C84CD375CCB3,ZC5916000FF791BE5,
     3        ZC7625A609CE8D7C4,ZC92C60E05CA9639D,ZCAC5BB4A0A3E8A4A,
     4        ZCC1DF5759EA4CD7F,ZCD1F4AE174B11F8F,ZCD7EF68A66DE76FA/
      DATA  Q/ZC3E8FC1718928CBE,Z46636C89756DA620,ZC91872A28FEA73B0,
     1        Z4B36B8A5D3E90971,ZCD3730269B3C3389/
C                                  COEFFICIENTS FOR 15.0 .LE. ABS(ARG)
C
      DATA PP/ZC066000000000003,Z412EBA8610DC8727,ZC1278895D7960507,
     1        Z407AA97ABC731DEF,ZBEF501441AC695A6,ZBEAFA571D0F30F33,
     2        Z3D67FC6278DA54F5,ZBC24B428E9FB4034/
      DATA QQ/ZC21F725A1EC2FD0C,Z42558A1FAB2B4E55,ZC23C3A5CDCBFF7E9,
     1        Z41DFB8ACD3C904B4,ZC111D7C19FFA40B9,Z3F8550B52893D8EF,
     2        ZBE242C13A56DB861/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      XCON = 170.0D0
      X = DABS(ARG)
      IF (X.LT.XSMALL) GO TO 30
      IF (IOPT.NE.1.AND.IOPT.NE.2) GO TO 40
      IF ((IOPT.EQ.1).AND.(X.GT.XMAX)) GO TO 35
      IF (X.GE.15.0D0) GO TO 10
C                                  XSMALL .LE. ABS(ARG) .LT. 15.0
      XX = X*X
      SUMP = P(1)
      DO 5 I=2,15
         SUMP = SUMP*XX+P(I)
    5 CONTINUE
      XX = XX-225.0D0
      SUMQ = ((((XX+Q(1))*XX+Q(2))*XX+Q(3))*XX+Q(4))*XX+Q(5)
      MMBSI0 = SUMP/SUMQ
      IF (IOPT.EQ.2) MMBSI0 = MMBSI0*DEXP(-X)
      GO TO 9005
C                                  15.0 .LE. ABS(ARG)
   10 XX = 1.0D0/X-CONST
      SUMP = ((((((PP(1)*XX+PP(2))*XX+PP(3))*XX+PP(4))*XX+PP(5))*XX+
     *PP(6))*XX+PP(7))*XX+PP(8)
      SUMQ = ((((((XX+QQ(1))*XX+QQ(2))*XX+QQ(3))*XX+QQ(4))*XX+QQ(5))*XX+
     *QQ(6))*XX+QQ(7)
      MMBSI0 = SUMP/SUMQ
      IF (IOPT.EQ.1) GO TO 15
      MMBSI0 = (MMBSI0-PP(1))/DSQRT(X)
      GO TO 9005
C                                  CALCULATION REFORMULATED TO PRESERVE
C                                    ACCURACY ON IBM EQUIPMENT AND TO
C                                    AVOID PREMATURE OVERFLOW
   15 IF (X.GT.XCON) GO TO 20
      A = DEXP(X)
      B = 1.0D0
      GO TO 25
   20 A = DEXP(X-40.0D0)
      B = DEXP40
   25 MMBSI0 = ((MMBSI0*A-PP(1)*A)/DSQRT(X))*B
      GO TO 9005
C                                  RETURN FOR ABS(ARG) .LT. XSMALL
   30 MMBSI0 = 1.0D0
      GO TO 9005
C                                  ERROR RETURN FOR ABS(ARG) .GT. XMAX
   35 MMBSI0 = XINF
      IER = 130
      GO TO 9000
   40 MMBSI0 = XINF
      IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMMBSI0)
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - MMBSK0
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MODIFIED BESSEL FUNCTION OF THE SECOND KIND
C                           OF ORDER ZERO
C
C   USAGE               - FUNCTION MMBSK0 (IOPT,ARG,IER)
C
C   ARGUMENTS    MMBSK0 - OUTPUT VALUE OF THE FUNCTION AT ARG. MMBSK0
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.)
C                IOPT   - INPUT OPTION PARAMETER.
C                         IF IOPT = 1, THE MODIFIED BESSEL FUNCTION OF
C                           THE SECOND KIND OF ORDER ZERO FOR ARGUMENT
C                           ARG IS EVALUATED. ARG MUST BE GREATER THAN
C                           XMIN AND LESS THAN OR EQUAL TO XMAX. XMIN IS
C                           OF THE ORDER OF 10**(-46) AND XMAX IS AT
C                           LEAST 87. THE EXACT VALUES OF XMIN AND XMAX
C                           MAY ALLOW LARGER RANGES FOR ARG ON SOME
C                           COMPUTERS. SEE THE PROGRAMMING NOTES IN THE
C                           MANUAL FOR THE EXACT VALUES.
C                         IF IOPT = 2, EXP(ARG)*THE MODIFIED BESSEL
C                           FUNCTION OF THE SECOND KIND OF ORDER ZERO
C                           FOR ARGUMENT ARG IS EVALUATED. ARG MUST BE
C                           GREATER THAN XMIN.
C                ARG    - INPUT ARGUMENT. SEE DESCRIPTION OF IOPT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ARG WAS LESS THAN
C                             XMIN OR IOPT WAS NOT 1 OR 2. MMBSK0 IS
C                             SET TO MACHINE INFINITY.
C                           IER = 130 INDICATES THAT ARG WAS OUT OF
C                             RANGE. MMBSK0 IS SET TO MACHINE INFINITY.
C                             HOWEVER, IF IOPT = 1 AND ARG IS GREATER
C                             THAN OR EQUAL TO XMAX, MMBSK0 IS SET
C                             TO ZERO.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMBSK0 (IOPT,ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      DOUBLE PRECISION   P(5),Q(3),PP(9),QQ(9),F(4),G(3),A,X,XX,
     *                   SUMF,SUMG,SUMP,SUMQ,TEMP,XINF,XMAX,CONST,
     *                   FUDGE,ONE,ZERO
      DATA               ONE/1.D0/,ZERO/0.D0/
      DATA               XINF/Z7FFFFFFFFFFFFFFF/,
     *                   FUDGE/Z4100000000000001/,
     *                   CONST/ZC01DADB014541EB2/,
     *                   XMAX/Z42B1DA5648BC5844/
C
C                                  COEFFICIENTS FOR ARG GREATER THAN 0
C                                  AND LESS THAN OR EQUAL TO 1.
C
      DATA               P(1)/ZC1345B4243F76BD1/,P(2)/ZC32018FE2BAD589D/
      DATA               P(3)/ZC4637EBF32BACD74/,P(4)/ZC556FB0F78607C57/
      DATA               P(5)/ZC52471589E70A64E/,Q(1)/ZC2E4AEFF10611FC4/
      DATA               Q(2)/Z4462EDB7CFFBE9BC/,Q(3)/ZC613A58E863D56DA/
      DATA               F(1)/ZC11A435C213CB1B9/,F(2)/ZC3128043E8448629/
      DATA               F(3)/ZC44545C8E11CF2BD/,F(4)/ZC5627036859396BF/
      DATA               G(1)/ZC2FAA6545795DF00/,G(2)/Z4474A9B691DA99C1/
      DATA               G(3)/ZC6189C0DA164E5B0/
C
C                                  COEFFICIENTS FOR ARG GREATER THAN 1.
C
      DATA               PP(1)/Z426706A7DF554543/,
     *                   PP(2)/Z43ADD8E945C002B9/,
     *                   PP(3)/Z444B49CD0D6B3B21/,
     *                   PP(4)/Z44CE173D5AC3D3EB/,
     *                   PP(5)/Z451091CFAA1AF762/,
     *                   PP(6)/Z44ACCA76AD4CF2C4/,
     *                   PP(7)/Z4439B53F6B0C57E7/,
     *                   PP(8)/Z43943238B5B339F6/,
     *                   PP(9)/Z428EC73B9E8B6D28/
      DATA               QQ(1)/Z42AC5AD379FF3E04/,
     *                   QQ(2)/Z43C94D887040D994/,
     *                   QQ(3)/Z4449F95E104C2DF5/,
     *                   QQ(4)/Z44B965D632B951C5/,
     *                   QQ(5)/Z44E2A1452C65B929/,
     *                   QQ(6)/Z448F3E83BCF7E10E/,
     *                   QQ(7)/Z442EF19A669A45AF/,
     *                   QQ(8)/Z437722298F7588BB/,
     *                   QQ(9)/Z4271EBAA1617ACF8/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      X=ARG
      IF(X .LE. ZERO) GO TO 40
      A = ONE
      IF (IOPT .NE. 1 .AND. IOPT .NE. 2) GO TO 50
      GO TO (5,10),IOPT
C                                  IOPT = 1
    5 IF(X.GT.XMAX) GO TO 45
      IF(X .LE. ONE) GO TO 15
      A=DEXP(-X)
      GO TO 25
C                                  IOPT = 2
   10 IF(X .GT. ONE) GO TO 25
      A=DEXP(X)
C                                  ARG IS GREATER THAN ZERO AND LESS
C                                  THAN 1
   15 TEMP=DLOG(X)
      IF(X .LT. 1.0D-39) GO TO 20
      XX = X*X
      SUMP=(((P(1)*XX+P(2))*XX+P(3))*XX+P(4))*XX+P(5)
      SUMQ=((XX+Q(1))*XX+Q(2))*XX+Q(3)
      SUMF=((F(1)*XX+F(2))*XX+F(3))*XX+F(4)
      SUMG=((XX+G(1))*XX+G(2))*XX+G(3)
      MMBSK0 = (SUMP/SUMQ - XX*SUMF*TEMP/SUMG - TEMP)*A
      IF (MMBSK0 .GE. ONE) MMBSK0 = MMBSK0+FUDGE
      GO TO 9005
C                                  SPECIAL CASE - SMALL ARG
   20 MMBSK0=CONST-TEMP
      GO TO 9005
C                                  ARG IS GREATER THAN 1.
   25 XX = ONE/X
      SUMP=PP(1)
      DO 30 I=2,9
         SUMP=SUMP*XX+PP(I)
   30 CONTINUE
      SUMQ = XX
      DO 35 I=1,8
         SUMQ=(SUMQ+QQ(I))*XX
   35 CONTINUE
      SUMQ=SUMQ+QQ(9)
      MMBSK0 = SUMP/SUMQ/DSQRT(X)*A
      GO TO 9005
C                                  TERMINAL ERROR - ARG IS OUT OF
C                                  RANGE
   40 IER = 130
      MMBSK0 = XINF
      GO TO 9000
   45 MMBSK0 = ZERO
      IER = 130
      GO TO 9000
C                                  TERMINAL ERROR - IOPT IS OUT OF
C                                  RANGE
   50 MMBSK0 = XINF
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HMMBSK0)
 9005 RETURN
      END

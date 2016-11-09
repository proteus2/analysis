C   IMSL ROUTINE NAME   - MMBSK1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MODIFIED BESSEL FUNCTION OF THE SECOND
C                           KIND OF ORDER ONE
C
C   USAGE               - FUNCTION MMBSK1 (IOPT,ARG,IER)
C
C   ARGUMENTS    MMBSK1 - OUTPUT VALUE OF THE FUNCTION AT ARG. MMBSK1
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE PRECISION/HARDWARE SECTION.)
C                IOPT   - INPUT OPTION PARAMETER.
C                         IF IOPT = 1, THE MODIFIED BESSEL FUNCTION OF
C                           THE SECOND KIND OF ORDER ONE FOR ARGUMENT
C                           ARG IS EVALUATED. ARG MUST BE GREATER THAN
C                           ARGMIN AND LESS THAN OR EQUAL TO ARGMAX.
C                           ARGMIN IS OF THE ORDER OF 10**(-75) AND
C                           ARGMAX IS AT LEAST 86. THE EXACT VALUES OF
C                           ARGMIN AND ARGMAX MAY ALLOW LARGER RANGES
C                           FOR ARG ON SOME COMPUTERS.
C                           SEE THE PROGRAMMING NOTES IN THE MANUAL
C                           FOR THE EXACT VALUES.
C                         IF IOPT = 2, EXP(ARG)*THE MODIFIED BESSEL
C                           FUNCTION OF THE SECOND KIND OF ORDER ONE
C                           FOR ARGUMENT ARG IS EVALUATED. ARG MUST
C                           BE GREATER THAN ARGMIN.
C                ARG    - INPUT ARGUMENT. ARG MUST BE TYPED APPROPRIATE-
C                           LY IN THE CALLING PROGRAM. (SEE THE PRE-
C                           CISION/HARDWARE SECTION.)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ARG WAS LESS THAN
C                             ARGMIN OR IOPT WAS NOT 1 OR 2. MMBSK1 IS
C                             SET TO MACHINE INFINITY.
C                           IER = 130 INDICATES THAT ARG WAS OUT OF
C                             RANGE. MMBSK1 IS SET TO MACHINE INFINITY.
C                             HOWEVER, IF IOPT = 1 AND ARG IS GREATER
C                             THAN OR EQUAL TO ARGMAX, MMBSK1 IS SET
C                             TO ZERO.
C
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMBSK1 (IOPT,ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      DOUBLE PRECISION   P(5),Q(4),PP(10),QQ(8),F(5),G(2)
      DOUBLE PRECISION   A,X,XX,SUMF,SUMG,SUMP,SUMQ,
     *                   XINF,XMAX,ARGMIN,ZERO,ONE
      DATA               ONE/1.D0/,ZERO/0.0D0/
      DATA               XINF/Z7FFFFFFFFFFFFFFF/,
     *                   ARGMIN/Z0210000000000001/,
     *                   XMAX/Z42B1DB0D7E674911/
C
C                                  COEFFICIENTS FOR ARG GREATER THAN
C                                  0. AND LESS THAN OR EQUAL TO 1.
C
      DATA               P(1)/Z401ECD23A42954E6/,P(2)/Z4218FF72AA1E8E93/
      DATA               P(3)/Z437052272DBC8C9C/,P(4)/Z44AD2D4F61EB21F8/
      DATA               P(5)/Z452BE874CD784BD2/
      DATA               Q(1)/Z4040000000000000/,Q(2)/ZC2465C1B283FCEB7/
      DATA               Q(3)/Z442464131D7172BF/,Q(4)/ZC5873065F38838A6/
      DATA               F(1)/Z3E8A1FF98FFA9100/,F(2)/Z40685AC008AA61B5/
      DATA               F(3)/Z421EEBDFC81F7AB5/,F(4)/Z433B7CAB71F3E3B0/
      DATA               F(5)/Z4421350CAA6B5B2A/
      DATA               G(1)/ZC2DDADBC5C5A3FCD/,G(2)/Z44426A1954D6B652/
C
C                                  COEFFICIENTS FOR ARG GREATER THAN 1.
C
      DATA               PP(1)/Z4011E0259AA8F01A/,
     *                   PP(2)/Z416F184985782E20/,
     *                   PP(3)/Z42659FE01E50B13E/,
     *                   PP(4)/Z4320705CAED6286D/,
     *                   PP(5)/Z434B71EB257D7C50/,
     *                   PP(6)/Z4359665F21891D19/,
     *                   PP(7)/Z4338495263DA8BE7/,
     *                   PP(8)/Z4312AA52E742DAA4/,
     *                   PP(9)/Z4230547FEC2848D4/,
     *                   PP(10)/Z412F54ABEF21C54C/
      DATA               QQ(1)/Z421E2BC91FAA1AB5/,
     *                   QQ(2)/Z42E41D72804713D9/,
     *                   QQ(3)/Z432A1032DF55B81B/,
     *                   QQ(4)/Z4339C9B14E56BA46/,
     *                   QQ(5)/Z4327E927469D8A3D/,
     *                   QQ(6)/Z42E06EA819073B1A/,
     *                   QQ(7)/Z4225AD3D5971818A/,
     *                   QQ(8)/Z4125C3B488C64098/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      X = ARG
      IF (X .LT. ARGMIN) GO TO 40
      A = ONE
      IF (IOPT .NE. 1 .AND. IOPT .NE. 2) GO TO 50
      GO TO (5,10),IOPT
C                                  IOPT = 1
    5 IF (X .GT. XMAX) GO TO 45
      IF (X .LE. ONE) GO TO 15
      A = DEXP(-X)
      GO TO 25
C                                  IOPT = 2
   10 IF (X .GT. ONE) GO TO 25
      A = DEXP(X)
C                                  ARG IS GREATER THAN 0 AND LESS THAN
C                                  OR EQUAL TO 1.
   15 IF (X .LT. 3.0D-39) GO TO 20
      XX = X * X
      SUMP = ((((P(1)*XX+P(2))*XX+P(3))*XX+P(4))*XX+P(5))*XX+Q(4)
      SUMQ = ((Q(1) * XX + Q(2)) * XX + Q(3)) * XX + Q(4)
      SUMF = (((F(1) * XX + F(2)) * XX + F(3)) * XX + F(4)) * XX + F(5)
      SUMG = (XX + G(1)) * XX + G(2)
      MMBSK1 = ((XX * DLOG(X) * SUMF / SUMG + SUMP / SUMQ) / X) * A
      GO TO 9005
C                                  RETURN FOR SMALL ARG
   20 MMBSK1 = ONE/X
      GO TO 9005
C                                  ARG GREATER THAN 1.
   25 XX = ONE/X
      SUMP = PP(1)
      DO 30 I = 2, 10
         SUMP = SUMP * XX + PP(I)
   30 CONTINUE
      SUMQ = XX
      DO 35 I = 1,7
         SUMQ = (SUMQ + QQ(I)) * XX
   35 CONTINUE
      SUMQ = SUMQ + QQ(8)
      MMBSK1 = A * SUMP / SUMQ / DSQRT(X)
      GO TO 9005
C                                  TERMINAL ERROR - ARG IS LESS THAN
C                                  ARGMIN
   40 IER = 130
      MMBSK1 = XINF
      GO TO 9000
C                                  TERMINAL ERROR - ARG IS GREATER THAN
C                                  177.855
   45 IER = 130
      MMBSK1 = ZERO
      GO TO 9000
C                                  TERMINAL ERROR - IOPT NOT IN RANGE
   50 MMBSK1 = XINF
      IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMMBSK1)
 9005 RETURN
      END

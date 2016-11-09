C   IMSL ROUTINE NAME   - MMLIND
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - INCOMPLETE ELLIPTIC INTEGRAL OF THE SECOND
C                           KIND
C
C   USAGE               - FUNCTION MMLIND (X,Y,Z,IER)
C
C   ARGUMENTS    MMLIND - OUTPUT VALUE OF THE INTEGRAL. MMLIND MUST BE
C                           TYPED APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE PRECISION/HARDWARE SECTION.)
C                X      - INPUT FIRST VARIABLE OF THE INCOMPLETE
C                           ELLIPTIC INTEGRAL. X MUST BE NONNEGATIVE.
C                           SEE THE REMARKS SECTION BELOW FOR FURTHER
C                           RESTRICTIONS ON X.
C                Y      - INPUT SECOND VARIABLE OF THE INCOMPLETE
C                           ELLIPTIC INTEGRAL. Y MUST BE NONNEGATIVE.
C                           SEE THE REMARKS SECTION BELOW FOR FURTHER
C                           RESTRICTIONS ON Y.
C                Z      - INPUT THIRD VARIABLE OF THE INCOMPLETE
C                           ELLIPTIC INTEGRAL. Z MUST BE POSITIVE.
C                           SEE THE REMARKS SECTION BELOW FOR FURTHER
C                           RESTRICTIONS ON Z.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER X, Y, OR Z
C                             IS NEGATIVE. MMLIND IS SET TO MACHINE
C                             INFINITY.
C                           IER = 130 INDICATES THAT EITHER X+Y OR Z
C                             IS LESS THAN ARGMIN. MMLIND IS SET TO
C                             MACHINE INFINITY.
C                           IER = 131 INDICATES THAT EITHER X, Y, OR Z
C                             IS GREATER THAN ARGMAX. MMLIND IS SET
C                             TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE SUM, X+Y, AND THE ARGUMENT, Z, MUST BE GREATER
C                THAN OR EQUAL TO ARGMIN. ARGMIN IS DEFINED AS FOLLOWS.
C
C                ARGMIN = THE MAXIMUM OF -
C                           3 * (MACHINE MINIMUM) ** (2/3)
C                                     AND
C                           3 / (MACHINE INFINITY) ** (2/3)
C
C                ALSO, EACH OF X, Y, AND Z MUST BE LESS THAN OR EQUAL
C                TO ARGMAX. ARGMAX IS DEFINED, APPROXIMATELY, AS
C
C                ARGMAX = ((0.085*(P)**(-1/6)) / MIN) ** (2/3)
C
C                WHERE P IS THE NUMBER OF DECIMAL DIGITS OF PRECISION
C                AVAILABLE AND MIN IS THE MACHINE MINIMUM.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMLIND(X,Y,Z,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   X,Y,Z
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   ARGMN1,ARGMN2,ARGMIN,ARGMAX,ERRTOL
      DOUBLE PRECISION   C1,C2,C3,C4,EA,EB,EC,ED,EF,EPSLON,LAMDA,SETA
      DOUBLE PRECISION   LOLIM,MU,POWER4,SIGMA,S1,S2,UPLIM,XN,XNDEV
      DOUBLE PRECISION   XNROOT,XINF,YN,YNDEV,YNROOT,ZN,ZNDEV,ZNROOT
      DATA               XINF/.723700557733226D+76/
      DATA               SETA/Z0010000000000000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      ERRTOL = 1.D-3
C                                  INITIALIZE UPPER AND LOWER LIMITS
      ARGMN1 = 3.D0 * (SETA**(2.D0/3.D0))
      ARGMN2 = 3.D0 / (XINF**(2.D0/3.D0))
      ARGMIN = DMAX1(ARGMN1,ARGMN2)
      ARGMAX = DEXP(2.D0/3.D0 * (DLOG(0.1D0*ERRTOL)-DLOG(SETA)))
C
      IF (DMIN1(X,Y,Z).GE.0.D0) GO TO 5
      IER = 129
      MMLIND = XINF
      GO TO 9000
    5 IF (DMIN1(X+Y,Z).GE.ARGMIN) GO TO 10
      IER = 130
      MMLIND = XINF
      GO TO 9000
   10 IF (DMAX1(X,Y,Z) .LE. ARGMAX) GO TO 15
      IER = 131
      MMLIND = XINF
      GO TO 9000
   15 XN = X
      YN = Y
      ZN = Z
C                                  BEGIN LOOP
      SIGMA = 0.D0
      POWER4 = 1.D0
C
   20 MU = (XN + YN + 3.D0 * ZN) * 0.2D0
      XNDEV = (MU - XN) / MU
      YNDEV = (MU - YN) / MU
      ZNDEV = (MU - ZN) / MU
      EPSLON = DMAX1(DABS(XNDEV),DABS(YNDEV),DABS(ZNDEV))
      IF (EPSLON .LT. ERRTOL) GO TO 25
      XNROOT = DSQRT(XN)
      YNROOT = DSQRT(YN)
      ZNROOT = DSQRT(ZN)
      LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
      SIGMA = SIGMA + POWER4 / (ZNROOT * (ZN + LAMDA))
      POWER4 = POWER4 * 0.25D0
      XN = (XN + LAMDA) * 0.25D0
      YN = (YN + LAMDA) * 0.25D0
      ZN = (ZN + LAMDA) * 0.25D0
      GO TO 20
C
   25 C1 = 3.D0 / 14.D0
      C2 = 1.D0 / 6.D0
      C3 = 9.D0 / 22.D0
      C4 = 3.D0 / 26.D0
      EA = XNDEV * YNDEV
      EB = ZNDEV * ZNDEV
      EC = EA - EB
      ED = EA - 6.D0 * EB
      EF = ED + EC + EC
      S1 = ED * (- C1 + 0.25D0 * C3 * ED - 1.5D0 * C4 * ZNDEV * EF)
      S2 = ZNDEV * (C2 * EF + ZNDEV * (- C3 * EC + ZNDEV * C4 * EA))
      MMLIND=3.D0 * SIGMA + POWER4 * (1.D0 + S1 + S2) / (MU * DSQRT(MU))
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMLIND)
 9005 RETURN
      END

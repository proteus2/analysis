C   IMSL ROUTINE NAME   - MMLINF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - INCOMPLETE ELLIPTIC INTEGRAL OF THE FIRST
C                           KIND
C
C   USAGE               - FUNCTION MMLINF (X,Y,Z,IER)
C
C   ARGUMENTS    MMLINF - OUTPUT VALUE OF THE INTEGRAL. MMLINF MUST BE
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
C                           ELLIPTIC INTEGRAL. Z MUST BE NONNEGATIVE.
C                           SEE THE REMARKS SECTION BELOW FOR FURTHER
C                           RESTRICTIONS ON Z.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER X, Y, OR Z
C                             IS NEGATIVE. MMLINF IS SET TO MACHINE
C                             INFINITY.
C                           IER = 130 INDICATES THAT ONE OF THE SUMS,
C                             X+Y, X+Z, OR Y+Z, IS LESS THAN ARGMIN.
C                             MMLINF IS SET TO MACHINE INFINITY.
C                           IER = 131 INDICATES THAT EITHER X, Y, OR Z
C                             IS GREATER THAN ARGMAX. MMLINF IS SET
C                             TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE SUMS, X+Y, X+Z, AND Y+Z MUST BE GREATER THAN OR
C                EQUAL TO ARGMIN. ARGMIN IS THE MACHINE MINIMUM * 5.
C                ALSO, EACH OF X, Y, AND Z MUST BE LESS THAN OR EQUAL
C                TO ARGMAX. ARGMAX IS THE MACHINE MAXIMUM DIVIDED BY 5.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMLINF(X,Y,Z,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   X,Y,Z
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   C1,C2,C3,E2,E3,EPSLON,LAMDA,SETA,XINF
      DOUBLE PRECISION   ARGMIN,MU,S,ARGMAX,XN,XNDEV,XNROOT
      DOUBLE PRECISION   YN,YNDEV,YNROOT,ERRTOL,ZN,ZNDEV,ZNROOT
      DATA               XINF/.723700557733226D+76/
      DATA               SETA/Z0010000000000000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      ERRTOL = 1.D-3
C                                  INITIALIZE UPPER AND LOWER LIMITS
      ARGMIN = 5.D0*SETA
      ARGMAX = XINF/5.D0
      IF (DMIN1(X,Y,Z).GE.0.D0) GO TO 5
      IER = 129
      MMLINF = XINF
      GO TO 9000
    5 IF (DMIN1(X+Y,X+Z,Y+Z).GE.ARGMIN) GO TO 10
      IER = 130
      MMLINF = XINF
      GO TO 9000
   10 IF (DMAX1(X,Y,Z) .LE. ARGMAX) GO TO 15
      IER = 131
      MMLINF = XINF
      GO TO 9000
   15 XN = X
      YN = Y
      ZN = Z
C                                  BEGIN LOOP
   20 MU = (XN + YN + ZN) / 3.D0
      XNDEV = 2.D0 - (MU + XN) / MU
      YNDEV = 2.D0 - (MU + YN) / MU
      ZNDEV = 2.D0 - (MU + ZN) / MU
      EPSLON = DMAX1(DABS(XNDEV),DABS(YNDEV),DABS(ZNDEV))
      IF (EPSLON .LT. ERRTOL) GO TO 25
      XNROOT = DSQRT(XN)
      YNROOT = DSQRT(YN)
      ZNROOT = DSQRT(ZN)
      LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
      XN = (XN + LAMDA) * 0.25D0
      YN = (YN + LAMDA) * 0.25D0
      ZN = (ZN + LAMDA) * 0.25D0
      GO TO 20
C                                  FINAL CALCULATION
   25 C1 = 1.D0 / 24.D0
      C2 = 3.D0 / 44.D0
      C3 = 1.D0 / 14.D0
      E2 = XNDEV * YNDEV - ZNDEV * ZNDEV
      E3 = XNDEV * YNDEV * ZNDEV
      S = 1.D0 + (C1 * E2 - 0.1D0 - C2 * E3) * E2 + C3 * E3
      MMLINF = S / DSQRT(MU)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMLINF)
 9005 RETURN
      END

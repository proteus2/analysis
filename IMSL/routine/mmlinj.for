C   IMSL ROUTINE NAME   - MMLINJ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - INCOMPLETE ELLIPTIC INTEGRAL OF THE THIRD
C                           KIND
C
C   USAGE               - FUNCTION MMLINJ (X,Y,Z,P,IER)
C
C   ARGUMENTS    MMLINJ - OUTPUT VALUE OF THE INTEGRAL. MMLINJ MUST BE
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
C                P      - INPUT FOURTH VARIABLE OF THE INCOMPLETE
C                           ELLIPTIC INTEGRAL. P MUST BE POSITIVE.
C                           SEE THE REMARKS SECTION BELOW FOR FURTHER
C                           RESTRICTIONS ON P.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER X, Y, OR Z
C                             IS NEGATIVE. MMLINJ IS SET TO MACHINE
C                             INFINITY.
C                           IER = 130 INDICATES THAT EITHER X+Y, X+Z,
C                             Y+Z, OR P IS LESS THAN ARGMIN. MMLINJ
C                             IS SET TO MACHINE INFINITY.
C                           IER = 131 INDICATES THAT EITHER X, Y, Z, OR
C                             P IS GREATER THAN ARGMAX. MMLINJ IS SET
C                             TO MACHINE INFINITY.
C                           IER = 132 INDICATES THAT AN ERROR OCCURRED
C                             IN IMSL ROUTINE MMLINC. MMLINJ IS SET TO
C                             MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - MMLINC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE SUMS, X+Y, X+Z, Y+Z, AND THE ARGUMENT, P, MUST BE
C                GREATER THAN OR EQUAL TO ARGMIN. ARGMIN IS DEFINED AS
C                THE CUBE ROOT OF (THE MACHINE MINIMUM * 5). ALSO, EACH
C                OF X, Y, Z, AND P MUST BE LESS THAN OR EQUAL TO ARGMAX.
C                ARGMAX IS DEFINED AS .3 * THE CUBE ROOT OF (THE MACHINE
C                MAXIMUM DIVIDED BY 5).
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMLINJ(X,Y,Z,P,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   X,Y,Z
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   ARGMIN,ARGMAX,ERRTOL
      DOUBLE PRECISION   ALFA,BETA,C1,C2,C3,C4,EA,EB,EC,E2,E3,MMLINC
      DOUBLE PRECISION   EPSLON,SETA,LAMDA,LOLIM,MU,P,PN,PNDEV
      DOUBLE PRECISION   POWER4,RC,SIGMA,S1,S2,S3,UPLIM,XN,XNDEV
      DOUBLE PRECISION   XNROOT,YN,YNDEV,YNROOT,XINF,ZN,ZNDEV,ZNROOT
      DATA               XINF/.723700557733226D+76/
      DATA               SETA/Z0010000000000000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      ERRTOL = 1.D-3
C                                  INITIALIZE UPPER AND LOWER LIMITS
      ARGMIN = (5.D0*SETA)**(1.D0/3.D0)
      ARGMAX = ((XINF/5.D0)**(1.D0/3.D0)) * .3D0
      IF (DMIN1(X,Y,Z).GE.0.D0) GO TO 5
      IER = 129
      MMLINJ = XINF
      GO TO 9000
    5 IF (DMIN1(X+Y,X+Z,Y+Z,P).GE.ARGMIN) GO TO 10
      IER = 130
      MMLINJ = XINF
      GO TO 9000
   10 IF (DMAX1(X,Y,Z,P) .LE. ARGMAX) GO TO 15
      IER = 131
      MMLINJ = XINF
      GO TO 9000
   15 XN = X
      YN = Y
      ZN = Z
C                                  BEGIN LOOP
C
      PN = P
      SIGMA = 0.D0
      POWER4 = 1.D0
C
   20 MU = (XN + YN + ZN + PN + PN) * 0.2D0
      XNDEV = (MU - XN) / MU
      YNDEV = (MU - YN) / MU
      ZNDEV = (MU - ZN) / MU
      PNDEV = (MU - PN) / MU
      EPSLON = DMAX1(DABS(XNDEV),DABS(YNDEV),DABS(ZNDEV),DABS(PNDEV))
      IF (EPSLON .LT. ERRTOL) GO TO 30
      XNROOT = DSQRT(XN)
      YNROOT = DSQRT(YN)
      ZNROOT = DSQRT(ZN)
      LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
      ALFA = PN * (XNROOT + YNROOT + ZNROOT) + XNROOT * YNROOT * ZNROOT
      ALFA = ALFA * ALFA
      BETA = PN * (PN + LAMDA) * (PN + LAMDA)
      RC = MMLINC(ALFA,BETA,JER)
      IF(JER.EQ.0) GO TO 25
      IER = 132
      MMLINJ = XINF
      GO TO 9000
   25 SIGMA = SIGMA + POWER4 * RC
      POWER4 = POWER4 * 0.25D0
      XN = (XN + LAMDA) * 0.25D0
      YN = (YN + LAMDA) * 0.25D0
      ZN = (ZN + LAMDA) * 0.25D0
      PN = (PN + LAMDA) * 0.25D0
      GO TO 20
C
   30 C1 = 3.D0 / 14.D0
      C2 = 1.D0 / 3.D0
      C3 = 3.D0 / 22.D0
      C4 = 3.D0 / 26.D0
      EA = XNDEV * (YNDEV + ZNDEV) + YNDEV * ZNDEV
      EB = XNDEV * YNDEV * ZNDEV
      EC = PNDEV * PNDEV
      E2 = EA - 3.D0 * EC
      E3 = EB + 2.D0 * PNDEV * (EA - EC)
      S1 = 1.D0 + E2 * (- C1 + 0.75D0 * C3 * E2 - 1.5D0 * C4 * E3)
      S2 = EB * (0.5D0 * C2 + PNDEV * (- C3 - C3 + PNDEV * C4))
      S3 = PNDEV * EA * (C2 - PNDEV * C3) - C2 * PNDEV * EC
      MMLINJ=3.D0 * SIGMA + POWER4 * (S1 + S2 + S3) / (MU * DSQRT(MU))
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMLINJ)
 9005 RETURN
      END

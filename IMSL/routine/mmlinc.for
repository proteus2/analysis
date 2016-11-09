C   IMSL ROUTINE NAME   - MMLINC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - ELEMENTARY INTEGRAL FROM WHICH INVERSE
C                           CIRCULAR FUNCTIONS, LOGARITHMS, OR INVERSE
C                           HYPERBOLIC FUNCTIONS MAY BE COMPUTED
C
C   USAGE               - FUNCTION MMLINC (X,Y,IER)
C
C   ARGUMENTS    MMLINC - OUTPUT VALUE OF THE INTEGRAL. MMLINC MUST BE
C                           TYPED APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE PRECISION/HARDWARE SECTION.)
C                X      - INPUT FIRST VARIABLE OF THE ELEMENTARY
C                           INTEGRAL. X MUST BE NONNEGATIVE. SEE THE
C                           REMARKS SECTION BELOW FOR FURTHER
C                           RESTRICTIONS ON X.
C                Y      - INPUT SECOND VARIABLE OF THE ELEMENTARY
C                           INTEGRAL. Y MUST BE POSITIVE. SEE THE
C                           REMARKS SECTION BELOW FOR FURTHER
C                           RESTRICTIONS ON Y.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER X IS
C                             NEGATIVE OR Y IS NOT POSITIVE. MMLINC
C                             IS SET TO MACHINE INFINITY.
C                           IER = 130 INDICATES THAT THE SUM, X + Y,
C                             IS LESS THAN ARGMIN. MMLINC IS SET TO
C                             MACHINE INFINITY.
C                           IER = 131 INDICATES THAT EITHER X OR Y IS
C                             GREATER THAN ARGMAX. MMLINC IS SET TO
C                             MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE SUM OF X AND Y MUST BE GREATER THAN OR EQUAL TO
C                ARGMIN. ARGMIN IS THE MACHINE MINIMUM * 5. ALSO, BOTH
C                X AND Y MUST BE LESS THAN OR EQUAL TO ARGMAX. ARGMAX
C                IS THE MACHINE MAXIMUM DIVIDED BY 5.
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMLINC(X,Y,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   X,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   C1,C2,LAMDA,ARGMIN,SETA,XINF,XYSUM
      DOUBLE PRECISION   MU,S,SN,ARGMAX,XN,YN,ERRTOL
      DATA               XINF/.723700557733226D+76/
      DATA               SETA/Z0010000000000000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      ERRTOL = 1.D-3
C                                  INITIALIZE UPPER AND LOWER LIMITS
      ARGMIN = 5.D0*SETA
      ARGMAX = XINF/5.D0
      IF (X .GE. 0.D0  .AND.  Y .GT. 0.D0) GO TO 5
      IER = 129
      MMLINC = XINF
      GO TO 9000
    5 XYSUM = X + Y
      IF (XYSUM .GE. ARGMIN) GO TO 10
      IER = 130
      MMLINC = XINF
      GO TO 9000
   10 IF (DMAX1(X,Y) .LE. ARGMAX) GO TO 15
      IER = 131
      MMLINC = XINF
      GO TO 9000
   15 XN = X
      YN = Y
C                                  BEGIN LOOP
   20 MU = (XN + YN + YN) / 3.D0
      SN = (YN + MU) / MU - 2.D0
      IF (DABS(SN) .LT. ERRTOL) GO TO 25
      LAMDA = 2.D0 * DSQRT(XN) * DSQRT(YN) + YN
      XN = (XN + LAMDA) * 0.25D0
      YN = (YN + LAMDA) * 0.25D0
      GO TO 20
C                                  FINAL CALCULATION
   25 C1 = 1.D0 / 7.D0
      C2 = 9.D0 / 22.D0
      S = SN * SN * (0.3D0 + SN * (C1 + SN * (0.375D0 + SN * C2)))
      MMLINC = (1.D0 + S) / DSQRT(MU)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMLINC)
 9005 RETURN
      END

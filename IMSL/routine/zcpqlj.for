C   IMSL ROUTINE NAME   - ZCPQLJ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION ZCPQLJ (NN,PT,REPSR1,RINFP,REPSP,RADIX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NN
      REAL               PT(NN),REPSR1,RINFP,REPSP,RADIX
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,L
      REAL               RHI,RLO,RMAX,RMIN,X,SC,ZERO,HALF,ONE
      DATA               ZERO,HALF,ONE/0.0,0.5,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      RHI = SQRT(RINFP)
      RLO = REPSP/REPSR1
      RMAX = ZERO
      RMIN = RINFP
C                                  RETURNS A SCALE FACTOR TO MULTIPLY
C                                    THE COEFFICIENTS OF THE POLYNOMIAL.
C                                    THE SCALING IS DONE TO AVOID
C                                    OVERFLOW AND TO AVOID UNDETECTED
C                                    UNDERFLOW INTERFERING WITH THE
C                                    CONVERGENCE CRITERION. THE FACTOR
C                                    IS A POWER OF THE BASE(RADIX).
C                                  PT - MODULUS OF COEFFICIENTS OF P
C                                  REPSR1,RINFP,REPSP,RADIX - CONSTANTS
C                                    DESCRIBING THE FLOATING POINT
C                                    ARITHMETIC.
C                                  FIND LARGEST AND SMALLEST MODULI OF
C                                    COEFFICIENTS.
      DO 5 I=1,NN
         X = PT(I)
         IF (X.GT.RMAX) RMAX = X
         IF (X.NE.ZERO.AND.X.LT.RMIN) RMIN = X
    5 CONTINUE
C                                  SCALE ONLY IF THERE ARE VERY LARGE
C                                    OR VERY SMALL COMPONENTS
      ZCPQLJ = ONE
      IF (RMIN.GE.RLO.AND.RMAX.LE.RHI) RETURN
      X = RLO/RMIN
      IF (X.GT.ONE) GO TO 10
      SC = ONE/(SQRT(RMAX)*SQRT(RMIN))
      GO TO 15
   10 SC = X
      IF (RINFP/SC.LT.RMAX) SC = ONE
   15 L = ALOG(SC)/ALOG(RADIX)+HALF
      ZCPQLJ = RADIX**L
      RETURN
      END

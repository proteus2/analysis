C   IMSL ROUTINE NAME   - MDNORD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NORMAL OR GAUSSIAN PROBABILITY DISTRIBUTION
C                           FUNCTION OF A DOUBLE PRECISION ARGUMENT
C
C   USAGE               - CALL MDNORD (Y,P)
C
C   ARGUMENTS    Y      - INPUT DOUBLE PRECISION VALUE AT WHICH THE
C                           FUNCTION IS TO BE EVALUATED.
C                P      - OUTPUT DOUBLE PRECISION PROBABILITY THAT A
C                           RANDOM VARIABLE HAVING A NORMAL (0,1)
C                           DISTRIBUTION WILL BE LESS THAN OR EQUAL
C                           TO Y.
C
C   REQD. IMSL ROUTINES - MERRCD=DERFC
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDNORD (Y,P)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   P,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   DERFC,SQR1D2
      DATA               SQR1D2/.7071067811865475D0/
C                                  FIRST EXECUTABLE STATEMENT
      P = -Y * SQR1D2
      IF (DABS(P) .LE. 13.2D0) GO TO 5
      P = 0.0D0
      IF (Y .LT. 0.0D0) RETURN
      P = 1.0D0
      RETURN
    5 P = .5D0 * DERFC(P)
      RETURN
      END

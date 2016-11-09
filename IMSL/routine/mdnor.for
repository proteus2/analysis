C   IMSL ROUTINE NAME   - MDNOR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NORMAL OR GAUSSIAN PROBABILITY DISTRIBUTION
C                           FUNCTION
C
C   USAGE               - CALL MDNOR (Y,P)
C
C   ARGUMENTS    Y      - INPUT VALUE AT WHICH FUNCTION IS TO BE
C                           EVALUATED.
C                P      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE
C                           HAVING A NORMAL (0,1) DISTRIBUTION WILL BE
C                           LESS THAN OR EQUAL TO Y.
C
C   REQD. IMSL ROUTINES - MERRC=ERFC
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDNOR  (Y,P)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               SQR1D2
      DATA               SQR1D2/.7071068/
C                                  FIRST EXECUTABLE STATEMENT
      P = -Y * SQR1D2
      IF (ABS(P) .LE. 13.2) GO TO 5
      P = 0.0
      IF (Y .LT. 0.0) RETURN
      P = 1.0
      RETURN
    5 P = .5 * ERFC(P)
      RETURN
      END

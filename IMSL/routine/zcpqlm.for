C   IMSL ROUTINE NAME   - ZCPQLM
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
      SUBROUTINE ZCPQLM (P1,P2,P3,P4)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P1,P2,P3,P4
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               REPSR1,RINFP,REPSP,RADIX
      DATA               REPSR1/Z3C100000/
      DATA               RINFP/Z7FFFFFFF/
      DATA               REPSP/Z00100000/
      DATA               RADIX/16.0/
C                                  ZCPQLM PROVIDES MACHINE CONSTANTS
C                                    USED IN VARIOUS PARTS OF THE
C                                    PROGRAM. THE USER MAY EITHER SET
C                                    THEM DIRECTLY OR USE THE STATEMENTS
C                                    BELOW TO COMPUTE THEM. THE MEANING
C                                    OF THE FOUR CONSTANTS ARE -
C                                  REPSR1 THE MAXIMUM RELATIVE
C                                    REPRESENTATION ERROR WHICH CAN BE
C                                    DESCRIBED AS THE SMALLEST POSITIVE
C                                    FLOATING-POINT NUMBER SUCH THAT
C                                    1.0D0 + ETA IS GREATER THAN 1.0D0
C                                  RINFP THE LARGEST FLOATING-POINT
C                                    NUMBER
C                                  REPSP THE SMALLEST POSITIVE FLOATING-
C                                    POINT NUMBER
C                                  RADIX THE BASE OF THE FLOATING-POINT
C                                    NUMBER SYSTEM USED
C                                  FIRST EXECUTABLE STATEMENT
      P1 = REPSR1
      P2 = RINFP
      P3 = REPSP
      P4 = RADIX
      RETURN
      END

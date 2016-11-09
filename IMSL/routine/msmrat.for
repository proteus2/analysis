C   IMSL ROUTINE NAME   - MSMRAT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - RATIO OF THE ORDINATE TO THE UPPER TAIL AREA
C                           OF THE STANDARDIZED NORMAL (GAUSSIAN)
C                           DISTRIBUTION
C
C   USAGE               - CALL MSMRAT (X,RM,IER)
C
C   ARGUMENTS    X      - INPUT VALUE AT WHICH THE RATIO IS CALCULATED.
C                           X MUST BE LESS THAN SUPRBD, WHICH IS AT
C                           LEAST 12.
C                RM     - OUTPUT RATIO. IF X IS LESS THAN SLRBD, RM IS
C                           SET TO 0.0. SLRBD IS NO LARGER THAN -13.
C                           THE EXACT VALUES OF SUPRBD AND SLRBD MAY
C                           ALLOW LARGER RANGES FOR X ON SOME COMPUTERS.
C                           SEE THE PROGRAMMING NOTES IN THE MANUAL FOR
C                           THE EXACT VALUES.
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES X IS GREATER THAN
C                             SUPRBD. RM IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - MERRC=ERFC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MSMRAT  (X,RM,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               X,RM
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               SPI,SQ1H,RINFP,SLRBD,SUPRBD
C                                  SPI IS 2.0/SQRT(2*PI)
      DATA               SPI/.7978846E0/
      DATA               SLRBD/-18.97366E0/
      DATA               SUPRBD/18.70000E0/
      DATA               SQ1H/.7071068/
      DATA               RINFP/Z7FFFFFFF/
C                                  INITIALIZE IER
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (X .GE. SLRBD) GO TO 5
C                                  X IS LESS THAN SLRBD. RM IS SET TO
C                                  ZERO AND RETURNED.
      RM = 0.0
      GO TO 9005
    5 IF (X .LE. SUPRBD) GO TO 10
C                                  TERMINAL ERROR - X IS GREATER THAN
C                                  SUPRBD. RM IS SET TO MACHINE
C                                  INFINITY.
      IER = 129
      RM = RINFP
      GO TO 9000
   10 RM = SPI*EXP(-X*X*.5)
      RM = RM/ERFC(SQ1H*X)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMSMRAT)
 9005 RETURN
      END

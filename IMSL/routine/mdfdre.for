C   IMSL ROUTINE NAME   - MDFDRE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - F PROBABILITY DISTRIBUTION FUNCTION (INTEGER
C                           OR FRACTIONAL DEGREES OF FREEDOM)
C
C   USAGE               - CALL MDFDRE (X,DFN,DFD,P,IER)
C
C   ARGUMENTS    X      - INPUT VALUE TO WHICH FUNCTION IS TO BE INTE-
C                           GRATED. X MUST BE GREATER THAN OR EQUAL TO
C                           ZERO.
C                DFN    - INPUT NUMERATOR DEGREES OF FREEDOM (MUST BE
C                           GREATER THAN ZERO)
C                DFD    - INPUT DENOMINATOR DEGREES OF FREEDOM (MUST BE
C                           GREATER THAN ZERO)
C                P      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE
C                           FROM AN F DISTRIBUTION HAVING DFN AND DFD
C                           DEGREES OF FREEDOM WILL BE LESS THAN OR
C                           EQUAL TO X.
C                IER    - ERROR INDICATOR. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES X IS LESS THAN ZERO
C                           IER = 130 INDICATES DFN AND/OR DFD IS LESS
C                             THAN OR EQUAL TO ZERO
C                           IER = 131 INDICATES THAT AN ERROR OCCURRED
C                             EITHER MLGAMD=DLGAMA OR MLGAMA=ALGAMA.
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/MDBETA,MLGAMA=ALGAMA,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDFDRE (X,DFN,DFD,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               X,DFN,DFD,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               FX,A,B
C                                  TRANSFORM VARIABLES FOR COMPATIBILITY
C                                  WITH BETA PDF
C                                  FIRST EXECUTABLE STATEMENT
      FX = DFD/(DFD + DFN*X)
      A = 0.5 * DFD
      B = 0.5 * DFN
C                                  INVOKE BETA PDF FOR PROBABILITY
      CALL MDBETA(FX,A,B,P,IER)
      IF (IER .NE. 0) GO TO 9000
      P = 1.0 - P
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMDFDRE)
 9005 RETURN
      END

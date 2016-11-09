C   IMSL ROUTINE NAME   - MDSMR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - KOLMOGOROV-SMIRNOV STATISTICS ASYMPTOTIC
C                           PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDSMR (X,P1,P2)
C
C   ARGUMENTS    X      - INPUT CONSTANT TO WHICH INTEGRATION IS
C                           PERFORMED. PROBABILITY THAT A RANDOM
C                           VARIABLE Z IS LESS THAN OR EQUAL TO X IS
C                           CALCULATED, WHERE Z IS DISTRIBUTED SMIRNOV
C                P1     - OUTPUT PROBABILITY (THAT Z IS LESS THAN OR
C                           EQUAL TO X) FOR A ONE SIDED TEST
C                P2     - OUTPUT PROBABILITY (THAT Z IS LESS THAN OR
C                           EQUAL TO X) FOR A TWO SIDED TEST
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDSMR  (X,P1,P2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P1,P2,X
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               CON1,CON2,SMXE,X1,X8
      DATA               SMXE/-180.2182/
      DATA               CON1/2.506628/
      DATA               CON2/9.869604/
C                                  FIRST EXECUTABLE STATEMENT
      P1 = 0.0
      X1 = -X*X
      P2 = X1+X1
      IF (P2.LT.SMXE) GO TO 5
      P1 = EXP(P2)
    5 IF (X.GE..22) GO TO 10
      P2 = 0.0
      GO TO 25
   10 X8 = 8.*X1
      IF (X.GT.0.8) GO TO 15
      P2 = (CON1/X)*EXP(CON2/X8)
      GO TO 25
   15 IF (X.LE.3.15) GO TO 20
      P2 = 1.0
      P1 = 1.0
      GO TO 30
   20 X8 = P1*P1
      X8 = X8*X8
      P2 = P1*X8*X8-X8+P1
      P2 = 1.0-P2-P2
   25 P1 = 1.0-P1
   30 RETURN
      END

C   IMSL ROUTINE NAME   - UERSET
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SET MESSAGE LEVEL FOR IMSL ROUTINE UERTST
C
C   USAGE               - CALL UERSET (LEVEL,LEVOLD)
C
C   ARGUMENTS    LEVEL  - NEW VALUE FOR MESSAGE LEVEL. (INPUT)
C                           OUTPUT FROM IMSL ROUTINE UERTST IS
C                           CONTROLLED SELECTIVELY AS FOLLOWS,
C                             LEVEL = 4 CAUSES ALL MESSAGES TO BE
C                                       PRINTED,
C                             LEVEL = 3 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 32,
C                             LEVEL = 2 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 64,
C                             LEVEL = 1 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 128,
C                             LEVEL = 0 ALL MESSAGE PRINTING IS
C                                       SUPPRESSED.
C                LEVOLD - PREVIOUS MESSAGE LEVEL. (OUTPUT)
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UERSET (LEVEL,LEVOLD)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LEVEL,LEVOLD
C                                  FIRST EXECUTABLE STATEMENT
      LEVOLD = LEVEL
      CALL UERTST (LEVOLD,6HUERSET)
      RETURN
      END

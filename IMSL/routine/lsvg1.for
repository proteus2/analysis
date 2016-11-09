C   IMSL ROUTINE NAME   - LSVG1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LSVG1  (A,B,CS,SN,SIG)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               A,B,CS,SN,SIG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               AA,BB
C                                  FIRST EXECUTABLE STATEMENT
      IF (ABS(A).LE.ABS(B)) GO TO 5
      AA = ABS(A+A)
      SIG = AA*SQRT(0.25+(B/AA)**2)
      CS = A/SIG
      SN = B/SIG
      RETURN
    5 IF (B.EQ.0.0) GO TO 10
      BB = ABS(B+B)
      SIG = BB*SQRT(0.25+(A/BB)**2)
      CS = A/SIG
      SN = B/SIG
      RETURN
   10 SIG = 0.0
      CS = 0.0
      SN = 1.0
      RETURN
      END

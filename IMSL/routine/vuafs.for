C   IMSL ROUTINE NAME   - VUAFS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MATRIX ADDITION (FULL + SYMMETRIC MATRICES)
C
C   USAGE               - CALL VUAFS (A,N,IA,B,C,IC)
C
C   ARGUMENTS    A      - N BY N MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                N      - ORDER OF MATRICES A,B, AND C. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           AN N BY N MATRIX IN SYMMETRIC STORAGE MODE.
C                C      - N BY N MATRIX STORED IN FULL STORAGE MODE
C                           CONTAINING THE SUM C = A+B. (OUTPUT)
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      MATRIX SUBTRACTION MAY BE DONE VIA THIS ROUTINE IF
C                (PRIOR TO ENTRY) THE USER MANIPULATES THE SIGNS OF THE
C                MATRICES TO GIVE THE DESIRED RESULT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VUAFS  (A,N,IA,B,C,IC)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IC
      REAL               A(IA,1),B(1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,L
C                                  FIRST EXECUTABLE STATEMENT
      L = 1
      DO 10 I = 1,N
         DO 5 J = 1,I
            C(I,J) = A(I,J)+B(L)
            IF (I .NE. J) C(J,I) = A(J,I)+B(L)
            L = L+1
    5    CONTINUE
   10 CONTINUE
      RETURN
      END

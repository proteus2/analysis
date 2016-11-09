C   IMSL ROUTINE NAME   - VCVTSF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION OF MATRICES
C                           (SYMMETRIC TO FULL)
C
C   USAGE               - CALL VCVTSF (A,N,B,IB)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N*(N+1)/2 CONTAINING
C                           AN N BY N SYMMETRIC MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE.
C                N      - ORDER OF MATRIX A. (INPUT)
C                B      - OUTPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           MATRIX A IN FULL STORAGE MODE.
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCVTSF (A,N,B,IB)
C
      REAL               A(1),B(IB,1)
C                                  FIRST EXECUTABLE STATEMENT
      I1 = (N*(N+1))/2
      J = N
    5 JP1 = J+1
      DO 10 K = 1,J
         I = JP1-K
         B(I,J) = A(I1)
         I1 = I1-1
   10 CONTINUE
      J = J-1
      IF (J .GE. 1) GO TO 5
      IF (N.LT.2) GO TO 25
      DO 20 I = 2,N
         I1 = I-1
         DO 15 J = 1,I1
            B(I,J) = B(J,I)
   15    CONTINUE
   20 CONTINUE
   25 RETURN
      END

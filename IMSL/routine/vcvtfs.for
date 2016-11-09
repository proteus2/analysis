C   IMSL ROUTINE NAME   - VCVTFS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION OF MATRICES (FULL
C                           TO SYMMETRIC)
C
C   USAGE               - CALL VCVTFS (A,N,IA,B)
C
C   ARGUMENTS    A      - INPUT MATRIX OF DIMENSION N BY N. A CONTAINS
C                           A SYMMETRIC MATRIX STORED IN FULL MODE.
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - OUTPUT VECTOR OF DIMENSION N*(N+1)/2
C                           CONTAINING MATRIX A IN SYMMETRIC STORAGE
C                           MODE.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCVTFS (A,N,IA,B)
C
      REAL               A(IA,1),B(1)
C                                  FIRST EXECUTABLE STATEMENT
      K = 1
      DO 10 I = 1,N
         DO 5 J = 1,I
            B(K) = A(J,I)
            K = K+1
    5    CONTINUE
   10 CONTINUE
      RETURN
      END

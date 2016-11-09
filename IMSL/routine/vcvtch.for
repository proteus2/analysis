C   IMSL ROUTINE NAME   - VCVTCH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION OF MATRICES (FULL
C                           COMPLEX TO HERMITIAN)
C
C   USAGE               - CALL VCVTCH (A,N,IA,H)
C
C   ARGUMENTS    A      - INPUT COMPLEX MATRIX OF DIMENSION N BY N. A
C                           CONTAINS A HERMITIAN MATRIX STORED IN FULL
C                           COMPLEX STORAGE MODE.
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                H      - OUTPUT COMPLEX VECTOR OF LENGTH N*(N+1)/2
C                           CONTAINING MATRIX A IN HERMITIAN STORAGE
C                           MODE. H CONTAINS THE ELEMENTS IN THE
C                           LOWER TRIANGULAR PART OF A STORED ROWWISE.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE DIAGONAL ELEMENTS OF THE OUTPUT MATRIX H WILL BE
C                THE REAL PARTS OF THE DIAGONAL ELEMENTS OF THE INPUT
C                MATRIX A.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCVTCH (A,N,IA,H)
      REAL               T(2)
      COMPLEX            A(IA,N),H(1),W
      EQUIVALENCE        (T(1),W)
C                                  FIRST EXECUTABLE STATEMENT
      K = 1
      DO 15 I = 1,N
         IM1 = I-1
         IF(I .EQ. 1) GO TO 10
         DO 5 J = 1,IM1
            H(K) = CONJG(A(J,I))
            K = K+1
    5    CONTINUE
   10    W = A(I,I)
         H(K) = T(1)
         K = K+1
   15 CONTINUE
      RETURN
      END

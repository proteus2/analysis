C   IMSL ROUTINE NAME   - VCVTHC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION OF MATRICES
C                           (HERMITIAN TO FULL COMPLEX)
C
C   USAGE               - CALL VCVTHC (H,N,B,IB)
C
C   ARGUMENTS    H      - INPUT COMPLEX VECTOR OF LENGTH N*(N+1)/2
C                           CONTAINING AN N BY N HERMITIAN MATRIX
C                           STORED IN HERMITIAN STORAGE MODE. H
C                           CONTAINS THE ELEMENTS IN THE LOWER TRIANGU-
C                           LAR PART OF THE HERMITIAN MATRIX STORED
C                           ROWWISE.
C                N      - ORDER OF MATRIX H. (INPUT)
C                B      - OUTPUT COMPLEX MATRIX OF DIMENSION N BY N
C                           CONTAINING MATRIX H IN FULL COMPLEX STORAGE
C                           MODE.
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
      SUBROUTINE VCVTHC (H,N,B,IB)
C
      REAL               T(2)
      COMPLEX            H(1),B(IB,N),W
      EQUIVALENCE        (T(1),W)
C                                  FIRST EXECUTABLE STATEMENT
      I1 = (N*(N+1))/2
      J = N
    5 JP1 = J+1
      DO 10 K=1,J
         I = JP1-K
         B(I,J) = CONJG(H(I1))
         I1 = I1-1
   10 CONTINUE
      W = B(J,J)
      B(J,J) = T(1)
      J = J-1
      IF(J .GE.1)GO TO 5
      IF(N .LT. 2) GO TO 25
      DO 20 I=2,N
         I1 = I-1
         DO 15 J=1,I1
            B(I,J) = CONJG(B(J,I))
   15    CONTINUE
   20 CONTINUE
   25 RETURN
      END

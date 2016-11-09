C   IMSL ROUTINE NAME   - VMULSF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (SYMMETRIC BY FULL
C                           MATRICES)
C
C   USAGE               - CALL VMULSF (A,N,B,M,IB,C,IC)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING AN
C                           N BY N SYMMETRIC MATRIX STORED IN SYMMETRIC
C                           STORAGE MODE.
C                N      - ORDER OF A AND NUMBER OF ROWS IN B.  (INPUT)
C                B      - N BY M MATRIX IN FULL STORAGE MODE.  (INPUT)
C                M      - NUMBER OF COLUMNS IN B AND NUMBER OF COLUMNS
C                           IN OUTPUT MATRIX C.  (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - N BY M MATRIX IN FULL STORAGE MODE CONTAINING
C                           THE PRODUCT C = A*B.  (OUTPUT)
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
C-----------------------------------------------------------------------
C
      SUBROUTINE VMULSF (A,N,B,M,IB,C,IC)
C
      REAL               A(1),B(IB,M),C(IC,M)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      DO 15 J = 1,M
         IBEG = 1
         DO 10 I = 1,N
            L = IBEG
            KK = 1
            SUM = 0.0D0
            DO 5 K = 1,N
               SUM = SUM + DBLE(A(L))*DBLE(B(K,J))
               IF (K .GE. I) KK = K
               L = L + KK
    5       CONTINUE
            C(I,J) = SUM
            IBEG = IBEG + I
   10    CONTINUE
   15 CONTINUE
      RETURN
      END

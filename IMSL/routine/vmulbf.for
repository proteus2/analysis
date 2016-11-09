C   IMSL ROUTINE NAME   - VMULBF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (BAND BY FULL MATRICES)
C
C   USAGE               - CALL VMULBF (A,IA,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) BAND MATRIX WITH N(2) LOWER
C                           CODIAGONALS AND N(3) UPPER CODIAGONALS. A
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(2)+N(3)+1). (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N(1) BY N(4) MATRIX IN FULL STORAGE MODE.
C                           (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 4 CONTAINING THE
C                           FOLLOWING INFORMATION:
C                           N(1) CONTAINS THE NUMBER OF ROWS IN A,B,
C                             AND C AND THE NUMBER OF COLUMNS IN A.
C                           N(2) CONTAINS THE NUMBER OF LOWER
C                             CODIAGONALS IN A.
C                           N(3) CONTAINS THE NUMBER OF UPPER
C                             CODIAGONALS IN A.
C                           N(4) CONTAINS THE NUMBER OF COLUMNS IN B
C                             AND C.
C                C      - N(1) BY N(4) MATRIX CONTAINING THE PRODUCT
C                           C = A*B. C IS STORED IN FULL STORAGE MODE.
C                           (OUTPUT)
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
      SUBROUTINE VMULBF (A,IA,B,IB,N,C,IC)
C
      INTEGER            N(4)
      REAL               A(IA,1),B(IB,1),C(IC,1)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      NLC = N(2)
      NUC = N(3)
      NCC = NLC+NUC+1
      M = N(4)
      NLCP1 = NLC+1
      NLCP2 = NLC+2
      NMNUC = NN-NUC
      DO 15 JC = 1,M
         NC = NCC
         KK = 1
         DO 10 I = 1,NN
            L = NLCP2-I
            IR = MAX0(L,1)
            IF (L .LE. 0) KK = KK+1
            K = KK
            SUM = 0.0D0
            DO 5 J = IR,NC
               SUM = SUM+DBLE(A(I,J))*DBLE(B(K,JC))
               K = K+1
    5       CONTINUE
            C(I,JC) = SUM
            IF (I .GE. NMNUC) NC = NC-1
   10    CONTINUE
   15 CONTINUE
      RETURN
      END

C   IMSL ROUTINE NAME   - VMULFB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (FULL BY BAND MATRICES)
C
C   USAGE               - CALL VMULFB (A,IA,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(2) FULL MATRIX STORED IN FULL
C                           STORAGE MODE. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N(2) BY N(2) BAND MATRIX WITH N(3) LOWER
C                           CODIAGONALS AND N(4) UPPER CODIAGONALS. B
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(2) BY
C                           (N(3)+N(4)+1). (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 4 CONTAINING THE
C                           FOLLOWING INFORMATION:
C                         N(1) CONTAINS THE NUMBER OF ROWS IN MATRICES
C                           A AND C.
C                         N(2) CONTAINS THE NUMBER OF COLUMNS IN MATRIX
C                           A,B AND C AND THE NUMBER OF ROWS IN MATRIX
C                           B.
C                         N(3) CONTAINS THE NUMBER OF LOWER CODIAGONALS
C                           IN MATRIX B.
C                         N(4) CONTAINS THE NUMBER OF UPPER CODIAGONALS
C                           IN MATRIX B.
C                C      - N(1) BY N(2) MATRIX CONTAINING THE PRODUCT
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
      SUBROUTINE VMULFB (A,IA,B,IB,N,C,IC)
C
      INTEGER            N(4)
      REAL               A(IA,1),B(IB,1),C(IC,1)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NRA = N(1)
      NCA = N(2)
      NLC = N(3)
      NUC = N(4)
      NLCP1 = NLC+1
      NC = NLCP1+NUC
      NCA1 = NCA-NLC
      DO 25 I = 1,NRA
         LCA = 1
         JCB1 = NLCP1
         IRB1 = 1
         ICNT = 0
         DO 20 J = 1,NCA
            IRB = IRB1
            JCB = JCB1
            LC = LCA
            SUM = 0.0D0
            L = MIN0(JCB,NC)
            IF (J .LE. NCA1) GO TO 5
            ICNT = ICNT+1
            L = L-ICNT
    5       DO 10 K = 1,L
               SUM = SUM+DBLE(A(I,LC))*DBLE(B(IRB,JCB))
               JCB = JCB-1
               IRB = IRB+1
               LC = LC+1
   10       CONTINUE
            C(I,J) = SUM
            IF (JCB1 .LT. NC) GO TO 15
            IRB1 = IRB1+1
            LCA = LCA+1
            GO TO 20
   15       JCB1 = JCB1+1
   20    CONTINUE
   25 CONTINUE
      RETURN
      END

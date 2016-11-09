C   IMSL ROUTINE NAME   - VMULSB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (SYMMETRIC BY BAND
C                            MATRICES)
C
C   USAGE               - CALL VMULSB (A,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) SYMMETRIC MATRIX . A IS STORED
C                           IN SYMMETRIC STORAGE MODE AND THEREFORE IS
C                           A VECTOR OF LENGTH N(1)*(N(1)+1)/2. (INPUT)
C                B      - N(1) BY N(1) BAND MATRIX WITH N(2) LOWER
C                           CODIAGONALS AND N(3) UPPER CODIAGONALS. B
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(2)+N(3)+1). (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 3 CONTAINING THE
C                           FOLLOWING INFORMATION:
C                         N(1) CONTAINS THE NUMBER OF ROWS AND COLUMNS
C                           IN MATRICES A,B,AND C.
C                         N(2) CONTAINS THE NUMBER OF LOWER CODIAGONALS
C                           IN MATRIX B.
C                         N(3) CONTAINS THE NUMBER OF UPPER CODIAGONALS
C                           IN MATRIX B.
C                C      - N(1) BY N(1) MATRIX CONTAINING THE PRODUCT
C                           C=A*B. C IS STORED IN FULL STORAGE MODE.
C                           (OUTPUT)
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VMULSB (A,B,IB,N,C,IC)
C
      INTEGER            N(3)
      REAL               A(1),B(IB,1),C(IC,1)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      N3P1 = N(3)+1
      NC = N3P1+N(2)
      NMN2 = NN-N(2)
      IBEG1 = 1
      DO 30 I = 1,NN
         KBEG = 1
         JZ = 1
         DO 25 J = 1,NN
            IF (J .EQ. 1) IBEGP = IBEG1
            SUM = 0.0D0
            KEND = N(2)+J
            IF (J .LE. N3P1) GO TO 15
            IF (JZ .GE. I) GO TO 5
            IBEGP = IBEGP+1
            GO TO 10
    5       IBEGP = IBEGP+JZ
   10       JZ = JZ+1
            KEND = NC
   15       IBEG = IBEGP
            LPD = JZ
            JJ = 1
            KENDP1 = KEND+KBEG
            ICNT = JZ-1
            DO 20 KK = KBEG,KEND
               ICNT = ICNT+1
               K = KENDP1-KK
               SUM = SUM+DBLE(A(IBEG))*DBLE(B(LPD,K))
               LPD = LPD+1
               IF (ICNT .GE. I) JJ = ICNT
               IBEG = IBEG+JJ
   20       CONTINUE
            C(I,J) = SUM
            IF (J .GE. NMN2) KBEG = KBEG+1
   25    CONTINUE
         IBEG1 = IBEG1+I
   30 CONTINUE
      RETURN
      END

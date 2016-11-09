C   IMSL ROUTINE NAME   - VMULBS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - MATRIX MULTIPLICATION (BAND BY SYMMETRIC
C                           MATRICES)
C
C   USAGE               - CALL VMULBS (A,IA,B,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) BAND MATRIX WITH N(2) LOWER
C                           CODIAGONALS AND N(3) UPPER CODIAGONALS. A
C                           IS STORED IN BAND MODE AND THEREFORE HAS
C                           DIMENSION N(1) BY (N(2)+N(3)+1). (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N(1) BY N(1) SYMMETRIC MATRIX. B IS STORED
C                           IN SYMMETRIC STORAGE MODE AND THEREFORE
C                           HAS DIMENSION N(1)*(N(1)+1)/2. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 3 CONTAINING THE
C                           FOLLOWING MATRIX INFORMATION -
C                         N(1) CONTAINS THE NUMBER OF ROWS AND COLUMNS
C                           IN MATRICES A,B,AND C.
C                         N(2) CONTAINS THE NUMBER OF LOWER CODIAGONALS
C                           IN MATRIX A.
C                         N(3) CONTAINS THE NUMBER OF UPPER CODIAGONALS
C                           IN MATRIX A.
C                C      - N(1) BY N(1) MATRIX CONTAINING THE PRODUCT
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
      SUBROUTINE VMULBS (A,IA,B,N,C,IC)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N(3),IA,IC
      REAL               A(IA,1),B(1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      KBEG = 1+N(2)
      KEND = KBEG+N(3)
      JCNT = 0
      N1MN3 = NN-N(3)
      DO 30 I = 1,NN
         IBEG1 = 1
         DO 25 J = 1,NN
            IBEG = IBEG1
            ICNT = 0
            SUM = 0.0D0
            IF(I .LE. N(2)+1) GO TO 15
            IBEG = IBEG+MIN0(JCNT,J-1)
            IF(JCNT .LT. J) GO TO 10
            DO 5 KKK = J,JCNT
               IBEG = IBEG+KKK
    5       CONTINUE
   10       ICNT = JCNT
   15       DO 20 K = KBEG,KEND
               ICNTO = ICNT
               ICNT = ICNT+1
               SUM = SUM+DBLE(A(I,K))*DBLE(B(IBEG))
               IBEG = IBEG+1
               IF(ICNT .GE. J) IBEG = IBEG+ICNTO
   20       CONTINUE
            C(I,J) = SUM
            IBEG1 = IBEG1+J
   25    CONTINUE
         IF(I .GT. N(2)) JCNT = JCNT+1
         IF(I .GE. N1MN3) KEND = KEND-1
         KBEG = MAX0(KBEG-1,1)
   30 CONTINUE
      RETURN
      END

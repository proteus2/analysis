C   IMSL ROUTINE NAME   - VUASB
C
C-----------------------------------------------------------------------
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MATRIX ADDITION (SYMMETRIC + BAND MATRICES)
C
C   USAGE               - CALL VUASB (A,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) SYMMETRIC MATRIX. A IS STORED
C                           IN SYMMETRIC STORAGE MODE AND THEREFORE
C                           IS STORED AS A VECTOR WITH LENGTH
C                           N(1)*(N(1)+1)/2. (INPUT)
C                B      - N(1) BY N(1) BAND MATRIX WITH N(2) LOWER
C                           CODIAGONALS AND N(3) UPPER CODIAGONALS.
C                           B IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(2)+N(3)+1). (INPUT)
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 3 CONTAINING THE
C                           FOLLOWING -
C                         N(1) CONTAINS THE NUMBER OF ROWS AND COLUMNS
C                           IN MATRICES A,B,AND C.
C                         N(2) CONTAINS THE NUMBER OF LOWER CODIAGONALS
C                           IN MATRIX B.
C                         N(3) CONTAINS THE NUMBER OF UPPER CODIAGONALS
C                           IN MATRIX B.
C                C      - N(1) BY N(1) MATRIX CONTAINING THE SUM
C                           C = A+B. C IS STORED IN FULL STORAGE MODE.
C                           (OUTPUT)
C                IC     - ROW DIMENSION OF C EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
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
      SUBROUTINE VUASB  (A,B,IB,N,C,IC)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IB,IC,N(1)
      REAL               A(1),B(IB,1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JC,JJC,J,KBEG,KEND,K,NLCP1,NLC,NMNUC,NN,NUC
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      NLC = N(2)
      NUC = N(3)
      NLCP1 = NLC+1
      KBEG = (NN*(NN+1))/2
      J = NN
    5 KEND = J+1
      DO 10 K = 1,J
         I = KEND-K
         C(I,J) = A(KBEG)
         KBEG = KBEG-1
   10 CONTINUE
      J = J-1
      IF (J .GE. 1) GO TO 5
      IF (NN .LT. 2) GO TO 25
      DO 20 I = 2,NN
         KBEG = I-1
         DO 15 J = 1,KBEG
            C(I,J) = C(J,I)
   15    CONTINUE
   20 CONTINUE
   25 KEND = NLCP1+NUC
      KBEG = NLCP1
      JJC = 1
      NMNUC = NN-NUC
      DO 35 I = 1,NN
         IF (I .GT. NLCP1) JJC = JJC+1
         JC = JJC
         DO 30 J = KBEG,KEND
            C(I,JC) = C(I,JC)+B(I,J)
            JC = JC+1
   30    CONTINUE
         KBEG = KBEG-1
         KBEG = MAX0(KBEG,1)
         IF (I .GE. NMNUC) KEND = KEND-1
   35 CONTINUE
      RETURN
      END

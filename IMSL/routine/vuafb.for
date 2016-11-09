C   IMSL ROUTINE NAME   - VUAFB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MATRIX ADDITION (FULL + BAND MATRICES)
C
C   USAGE               - CALL VUAFB (A,IA,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) MATRIX STORED IN FULL STORAGE
C                           MODE. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N(1) BY N(1) BAND MATRIX WITH N(2) LOWER
C                           CODIAGONALS AND N(3) UPPER CODIAGONALS. B
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(2)+N(3)+1). (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 3 CONTAINING THE
C                           FOLLOWING -
C                         N(1) CONTAINS THE NUMBER OF ROWS AND COLUMNS
C                           IN MATRICES A, B, AND C.
C                         N(2) CONTAINS THE NUMBER OF LOWER CODIAGONALS
C                           IN MATRIX B.
C                         N(3) CONTAINS THE NUMBER OF UPPER CODIAGONALS
C                           IN MATRIX B.
C                C      - N(1) BY N(1) MATRIX CONTAINING THE SUM
C                           C = A+B. C IS STORED IN FULL STORAGE MODE.
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
C   REMARKS      MATRIX SUBTRACTION MAY BE DONE VIA THIS ROUTINE IF
C                (PRIOR TO ENTRY) THE USER MANIPULATES THE SIGNS OF THE
C                MATRICES TO GIVE THE DESIRED RESULT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VUAFB  (A,IA,B,IB,N,C,IC)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,IC,N(1)
      REAL               A(IA,1),B(IB,1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JC,JJC,J,KBEG,KEND,NLCP1,NLC,NMNUC,NN,NUC
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      NLC = N(2)
      NUC = N(3)
      NLCP1 = NLC+1
      KEND = NLCP1+NUC
      DO 10 I = 1,NN
         DO 5 J = 1,NN
            C(I,J) = A(I,J)
    5    CONTINUE
   10 CONTINUE
      KBEG = NLCP1
      JJC = 1
      NMNUC = NN-NUC
      DO 20 I = 1,NN
         IF (I .GT. NLCP1) JJC = JJC+1
         JC = JJC
         DO 15 J = KBEG,KEND
            C(I,JC) = A(I,JC)+B(I,J)
            JC = JC+1
   15    CONTINUE
         KBEG = KBEG-1
         KBEG = MAX0(KBEG,1)
         IF (I .GE. NMNUC) KEND = KEND-1
   20 CONTINUE
      RETURN
      END

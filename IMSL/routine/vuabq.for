C   IMSL ROUTINE NAME   - VUABQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX ADDITION (BAND + BAND SYMMETRIC
C                           MATRICES)
C
C   USAGE               - CALL VUABQ (A,IA,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) BAND MATRIX WITH N(2) LOWER
C                           CODIAGONALS AND N(3) UPPER CODIAGONALS. A
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(2)+N(3)+1). (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N(1) BY N(1) SYMMETRIC BAND MATRIX WITH N(4)
C                           UPPER OR LOWER CODIAGONALS. B IS STORED IN
C                           BAND SYMMETRIC STORAGE MODE AND THEREFORE
C                           HAS DIMENSION N(1) BY (N(4)+1). (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - INPUT/OUTPUT VECTOR OF LENGTH 6 CONTAINING
C                           THE FOLLOWING:
C                           N(1) CONTAINS THE NUMBER OF ROWS AND
C                             COLUMNS IN A,B AND C. (INPUT)
C                           N(2) CONTAINS THE NUMBER OF LOWER
C                             CODIAGONALS IN MATRIX A. (INPUT)
C                           N(3) CONTAINS THE NUMBER OF UPPER
C                             CODIAGONALS IN MATRIX A .(INPUT)
C                           N(4) CONTAINS THE NUMBER OF UPPER OR LOWER
C                             CODIAGONALS IN MATRIX B. (INPUT)
C                           N(5) CONTAINS THE NUMBER OF LOWER
C                             CODIAGONALS IN MATRIX C. (OUTPUT)
C                           N(6) CONTAINS THE NUMBER OF UPPER
C                             CODIAGONALS IN MATRIX C. (OUTPUT)
C                C      - N(1) BY N(1) BAND MATRIX WITH N(5) LOWER
C                           CODIAGONALS AND N(6) UPPER CODIAGONALS. C
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(5)+N(6)+1). (OUTPUT)
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
C                (PRIOR TO ENTRY) THE USER MANIPULATES THE SIGNS OF
C                THE MATRICES TO GIVE THE DESIRED RESULT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VUABQ  (A,IA,B,IB,N,C,IC)
C
      DIMENSION          A(IA,1),B(IB,1),C(IC,1),N(1)
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      NLA = N(2)
      NUA = N(3)
      NCS = N(4)
      N(5) = MAX0(NLA,NCS)
      N(6) = MAX0(NUA,NCS)
      NLC = N(5)
      NUC = N(6)
      NLCP1 = NLC+1
      NCSP1 = NCS+1
      NN1 = NLCP1+NUC
      DO 10 I = 1,NN
         DO 5 J = 1,NN1
            C(I,J) = ZERO
    5    CONTINUE
   10 CONTINUE
C                                  MOVE BAND SYMMETRIC MATRIX TO
C                                  MATRIX C, IN BAND FORM
      DO 15 I = 1,NN
         C(I,NLCP1) = B(I,NCSP1)
   15 CONTINUE
      NN1 = NN+1
      IF (NCS .LE. 0) GO TO 30
      DO 25 I = 1,NCS
         I1 = I+1
         K = 1
         DO 20 J = I1,NN
            C(J,NLCP1-I) = B(J,NCSP1-I)
            C(K,NLCP1+I) = B(J,NCSP1-I)
            K = K+1
   20    CONTINUE
   25 CONTINUE
   30 J = NLA+1
C                                  INSERT ZEROES INTO THE UNUSED
C                                  CORNERS OF THE BAND MATRIX A
      IF (NLA .LE. 0) GO TO 45
      DO 40 I = 1,NLA
         J = J-1
         DO 35  JJ = 1,J
            A(I,JJ) = ZERO
   35    CONTINUE
   40 CONTINUE
   45 IF (NUA .LT. 1) GO TO 60
      K = NN-NUA+1
      NCA = NLA+NUA+1
      J = NCA+1
      DO 55 I = K,NN
         J = J-1
         DO 50 JJ = J,NCA
            A(I,JJ) = ZERO
   50    CONTINUE
   55 CONTINUE
C                                  ADD THE APPROPRIATE PORTIONS OF THE
C                                  MODIFIED MATRICES
   60 K = NLCP1-NLA
      J = NLCP1+NUA
      DO 65  I = 1,NN
         NN1 = 1
         DO 65  I1 = K,J
            C(I,I1) = C(I,I1)+A(I,NN1)
            NN1 = NN1+1
   65    CONTINUE
   70 CONTINUE
      RETURN
      END

C   IMSL ROUTINE NAME   - LEQT2B
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - BAND STORAGE MODE
C                           - HIGH ACCURACY SOLUTION
C
C   USAGE               - CALL LEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,
C                             XL,IER)
C
C   ARGUMENTS    A      - THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B, WHERE A IS ASSUMED TO BE AN N BY N
C                           BAND MATRIX. A IS STORED IN BAND STORAGE
C                           MODE AND THEREFORE HAS DIMENSION N BY
C                           (NLC+NUC+1). (INPUT)
C                N      - ORDER OF A AND THE NUMBER OF ROWS IN B.
C                           (INPUT)
C                NLC    - NUMBER OF LOWER CODIAGONALS IN MATRIX A.
C                           (INPUT)
C                NUC    - NUMBER OF UPPER CODIAGONALS IN MATRIX A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT/OUTPUT MATRIX OF DIMENSION N BY M.
C                         ON INPUT, B CONTAINS THE M RIGHT HAND SIDES
C                           OF THE EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION MATRIX X REPLACES B.
C                           IF IJOB = 1, B IS NOT USED.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES:
C                           I = 0, FACTOR THE MATRIX A AND SOLVE THE
C                             EQUATION AX = B.
C                           I = 1, FACTOR THE MATRIX A.
C                           I = 2, SOLVE THE EQUATION AX = B. THIS
C                             OPTION IMPLIES THAT MATRIX A HAS ALREADY
C                             BEEN FACTORED BY LEQT2B USING
C                             IJOB = 0 OR 1. IN THIS CASE THE OUTPUT
C                             MATRICES U AND XL MUST HAVE BEEN SAVED
C                             FOR REUSE IN THE CALL TO LEQT2B.
C                U      - OUTPUT MATRIX OF DIMENSION N BY (NUC+NLC+3)
C                           CONTAINING MATRIX U OF THE L-U DECOMPOSITION
C                           OF A ROWWISE PERMUTATION OF MATRIX A. U IS
C                           STORED IN BAND STORAGE MODE OCCUPYING THE
C                           FIRST N*(NUC+NLC+1) LOCATIONS. THE REMAINING
C                           2*N LOCATIONS ARE USED AS WORKING STORAGE.
C                IU     - ROW DIMENSION OF MATRIX U EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                XL     - WORK AREA OF DIMENSION N*(NLC+1). THE FIRST
C                           NLC*N LOCATIONS OF XL CONTAIN COMPONENTS OF
C                           THE L-U DECOMPOSITION OF A ROWWISE
C                           PERMUTATION OF A. THE LAST N LOCATIONS
C                           CONTAIN THE PIVOT INDICES. (OUTPUT)
C                IER    - ERROR INDICATOR. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE).
C                           IER=130 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE. THE
C                             MATRIX A IS TOO ILL-CONDITIONED.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1B,UERTST,UGETIO
C                       - DOUBLE/LEQT1B,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
C
      DIMENSION          A(IA,1),U(IU,1),XL(N,1),B(IB,1)
      DOUBLE PRECISION   SUM
      DATA               ZERO/0.0/,ITMAX/50/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NLCP1 = NLC+1
      NCC = NUC+NLCP1
      NLCP2 = NLC+2
      NMNUC = N-NUC
      NU2 = NCC+1
      NU3 = NCC+2
      IF (IJOB .EQ. 2) GO TO 15
C                                  SAVE MATRIX A
      DO 10 I = 1,N
         DO 5 J = 1,NCC
            U(I,J) = A(I,J)
    5    CONTINUE
   10 CONTINUE
C                                  FACTOR MATRIX A
      CALL LEQT1B(U,N,NLC,NUC,IU,B,M,IB,1,XL,IER)
      IF (IER .NE. 0) GO TO 9000
      IF (IJOB .EQ. 1) GO TO 9005
C                                  SAVE THE RIGHT HAND SIDES
   15 DO 60 J = 1,M
         DO 20 I = 1,N
            U(I,NU2) = B(I,J)
   20    CONTINUE
C                                  OBTAIN A SOLUTION
         CALL LEQT1B(U,N,NLC,NUC,IU,U(1,NU2),1,IU,2,XL,IER)
C                                  COMPUTE THE NORM OF THE SOLUTION
         XNORM = ZERO
         DO 25 I = 1,N
            XNORM = AMAX1(XNORM, ABS(U(I,NU2)))
   25    CONTINUE
         IF (XNORM .EQ. ZERO) GO TO 60
C                                  COMPUTE THE RESIDUALS
         DO 45 ITER = 1,ITMAX
            NC = NCC
            KK = 1
            DO 35 I = 1,N
               SUM = B(I,J)
               L = NLCP2-I
               IR = MAX0(L,1)
               IF (L .LE. 0) KK = KK+1
               K = KK
               DO 30 JJ = IR,NC
                  SUM = SUM-DBLE(A(I,JJ))*DBLE(U(K,NU2))
                  K = K+1
   30          CONTINUE
               U(I,NU3) = SUM
               IF (I .GE. NMNUC) NC = NC-1
   35       CONTINUE
            CALL LEQT1B(U,N,NLC,NUC,IU,U(1,NU3),1,IU,2,XL,IER)
            DXNORM = ZERO
C                                  UPDATE THE SOLUTION
            DO 40 I = 1,N
               U(I,NU2) = U(I,NU2)+U(I,NU3)
               DXNORM = AMAX1(DXNORM, ABS(U(I,NU3)))
   40       CONTINUE
            IF (XNORM+DXNORM .EQ. XNORM) GO TO 50
   45    CONTINUE
         IER = 130
C                                  STORE THE SOLUTION
   50    DO 55 JK = 1,N
            B(JK,J) = U(JK,NU2)
   55    CONTINUE
         IF (IER .NE. 0) GO TO 9000
   60 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HLEQT2B)
 9005 RETURN
      END

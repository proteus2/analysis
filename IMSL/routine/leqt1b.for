C   IMSL ROUTINE NAME   - LEQT1B
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - BAND STORAGE
C                           MODE - SPACE ECONOMIZER SOLUTION
C
C   USAGE               - CALL LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,
C                           IER)
C
C   ARGUMENTS    A      - INPUT/OUTPUT MATRIX OF DIMENSION N BY
C                           (NUC+NLC+1). SEE PARAMETER IJOB.
C                N      - ORDER OF MATRIX A AND THE NUMBER OF ROWS IN
C                           B. (INPUT)
C                NLC    - NUMBER OF LOWER CODIAGONALS IN MATRIX A.
C                           (INPUT)
C                NUC    - NUMBER OF UPPER CODIAGONALS IN MATRIX A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT/OUTPUT MATRIX OF DIMENSION N BY M.
C                           ON INPUT, B CONTAINS THE M RIGHT-HAND SIDES
C                           OF THE EQUATION AX = B. ON OUTPUT, THE
C                           SOLUTION MATRIX X REPLACES B. IF IJOB = 1,
C                           B IS NOT USED.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES WHEN
C                           I = 0, FACTOR THE MATRIX A AND SOLVE THE
C                             EQUATION AX = B. ON INPUT, A CONTAINS THE
C                             COEFFICIENT MATRIX OF THE EQUATION AX = B,
C                             WHERE A IS ASSUMED TO BE AN N BY N BAND
C                             MATRIX. A IS STORED IN BAND STORAGE MODE
C                             AND THEREFORE HAS DIMENSION N BY
C                             (NLC+NUC+1). ON OUTPUT, A IS REPLACED
C                             BY THE U MATRIX OF THE L-U DECOMPOSITION
C                             OF A ROWWISE PERMUTATION OF MATRIX A. U
C                             IS STORED IN BAND STORAGE MODE.
C                           I = 1, FACTOR THE MATRIX A. A CONTAINS THE
C                             SAME INPUT/OUTPUT INFORMATION AS IF
C                             IJOB = 0.
C                           I = 2, SOLVE THE EQUATION AX = B. THIS
C                             OPTION IMPLIES THAT LEQT1B HAS ALREADY
C                             BEEN CALLED USING IJOB = 0 OR 1 SO THAT
C                             THE MATRIX A HAS ALREADY BEEN FACTORED.
C                             IN THIS CASE, OUTPUT MATRICES A AND XL
C                             MUST HAVE BEEN SAVED FOR REUSE IN THE
C                             CALL TO LEQT1B.
C                XL     - WORK AREA OF DIMENSION N*(NLC+1). THE FIRST
C                           NLC*N LOCATIONS OF XL CONTAIN COMPONENTS OF
C                           THE L MATRIX OF THE L-U DECOMPOSITION OF A
C                           ROWWISE PERMUTATION OF A. THE LAST N
C                           LOCATIONS CONTAIN THE PIVOT INDICES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE).
C
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
C
      DIMENSION          A(IA,1),XL(N,1),B(IB,1)
      DATA               ZERO/0./,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JBEG = NLC+1
      NLC1 = JBEG
      IF (IJOB .EQ. 2) GO TO 80
      RN = N
C                                  RESTRUCTURE THE MATRIX
C                                  FIND RECIPROCAL OF THE LARGEST
C                                  ABSOLUTE VALUE IN ROW I
      I = 1
      NC = JBEG+NUC
      NN = NC
      JEND = NC
      IF (N .EQ. 1 .OR. NLC .EQ. 0) GO TO 25
    5 K = 1
      P = ZERO
      DO 10 J = JBEG,JEND
         A(I,K) = A(I,J)
         Q =  ABS(A(I,K))
         IF (Q .GT. P) P = Q
         K = K+1
   10 CONTINUE
      IF (P .EQ. ZERO) GO TO 135
      XL(I,NLC1) = ONE/P
      IF (K .GT. NC) GO TO 20
      DO 15 J = K,NC
         A(I,J) = ZERO
   15 CONTINUE
   20 I = I+1
      JBEG = JBEG-1
      IF (JEND-JBEG .EQ. N) JEND = JEND-1
      IF (I .LE. NLC) GO TO 5
      JBEG = I
      NN = JEND
   25 JEND = N-NUC
      DO 40 I = JBEG,N
         P = ZERO
         DO 30 J = 1,NN
            Q =  ABS(A(I,J))
            IF (Q .GT. P) P = Q
   30    CONTINUE
         IF (P .EQ. ZERO) GO TO 135
         XL(I,NLC1) = ONE/P
         IF (I .EQ. JEND) GO TO 37
         IF (I .LT. JEND) GO TO 40
         K = NN+1
         DO 35 J = K,NC
            A(I,J) = ZERO
   35    CONTINUE
   37    NN = NN-1
   40 CONTINUE
      L = NLC
C                                  L-U DECOMPOSITION
      DO 75 K = 1,N
         P =  ABS(A(K,1))*XL(K,NLC1)
         I = K
         IF (L .LT. N) L = L+1
         K1 = K+1
         IF (K1 .GT. L) GO TO 50
         DO 45 J = K1,L
            Q = ABS(A(J,1))*XL(J,NLC1)
            IF (Q .LE. P) GO TO 45
            P = Q
            I = J
   45    CONTINUE
   50    XL(I,NLC1) = XL(K,NLC1)
         XL(K,NLC1) = I
C                                  SINGULARITY FOUND
         Q = RN+P
         IF (Q .EQ. RN) GO TO 135
C                                  INTERCHANGE ROWS I AND K
         IF (K .EQ. I) GO TO 60
         DO 55 J = 1,NC
            P = A(K,J)
            A(K,J) = A(I,J)
            A(I,J) = P
   55    CONTINUE
   60    IF (K1 .GT. L) GO TO 75
         DO 70 I = K1,L
            P = A(I,1)/A(K,1)
            IK = I-K
            XL(K1,IK) = P
            DO 65 J = 2,NC
               A(I,J-1) = A(I,J)-P*A(K,J)
   65    CONTINUE
         A(I,NC) = ZERO
   70    CONTINUE
   75 CONTINUE
      IF (IJOB .EQ. 1) GO TO 9005
C                                  FORWARD SUBSTITUTION
   80 L = NLC
      DO 105 K = 1,N
         I = XL(K,NLC1)
         IF (I .EQ. K) GO TO 90
         DO 85 J = 1,M
            P = B(K,J)
            B(K,J) = B(I,J)
            B(I,J) = P
   85    CONTINUE
   90    IF (L .LT. N) L = L+1
         K1 = K+1
         IF (K1 .GT. L) GO TO 105
         DO 100 I = K1,L
            IK = I-K
            P = XL(K1,IK)
            DO 95 J = 1,M
               B(I,J) = B(I,J)-P*B(K,J)
   95       CONTINUE
  100    CONTINUE
  105 CONTINUE
C                                  BACKWARD SUBSTITUTION
      JBEG = NUC+NLC
      DO 125 J = 1,M
         L = 1
         K1 = N+1
         DO 120 I = 1,N
            K = K1-I
            P = B(K,J)
            IF (L .EQ. 1) GO TO 115
            DO 110 KK = 2,L
               IK = KK+K
               P = P-A(K,KK)*B(IK-1,J)
  110       CONTINUE
  115       B(K,J) = P/A(K,1)
            IF (L .LE. JBEG) L = L+1
  120    CONTINUE
  125 CONTINUE
      GO TO 9005
  135 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HLEQT1B)
 9005 RETURN
      END

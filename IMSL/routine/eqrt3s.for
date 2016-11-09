C   IMSL ROUTINE NAME   - EQRT3S
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - THE SMALLEST (OR LARGEST) EIGENVALUES OF A
C                           TRIDIAGONAL MATRIX IN ALGEBRAIC VALUE WHOSE
C                           SUM EXCEEDS A GIVEN VALUE
C
C   USAGE               - CALL EQRT3S (D,E2,N,VALUE,M,ISW,INFER,IER)
C
C   ARGUMENTS    D      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           DIAGONAL OF THE MATRIX. THE COMPUTED
C                           EIGENVALUES REPLACE THE FIRST M COMPONENTS
C                           OF THE VECTOR D IN NON-DECREASING ORDER,
C                           WHILE THE REMAINING COMPONENTS ARE DES-
C                           TROYED.
C                E2     - INPUT VECTOR OF LENGTH N CONTAINING INFOR-
C                           MATION ABOUT THE OFF-DIAGONAL ELEMENTS
C                           OF THE MATRIX. SEE REMARKS.
C                N      - INPUT ORDER OF THE TRIDIAGONAL MATRIX.
C                VALUE  - INPUT. THE SUM OF THE M EIGENVALUES
C                           FOUND WILL EXCEED ABS(VALUE).
C                M      - OUTPUT NUMBER OF EIGENVALUES FOUND.
C                ISW    - INPUT SCALAR CONTAINING
C                           1 IF THE MATRIX IS KNOWN TO BE POSITIVE
C                             DEFINITE.
C                           0 IF NO SUCH KNOWLEDGE IS AVAILABLE.
C                INFER  - OUTPUT SCALAR CONTAINING ZERO OR THE
C                           VALUE K, INDICATING THAT SUCCESSIVE
C                           ITERATES TO THE K-TH EIGENVALUE WERE
C                           NOT MONOTONE INCREASING.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           TERMINAL ERROR
C                             IER = 129 INDICATES THAT ISW=1 BUT MATRIX
C                               IS NOT POSITIVE DEFINITE
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  AS WRITTEN, THE ROUTINE COMPUTES THE M SMALLEST
C                EIGENVALUES. TO COMPUTE THE LARGEST M EIGENVALUES,
C                REVERSE THE SIGN OF EACH ELEMENT OF D BEFORE AND
C                AFTER CALLING THE ROUTINE. IN THIS CASE, ISW MUST
C                BE SET TO ZERO.
C            2.  LET T BE THE MATRIX FOR WHICH EIGENVALUES ARE BEING
C                COMPUTED. IF T IS SYMMETRIC, ELEMENT I OF THE INPUT
C                VECTOR E2 IS THE SQUARE OF THE SUB-DIAGONAL ELEMENT
C                APPEARING IN ROW I OF THE MATRIX (T(I,I-1)**2). THIS
C                SUBROUTINE CAN ALSO HANDLE NON-SYMMETRIC MATRICES, IN
C                WHICH CASE ELEMENT I OF THE INPUT VECTOR E2 SHOULD BE
C                THE PRODUCT OF THE SUB-DIAGONAL ELEMENT IN ROW I WITH
C                THE SUPER-DIAGONAL ELEMENT IN ROW (I-1) OF THE MATRIX
C                (T(I,I-1)*T(I-1,I)). E2(1) IS NOT USED IN EITHER CASE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQRT3S  (D,E2,N,VALUE,M,ISW,INFER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,ISW,INFER,IER
      REAL               D(1),E2(1),VALUE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NP1,II,I,NM1,NMPK,MP1,J,IMM,JJ
      REAL               DLAM,SUM,RLDLP,DELTA,E,EP,ERR,P,Q,QP,RS,ETA,
     1                   TOT,ZERO,ONE,R,S
      DATA               RLDLP/Z3C100000/
      DATA               ETA/Z00100000/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  RLDLP REPRESENTS MACHINE PRECISION
C                                  FIRST EXECUTABLE STATEMENT
      INFER = 0
      IER = 0
      M = 0
      SUM = ZERO
      S = ZERO
      Q = ZERO
      ERR = ZERO
      E2(1) = ZERO
      DLAM = ZERO
C                                  LOWER BOUND FOR EIGENVALUES FROM
C                                  GERSHGORIN, INITIAL SHIFT
      TOT = D(1)
      NP1 = N+1
      DO 5 II=1,N
         I = NP1-II
         P = Q
         Q = SQRT(E2(I))
         E = D(I)-P-Q
         IF (E.LT.TOT) TOT = E
    5 CONTINUE
      IF (ISW.EQ.1.AND.TOT.LT.ZERO) GO TO 15
      DO 10 I=1,N
         D(I) = D(I)-TOT
   10 CONTINUE
      GO TO 20
   15 TOT = ZERO
   20 M = M+1
C                                  NEXT QR TRANSFORMATION
   25 TOT = TOT+S
      DELTA = D(N)-S
      I = N
      E = ABS(RLDLP*TOT)
      IF (DLAM.LT.E) DLAM = E
      IF (DELTA.GT.DLAM) GO TO 30
      IF (-DELTA.LE.DLAM) GO TO 60
      IER = 129
      GO TO 9000
   30 E = E2(N)/DELTA
      QP = DELTA+E
      P = ONE
      NM1 = N-1
      IF (M.GT.NM1) GO TO 50
      NMPK = NM1+M
      DO 45 II=M,NM1
         I = NMPK-II
         Q = D(I)-S-E
         R = Q/QP
         P = P*R+ONE
C                                  AVOID UNDERFLOW
         IF ((E.GE.1.0).OR.(E.EQ.0.0)) GO TO 35
         IF (R.LT.ETA/E) E = ZERO
   35    EP = E*R
         D(I+1) = QP+EP
         DELTA = Q-EP
         IF (DELTA.GT.DLAM) GO TO 40
         IF (-DELTA.LE.DLAM) GO TO 60
         IER = 129
         GO TO 9000
   40    E = E2(I)/Q
         QP = DELTA+E
         E2(I+1) = QP*EP
   45 CONTINUE
   50 D(M) = QP
      S = QP/P
      IF ((TOT+S).GT.TOT) GO TO 25
      INFER = M
C                                  IRREGULAR END OF ITERATION
C                                  DEFLATE MINIMUM DIAGONAL ELEMENT
      S = ZERO
      I = M
      DELTA = QP
      IF (M.GE.N) GO TO 60
      MP1 = M+1
      DO 55 J=MP1,N
         IF (D(J).GE.DELTA) GO TO 55
         I = J
         DELTA = D(J)
   55 CONTINUE
C                                  CONVERGENCE
   60 IF (I.LT.N) E2(I+1) = E2(I)*E/QP
      IF (I.LE.M) GO TO 70
      IMM = I-M
      DO 65 JJ=1,IMM
         J = I-JJ
         D(J+1) = D(J)-S
         E2(J+1) = E2(J)
   65 CONTINUE
   70 D(M) = TOT
      E2(M) = ERR+ABS(DELTA)
      SUM = SUM+ABS(D(M))
      IF (SUM.GT.ABS(VALUE)) GO TO 9005
      ERR = E2(M)
      IF (M.NE.N) GO TO 20
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HEQRT3S)
 9005 RETURN
      END

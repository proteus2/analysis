C   IMSL ROUTINE NAME   - EQRT1S
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - SMALLEST OR LARGEST M EIGENVALUES OF A
C                           SYMMETRIC TRIDIAGONAL MATRIX
C
C   USAGE               - CALL EQRT1S (D,E2,N,M,ISW,IER)
C
C   ARGUMENTS    D      - INPUT VECTOR OF LENGTH N CONTAINING
C                           THE DIAGONAL ELEMENTS OF THE MATRIX.  THE
C                           COMPUTED EIGENVALUES REPLACE THE FIRST M
C                           COMPONENTS OF THE VECTOR D IN NON-
C                           DECREASING SEQUENCE, WHILE THE REMAINING
C                           COMPONENTS ARE LOST.
C                E2     - INPUT VECTOR OF LENGTH N CONTAINING
C                           THE SQUARES OF THE OFF-DIAGONAL ELEMENTS
C                           OF THE MATRIX.  INPUT E2 IS DESTROYED.
C                N      - INPUT SCALAR CONTAINING THE ORDER OF THE
C                           MATRIX.
C                M      - INPUT SCALAR CONTAINING THE NUMBER OF
C                           SMALLEST EIGENVALUES DESIRED (M IS
C                           LESS THAN OR EQUAL TO N).
C                ISW    - INPUT SCALAR MEANING AS FOLLOWS -
C                           ISW=1 MEANS THAT THE MATRIX IS KNOWN TO BE
C                             POSITIVE DEFINITE.
C                           ISW=0 MEANS THAT THE MATRIX IS NOT KNOWN
C                             TO BE POSITIVE DEFINITE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           WARNING ERROR
C                             IER = 33 INDICATES THAT SUCCESSIVE
C                               ITERATES TO THE K-TH EIGENVALUE WERE NOT
C                               MONOTONE INCREASING. THE VALUE K IS
C                               STORED IN E2(1).
C                           TERMINAL ERROR
C                             IER = 130 INDICATES THAT ISW=1 BUT MATRIX
C                               IS NOT POSITIVE DEFINITE
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      AS WRITTEN, THE ROUTINE COMPUTES THE M SMALLEST
C                EIGENVALUES. TO COMPUTE THE M LARGEST EIGENVALUES,
C                REVERSE THE SIGN OF EACH ELEMENT OF D BEFORE AND
C                AFTER CALLING THE ROUTINE. IN THIS CASE, ISW MUST
C                EQUAL ZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQRT1S (D,E2,N,M,ISW,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,ISW,IER
      REAL               D(N),E2(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II,I,JJ,J,K1,K
      REAL               DELTA,DLAM,EP,ERR,F,P,QP,Q,RLDLP,R,S,TOT
C                                  RLDLP = MACHINE PRECISION
!     DATA               RLDLP/Z3C100000/
      DATA               RLDLP/Z'34000000'/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DLAM = 0.0E0
      ERR = 0.0E0
      S = 0.0E0
C                                  LOOK FOR SMALL SUB-DIAGONAL ENTRIES
C                                  DEFINE INITIAL SHIFT FROM LOWER
C                                  GERSCHGORIN BOUND.
      TOT = D(1)
      Q = 0.0E0
      J = 0
      DO 15 I=1,N
         P = Q
         IF (I.EQ.1) GO TO 5
         IF (P.GT.RLDLP*(ABS(D(I))+ABS(D(I-1)))) GO TO 10
    5    E2(I) = 0.0E0
C                                  COUNT IF E2(I) HAS UNDERFLOWED
   10    IF (E2(I).EQ.0.0E0) J = J+1
         Q = 0.0E0
         IF (I.NE.N) Q = SQRT(E2(I+1))
         TOT = AMIN1(D(I)-P-Q,TOT)
   15 CONTINUE
      IF (ISW.EQ.1.AND.TOT.LT.0.0E0) GO TO 25
      DO 20 I=1,N
   20 D(I) = D(I)-TOT
      GO TO 30
   25 TOT = 0.0E0
   30 DO 90 K=1,M
C                                  NEXT QR TRANSFORMATION
   35    TOT = TOT+S
         DELTA = D(N)-S
         I = N
         F = ABS(RLDLP*TOT)
         IF (DLAM.LT.F) DLAM = F
         IF (DELTA.GT.DLAM) GO TO 40
         IF (DELTA.GE.(-DLAM)) GO TO 75
         IER = 130
         GO TO 9000
C                                  REPLACE SMALL SUB-DIAGONAL SQUARES
C                                  BY ZERO TO REDUCE THE INCIDENCE OF
C                                  UNDERFLOWS
   40    IF (K.EQ.N) GO TO 50
         K1 = K+1
         DO 45 J=K1,N
            IF (E2(J).LE.(RLDLP*(D(J)+D(J-1)))**2) E2(J) = 0.0E0
   45    CONTINUE
   50    F = E2(N)/DELTA
         QP = DELTA+F
         P = 1.0E0
         IF (K.EQ.N) GO TO 65
         K1 = N-K
         DO 60 II=1,K1
            I = N-II
            Q = D(I)-S-F
            R = Q/QP
            P = P*R+1.0E0
            EP = F*R
            D(I+1) = QP+EP
            DELTA = Q-EP
            IF (DELTA.GT.DLAM) GO TO 55
            IF (DELTA.GE.(-DLAM)) GO TO 75
            IER = 130
            GO TO 9000
   55       F = E2(I)/Q
            QP = DELTA+F
            E2(I+1) = QP*EP
   60    CONTINUE
   65    D(K) = QP
         S = QP/P
         IF (TOT+S.GT.TOT) GO TO 35
         IER = 33
         E2(1) = K
C                                  SET ERROR -- IRREGULAR END
C                                  DEFLATE MINIMUM DIAGONAL ELEMENT
         S = 0.0E0
         DELTA = QP
         DO 70 J=K,N
            IF (D(J).GT.DELTA) GO TO 70
            I = J
            DELTA = D(J)
   70    CONTINUE
C                                  CONVERGENCE
   75    IF (I.LT.N) E2(I+1) = E2(I)*F/QP
         IF (I.EQ.K) GO TO 85
         K1 = I-K
         DO 80 JJ=1,K1
            J = I-JJ
            D(J+1) = D(J)-S
            E2(J+1) = E2(J)
   80    CONTINUE
   85    D(K) = TOT
         ERR = ERR+ABS(DELTA)
         E2(K) = ERR
   90 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HEQRT1S)
 9005 RETURN
      END


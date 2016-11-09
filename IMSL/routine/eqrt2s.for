C   IMSL ROUTINE NAME   - EQRT2S
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A SYMMETRIC TRIDIAGONAL MATRIX USING THE
C                           QL METHOD.
C
C   USAGE               - CALL EQRT2S (D,E,N,Z,IZ,IER)
C
C   ARGUMENTS    D      - ON INPUT, THE VECTOR D OF LENGTH N CONTAINS
C                           THE DIAGONAL ELEMENTS OF THE SYMMETRIC
C                           TRIDIAGONAL MATRIX T.
C                           ON OUTPUT, D CONTAINS THE EIGENVALUES OF
C                           T IN ASCENDING ORDER.
C                E      - ON INPUT, THE VECTOR E OF LENGTH N CONTAINS
C                           THE SUB-DIAGONAL ELEMENTS OF T IN POSITION
C                           2,...,N. ON OUTPUT, E IS DESTROYED.
C                N      - ORDER OF TRIDIAGONAL MATRIX T.(INPUT)
C                Z      - ON INPUT, Z CONTAINS THE IDENTITY MATRIX OF
C                           ORDER N.
C                           ON OUTPUT, Z CONTAINS THE EIGENVECTORS
C                           OF T. THE EIGENVECTOR IN COLUMN J OF Z
C                           CORRESPONDS TO THE EIGENVALUE D(J).
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IF IZ IS LESS THAN N, THE
C                           EIGENVECTORS ARE NOT COMPUTED. IN THIS CASE
C                           Z IS NOT USED.
C                IER    - ERROR PARAMETER
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQRT2S FAILED
C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                             AND EIGENVECTORS 1,...,J-1 HAVE BEEN
C                             COMPUTED CORRECTLY, BUT THE EIGENVALUES
C                             ARE UNORDERED.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQRT2S (D,E,N,Z,IZ,IER)
C
      DIMENSION          D(1),E(1),Z(IZ,1)
!     DATA               RDELP/Z3C100000/
      DATA               RDELP/Z'34000000'/
      DATA               ZERO,ONE/0.0,1.0/
C                                  MOVE THE LAST N-1 ELEMENTS
C                                  OF E INTO THE FIRST N-1 LOCATIONS
C                                  FIRST EXECUTABLE STATEMENT
      IER  = 0
      IF (N .EQ. 1) GO TO 9005
      DO 5  I=2,N
         E(I-1) = E(I)
    5 CONTINUE
      E(N) = ZERO
      B = ZERO
      F = ZERO
      DO  60  L=1,N
         J = 0
         H = RDELP*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C                                  LOOK FOR SMALL SUB-DIAGONAL ELEMENT
         DO 10  M=L,N
            K=M
            IF (ABS(E(K)) .LE. B) GO TO 15
   10    CONTINUE
   15    M = K
         IF (M.EQ.L) GO TO 55
   20    IF (J .EQ. 30) GO TO 85
         J = J+1
         L1 = L+1
         G = D(L)
         P = (D(L1)-G)/(E(L)+E(L))
         R = ABS(P)
         IF (RDELP*ABS(P) .LT. 1.0) R = SQRT(P*P+ONE)
         D(L) = E(L)/(P+SIGN(R,P))
         H = G-D(L)
         DO 25 I = L1,N
            D(I) = D(I)-H
   25    CONTINUE
         F = F+H
C                                  QL TRANSFORMATION
         P = D(M)
         C = ONE
         S = ZERO
         MM1 = M-1
         MM1PL = MM1+L
         IF (L.GT.MM1) GO TO 50
         DO 45 II=L,MM1
            I = MM1PL-II
            G = C*E(I)
            H = C*P
            IF (ABS(P).LT.ABS(E(I))) GO TO 30
            C = E(I)/P
            R = SQRT(C*C+ONE)
            E(I+1) = S*P*R
            S = C/R
            C = ONE/R
            GO TO 35
   30       C = P/E(I)
            R = SQRT(C*C+ONE)
            E(I+1) = S*E(I)*R
            S = ONE/R
            C = C*S
   35       P = C*D(I)-S*G
            D(I+1) = H+S*(C*G+S*D(I))
            IF (IZ .LT. N) GO TO 45
C                                  FORM VECTOR
            DO 40 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I)+C*H
               Z(K,I) = C*Z(K,I)-S*H
   40       CONTINUE
   45    CONTINUE
   50    E(L) = S*P
         D(L) = C*P
         IF ( ABS(E(L)) .GT.B) GO TO 20
   55    D(L) = D(L) + F
   60 CONTINUE
C                                  ORDER EIGENVALUES AND EIGENVECTORS
      DO  80  I=1,N
         K = I
         P = D(I)
         IP1 = I+1
         IF (IP1.GT.N) GO TO 70
         DO 65  J=IP1,N
            IF (D(J) .GE. P) GO TO 65
            K = J
            P = D(J)
   65    CONTINUE
   70    IF (K.EQ.I) GO TO 80
         D(K) = D(I)
         D(I) = P
         IF (IZ .LT. N) GO TO 80
         DO 75 J = 1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
   75    CONTINUE
   80 CONTINUE
      GO TO 9005
   85 IER = 128+L
 9000 CONTINUE
      CALL UERTST(IER,6HEQRT2S)
 9005 RETURN
      END


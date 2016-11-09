C   IMSL ROUTINE NAME   - EBALAF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBALAF (A,N,IA,D,K,L)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,K,L
      REAL               A(IA,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L1,K1,K1P1,K11,JJ,J,I,LL,NOCONV
      REAL               R,C,F,G,B,S,B2,ONE,ZERO,P95
      DATA               B/16.0/,B2/256.0/
      DATA               ZERO/0.0/,ONE/1.0/,P95/.95/
C                                  REDUCE NORM A BY DIAGONAL SIMILARITY
C                                  TRANSFORMATION STORED IN D
C                                  FIRST EXECUTABLE STATEMENT
      L1 = 1
      K1 = N
C                                  SEARCH FOR ROWS ISOLATING AN EIGEN-
C                                    VALUE AND PUSH THEM DOWN
    5 K1P1 = K1+1
      IF (K1.LT.1) GO TO 35
      K11=K1
      DO 30 JJ=1,K11
         J = K1P1-JJ
         R = ZERO
         DO 10 I=1,K1
            IF (I.EQ.J) GO TO 10
            R=R+ABS(A(J,I))
   10    CONTINUE
         IF (R.NE.ZERO) GO TO 30
         D(K1) = J
         IF (J.EQ.K1) GO TO 25
         DO 15 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,K1)
            A(I,K1) = F
   15    CONTINUE
         DO 20 I=L1,N
            F = A(J,I)
            A(J,I) = A(K1,I)
            A(K1,I) = F
   20    CONTINUE
   25    K1 = K1-1
         GO TO 5
   30 CONTINUE
C                                  SEARCH FOR COLUMNS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM LEFT
   35 IF (K1.LT.L1) GO TO 65
      LL = L1
      DO 60 J=LL,K1
         C = ZERO
         DO 40 I=L1,K1
            IF (I.EQ.J) GO TO 40
            C = C+ABS(A(I,J))
   40    CONTINUE
         IF (C.NE.ZERO) GO TO 60
         D(L1) = J
         IF (J.EQ.L1) GO TO 55
         DO 45 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,L1)
            A(I,L1) = F
   45    CONTINUE
         DO 50  I=L1,N
            F = A(J,I)
            A(J,I) = A(L1,I)
            A(L1,I) = F
   50    CONTINUE
   55    L1 = L1+1
         GO TO 35
   60 CONTINUE
C                                  NOW BALANCE THE SUBMATRIX IN ROWS
C                                    L1 THROUGH K1
   65 K = L1
      L = K1
      IF (K1.LT.L1) GO TO 75
      DO 70  I=L1,K1
         D(I) = ONE
   70 CONTINUE
   75 NOCONV = 0
      IF (K1.LT.L1) GO TO 120
      DO 115 I=L1,K1
         C = ZERO
         R = ZERO
         DO 80 J=L1,K1
            IF (J.EQ.I) GO TO 80
            C = C+ABS(A(J,I))
            R = R+ABS(A(I,J))
   80    CONTINUE
         G = R/B
         F = ONE
         S = C+R
   85    IF (C.GE.G) GO TO 90
         F = F * B
         C = C*B2
         GO TO 85
   90    G = R*B
   95    IF (C.LT.G) GO TO 100
         F = F/B
         C = C/B2
         GO TO 95
C                                  NOW BALANCE
  100    IF ((C+R)/F.GE.P95*S) GO TO 115
         G = ONE/F
         D(I) = D(I)*F
         NOCONV = 1
         DO 105 J=L1,N
            A(I,J) = A(I,J)*G
  105    CONTINUE
         DO 110 J=1,K1
            A(J,I) = A(J,I)*F
  110    CONTINUE
  115 CONTINUE
  120 IF (NOCONV.EQ.1) GO TO 75
      RETURN
      END

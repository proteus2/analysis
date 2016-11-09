C   IMSL ROUTINE NAME   - EHOUSS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHOUSS (A,N,D,E,E2)
C
      DIMENSION          A(1),D(N),E(N),E2(N)
      REAL               A,D,E,E2,ZERO,H,SCALE,F,G,HH
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      NP1 = N+1
      NN = (N*NP1)/2-1
      NBEG = NN+1-N
      DO 70 II = 1,N
         I = NP1-II
         L = I-1
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 10
C                                  SCALE ROW (ALGOL TOL THEN NOT NEEDED)
         NK = NN
         DO 5 K = 1,L
            SCALE = SCALE+ABS(A(NK))
            NK = NK-1
    5    CONTINUE
         IF (SCALE .NE. ZERO) GO TO 15
   10    E(I) = ZERO
         E2(I) = ZERO
         GO TO 65
   15    NK = NN
         DO 20 K = 1,L
            A(NK) = A(NK)/SCALE
            H = H+A(NK)*A(NK)
            NK = NK-1
   20    CONTINUE
         E2(I) = SCALE*SCALE*H
         F = A(NN)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE*G
         H = H-F*G
         A(NN) = F-G
         IF (L .EQ. 1) GO TO 55
         F = ZERO
         JK1 = 1
         DO 40 J = 1,L
            G = ZERO
            IK = NBEG+1
            JK = JK1
C                                  FORM ELEMENT OF A*U
            DO 25 K = 1,J
               G = G+A(JK)*A(IK)
               JK = JK+1
               IK = IK+1
   25       CONTINUE
            JP1 = J+1
            IF (L .LT. JP1) GO TO 35
            JK = JK+J-1
            DO 30 K = JP1,L
               G = G+A(JK)*A(IK)
               JK = JK+K
               IK = IK+1
   30       CONTINUE
C                                  FORM ELEMENT OF P
   35       E(J) = G/H
            F = F+E(J)*A(NBEG+J)
            JK1 = JK1+J
   40    CONTINUE
         HH = F/(H+H)
C                                  FORM REDUCED A
         JK = 1
         DO 50 J = 1,L
            F = A(NBEG+J)
            G = E(J)-HH*F
            E(J) = G
            DO 45 K = 1,J
               A(JK) = A(JK)-F*E(K)-G*A(NBEG+K)
               JK = JK+1
   45       CONTINUE
   50    CONTINUE
   55    DO 60 K = 1,L
            A(NBEG+K) = SCALE*A(NBEG+K)
   60    CONTINUE
   65    D(I) = A(NBEG+I)
         A(NBEG+I) = H*SCALE*SCALE
         NBEG = NBEG-I+1
         NN = NN-I
   70 CONTINUE
      RETURN
      END


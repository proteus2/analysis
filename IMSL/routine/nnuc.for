C   IMSL ROUTINE NAME   - NNUC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINES
C                           NAK1, NMCC, AND NMKN
C
C   REQD. IMSL ROUTINES - VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NNUC   (X,N,EPS,IR,R,RANK,S,T)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IR(1)
      REAL               X(1),EPS,R(1),RANK(1),S,T
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JDON,JJ,JJJ,J2,K,K1,L,N1
      REAL               Y
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I = 1,N
         R(I) = X(I)
    5 IR(I) = I
C                                  SORT ELEMENTS OF VECTOR R INTO
C                                  ASCENDING SEQUENCE SAVING
C                                  PERMUTATIONS
      CALL VSRTR (R,N,IR)
      S = 0.0
      T = 0.0
      N1 = N-1
      L = 1
   10 DO 30 J = L,N1
         JJ = J
         Y = R(J)
         IF (ABS(Y-R(J+1)) .GT. EPS) GO TO 28
C                                  COUNT THE NUMBER OF TIES
         K = 1
         J2 = J+2
         IF (J2 .GT. N) GO TO 20
         DO 15 I = J2,N
            IF (ABS(Y-R(I)) .GT. EPS) GO TO 20
   15       K = K+1
   20    Y = J+.5*K
      K1 = K+1
      DO 25 I = 1,K1
         JJ = J+I-1
         JJJ = IR(JJ)
   25    RANK(JJJ) = Y
         I = K*(K+1)
         S = S+I
         T = T+I*(K+2)
         GO TO 35
   28 JDON = IR(J)
   30 RANK(JDON) = J
   35 L = JJ+1
      IF (L .LE. N1) GO TO 10
      IF (L .NE. N) GO TO 40
      JDON = IR(N)
      RANK(JDON) = N
   40 CONTINUE
      RETURN
      END

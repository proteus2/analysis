C   IMSL ROUTINE NAME   - LUELPB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ELIMINATION PART OF SOLUTION OF AX=B -
C                           POSITIVE DEFINITE BAND SYMMETRIC MATRIX -
C                           BAND SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL LUELPB (UL,B,N,NC,IA,X)
C
C   ARGUMENTS    UL     - THE RESULT L COMPUTED IN THE ROUTINE LUDAPB
C                           WHERE A = L*L-TRANSPOSE. L IS A LOWER BAND
C                           MATRIX STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N X (NC+1). THE
C                           MAIN DIAGONAL ELEMENTS OF L ARE STORED IN
C                           RECIPROCAL FORM. (INPUT)
C                B      - VECTOR OF LENGTH N CONTAINING THE RIGHT HAND
C                           SIDE OF THE EQUATION AX = B. (INPUT)
C                N      - ORDER OF A AND THE LENGTH OF B AND X. (INPUT)
C                NC     - NUMBER OF LOWER CODIAGONALS OF A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX UL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                X      - VECTOR OF LENGTH N CONTAINING THE SOLUTION TO
C                           THE EQUATION AX = B. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUELPB (UL,B,N,NC,IA,X)
C
      REAL               UL(IA,1),B(1),X(1),ZERO,SUM
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  SOLUTION LY = B
      NC1 = NC+1
      IW = 0
      L = 0
      DO 15 I = 1,N
         SUM = B(I)
         IF (NC .LE. 0) GO TO 10
         IF (IW .EQ. 0) GO TO 9
         L = L+1
         IF (L .GT. NC) L = NC
         K = NC1-L
         KL = I-L
         DO 5 J = K,NC
            SUM = SUM -X(KL) * UL(I,J)
            KL = KL+1
    5    CONTINUE
         GO TO 10
    9    IF (SUM .NE. ZERO) IW = 1
   10    X(I) = SUM*UL(I,NC1)
   15 CONTINUE
C                                  SOLUTION UX = Y
   20 X(N) = X(N)*UL(N,NC1)
      IF (N .LE. 1) GO TO 40
      N1 = N+1
      DO 35 I = 2,N
         K = N1-I
         SUM = X(K)
         IF (NC .LE. 0) GO TO 30
         KL = K+1
         K1 = MIN0(N,K+NC)
         L = 1
         DO 25 J = KL,K1
            SUM = SUM -X(J) * UL(J,NC1-L)
            L = L+1
   25    CONTINUE
   30    X(K) = SUM*UL(K,NC1)
   35 CONTINUE
   40 RETURN
      END

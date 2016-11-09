C   IMSL ROUTINE NAME   - VCVTFQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION (FULL TO BAND
C                           SYMMETRIC STORAGE MODE)
C
C   USAGE               - CALL VCVTFQ (A,N,NC,IA,AA,IAA)
C
C   ARGUMENTS    A      - INPUT. N BY N BAND SYMMETRIC MATRIX STORED
C                           IN FULL STORAGE MODE.
C                N      - INPUT. ORDER OF MATRIX A.
C                NC     - INPUT. NUMBER OF UPPER OR LOWER
C                           CODIAGONALS IN MATRICES A AND AA.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                AA     - OUTPUT. N BY (NC+1) MATRIX CONTAINING A IN
C                           BAND SYMMETRIC STORAGE MODE.
C                IAA    - ROW DIMENSION OF MATRIX AA EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCVTFQ (A,N,NC,IA,AA,IAA)
C
      REAL               A(IA,N),AA(IAA,1),ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  MOVE FIRST NC+1 ROWS
      K = 1
      NCP1 = NC+1
      DO 10 I=1,NCP1
         KK = NCP1
         KP1 = K+1
         DO 5 J=1,K
            LL = KP1-J
            AA(I,KK) = A(I,LL)
            KK = KK-1
    5    CONTINUE
         K = K+1
   10 CONTINUE
C                                  MOVE REMAINING ROWS
      IF(N .EQ. 1) GO TO 35
      KK = 2
      IF (K .GT. N) GO TO 23
      DO 20 I=K,N
         LL = KK
         DO 15 J = 1,NCP1
            AA(I,J) = A(I,LL)
            LL = LL+1
   15    CONTINUE
         KK = KK+1
   20 CONTINUE
   23 IF(NC .EQ. 0) GO TO 35
C                                  ZERO CORNERS OF THE MATRIX
      K = NC
      L = K
      DO 30 I = 1,K
         DO 25 J = 1,L
            AA(I,J) = ZERO
   25    CONTINUE
         L = L-1
   30 CONTINUE
   35 RETURN
      END

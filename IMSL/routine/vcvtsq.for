C   IMSL ROUTINE NAME   - VCVTSQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION (SYMMETRIC TO BAND
C                           SYMMETRIC STORAGE MODE)
C
C   USAGE               - CALL VCVTSQ (A,N,NC,B,IB)
C
C   ARGUMENTS    A      - VECTOR OF LENGTH (N*(N+1))/2 CONTAINING AN
C                           N BY N SYMMETRIC MATRIX IN SYMMETRIC
C                           STORAGE MODE. (INPUT)
C                N      - THE ORDER OF MATRICES A AND B. (INPUT)
C                NC     - THE NUMBER OF UPPER OR LOWER CODIAGONALS IN
C                           MATRICES A AND B. (INPUT)
C                B      - N BY (NC+1) MATRIX CONTAINING A IN BAND
C                           SYMMETRIC STORAGE MODE. (OUTPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
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
      SUBROUTINE VCVTSQ (A,N,NC,B,IB)
C
      REAL               A(1),B(IB,1),ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      NP1 = N+1
      NCP1 = NC+1
      NCP2 = NCP1+1
      NMNC = N-NC
      JJ = NCP1
      KK = (N*(N+1))/2
      K = KK
      DO 10 I=1,N
         II = NP1-I
         DO 5 J=1,JJ
            IJ = NCP2-J
            B(II,IJ) = A(K)
            K = K-1
    5    CONTINUE
         K = KK-II
         KK = K
         IF (I .GE. NMNC) JJ = JJ-1
   10 CONTINUE
      IF (NCP1 .EQ. 1) GO TO 25
C                                  ZERO UPPER LEFT CORNER
      JJ = NC
      DO 20 I = 1,NC
         DO 15 J = 1,JJ
            B(I,J) = ZERO
   15    CONTINUE
         JJ = JJ-1
   20 CONTINUE
   25 RETURN
      END

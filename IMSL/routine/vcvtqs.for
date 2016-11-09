C   IMSL ROUTINE NAME   - VCVTQS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION (BAND SYMMETRIC TO
C                           SYMMETRIC STORAGE MODE)
C
C   USAGE               - CALL VCVTQS (A,N,NC,IA,B)
C
C   ARGUMENTS    A      - INPUT. N BY N BAND SYMMETRIC MATRIX STORED
C                           IN BAND SYMMETRIC STORAGE MODE. MATRIX A
C                           IS DIMENSIONED N BY (NC+1).
C                N      - INPUT. ORDER OF MATRICES A AND B.
C                NC     - INPUT. THE NUMBER OF UPPER OR LOWER
C                           CODIAGONALS IN MATRICES A AND B.
C                IA     - INPUT. ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                B      - OUTPUT. VECTOR OF LENGTH (N*(N+1))/2
C                           CONTAINING A IN SYMMETRIC STORAGE MODE.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCVTQS (A,N,NC,IA,B)
C
      REAL               A(IA,1),B(1),ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      NP1 = N+1
      NCP1 = NC+1
      NCP2 = NCP1+1
      NMNC = N-NC
      KK = (N*(N+1))/2
      K = KK
      JJ = NCP1
      DO 10 I=1,N
         II = NP1-I
         DO 5 J=1,JJ
            IJ = NCP2-J
            B(K) = A(II,IJ)
            K = K-1
    5    CONTINUE
         IF (I .GE. NMNC) JJ = JJ-1
         K = KK-II
         KK = K
   10 CONTINUE
      IF (NCP1 .EQ. N) GO TO 25
      II = (NCP1*NCP2)/2+1
      IJ = NCP2
      JJ = 1
      NMNC = N-NCP1
      DO 20 I=1,NMNC
         K = II
         DO 15 J = 1,JJ
            B(K) = ZERO
            K = K+1
   15    CONTINUE
         JJ = JJ+1
         II = II+IJ
         IJ = IJ+1
   20 CONTINUE
   25 RETURN
      END

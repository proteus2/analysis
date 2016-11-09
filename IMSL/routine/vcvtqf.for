C   IMSL ROUTINE NAME   - VCVTQF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION (BAND SYMMETRIC
C                           TO FULL STORAGE MODE)
C
C   USAGE               - CALL VCVTQF (A,N,NC,IA,B,IB)
C
C   ARGUMENTS    A      - INPUT. N BY N BAND SYMMETRIC MATRIX TO BE
C                           CONVERTED TO FULL STORAGE MODE. A IS
C                           STORED IN BAND SYMMETRIC STORAGE MODE
C                           AND THEREFORE HAS DIMENSION N BY (NC+1).
C                N      - INPUT. ORDER OF MATRICES A AND B.
C                NC     - INPUT. NUMBER OF UPPER OR LOWER
C                           CODIAGONALS IN MATRICES A AND B.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - OUTPUT. N BY N MATRIX CONTAINING A BAND
C                           SYMMETRIC MATRIX IN FULL STORAGE MODE.
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCVTQF (A,N,NC,IA,B,IB)
C
      REAL               A(IA,1),B(IB,N),ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  MOVE FIRST NC+1 ROWS
      NCP1 = NC+1
      L = NCP1
      DO 10 I=1,NCP1
         K = 1
         DO 5 J=L,NCP1
            B(I,K) = A(I,J)
            B(K,I) = B(I,K)
            K = K+1
    5    CONTINUE
         L = L-1
   10 CONTINUE
      IF(NCP1 .EQ. N) GO TO 35
C                                  MOVE REMAINING ROWS
      L = K
      IG = 2
      DO 20 I=K,N
         JJ = NCP1
         LL = L+IG
         DO 15 J=IG,L
            KK = LL-J
            B(I,KK) = A(I,JJ)
            B(KK,I) = B(I,KK)
            JJ = JJ-1
   15    CONTINUE
         L = L+1
         IG = IG+1
   20 CONTINUE
C                                  ZERO CORNERS OF THE OUTPUT MATRIX
      L = NCP1+1
      K = 1
      DO 30 I=L,N
         DO 25 J=1,K
            B(I,J) = ZERO
            B(J,I) = B(I,J)
   25    CONTINUE
         K = K+1
   30 CONTINUE
   35 RETURN
      END

C   IMSL ROUTINE NAME   - VCVTBF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION OF MATRICES (BAND
C                           TO FULL STORAGE MODE)
C
C   USAGE               - CALL VCVTBF (A,N,NUC,NLC,IA,B,IB)
C
C   ARGUMENTS    A      - INPUT. N BY N BAND MATRIX TO BE CONVERTED.
C                           A IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N BY (NLC+1+NUC).
C                N      - INPUT. ORDER OF MATRIX B.
C                NUC    - INPUT. NUMBER OF UPPER CODIAGONALS IN
C                           MATRICES A AND B.
C                NLC    - INPUT. NUMBER OF LOWER CODIAGONALS IN
C                           MATRICES A AND B.
C                IA     - INPUT. ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                B      - OUTPUT. N BY N MATRIX CONTAINING A IN FULL
C                           STORAGE MODE.
C                IB     - INPUT. ROW DIMENSION OF MATRIX B EXACTLY AS
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
      SUBROUTINE VCVTBF (A,N,NUC,NLC,IA,B,IB)
C
      REAL               A(IA,1),B(IB,1),ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  MOVE FIRST NLC+1 ROWS
      NMNUC = N-NUC
      JJ = NLC
      NLCP1 = NLC+1
      NEPT = NMNUC
      L = NUC+1
      DO 10 I=1,NLCP1
         DO 5 J=1,L
            K = JJ+J
            B(I,J) = A(I,K)
    5    CONTINUE
         IF(I .LT. NEPT) L = L+1
         JJ = JJ-1
   10 CONTINUE
      IF(NLCP1 .EQ. N) GO TO 35
C                                  MOVE REMAINING ROWS
      NLCP1 = NLCP1+1
      L = 2
      NEPT = MIN0(NLC+2+NUC,N)
      K = NEPT-1
      DO 20 I=NLCP1,N
         KK = K
         LPNEPT = L+NEPT
         DO 15 J=L,NEPT
            JJ=LPNEPT-J
            B(I,JJ) = A(I,KK)
            KK = KK-1
   15    CONTINUE
         IF(I.GE.NMNUC) K = K-1
         L = L+1
         NEPT = MIN0(NEPT+1,N)
   20 CONTINUE
C                                  ZERO LOWER LEFT CORNER
      K = 1
      DO 30 I=NLCP1,N
         DO 25 J=1,K
            B(I,J) = ZERO
   25    CONTINUE
         K = K+1
   30 CONTINUE
   35 L = NUC+1
      IF(L .EQ. N) GO TO 50
C                                  ZERO UPPER RIGHT CORNER
      K = N-L
      L = L+1
      DO 45 I=1,K
         DO 40 J=L,N
            B(I,J) = ZERO
   40    CONTINUE
         L = L+1
   45 CONTINUE
   50 RETURN
      END

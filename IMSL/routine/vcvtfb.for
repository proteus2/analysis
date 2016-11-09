C   IMSL ROUTINE NAME   - VCVTFB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - STORAGE MODE CONVERSION OF MATRICES (FULL
C                           TO BAND STORAGE MODE)
C
C   USAGE               - CALL VCVTFB (A,N,NUC,NLC,IA,B,IB)
C
C   ARGUMENTS    A      - INPUT. N BY N BAND MATRIX STORED IN FULL
C                           STORAGE MODE.
C                N      - INPUT. ORDER OF MATRIX A.
C                NUC    - INPUT. NUMBER OF UPPER CODIAGONALS IN
C                           MATRICES A AND B.
C                NLC    - INPUT. NUMBER OF LOWER CODIAGONALS IN
C                           MATRICES A AND B.
C                IA     - INPUT. ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                B      - OUTPUT. N BY (NUC+1+NLC) MATRIX CONTAINING
C                           A IN BAND STORAGE MODE.
C                           THE USER SHOULD NOTE THAT NO STORAGE WILL
C                           BE SAVED IF (NLC+1+NUC) IS GREATER THAN OR
C                           EQUAL TO N. IN THE CASE WHERE (NUC+1+NLC)
C                           IS GREATER THAN N, MORE STORAGE WILL BE
C                           REQUIRED TO STORE B THAN WAS REQUIRED TO
C                           STORE A.
C                IB     - INPUT. ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCVTFB (A,N,NUC,NLC,IA,B,IB)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NUC,NLC,IA,IB
      REAL               A(IA,N),B(IB,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NMNUC,NUCP1,NLCP1,KSTOP,KKK,KK,I,LL,J,JJ
      REAL               ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  MOVE FIRST NLC+1 ROWS
      NMNUC = N-NUC
      NUCP1 = NUC+1
      NLCP1 = NLC+1
      KSTOP = NUCP1
      KKK = NUCP1+NLC
      KK = KKK
      DO 10 I=1,NLCP1
         LL = KK
         DO 5 J=1,KSTOP
            JJ = KSTOP+1-J
            B(I,LL) = A(I,JJ)
            LL = LL-1
    5    CONTINUE
         IF(I.GE.NMNUC) KK = KK-1
         KSTOP = MIN0(KSTOP+1,N)
   10 CONTINUE
      IF(NLCP1 .EQ. N) GO TO 25
C                                  MOVE REMAINING ROWS
      JJ = 2
      NLCP1 = NLCP1+1
      KSTOP = MIN0(N,KKK+1)
      DO 20 I=NLCP1,N
         LL = 1
         DO 15 J=JJ,KSTOP
            B(I,LL) = A(I,J)
            LL = LL+1
   15    CONTINUE
         JJ = JJ+1
         KSTOP = MIN0(KSTOP+1,N)
   20 CONTINUE
   25 IF(NLC .EQ. 0) GO TO 45
C                                  ZERO UPPER LEFT CORNER
      KSTOP = NLC
      DO 40 I=1,NLC
         DO 35 J=1,KSTOP
            B(I,J) = ZERO
   35    CONTINUE
         KSTOP = KSTOP-1
   40 CONTINUE
   45 IF(NUC .EQ. 0) GO TO 60
C                                  ZERO LOWER RIGHT CORNER
      KSTOP = KKK
      LL = NMNUC+1
      DO 55 I=LL,N
         DO 50 J=KKK,KSTOP
            B(I,J) = ZERO
   50    CONTINUE
         KKK = KKK-1
   55 CONTINUE
   60 RETURN
      END

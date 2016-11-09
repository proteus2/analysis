C   IMSL ROUTINE NAME   - VMULQQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (BAND SYMMETRIC STORAGE
C                           MODE)
C
C   USAGE               - CALL VMULQQ (A,N,NCA,IA,B,NCB,IB,C,IC)
C
C   ARGUMENTS    A      - N BY N BAND SYMMETRIC MATRIX STORED IN BAND
C                           SYMMETRIC STORAGE MODE.  MATRIX A IS
C                           DIMENSIONED N BY (NCA + 1).  (INPUT)
C                N      - ORDER OF MATRICES A, B, AND C.  (INPUT)
C                NCA    - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX A.  (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N BY N BAND SYMMETRIC MATRIX STORED IN BAND
C                           SYMMETRIC STORAGE MODE.  MATRIX B IS
C                           DIMENSIONED N BY (NCB + 1).  (INPUT)
C                NCB    - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX B.  (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - OUTPUT MATRIX CONTAINING THE PRODUCT C = A*B.
C                           C IS A BAND MATRIX WITH M UPPER AND M LOWER
C                           CO-DIAGONALS, WHERE M = MIN(N-1,NCA+NCB).
C                           C IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N BY
C                           MIN(2*N-1,2*(NCA+NCB)+1).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
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
      SUBROUTINE VMULQQ (A,N,NCA,IA,B,NCB,IB,C,IC)
C
      REAL               A(IA,1),B(IB,1),C(IC,1)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NCAP1 = NCA + 1
      NCAP2 = NCA + 2
      NCBP1 = NCB + 1
      NCBP2 = NCB + 2
      NM1 = N+N-1
      NCC = 2*(NCA+NCB) +1
      NCC = MIN0(NCC,NM1)
      ND = (NCC+3)/2
      DO 65 I = 1,N
         DO 2 K1 = 1,NCC
            C(I,K1) = 0.
    2    CONTINUE
         NDMI = ND-I
         M = MAX0(1,NDMI)
         DO 60 J = 1,N
            JPTA = 0
            JPTB = 0
            LEDZRA = NCAP2
            LEDZRB = NCBP2
C                                  SET INDICES OF A AND B
            IRA = I
            IRB = J
            JCA = MAX0(NCAP2-I,1)
            JCB = MAX0(NCBP2-J,1)
            SUM = 0.D0
C                                  MULTIPLY ACROSS ROW OF A AND DOWN
C                                  COLUMN OF B
            DO 50 K = 1,N
               NB0 = 0
C                                  TEST FOR BEGINNING AND ENDING ZEROES
C                                  IN ROW OF A AND COLUMN OF B
               IF (JCA .EQ. 0) GO TO 52
               IF (I .LT. LEDZRA) GO TO 5
               LEDZRA = LEDZRA + 1
               IF (J .LT. LEDZRB) GO TO 35
               LEDZRB = LEDZRB + 1
               GO TO 50
    5          IF (JCB .EQ. 0) GO TO 52
               IF (J .LT. LEDZRB) GO TO 10
               LEDZRB = LEDZRB + 1
               NB0 = 1
               GO TO 15
   10          SUM = SUM + DBLE(A(IRA,JCA)) * DBLE(B(IRB,JCB))
C                                  INCREMENT OR DECREMENT INDICES OF A
   15          IF (JPTA .NE. 0) GO TO 20
               IF (JCA .LT. NCAP1) GO TO 25
               JPTA = 1
   20          JCA = JCA - 1
               IRA = IRA + 1
               GO TO 30
   25          JCA = JCA + 1
   30          IF (NB0 .EQ. 1) GO TO 50
C                                  INCREMENT OR DECREMENT INDICES OF B
   35          IF (JPTB .NE. 0) GO TO 40
               IF (JCB .LT. NCBP1) GO TO 45
               JPTB = 1
   40          JCB = JCB - 1
               JCB = MAX0(JCB,0)
               IRB = IRB + 1
               GO TO 50
   45          JCB = JCB + 1
   50       CONTINUE
   52       IF (NDMI .GE. 1) GO TO 55
            NDMI = NDMI + 1
            GO TO 60
   55       C(I,M) = SUM
            M = M + 1
            IF (M .GT. NCC) GO TO 65
   60    CONTINUE
   65 CONTINUE
      RETURN
      END

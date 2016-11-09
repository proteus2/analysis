C   IMSL ROUTINE NAME   - VMULSQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (SYMMETRIC BY BAND
C                           SYMMETRIC MATRICES)
C
C   USAGE               - CALL VMULSQ (A,N,B,NC,IB,C,IC)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING AN
C                           N BY N SYMMETRIC MATRIX STORED IN SYMMETRIC
C                           STORAGE MODE.
C                N      - ORDER OF MATRICES A, B, AND C.  (INPUT)
C                B      - N BY N BAND SYMMETRIC MATRIX STORED IN BAND
C                           SYMMETRIC STORAGE MODE.  MATRIX B IS
C                           DIMENSIONED N BY (NC + 1).  (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX B.  (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - N X N MATRIX CONTAINING THE PRODUCT C = A*B
C                           STORED IN FULL STORAGE MODE.  (OUTPUT)
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
      SUBROUTINE VMULSQ (A,N,B,NC,IB,C,IC)
C
      REAL               A(1),B(IB,1),C(IC,N)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NCP1 = NC + 1
      NCP2 = NC + 2
      IRA = 1
      DO 40 J = 1,N
         DO 35 I = 1,N
C                                  SET INDEX OF VECTOR A AND INDICES
C                                  OF MATRIX B
            INDXA = IRA
            IRB = I
            JCB = MAX0(NCP2-I,1)
            JPTB = 0
            KAD = 1
            LEDZRO = NCP2
            SUM = 0.0D0
            DO 25 K = 1,N
C                                  TEST FOR ENDING ZERO IN COLUMN OF
C                                  MATRIX B, AND IF FOUND, STORE PRODUCT
               IF (JCB .EQ. 0) GO TO 30
C                                  TEST FOR LEADING ZERO IN COLUMN OF
C                                  MATRIX B, AND IF FOUND, INCREMENT
C                                  INDEX OF VECTOR A
               IF (K .GE. J) KAD = K
               IF (I .LT. LEDZRO) GO TO 5
               LEDZRO = LEDZRO + 1
               GO TO 20
    5          SUM = SUM + DBLE(A(INDXA)) * DBLE(B(IRB,JCB))
C                                  INCREMENT OR DECREMENT COLUMN INDEX
C                                  OF MATRIX B AND INCREMENT ROW INDEX
C                                  OF MATRIX B
               IF (JPTB .NE. 0) GO TO 10
               IF (JCB .LT. NCP1) GO TO 15
               JPTB = 1
   10          JCB = JCB - 1
               IRB = IRB + 1
               GO TO 20
   15          JCB = JCB + 1
C                                  INCREMENT INDEX OF VECTOR A
   20          INDXA = INDXA + KAD
   25       CONTINUE
   30       C(J,I) = SUM
   35    CONTINUE
C                                  SET INDEX OF VECTOR A TO FIRST
C                                  ELEMENT OF NEW ROW
         IRA = IRA + J
   40 CONTINUE
      RETURN
      END

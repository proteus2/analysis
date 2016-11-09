C   IMSL ROUTINE NAME   - VMULQS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (BAND SYMMETRIC BY
C                           SYMMETRIC MATRICES)
C
C   USAGE               - CALL VMULQS (A,N,NC,IA,B,C,IC)
C
C   ARGUMENTS    A      - N BY N BAND SYMMETRIC MATRIX STORED IN BAND
C                           SYMMETRIC STORAGE MODE.  MATRIX A IS
C                           DIMENSIONED N BY (NC + 1).  (INPUT)
C                N      - ORDER OF MATRICES A,B, AND C.  (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX A.  (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING AN
C                           N BY N SYMMETRIC MATRIX STORED IN SYMMETRIC
C                           STORAGE MODE.
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
      SUBROUTINE VMULQS (A,N,NC,IA,B,C,IC)
C
      REAL               A(IA,1),B(1),C(IC,N)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NCP1 = NC + 1
      NCP2 = NC + 2
      IRB = 1
      DO 40 J = 1,N
         DO 35 I = 1,N
C                                  SET INDEX OF VECTOR B AND INDICES
C                                  OF MATRIX A
            INDXB = IRB
            IRA = I
            JCA = MAX0(NCP2-I,1)
            JPTA = 0
            KAD = 1
            LEDZRO = NCP2
            SUM = 0.0D0
            DO 25 K = 1,N
C                                  TEST FOR ENDING ZERO IN COLUMN OF
C                                  MATRIX A, AND IF FOUND, STORE PRODUCT
               IF (JCA .EQ. 0) GO TO 30
C                                  TEST FOR LEADING ZERO IN COLUMN OF
C                                  MATRIX A, AND IF FOUND, INCREMENT
C                                  INDEX OF VECTOR B
               IF (K .GE. J) KAD = K
               IF (I .LT. LEDZRO) GO TO 5
               LEDZRO = LEDZRO + 1
               GO TO 20
    5          SUM = SUM + DBLE(A(IRA,JCA)) * DBLE(B(INDXB))
C                                  INCREMENT OR DECREMENT COLUMN INDEX
C                                  OF MATRIX A AND INCREMENT ROW INDEX
C                                  OF MATRIX A
               IF (JPTA .NE. 0) GO TO 10
               IF (JCA .LT.NCP1) GO TO 15
               JPTA = 1
   10          JCA = JCA - 1
               IRA = IRA + 1
               GO TO 20
   15          JCA = JCA + 1
C                                  INCREMENT INDEX OF VECTOR B
   20          INDXB = INDXB + KAD
   25       CONTINUE
   30       C(I,J) = SUM
   35    CONTINUE
C                                  SET INDEX OF VECTOR B TO FIRST
C                                  ELEMENT OF NEW ROW
         IRB = IRB + J
   40 CONTINUE
      RETURN
      END

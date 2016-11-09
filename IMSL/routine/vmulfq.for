C   IMSL ROUTINE NAME   - VMULFQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (FULL BY BAND SYMMETRIC
C                           MATRICES)
C
C   USAGE               - CALL VMULFQ (A,M,N,IA,B,NC,IB,C,IC)
C
C   ARGUMENTS    A      - M BY N MATRIX IN FULL STORAGE MODE. (INPUT)
C                M      - NUMBER OF ROWS IN A AND C.  (INPUT)
C                N      - NUMBER OF COLUMNS IN A AND C AND ORDER OF B.
C                           (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N BY N BAND SYMMETRIC MATRIX.  B IS STORED IN
C                           BAND SYMMETRIC STORAGE MODE AND THEREFORE
C                           HAS DIMENSION N BY (NC + 1). (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX B. (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - M BY N MATRIX CONTAINING THE PRODUCT C = A*B.
C                           C IS STORED IN FULL STORAGE MODE.  (OUTPUT)
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
      SUBROUTINE VMULFQ (A,M,N,IA,B,NC,IB,C,IC)
C
      REAL               A(IA,N),B(IB,1),C(IC,N)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NCP1 = NC + 1
      NCP2 = NC + 2
      DO 35 ICC = 1,N
         DO 30 IRA = 1,M
            SUM = 0.0D0
C                                  SET INDICES OF MATRIX B
            IRB = ICC
            JCB = MAX0(NCP2-ICC,1)
            LEDZRO = NCP2
            K = 0
            DO 20 JCA = 1,N
C                                  CHECK FOR ENDING ZERO IN COLUMN OF
C                                  MATRIX B AND, IF FOUND, STORE THE
C                                  ACCUMULATED PRODUCT IN MATRIX C
               IF (JCB .EQ. 0) GO TO 25
C                                  CHECK FOR LEADING ZERO IN COLUMN OF
C                                  MATRIX B AND, IF FOUND, INCREMENT
C                                  COLUMN INDEX OF MATRIX A
               IF (ICC .LT. LEDZRO) GO TO 5
               LEDZRO = LEDZRO + 1
               GO TO 20
    5          SUM = SUM + DBLE(A(IRA,JCA)) * DBLE(B(IRB,JCB))
C
C                                  CHECK FOR INCREMENTING OR DECRE-
C                                  MENTING COLUMN INDEX OF MATRIX B
               IF (K .NE. 0) GO TO 15
               IF (JCB .GE. NCP1) GO TO 10
               JCB = JCB + 1
               GO TO 20
   10          K = 1
   15          JCB = JCB - 1
               IRB = IRB + 1
   20       CONTINUE
   25       C(IRA,ICC) = SUM
   30    CONTINUE
   35 CONTINUE
      RETURN
      END

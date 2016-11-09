C   IMSL ROUTINE NAME   - VMULQF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (BAND SYMMETRIC BY
C                           FULL MATRICES)
C
C   USAGE               - CALL VMULQF (A,M,NC,IA,B,N,IB,C,IC)
C
C   ARGUMENTS    A      - M BY M BAND SYMMETRIC MATRIX.  A IS STORED IN
C                           BAND SYMMETRIC STORAGE MODE AND THERE-
C                           FORE HAS DIMENSION M BY (NC + 1). (INPUT)
C                M      - ORDER OF MATRIX A AND NUMBER OF ROWS IN
C                           MATRICES B AND C.  (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - M BY N MATRIX IN FULL STORAGE MODE. (INPUT)
C                N      - NUMBER OF COLUMNS IN B AND C.  (INPUT)
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
      SUBROUTINE VMULQF (A,M,NC,IA,B,N,IB,C,IC)
C
      REAL               A(IA,1),B(IB,1),C(IC,1)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
      NCP1 = NC + 1
      NCP2 = NC + 2
      DO 35 IRC = 1,M
         DO 30 JCB = 1,N
            SUM = 0.0D0
C                                  SET INDICES OF MATRIX A
            IRA = IRC
            JCA = MAX0(NCP2-IRC,1)
            LEDZRO = NCP2
            K = 0
            DO 20 IRB = 1,M
C                                  CHECK FOR ENDING ZERO IN ROW OF
C                                  MATRIX A AND, IF FOUND, STORE THE
C                                  ACCUMULATED PRODUCT IN MATRIX C
               IF ( JCA .EQ. 0) GO TO 25
C                                  CHECK FOR LEADING ZERO IN ROW OF
C                                  MATRIX A AND, IF FOUND, INCREMENT
C                                  ROW INDEX OF MATRIX B
               IF (IRC .LT. LEDZRO) GO TO 5
               LEDZRO = LEDZRO + 1
               GO TO 20
    5          SUM = SUM + DBLE(A(IRA,JCA)) * DBLE(B(IRB,JCB))
C
C                                  CHECK FOR INCREMENTING OR DECREMENTIN
C                                  COLUMN INDEX OF MATRIX A
               IF (K .NE. 0) GO TO 15
               IF (JCA .GE. NCP1) GO TO 10
               JCA = JCA + 1
               GO TO 20
   10          K = 1
   15          JCA = JCA - 1
               IRA = IRA + 1
   20       CONTINUE
   25       C(IRC,JCB) = SUM
   30    CONTINUE
   35 CONTINUE
      RETURN
      END


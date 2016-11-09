C   IMSL ROUTINE NAME   - LEQ2C
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - COMPLEX MATRIX -
C                           HIGH ACCURACY SOLUTION
C
C   USAGE               - CALL LEQ2C (A,N,IA,B,M,IB,IJOB,WA,WK,IER)
C
C   ARGUMENTS    A      - INPUT COMPLEX MATRIX OF DIMENSION N BY N
C                           CONTAINING THE COMPLEX COEFFICIENTS OF THE
C                           EQUATION AX = B.
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT COMPLEX MATRIX OF DIMENSION N BY M
C                           CONTAINING THE M COMPLEX VALUED RIGHT HAND
C                           SIDES OF THE EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION MATRIX X REPLACES B.
C                           IF IJOB = 1, B IS NOT USED.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES,
C                           I = 0, FACTOR THE MATRIX AND SOLVE THE
C                             EQUATION AX = B.
C                           I = 1, FACTOR THE MATRIX A.
C                           I = 2, SOLVE THE EQUATION AX = B. THIS
C                             OPTION IMPLIES THAT LEQ2C HAS ALREADY
C                             BEEN CALLED USING IJOB = 0 OR 1 SO THAT
C                             THE MATRIX HAS ALREADY BEEN FACTORED. IN
C                             THIS CASE WORK AREAS WA AND WK MUST HAVE
C                             BEEN SAVED FOR REUSE IN THE CALL TO
C                             LEQ2C.
C                WA     - COMPLEX WORK AREA OF LENGTH N*(N+2). THE
C                           FIRST N*N LOCATIONS CONTAIN THE LU
C                           DECOMPOSITION OF A ROWWISE PERMUTATION OF
C                           A. THE REMAINDER OF WA IS USED AS WORK
C                           SPACE.
C                WK     - WORK AREA OF LENGTH N. ON OUTPUT WK CONTAINS
C                           THE PIVOT INDICES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE).
C                           IER=130 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE. THE
C                             MATRIX IS TOO ILL CONDITIONED.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1C,UERTST,UGETIO
C                       - DOUBLE/LEQT1C,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  WHEN IJOB=1, ARGUMENTS B, M AND IB ARE NOT USED BY
C                LEQ2C.
C            2.  WHEN IJOB=0 OR 2, B IS REPLACED WITH THE SOLUTION X.
C            3.  LEQ2C CAN BE USED TO COMPUTE THE INVERSE OF A COMPLEX
C                MATRIX. THIS IS DONE BY CALLING LEQ2C WITH M=N,
C                B=THE N BY N IDENTITY MATRIX AND IJOB=0. WHEN N IS
C                LARGE, IT MAY BE MORE PRACTICAL TO COMPUTE THE INVERSE
C                A COLUMN AT A TIME. TO DO THIS, FIRST CALL LEQ2C WITH
C                IJOB=1 TO FACTOR A. MAKE SUCCEEDING CALLS WITH M=1, B
C                =A COLUMN OF THE IDENTITY MATRIX AND IJOB=2. B WILL BE
C                REPLACED BY THE CORRESPONDING COLUMN OF A INVERSE.
C            4.  THE DETERMINANT OF A CAN BE COMPUTED AFTER LEQ2C HAS
C                BEEN CALLED AS FOLLOWS
C
C                  DET = (1.0,0.0)
C                  DO 5 I = 1,N
C                     IPVT = WK(I)
C                     IF (IPVT .NE. I) DET = -DET
C                     INDX = I + (I-1)*N
C                     DET = DET*WA(INDX)
C                5 CONTINUE
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQ2C  (A,N,IA,B,M,IB,IJOB,WA,WK,IER)
C
      COMPLEX            A(IA,1),B(IB,1),WA(N,1),TEMPA,TEMPB,TEMPC
      REAL               WK(N),TA(2),TB(2),TC(2)
      REAL               AR,AI,BR,BI,CR,CI,DXNORM,XNORM,ZERO
      DOUBLE PRECISION   SUM
      EQUIVALENCE        (TA(1),TEMPA),(TB(1),TEMPB),(TC(1),TEMPC),
     *                   (TA(1),AR),(TA(2),AI),(TB(1),BR),(TB(2),BI),
     *                   (TC(1),CR),(TC(2),CI)
      DATA               ZERO/0.0/
      DATA               ITMAX/50/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      N1 = N+1
      N2 = N+2
      IF (IJOB .EQ. 2) GO TO 15
C                                  SAVE MATRIX A
      DO 10 I = 1,N
         DO 5 J = 1,N
            WA(I,J) = A(I,J)
    5    CONTINUE
   10 CONTINUE
C                                  FACTOR MATRIX A
      CALL LEQT1C (WA,N,N,B,M,IB,1,WK,IER)
      IF (IER .NE. 0) GO TO 9000
      IF (IJOB .EQ. 1) GO TO 9005
C                                  SAVE THE RIGHT HAND SIDES
   15 DO 65 J = 1,M
         DO 20 I = 1,N
            WA(I,N1) = B(I,J)
   20    CONTINUE
C                                  OBTAIN A SOLUTION
         CALL LEQT1C(WA,N,N,WA(1,N1),1,N,2,WK,IER)
C                                  COMPUTE THE NORM OF THE SOLUTION
         XNORM = ZERO
         DO 25 I = 1,N
            TEMPA = WA(I,N1)
            XNORM = AMAX1(XNORM,ABS(AR),ABS(AI))
   25    CONTINUE
         IF (XNORM .EQ. ZERO) GO TO 65
C                                  COMPUTE RESIDUALS
         DO 50 ITER = 1,ITMAX
            DO 40 I = 1,N
               TEMPB = B(I,J)
               SUM = BR
               DO 30 JJ = 1,N
                  TEMPA = A(I,JJ)
                  TEMPB = WA(JJ,N1)
                  SUM = SUM-DBLE(AR)*DBLE(BR)
                  SUM = SUM+DBLE(AI)*DBLE(BI)
   30          CONTINUE
               CR = SUM
               TEMPB = B(I,J)
               SUM = BI
               DO 35 JJ = 1,N
                  TEMPA = A(I,JJ)
                  TEMPB = WA(JJ,N1)
                  SUM = SUM-DBLE(AR)*DBLE(BI)
                  SUM = SUM-DBLE(BR)*DBLE(AI)
   35          CONTINUE
               CI = SUM
               WA(I,N2) = TEMPC
   40       CONTINUE
            CALL LEQT1C(WA,N,N,WA(1,N2),1,N,2,WK,IER)
            DXNORM = ZERO
C                                  UPDATE THE SOLUTION
            DO 45 I = 1,N
               WA(I,N1) = WA(I,N1)+WA(I,N2)
               TEMPA = WA(I,N2)
               DXNORM = AMAX1(DXNORM,ABS(AR),ABS(AI))
   45       CONTINUE
            IF (XNORM+DXNORM .EQ. XNORM) GO TO 55
   50    CONTINUE
         IER = 130
C                                  STORE THE SOLUTION
   55    DO 60 JK = 1,N
            B(JK,J) = WA(JK,N1)
   60    CONTINUE
         IF (IER .NE. 0) GO TO 9000
   65 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HLEQ2C )
 9005 RETURN
      END


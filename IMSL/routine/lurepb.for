C   IMSL ROUTINE NAME   - LUREPB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REFINEMENT OF SOLUTION TO LINEAR EQUATIONS -
C                           POSITIVE DEFINITE BAND SYMMETRIC MATRIX -
C                           BAND SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL LUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,
C                           IER)
C
C   ARGUMENTS    A      - THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B, WHERE A IS AN N X N POSITIVE
C                           DEFINITE BAND SYMMETRIC MATRIX. A IS STORED
C                           IN BAND SYMMETRIC STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N BY (NC+1).(INPUT)
C                N      - ORDER OF THE MATRICES A AND UL. ALSO THE
C                           LENGTH OF VECTORS B,X, AND RES. (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN THE
C                           COEFFICIENT MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                UL     - THE LU DECOMPOSITION AS SUPPLIED BY THE
C                           IMSL ROUTINE LUDAPB. (INPUT)
C                IU     - ROW DIMENSION OF MATRIX UL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT VECTOR OF LENGTH N CONTAINING THE RIGHT
C                           HAND SIDE OF THE EQUATION AX = B.
C                X      - ON INPUT, AN ESTIMATE TO THE SOLUTION OF
C                           AX = B. ON OUTPUT, THE REFINED RESULT.
C                IDGT   - APPROXIMATE NUMBER OF DIGITS IN THE ANSWER
C                           WHICH WERE UNCHANGED AFTER IMPROVEMENT.
C                           (OUTPUT)
C                RES    - VECTOR OF LENGTH N USED AS WORK STORAGE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE.  THE
C                             MATRIX A IS TOO ILL-CONDITIONED.
C
C   REQD. IMSL ROUTINES - SINGLE/LUELPB,UERTST,UGETIO
C                       - DOUBLE/LUELPB,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
C
      REAL               A(IA,1),B(1),UL(IU,1),RES(1),X(1),XNORM,
     *                   RDELP,ZERO,DXNORM
      DOUBLE PRECISION   SUM
      DATA               ITMAX/50/
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NC1 = NC+1
      XNORM = ZERO
      DO 5 I = 1,N
         XNORM = AMAX1(XNORM,ABS(X(I)))
    5 CONTINUE
      IF (XNORM .NE. ZERO) GO TO 10
      IDGT = 50
      GO TO 9005
C                                  ITERATION LOOP
   10 DO 45 ITER = 1,ITMAX
         L = 0
         DO 30 I = 1,N
            K = MAX0(1,NC1-L)
            LK = MAX0(1,I-NC)
            L = I
            SUM = B(I)
            DO 15 J = K,NC1
               SUM = SUM+DBLE(-A(I,J))*DBLE(X(LK))
               LK = LK+1
   15       CONTINUE
            LL = MIN0(NC,N-I)
            IF (LL .LT. 1) GO TO 25
            DO 20 J = 1,LL
               SUM = SUM+DBLE(-A(I+J,NC1-J))*DBLE(X(LK))
               LK = LK+1
   20       CONTINUE
   25       RES(I) = SUM
   30    CONTINUE
         CALL LUELPB(UL,RES,N,NC,IU,RES)
         DXNORM = ZERO
         XNORM = ZERO
         DO 35 I = 1,N
            X(I) = X(I)+RES(I)
            DXNORM = AMAX1(DXNORM,ABS(RES(I)))
            XNORM = AMAX1(XNORM,ABS(X(I)))
   35    CONTINUE
         IF (ITER .NE. 1) GO TO 40
         IF (DXNORM .GT. ZERO) IDGT = -ALOG10(DXNORM/XNORM)
   40    IF (XNORM + DXNORM .EQ. XNORM)GO TO 9005
   45 CONTINUE
C                                  ITERATIVE IMPROVEMENT FAILED
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HLUREPB)
 9005 RETURN
      END

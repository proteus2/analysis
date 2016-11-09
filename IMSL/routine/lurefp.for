
 
C   IMSL ROUTINE NAME   - LUREFP
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM/SINGLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REFINEMENT OF SOLUTION TO LINEAR EQUATIONS -
C                           POSITIVE DEFINITE MATRIX - SYMMETRIC
C                           STORAGE MODE
C
C   USAGE               - CALL LUREFP (A,B,UL,N,X,IDGT,RES,IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE POSITIVE DEFINITE SYMMETRIC COEFFICIENT
C                           MATRIX OF THE EQUATION AX = B. A IS AN N BY
C                           N MATRIX STORED IN SYMMETRIC STORAGE MODE.
C                B      - INPUT VECTOR OF LENGTH N CONTAINING THE RIGHT
C                           HAND SIDE OF THE EQUATION AX = B.
C                UL     - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE LU DECOMPOSITION OF A AS SUPPLIED BY
C                           IMSL ROUTINE LUDECP.
C                N      - ORDER OF MATRICES A AND UL. (INPUT)
C                X      - ON INPUT, VECTOR OF LENGTH N CONTAINING AN
C                           ESTIMATE TO THE SOLUTION OF AX = B.
C                         ON OUTPUT, THE REFINED RESULT REPLACES X.
C                IDGT   - APPROXIMATE NUMBER OF DIGITS IN THE ANSWER
C                           WHICH WERE UNCHANGED AFTER IMPROVEMENT.
C                           (OUTPUT)
C                RES    - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           RESIDUALS FROM THE LAST ITERATION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE. MATRIX
C                             A IS TOO ILL-CONDITIONED.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - SINGLE/LUELMP,UERTST,UGETIO
C                       - DOUBLE/LUELMP,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUREFP (A,B,UL,N,X,IDGT,RES,IER)
C
      REAL               A(1),B(N),UL(1),X(N),RES(N),XNORM,DXNORM,ZERO
      DOUBLE PRECISION   SUM
      DATA               ITMAX/50/
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      XNORM = ZERO
      DO 5 I = 1,N
         XNORM = AMAX1(XNORM,ABS(X(I)))
    5 CONTINUE
      IF (XNORM .NE. ZERO) GO TO 10
      IDGT = 50
      GO TO 9005
C                                  ITERATION LOOP
   10 DO 40 ITER = 1,ITMAX
         L = 1
         DO 25 I = 1,N
            K = L
C                                  IT IS ESSENTIAL THAT A(K)*X(J)
C                                  YIELD A DOUBLE PRECISION RESULT
            SUM = B(I)
            DO 20 J = 1,N
               SUM = SUM-DBLE(A(K))*DBLE(X(J))
               IF (J .GE. I) GO TO 15
               K = K+1
               GO TO 20
   15          K = K+J
   20       CONTINUE
            RES(I) = SUM
            L = L+I
   25    CONTINUE
         CALL LUELMP (UL,RES,N,RES)
         DXNORM = ZERO
         XNORM = ZERO
         DO 30 I = 1,N
            X(I) = X(I)+RES(I)
            DXNORM = AMAX1(DXNORM,ABS(RES(I)))
            XNORM = AMAX1(XNORM,ABS(X(I)))
   30    CONTINUE
         IF (ITER .NE. 1) GO TO 35
         IDGT = 50
         IF (DXNORM .NE. ZERO) IDGT = -ALOG10(DXNORM/XNORM)
   35    IF (XNORM + DXNORM .EQ. XNORM) GO TO 9005
   40 CONTINUE
C                                  ITERATION DID NOT CONVERGE
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HLUREFP)
 9005 RETURN
      END
 
R; T=0.03/0.25 22:59:05
NVERGE
      IER = 129
 9000 CONTINUE
      CALL UE
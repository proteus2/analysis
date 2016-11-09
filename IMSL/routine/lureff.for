C   IMSL ROUTINE NAME   - LUREFF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REFINEMENT OF SOLUTION TO LINEAR EQUATIONS -
C                           FULL STORAGE MODE
C
C   USAGE               - CALL LUREFF (A,B,UL,IPVT,N,IA,X,IDGT,RES,DX,
C                           IER)
C
C   ARGUMENTS    A      - THE COEFFICIENT MATRIX, AX=B, WHERE A
C                           IS N X N. (INPUT)
C                B      - THE RIGHT HAND SIDE, A VECTOR OF SIZE N.
C                           (INPUT)
C                UL     - A GIVEN N X N MATRIX, UL IS THE LU
C                           DECOMPOSITION OF A AS SUPPLIED BY IMSL
C                           ROUTINE LUDATF. (INPUT)
C                IPVT   - A GIVEN VECTOR OF PIVOT INDICES OF SIZE N AS
C                           SUPPLIED BY IMSL ROUTINE LUDATF. (INPUT)
C                N      - ORDER OF A AND UL, AND ALSO THE LENGTH OF
C                           B, IPVT, X, RES, AND DX. (INPUT)
C                IA     - ROW DIMENSION OF A AND UL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                X      - AN INPUT VECTOR OF SIZE N, X IS AN ESTIMATE
C                           TO THE SOLUTION OF AX=B. ON OUTPUT, THE
C                           IMPROVED RESULT OVERWRITES THE INPUT VECTOR
C                           X.
C                IDGT   - APPROXIMATE NUMBER OF DIGITS IN THE ANSWER
C                           WHICH WERE UNCHANGED AFTER IMPROVEMENT.
C                           (OUTPUT)
C                RES    - THE RESIDUAL VECTOR OF SIZE N, USED AS A WORK
C                           VECTOR.
C                DX     - A WORK VECTOR OF SIZE N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES ITERATIVE IMPROVEMENT
C                             FAILED. MATRIX IS TOO ILL CONDITIONED.
C
C   REQD. IMSL ROUTINES - SINGLE/LUELMF,UERTST,UGETIO
C                       - DOUBLE/LUELMF,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUREFF (A,B,UL,IPVT,N,IA,X,IDGT,RES,DX,IER)
C
      DIMENSION          A(IA,1),UL(IA,1),B(1),X(1),RES(1),DX(1),IPVT(1)
      DOUBLE PRECISION   SUM
      DATA               ITMAX/50/,ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      XNORM = ZERO
      DO 10 I=1,N
         XNORM = AMAX1(XNORM,ABS(X(I)))
   10 CONTINUE
      IF (XNORM .NE. ZERO) GO TO 20
      IDGT = 50
      GO TO 9005
   20 DO 45 ITER=1,ITMAX
         DO 30 I=1,N
            SUM = DBLE(B(I))
            DO 25 J=1,N
               SUM = SUM - DBLE(A(I,J)) * DBLE(X(J))
   25       CONTINUE
            RES(I) = SUM
   30    CONTINUE
         CALL LUELMF (UL,RES,IPVT,N,IA,DX)
         DXNORM = ZERO
         XNORM = ZERO
         DO 35 I=1,N
            X(I) = X(I) + DX(I)
            DXNORM = AMAX1(DXNORM,ABS(DX(I)))
            XNORM = AMAX1(XNORM,ABS(X(I)))
   35    CONTINUE
         IF (ITER .NE. 1) GO TO 40
         IDGT = 50
         IF (DXNORM .NE. ZERO) IDGT = -ALOG10(DXNORM/XNORM)
   40    IF (XNORM+DXNORM .EQ. XNORM) GO TO 9005
   45 CONTINUE
C                                  ITERATION DID NOT CONVERGE
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HLUREFF)
 9005 RETURN
      END

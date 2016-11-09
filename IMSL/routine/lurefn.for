C   IMSL ROUTINE NAME   - LUREFN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE LEQT2F
C
C   REQD. IMSL ROUTINES - SINGLE/LUELMN,UERTST,UGETIO
C                       - DOUBLE/LUELMN,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUREFN (A,IA,N,UL,IUL,B,IDGT,APVT,X,RES,DX,IER)
C
      DIMENSION          A(IA,1),UL(IUL,1),B(1),X(1),RES(1),DX(1)
      DIMENSION          APVT(1)
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
         CALL LUELMN (UL,IUL,N,RES,APVT,DX)
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
      CALL UERTST(IER,6HLUREFN)
 9005 RETURN
      END

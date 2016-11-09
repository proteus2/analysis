C   IMSL ROUTINE NAME   - ZSPWC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2
C                       - DOUBLE/VBLA=DNRM2
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSPWC (N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,LR
      REAL               R(LR),DIAG(N),QTB(N),DELTA,X(N),WA1(N),WA2(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JJ,JP1,J,K,L
      REAL               ALPHA,BNORM,EPSMCH,GNORM,ONE,QNORM,SGNORM,
     *                   SPMPAR,SUM,TEMP,ZERO
      REAL               SNRM2
      DATA               SPMPAR /Z3C100000/
      DATA               ONE,ZERO /1.0E0,0.0E0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
C                                  FIRST, CALCULATE THE GAUSS-NEWTON
C                                  DIRECTION.
      JJ = (N*(N+1))/2+1
      DO 25 K=1,N
         J = N-K+1
         JP1 = J+1
         JJ = JJ-K
         L = JJ+1
         SUM = ZERO
         IF (N.LT.JP1) GO TO 10
         DO 5 I=JP1,N
            SUM = SUM+R(L)*X(I)
            L = L+1
    5    CONTINUE
   10    CONTINUE
         TEMP = R(JJ)
         IF (TEMP.NE.ZERO) GO TO 20
         L = J
         DO 15 I=1,J
            TEMP = AMAX1(TEMP,ABS(R(L)))
            L = L+N-I
   15    CONTINUE
         TEMP = EPSMCH*TEMP
         IF (TEMP.EQ.ZERO) TEMP = EPSMCH
   20    CONTINUE
         X(J) = (QTB(J)-SUM)/TEMP
   25 CONTINUE
C                                  TEST WHETHER THE GAUSS-NEWTON
C                                  DIRECTION IS ACCEPTABLE.
      DO 30 J=1,N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   30 CONTINUE
      QNORM = SNRM2(N,WA2,1)
      IF (QNORM.LE.DELTA) GO TO 70
C                                  THE GAUSS-NEWTON DIRECTION IS NOT
C                                  ACCEPTABLE. NEXT, CALCULATE THE
C                                  SCALED GRADIENT DIRECTION.
      L = 1
      DO 40 J=1,N
         TEMP = QTB(J)
         DO 35 I=J,N
            WA1(I) = WA1(I)+R(L)*TEMP
            L = L+1
   35    CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   40 CONTINUE
C                                  CALCULATE THE NORM OF THE SCALED
C                                  GRADIENT AND TEST FOR THE SPECIAL
C                                  CASE IN WHICH THE SCALED GRADIENT IS
C                                  ZERO.
      GNORM = SNRM2(N,WA1,1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (GNORM.EQ.ZERO) GO TO 60
C                                  CALCULATE THE POINT ALONG THE SCALED
C                                  GRADIENT AT WHICH THE QUADRATIC IS
C                                  MINIMIZED.
      DO 45 J=1,N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   45 CONTINUE
      L = 1
      DO 55 J=1,N
         SUM = ZERO
         DO 50 I=J,N
            SUM = SUM+R(L)*WA1(I)
            L = L+1
   50    CONTINUE
         WA2(J) = SUM
   55 CONTINUE
      TEMP = SNRM2(N,WA2,1)
      SGNORM = (GNORM/TEMP)/TEMP
C                                  TEST WHETHER THE SCALED GRADIENT
C                                  DIRECTION IS ACCEPTABLE.
      ALPHA = ZERO
      IF (SGNORM.GE.DELTA) GO TO 60
C                                  THE SCALED GRADIENT DIRECTION IS NOT
C                                  ACCEPTABLE. FINALLY, CALCULATE THE
C                                  POINT ALONG THE DOGLEG AT WHICH THE
C                                  QUADRATIC IS MINIMIZED.
      BNORM = SNRM2(N,QTB,1)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP-(DELTA/QNORM)*(SGNORM/DELTA)**2+SQRT((TEMP-(DELTA
     */QNORM))**2+(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE-(SGNORM/DELTA)**2))/TEMP
   60 CONTINUE
C                                  FORM APPROPRIATE CONVEX COMBINATION
C                                  OF THE GAUSS-NEWTON DIRECTION AND THE
C                                  SCALED GRADIENT DIRECTION.
      TEMP = (ONE-ALPHA)*AMIN1(SGNORM,DELTA)
      DO 65 J=1,N
         X(J) = TEMP*WA1(J)+ALPHA*X(J)
   65 CONTINUE
   70 CONTINUE
      RETURN
      END

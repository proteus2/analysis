C   IMSL ROUTINE NAME   - ZSPWB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSPWB (FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,
     *                   WA2,PAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,LDFJAC,IFLAG,ML,MU
      REAL               X(N),FVEC(N),FJAC(LDFJAC,N),EPSFCN,WA1(N),
     *                   WA2(N),PAR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,MSUM
      REAL               EPSMCH,EPS,H,SPMPAR,TEMP,ZERO
      DATA               SPMPAR /Z3C100000/
      DATA               ZERO /0.0E0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
      EPS = SQRT(AMAX1(EPSFCN,EPSMCH))
      MSUM = ML+MU+1
      IF (MSUM.LT.N) GO TO 20
C                                  COMPUTATION OF DENSE APPROXIMATE
C                                  JACOBIAN.
      DO 10 J=1,N
         TEMP = X(J)
         H = EPS*ABS(TEMP)
         IF (H.EQ.ZERO) H = EPS
         X(J) = TEMP+H
         CALL FCN(X,WA1,N,PAR)
         IF (IFLAG.LT.0) GO TO 15
         X(J) = TEMP
         DO 5 I=1,N
            FJAC(I,J) = (WA1(I)-FVEC(I))/H
    5    CONTINUE
   10 CONTINUE
   15 CONTINUE
      GO TO 50
   20 CONTINUE
C                                  COMPUTATION OF BANDED APPROXIMATE
C                                  JACOBIAN.
      DO 40 K=1,MSUM
         DO 25 J=K,N,MSUM
            WA2(J) = X(J)
            H = EPS*ABS(WA2(J))
            IF (H.EQ.ZERO) H = EPS
            X(J) = WA2(J)+H
   25    CONTINUE
         CALL FCN(X,WA1,N,PAR)
         IF (IFLAG.LT.0) GO TO 45
         DO 35 J=K,N,MSUM
            X(J) = WA2(J)
            H = EPS*ABS(WA2(J))
            IF (H.EQ.ZERO) H = EPS
            DO 30 I=1,N
               FJAC(I,J) = ZERO
               IF (I.GE.J-MU .AND. I.LE.J+ML) FJAC(I,J) =
     *         (WA1(I)-FVEC(I))/H
   30       CONTINUE
   35    CONTINUE
   40 CONTINUE
   45 CONTINUE
   50 CONTINUE
      RETURN
      END

C   IMSL ROUTINE NAME   - ZSCNU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSCNT
C
C   REQD. IMSL ROUTINES - SINGLE/GGUBFS,LEQT2F,LUDATN,LUELMN,LUREFN,
C                           UERTST,UGETIO
C                       - DOUBLE/GGUBFS,LEQT2F,LUDATN,LUELMN,LUREFN,
C                           UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSCNU  (X,N,FCN,NDIGIT,N1,A,Z,Y,XNORM,B,WK,
     1 MAXIT,PAR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,MAXIT,N,NDIGIT,N1
      REAL               A(N1,1),B(1),PAR(1),WK(1),X(1),XNORM(1),
     1                   Y(1),Z(N,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBMAX,IBNORM,IDGT,IEVAL,ITER,J,JER,JI,JS,
     1                   MI,NRS,NSTART
      REAL               BIG,BNORM,CFACT,DX,EPS,HALF,HLMAX,RACC,REPS,
     1                   RRX,RX,SFACT,SMALL,TEST,TN,TR
      DOUBLE PRECISION   DSEED
      DATA               SMALL/Z3C100000/
C                                  CFACT - CONVERGENCE FACTOR FOR
C                                    SETTING NSTART=0
C                                  BIG - LARGEST VALUE FOR EPS
C                                  SMALL - MACHINE EPS
C                                  SFACT - SIGNIFICANT IMPROVEMENT
C                                    FACTOR
C                                  IBMAX - MAX NUMBER OF ITERATIONS
C                                    WITHOUT SIGNIFICANT IMPROVEMENT
C                                  NRS - MAX NUMBER OF RESTARTS WITH
C                                    SAME EPS
C                                  RACC - RELATIVE ACCURACY REQUIREMENT
C                                    (SMALL .LE. RACC .LE. 0.1)
C                                  REPS - NOMINAL EPS FOR STARTING
C                                    PROCEDURE
C                                  HLMAX - MAX VALUE FOR HALF
C
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DSEED = 12345.0D0
      CFACT = 0.99
      BIG = 5.0E5
      SFACT = 0.1
      IBMAX = 50
      NRS = 2
      RACC = AMIN1(AMAX1(SMALL,10.0**(-NDIGIT)),0.1)
      REPS = SQRT(SMALL)
      HLMAX = 3.0
C                                  INITIALIZATION OF OTHER VARIABLES
      ITER = 0
      IEVAL = 0
      NSTART = 0
      IBNORM = 0
      RX = 1.0
      RRX = 0.0
      EPS = 0.0
      DO 5 I=1,N
         B(I) = 0.
         A(N1,I) = 1.0
         Z(I,N1) = X(I)
    5 CONTINUE
      B(N1) = 1.
      A(N1,N1) = 1.0
      JI = N1
      IEVAL = IEVAL+1
      CALL FCN (Z(1,N1),A(1,N1),N,PAR)
      BNORM = 0.0
      DO 10 I=1,N
         BNORM = BNORM+A(I,N1)*A(I,N1)
   10 CONTINUE
      XNORM(N1) = BNORM
C                                  STARTING PROCEDURE TO GENERATE N
C                                    VECTORS IN AN EPS NEIGHBORHOOD OF
C                                    X
   15 IF (NSTART.EQ.NRS) EPS = EPS*10.0
      IF (NSTART.EQ.0) EPS = AMIN1(RX,REPS)
      IF (EPS.GT.BIG) GO TO 120
      NSTART = MOD(NSTART,NRS)+1
      DO 30 J=1,N
         DO 25 I=1,N
   20       TR = (GGUBFS(DSEED)-0.5)*2.0
            IF (ABS(TR).LT.0.1) GO TO 20
            Z(I,J) = X(I)+AMAX1(ABS(X(I)),0.1)*TR*EPS
   25    CONTINUE
         IEVAL = IEVAL+1
         CALL FCN (Z(1,J),A(1,J),N,PAR)
   30 CONTINUE
C                                  CALCULATE N NORMS
      DO 35 J=1,N
         XNORM(J) = 0.
      DO 35 I=1,N
         XNORM(J) = XNORM(J)+A(I,J)*A(I,J)
   35 CONTINUE
C                                  FIND JI SO THAT
C                                    XNORM(JI).LE.XNORM(J) FOR J=1,N1
C                                    FIND JS SO THAT
C                                    XNORM(JS).GE.XNORM(J) FOR J=1,N1
   40 JI = N1
      JS = JI
      DO 45 J=1,N
         IF (XNORM(J).GT.XNORM(JS)) JS = J
         IF (XNORM(J).LT.XNORM(JI)) JI = J
   45 CONTINUE
C                                  TEST FOR F(X)=0
      IF (XNORM(JI).EQ.0.) GO TO 125
C                                  TEST FOR NON-CONVERGENCE
      IF (XNORM(JI).GT.SFACT*BNORM) GO TO 50
      BNORM = XNORM(JI)
      IBNORM = ITER
   50 IF ((ITER-IBNORM).GT.IBMAX) GO TO 120
C                                  INCREMENT AND CHECK ITERATION COUNT
      ITER = ITER+1
      IF (ITER.GE.MAXIT) GO TO 115
C                                  SOLVE FOR MULTIPLIERS
      DO 55 MI=1,N1
         Y(MI) = B(MI)
   55 CONTINUE
      IDGT = 0
      CALL LEQT2F (A,1,N1,N1,Y,IDGT,WK,JER)
C                                  TEST FOR SINGULARITY OF A
      IF (JER.NE.0) GO TO 85
C                                  COMPUTE NEXT APPROXIMATION X
      DO 65 I=1,N
         DX = 0.0
         DO 60 J=1,N1
            DX = DX + Y(J)*Z(I,J)
   60    CONTINUE
         X(I) = DX
   65 CONTINUE
      HALF = 0.
C                                  EVALUATE FUNCTION AT X
   70 IEVAL = IEVAL+1
      CALL FCN (X,Y,N,PAR)
C                                  CALCULATE NORM OF F(X)
      TN = 0.
      DO 75 I=1,N
         TN = TN+Y(I)*Y(I)
   75 CONTINUE
C                                  TEST FOR IMPROVEMENT
      IF (TN.LT.XNORM(JS)) GO TO 95
C                                  EMPLOY HALVING TECHNIQUE
      HALF = HALF+1.
      IF (HALF.GT.HLMAX) GO TO 85
      DO 80 I=1,N
         X(I) = (X(I)+HALF*Z(I,JI))/(HALF+1.0)
   80 CONTINUE
      GO TO 70
C                                  START AGAIN AT BEST APPROXIMATION
   85 IF (JI.EQ.N1) GO TO 15
      XNORM(N1) = XNORM(JI)
      DO 90 I=1,N
         Z(I,N1) = Z(I,JI)
         A(I,N1) = A(I,JI)
   90 CONTINUE
      GO TO 15
C                                  TEST FOR RELATIVE(COMPONENT)
C                                    CONVERGENCE OF APPROXIMATIONS
   95 IF ((HALF.NE.0.).OR.(ITER.EQ.1)) GO TO 105
      RX = SMALL
      DO 100 I=1,N
         RX = AMAX1(RX,ABS(X(I)-Z(I,JI))/AMAX1(ABS(X(I)),0.1))
  100 CONTINUE
      RRX = AMAX1(-ALOG10(RX),0.0)
      IF (RX.LE.RACC) GO TO 125
C                                  REPLACE WORST APPROXIMATION WITH X
  105 IF (TN.LT.CFACT*XNORM(JI)) NSTART = 0
      XNORM(JS) = TN
      DO 110 I=1,N
         Z(I,JS) = X(I)
         A(I,JS) = Y(I)
  110 CONTINUE
      GO TO 40
C                                  SET IER AND RETURN
  115 IER = 129
      GO TO 125
  120 IER = 130
  125 DO 130 I=1,N
         X(I) = Z(I,JI)
  130 CONTINUE
      RETURN
      END

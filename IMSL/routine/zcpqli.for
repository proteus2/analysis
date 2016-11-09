C   IMSL ROUTINE NAME   - ZCPQLI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION ZCPQLI (NN,PT,Q)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NN
      REAL               Q(NN),PT(NN)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,I
      REAL               X
      REAL               XM,F,DX,DF,ZERO,PT1,PTOO5
      DATA               ZERO,PT1,PTOO5/0.0,0.1,0.005/
C                                  FIRST EXECUTABLE STATEMENT
      N = NN-1
C                                  CAUCHY COMPUTES A LOWER BOUND ON THE
C                                    MODULI OF THE ZEROS OF A
C                                    POLYNOMIAL - PT IS THE MODULUS OF
C                                    THE COEFFICIENTS
      PT(NN) = -PT(NN)
C                                  COMPUTE UPPER ESTIMATE OF BOUND
      X = EXP((ALOG(-PT(NN))-ALOG(PT(1)))/N)
      IF (PT(N).EQ.ZERO) GO TO 5
C                                  IF NEWTON STEP AT THE ORIGIN IS
C                                    BETTER, USE IT.
      XM = -PT(NN)/PT(N)
      IF (XM.LT.X) X = XM
C                                  CHOP THE INTERVAL (0,X) UNITL F.LE.0
    5 XM = X*PT1
      F = PT(1)
      DO 10 I=2,NN
         F = F*XM+PT(I)
   10 CONTINUE
      IF (F.LE.ZERO) GO TO 15
      X = XM
      GO TO 5
   15 DX = X
C                                  DO NEWTON ITERATION UNTIL X
C                                    CONVERGES TO TWO DECIMAL PLACES
   20 IF (X.EQ.0.0) GO TO 35
      IF (ABS(DX/X).LE.PTOO5) GO TO 35
      Q(1) = PT(1)
      DO 25 I=2,NN
         Q(I) = Q(I-1)*X+PT(I)
   25 CONTINUE
      F = Q(NN)
      DF = Q(1)
      DO 30 I=2,N
         DF = DF*X+Q(I)
   30 CONTINUE
      DX = F/DF
      X = X-DX
      GO TO 20
   35 ZCPQLI = X
      RETURN
      END

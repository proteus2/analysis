C   IMSL ROUTINE NAME   - ZRPQLD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLD (SSS,NZ,IFLAG)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NZ,IFLAG
      REAL               SSS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,I
      REAL               ARE,EE,ETA,OMP,RMP,RMS,RMRE
      REAL               P(101),QP(101),RK(101),QK(101),SVK(101)
      REAL               SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   PV,RKV,T,S,ZERO,PT001
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
      DATA               ZERO/0.0/,PT001/0.001/
C                                  VARIABLE-SHIFT H POLYNOMIAL
C                                    ITERATION FOR A REAL ZERO SSS -
C                                    STARTING ITERATE
C                                  NZ - NUMBER OF ZERO FOUND
C                                  IFLAG - FLAG TO INDICATE A PAIR OF
C                                    ZEROS NEAR REAL AXIS
C                                  FIRST EXECUTABLE STATEMENT
      NZ = 0
      S = SSS
      IFLAG = 0
      J = 0
C                                  MAIN LOOP
    5 PV = P(1)
C                                  EVALUATE P AT S
      QP(1) = PV
      DO 10 I=2,NN
         PV = PV*S+P(I)
         QP(I) = PV
   10 CONTINUE
      RMP = ABS(PV)
C                                  COMPUTE A RIGOROUS BOUND ON THE
C                                    ERROR IN EVALUATING P
      RMS = ABS(S)
      EE = (RMRE/(ARE+RMRE))*ABS(QP(1))
      DO 15 I=2,NN
   15 EE = EE*RMS+ABS(QP(I))
C                                  ITERATION HAS CONVERGED SUFFICIENTLY
C                                    IF THE POLYNOMIAL VALUE IS LESS
C                                    THAN 20 TIMES THIS BOUND
      IF (RMP.GT.20.*((ARE+RMRE)*EE-RMRE*RMP)) GO TO 20
      NZ = 1
      SZR = S
      SZI = ZERO
      RETURN
   20 J = J+1
C                                  STOP ITERATION AFTER 10 STEPS
      IF (J.GT.10) RETURN
      IF (J.LT.2) GO TO 25
      IF (ABS(T).GT.PT001*ABS(S-T).OR.RMP.LE.OMP) GO TO 25
C                                  A CLUSTER OF ZEROS NEAR THE REAL
C                                    AXIS HAS BEEN ENCOUNTERED RETURN
C                                    WITH IFLAG SET TO INITIATE A
C                                    QUADRATIC ITERATION
      IFLAG = 1
      SSS = S
      RETURN
C                                  RETURN IF THE POLYNOMIAL VALUE HAS
C                                    INCREASED SIGNIFICANTLY
   25 OMP = RMP
C                                  COMPUTE T, THE NEXT POLYNOMIAL, AND
C                                    THE NEW ITERATE
      RKV = RK(1)
      QK(1) = RKV
      DO 30 I=2,N
         RKV = RKV*S+RK(I)
         QK(I) = RKV
   30 CONTINUE
      IF (ABS(RKV).LE.ABS(RK(N))*10.*ETA) GO TO 40
C                                  USE THE SCALED FORM OF THE
C                                    RECURRENCE IF THE VALUE OF K AT S
C                                    IS NONZERO
      T = -PV/RKV
      RK(1) = QP(1)
      DO 35 I=2,N
   35 RK(I) = T*QK(I-1)+QP(I)
      GO TO 50
C                                  USE UNSCALED FORM
   40 RK(1) = ZERO
      DO 45 I=2,N
   45 RK(I) = QK(I-1)
   50 RKV = RK(1)
      DO 55 I=2,N
   55 RKV = RKV*S+RK(I)
      T = ZERO
      IF (ABS(RKV).GT.ABS(RK(N))*10.*ETA) T = -PV/RKV
      S = S+T
      GO TO 5
      END

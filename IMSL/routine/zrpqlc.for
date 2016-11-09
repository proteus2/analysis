C   IMSL ROUTINE NAME   - ZRPQLC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   REQD. IMSL ROUTINES - ZRPQLE,ZRPQLF,ZRPQLG,ZRPQLH,ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLC (UU,VV,NZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NZ
      REAL               UU,VV
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,I,ITYPE
      REAL               ARE,EE,ETA,OMP,RELSTP,RMP,RMRE,T,ZM
      REAL               P(101),QP(101),RK(101),QK(101),SVK(101)
      REAL               SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   UI,VI,ZERO,PT01,ONE
      LOGICAL            TRIED
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
      DATA               ZERO,PT01,ONE/0.0,0.01,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      NZ = 0
C                                  VARIABLE-SHIFT K-POLYNOMIAL
C                                    ITERATION FOR A QUADRATIC FACTOR
C                                    CONVERGES ONLY IF THE ZEROS ARE
C                                    EQUIMODULAR OR NEARLY SO
C                                  UU,VV - COEFFICIENTS OF STARTING
C                                    QUADRATIC
C                                  NZ - NUMBER OF ZERO FOUND
      TRIED = .FALSE.
      U = UU
      V = VV
      J = 0
C                                  MAIN LOOP
    5 CALL ZRPQLI (ONE,U,V,SZR,SZI,RLZR,RLZI)
C                                  RETURN IF ROOTS OF THE QUADRATIC ARE
C                                    REAL AND NOT CLOSE TO MULTIPLE OR
C                                    NEARLY EQUAL AND OF OPPOSITE SIGN
      IF ( ABS(ABS(SZR)-ABS(RLZR)).GT.PT01*ABS(RLZR)) RETURN
C                                  EVALUATE POLYNOMIAL BY QUADRATIC
C                                    SYNTHETIC DIVISION
      CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
      RMP = ABS(RA-SZR*RB)+ABS(SZI*RB)
C                                  COMPUTE A RIGOROUS BOUND ON THE
C                                    ROUNDING ERROR IN EVALUTING P
      ZM = SQRT(ABS(V))
      EE = 2.*ABS(QP(1))
      T = -SZR*RB
      DO 10 I=2,N
   10 EE = EE*ZM+ABS(QP(I))
      EE = EE*ZM+ABS(RA+T)
      EE = (5.*RMRE+4.*ARE)*EE-(5.*RMRE+2.*ARE)*(ABS(RA+T)+
     1     ABS(RB)*ZM)+2.*ARE*ABS(T)
C                                  ITERATION HAS CONVERGED SUFFICIENTLY
C                                    IF THE POLYNOMIAL VALUE IS LESS
C                                    THAN 20 TIMES THIS BOUND
      IF (RMP.GT.20.*EE) GO TO 15
      NZ = 2
      RETURN
   15 J = J+1
C                                  STOP ITERATION AFTER 20 STEPS
      IF (J.GT.20) RETURN
      IF (J.LT.2) GO TO 25
      IF (RELSTP.GT..01.OR.RMP.LT.OMP.OR.TRIED) GO TO 25
C                                  A CLUSTER APPEARS TO BE STALLING THE
C                                    CONVERGENCE. FIVE FIXED SHIFT
C                                    STEPS ARE TAKEN WITH A U,V CLOSE
C                                    TO THE CLUSTER
      IF (RELSTP.LT.ETA) RELSTP = ETA
      RELSTP = SQRT(RELSTP)
      U = U-U*RELSTP
      V = V+V*RELSTP
      CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
      DO 20 I=1,5
         CALL ZRPQLE (ITYPE)
         CALL ZRPQLF (ITYPE)
   20 CONTINUE
      TRIED = .TRUE.
      J = 0
   25 OMP = RMP
C                                  CALCULATE NEXT K POLYNOMIAL AND NEW
C                                    U AND V
      CALL ZRPQLE (ITYPE)
      CALL ZRPQLF (ITYPE)
      CALL ZRPQLE (ITYPE)
      CALL ZRPQLG (ITYPE,UI,VI)
C                                  IF VI IS ZERO THE ITERATION IS NOT
C                                    CONVERGING
      IF (VI.EQ.ZERO) RETURN
      RELSTP = ABS((VI-V)/VI)
      U = UI
      V = VI
      GO TO 5
      END

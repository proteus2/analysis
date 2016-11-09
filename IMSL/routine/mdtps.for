C   IMSL ROUTINE NAME   - MDTPS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - CUMULATIVE PROBABILITY AND, OPTIONALLY,
C                           INDIVIDUAL TERMS OF THE POISSON PROBABILITY
C                           DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDTPS (K,RLAM,IOPT,T,P)
C
C   ARGUMENTS    K      - INPUT VALUE. THE SUM OF THE POISSON IS TO AND
C                           INCLUDING K.
C                RLAM   - INPUT PARAMETER OF THE POISSON. RLAM MUST BE
C                           POSITIVE.
C                IOPT   - INPUT OPTION PARAMETER. IOPT = 1 IMPLIES
C                           INDIVIDUAL TERMS ARE DESIRED.
C                T      - OUTPUT VECTOR OF LENGTH K+1 CONTAINING THE
C                           INDIVIDUAL TERMS (PROBABILITIES) FROM THE
C                           POISSON DISTRIBUTION CORRESPONDING TO
C                           POISSON RANDOM VARIABLE VALUES 0,1,...,K.
C                P      - OUTPUT CUMULATIVE PROBABILITY FROM THE
C                           POISSON PROBABILITY DISTRIBUTION. SUM OF
C                           THE T VECTOR.
C
C   REQD. IMSL ROUTINES - H32/MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDTPS (K,RLAM,IOPT,T,P)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,IOPT
      REAL               T(1),RLAM,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K1,I,JJ,ICNT,KCNT,J
      DOUBLE PRECISION   X,Y,P1,EPS,ALEPS,X2,Y2,G,H,P2,PE,TEMP,DEPS
      DOUBLE PRECISION   DLGAMA
      DATA               EPS/1.D-78/,ALEPS/-179.6016372535356D0/
      DATA               DEPS/Z3410000000000000/
C                                  FIRST EXECUTABLE STATEMENT
      K1 = K+1
      IF (IOPT.NE.1) GO TO 10
      DO 5 I=1,K1
         T(I) = 0.
    5 CONTINUE
C                                  SPECIAL CASE - LAMBDA = 0.
   10 IF (RLAM.GT.DEPS) GO TO 15
      PE = 1.D0
      IF (IOPT.EQ.1) T(1) = 1.0
      GO TO 60
C                                  INITIALIZATION FOR FORWARD
C                                  CALCULATION
   15 X = RLAM
      Y = 1.D0
      JJ = 1
      P1 = -RLAM
      ICNT = P1/ALEPS
      P1 = P1-ICNT*ALEPS
      P1 = DEXP(P1)
C                                  INITIALIZATION FOR BACKWARD
C                                  CALCULATION
      X2 = K
      Y2 = RLAM
      G = X2*DLOG(Y2)
      H = K1
      H = DLGAMA(H)
      P2 = -Y2+G-H
      KCNT = P2/ALEPS
      P2 = P2-KCNT*ALEPS
      P2 = DEXP(P2)
      G = 1.D0
      H = 1.D0
      IF (ICNT.EQ.0) G = 1.D0-P1
      IF (KCNT.EQ.0) H = 1.D0-P2
      PE = 0.D0
C                                  DETERMINE AT WHICH END CALCULATION
C                                  IS OCCURRING
   20 J = ICNT-KCNT
      IF (J) 45,25,30
   25 IF (P1.GT.P2) GO TO 45
C                                  FORWARD CALCULATION
   30 IF (ICNT.NE.0) GO TO 35
C                                  SCALING IS UNNECESSARY
C                                  STORE INDIVIDUAL TERM
      IF (IOPT.EQ.1) T(JJ) = P1
      PE = PE+P1
C                                  ALL TERMS ARE ACCOUNTED FOR
   35 IF (JJ.EQ.K1) GO TO 60
C                                  CALCULATE NEXT TERM RECURSIVELY
      P1 = P1*X/Y
      IF (P1.LT.H) GO TO 40
C                                  RESCALE
      TEMP = P1*EPS
      IF(TEMP.EQ.0.0D0) GO TO 40
      P1 = TEMP
      ICNT = ICNT-1
   40 JJ = JJ+1
      Y = Y+1.D0
      GO TO 20
C                                  BACKWARD CALCULATION
   45 IF (KCNT.NE.0) GO TO 50
C                                  SCALING IS UNNECESSARY
C                                  STORE INDIVIDUAL TERM
      IF (IOPT.EQ.1) T(K1) = P2
      PE = PE+P2
C                                  ALL TERMS ARE ACCOUNTED FOR
   50 IF (JJ.EQ.K1) GO TO 60
C                                  CALCULATE NEXT TERM RECURSIVELY
      P2 = P2*X2/Y2
      IF (P2.LT.G) GO TO 55
C                                  RESCALE
      TEMP = P2*EPS
      IF(TEMP.EQ.0.0D0) GO TO 55
      P2 = TEMP
      KCNT = KCNT-1
   55 K1 = K1-1
      X2 = X2-1.D0
      GO TO 20
   60 P = PE
      IF(P.GT.1.0) P = 1.0E0
      RETURN
      END

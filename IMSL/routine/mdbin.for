C   IMSL ROUTINE NAME   - MDBIN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - BINOMIAL PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDBIN (K,N,P,PS,PK,IER)
C
C   ARGUMENTS    K      - NUMBER OF SUCCESSES (INPUT)
C                N      - NUMBER OF BERNOULLI TRIALS (INPUT)
C                P      - PROBABILITY OF SUCCESS ON EACH TRIAL (INPUT)
C                PS     - PROBABILITY THAT THE NUMBER OF SUCCESSES IS
C                           K OR LESS (OUTPUT)
C                PK     - PROBABILITY THAT THE NUMBER OF SUCCESSES IS
C                           EXACTLY K (OUTPUT)
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT K WAS SPECIFIED
C                             LESS THAN ZERO OR GREATER THAN N.
C                           IER = 130 INDICATES THAT P WAS SPECIFIED
C                             GREATER THAN 1 OR LESS THAN 0.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDBIN  (K,N,P,PS,PK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,N,IER
      REAL               P,PS,PK
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K1,ICNT,J
      REAL               P1,Q1,XN,XX,ALQN,QP,XJ,SML,ALSML
      DATA               SML/Z00100000/
C                                  FIRST EXECUTABLE STATEMENT
      ALSML = ALOG(SML)
      IER = 0
      PK = 0.
      PS = 0.
      IF (K .LE. N .AND. K .GE. 0) GO TO 5
C                                  TERMINAL ERROR - K SPECIFIED
C                                  INCORRECTLY
      IER = 129
      GO TO 9000
    5 IF (P .LE. 1. .AND. P .GE. 0.) GO TO 10
C                                  TERMINAL ERROR - P SPECIFIED
C                                  INCORRECTLY
      IER = 130
      GO TO 9000
C                                  SPECIAL CASE FOR P = 0.
   10 IF (P .NE. 0.) GO TO 15
      PS = 1.
      IF (K .NE. 0) GO TO 9005
      PK = 1.
      GO TO 9005
   15 IF (P .NE. 1.) GO TO 20
C                                  SPECIAL CASE FOR P = 1.
      IF (K .NE. N) GO TO 9005
      PK = 1.
      PS = 1.
      GO TO 9005
   20 P1 = P
      Q1 = 1.0-P
      K1 = K
      XN = N
      XX = XN*P
      IF (K .LE. XX) GO TO 25
      P1 = Q1
      Q1 = P
      K1 = N-K
   25 ALQN = XN*ALOG(Q1)
      ICNT = ALQN/ALSML
      ALQN = ALQN-ICNT*ALSML
      PK = EXP(ALQN)
      IF (K1 .EQ. 0) GO TO 35
      QP = P1/Q1
      XJ = 0.0
      XN = XN+1.0
      DO 30 J = 1,K1
         IF (ICNT .EQ. 0) PS = PS+PK
         XJ = XJ+1.0
         PK = PK*(QP*(XN-XJ))
         IF (PK .LT. XJ) GO TO 30
         PK = PK*SML
         ICNT = ICNT-1
   30    PK = PK/XJ
   35 IF (ICNT .NE. 0) PK = 0.0
      IF (K .GT. XX) GO TO 40
      PS = PS+PK
      GO TO 9005
   40 PS = 1.0-PS
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMDBIN )
 9005 RETURN
      END

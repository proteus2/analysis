C   IMSL ROUTINE NAME   - GGPOS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - POISSON RANDOM DEVIATE GENERATOR WHERE THE
C                           POISSON PARAMETER DOES NOT CHANGE OFTEN
C
C   USAGE               - CALL GGPOS (RLAM,DSEED,NR,K,IER)
C
C   ARGUMENTS    RLAM   - INPUT POISSON PARAMETER. RLAM MUST BE
C                           GREATER THAN ZERO.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT.  THE NUMBER OF RANDOM NUMBERS TO BE
C                           GENERATED.
C                K      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           RANDOM DEVIATES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT RLAM WAS SPECIFIED
C                             LESS THAN OR EQUAL TO ZERO.
C
C   REQD. IMSL ROUTINES - GGNML,GGUBS,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGPOS   (RLAM,DSEED,NR,K,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,K(NR),IER
      REAL               RLAM
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NN,J,NMIN,I
      REAL               RR(1),RLAM1,P,PP(85),PS,PN,X,RS,Y,EPS,R
      DATA               EPS/Z3C100000/
      DATA               R/0./,RLAM1/-1.0/,NN/0/,J/0/,P/0.0/,PP/85*0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (RLAM .GT. 0.0) GO TO 1
C                                  TERMINAL ERROR - RLAM WAS SPECIFIED
C                                  LESS THAN OR EQUAL TO ZERO
      IER = 129
      GO TO 9000
    1 IF(RLAM .LE. 50) GO TO 5
C                                  FOR RLAM GREATER THAN 50 USE
C                                  NORMAL APPROXIMATION
      RS = SQRT(RLAM)
      DO 3 I=1,NR
         CALL GGNML(DSEED,1,RR)
         R = RS*RR(1)+RLAM
         Y = R-AINT(R)
         K(I) = R
         IF(Y.GT..5) K(I) = K(I) + 1
         IF(K(I).LT.0) K(I) = 0
    3 CONTINUE
      GO TO 9005
    5 IF (RLAM .EQ. RLAM1) GO TO 10
C                                  THIS RLAM IS DIFFERENT FROM THE
C                                  PREVIOUS ONE
      NMIN = 1
      J = 1
      P = EXP(-RLAM)
      PP(1) = P
      RLAM1 = RLAM
   10 DO 55 I = 1,NR
         PS =PP(J)
         K(I) = 0
         CALL GGUBS(DSEED,1,RR)
         R = RR(1)
         IF (R .GT. PS) GO TO 30
C                                  FIND PROBABILITY BY SEARCHING THE
C                                  STORED VALUES
         NMIN = 1
         NMAX = J+1
   15    NN = (NMAX+NMIN-1)/2
         IF (NMAX-NMIN .LE. 1) GO TO 25
         IF (R .LE. PP(NN)) GO TO 20
C                                  UPPER HALF
         NMIN = NN+1
         GO TO 15
C                                  LOWER HALF
   20    NMAX = NN+1
         GO TO 15
   25    IF (NMIN .EQ. 85) GO TO 30
         K(I) = NMIN-1
         GO TO 55
C                                  CALCULATE ADDITIONAL TERMS
   30    PN = P
         NN = J+1
         X = J-1
   35    X = X+1.0
         IF (R .LE. PS) GO TO 50
         PN = PN*RLAM/X
         PS = PS+PN
         IF (PN .LE. EPS*PS) GO TO 45
         IF (NN .GT. 85) GO TO 40
C                                  SHOULD ADDITIONAL TERMS BE SAVED
         PP(NN) = PS
         P = PN
         J = NN
   40    NN = NN+1
         GO TO 35
   45    CALL GGUBS(DSEED,1,RR)
         R = RR(1)
         K(I) = (.5+R)*RLAM+NN
         GO TO 55
   50    K(I) = NN-2
   55 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GGPOS ')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - GGPON
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - POISSON RANDOM DEVIATE GENERATOR WHERE THE
C                           POISSON PARAMETER CHANGES FREQUENTLY
C
C   USAGE               - CALL GGPON (RLAM,DSEED,NR,K,IER)
C
C   ARGUMENTS    RLAM   - INPUT POISSON PARAMETER. RLAM MUST BE
C                           GREATER THAN 0.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - THE NUMBER OF RANDOM NUMBERS TO BE GENERATED.
C                           (INPUT)
C                K      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           DEVIATES.
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
      SUBROUTINE GGPON  (RLAM,DSEED,NR,K,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,K(NR),IER
      REAL               RLAM
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               Z,T,RS,R(1),Y,RR
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (RLAM .GT. 0.0) GO TO 5
C                                  TERMINAL ERROR - RLAM WAS SPECIFIED
C                                  LESS THAN OR EQUAL TO ZERO
      IER = 129
      GO TO 9000
    5 IF(RLAM .LE. 50) GO TO 15
C                                  FOR RLAM GREATER THAN 50 USE NORMAL
C                                  APPROXIMATION
      RS = SQRT(RLAM)
      DO 10 I=1,NR
         CALL GGNML(DSEED,1,R)
         RR = RS*R(1)+RLAM
         Y = RR-AINT(RR)
         K(I) = RR
         IF(Y.GT..5) K(I) = K(I)+1
         IF(K(I).LT.0) K(I) = 0
   10 CONTINUE
      GO TO 9005
   15 Z = EXP(-RLAM)
      DO 25 I = 1,NR
         T = 1.0
         K(I) = 0
C                                  OBTAIN UNIFORM RANDOM DEVIATE
   20    CALL GGUBS(DSEED,1,R)
         T = T*R(1)
         IF (T .LE. Z) GO TO 25
C                                  BUILD POISSON RANDOM DEVIATE
         K(I) = K(I)+1
      GO TO 20
   25 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GGPON ')
 9005 RETURN
      END

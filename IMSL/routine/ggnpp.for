C   IMSL ROUTINE NAME   - GGNPP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NONHOMOGENEOUS POISSON PROCESS GENERATOR
C                           WITH RATE FUNCTION LAMBDA(T) - FIXED
C                           INTERVAL, FIXED NUMBER, OR ONE AT A TIME.
C
C   USAGE               - CALL GGNPP  (DSEED,TL,TU,NUB,FUNLAM,RLAMAX,
C                           RLAMIN,IOPT,N,R,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0,2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                TL     - INPUT. LOWER ENDPOINT OF THE TIME INTERVAL
C                           (TL,TU). USUALLY TL=0.
C                TU     - INPUT. UPPER ENDPOINT OF THE TIME INTERVAL
C                           (TL,TU).
C                NUB    - INPUT. UPPER BOUND FOR THE NUMBER OF EVENTS
C                           TO BE GENERATED.  TO HAVE REASONABLE
C                           ASSURANCE THAT THE FULL PROCESS THROUGH
C                           TIME TU IS GENERATED, CALCULATE NUB
C                           AS FOLLOWS..
C                             X = RLAMAX * (TU-TL)
C                             NUB = X + 10. * (SQRT(X))
C                FUNLAM - USER-SUPPLIED FUNCTION SUBPROGRAM WHICH
C                           PROVIDES THE VALUE OF THE RATE FUNCTION
C                           LAMBDA(T) FOR EACH GIVEN VALUE OF THE
C                           TIMES-TO-EVENT IN A FIXED INTERVAL (TL,TU).
C                           THE USAGE IS AS FOLLOWS..
C                                    A = FUNLAM(T)
C                           FUNLAM MUST APPEAR IN AN EXTERNAL SPECIFICA-
C                           TION STATEMENT IN THE CALLING PROGRAM.
C                RLAMAX - INPUT. THE MAXIMUM VALUE OF THE RATE FUNCTION
C                           LAMBDA(T) IN A GIVEN INTERVAL (TL,TU).
C                RLAMIN - INPUT. THE MINIMUM VALUE OF THE RATE FUNCTION
C                           LAMBDA(T) IN A GIVEN INTERVAL (TL,TU).
C                           IF RLAMIN IS NOT KNOWN, SET RLAMIN = 0.
C                IOPT   - INPUT. OPTION SWITCH
C                           IF IOPT = 0, OUTPUT VECTOR R WILL CONTAIN
C                             THE TIMES TO EVENTS.
C                           IF IOPT = 1, OUTPUT VECTOR R WILL CONTAIN
C                             THE TIMES BETWEEN SUCCESSIVE EVENTS.
C                N      - OUTPUT. THE NUMBER OF EVENTS ACTUALLY
C                           GENERATED.  IF N IS LESS THAN NUB,
C                           THE TIME TU WAS EXCEEDED BEFORE NUB
C                           EVENTS WERE REALIZED.
C                R      - OUTPUT. VECTOR OF LENGTH N CONTAINING
C                           THE TIMES TO EVENTS IF IOPT = 0
C                           OR THE TIMES BETWEEN SUCCESSIVE EVENTS
C                           IF IOPT = 1.  R SHOULD BE DIMENSIONED TO
C                           BE OF LENGTH NUB.
C                IER    - ERROR PARAMETER.   (OUTPUT)
C                         WARNING ERROR (WITH FIX)
C                           IER=66 IMPLIES THAT RLAMAX IS POSITIVE BUT
C                             RLAMIN IS LESS THAN ZERO.  RLAMIN WILL BE
C                             SET TO ZERO.
C                         TERMINAL ERROR
C                           IER=130 IMPLIES THAT TL IS GREATER THAN OR
C                             EQUAL TO TU.
C                           IER=131 IMPLIES THAT RLAMAX IS NEGATIVE.
C                           IER=132 IMPLIES THAT BOTH RLAMAX AND RLAMIN
C                             ARE POSITIVE BUT RLAMIN IS GREATER THAN OR
C                             EQUAL TO RLAMAX.
C                           IER=133 IMPLIES THAT A COMPUTED VALUE OF
C                             FUNLAM IS NEGATIVE.
C                           IER=134 IMPLIES THAT FUNLAM COMPUTED
C                             AT A POINT IS GREATER THAN RLAMAX.
C
C   REQD. IMSL ROUTINES - GGEXN,GGUBS,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  TU MUST BE POSITIVE, TL MUST BE NONNEGATIVE.
C            2.  FUNLAM MUST BE LEFT-CONTINUOUS AND NONNEGATIVE IN THE
C                INTERVAL (TL,TU). GENERALLY FUNLAM WILL BE CONTINUOUS
C                BUT IT CAN ALSO BE PIECEWISE CONSTANT FUNCTION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGNPP  (DSEED,TL,TU,NUB,FUNLAM,RLAMAX,RLAMIN,
     1                   IOPT,N,R,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NUB,IOPT,N,IER
      REAL               TL,TU,RLAMAX,RLAMIN,R(1)
      DOUBLE PRECISION   DSEED
      EXTERNAL           FUNLAM
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      REAL               EXPD,Y,E(1),U(1),FLMT,VLMT,TEMP,RLAM
C
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CHECK VALUES OF TL AND TU
      IF (TL.LT.TU) GO TO 5
      IER = 130
      GO TO 9000
C                                  CHECK VALUES OF RLAMAX AND RLAMIN
    5 IF (RLAMAX.GT.0.) GO TO 10
      IER = 131
      GO TO 9000
   10 IF (RLAMIN.GE.0.) GO TO 15
      IER = 66
      CALL UERTST (IER,6HGGNPP )
      RLAMIN = 0.
      GO TO 20
   15 IF (RLAMIN.LT.RLAMAX) GO TO 20
      IER = 132
      GO TO 9000
C                                  ASSIGN INITIAL VALUES TO N
   20 FLMT = RLAMIN/RLAMAX
      N = 0
      Y = TL
   25 CALL GGEXN (DSEED,1.,1,E)
      EXPD = E(1)/RLAMAX
      Y = Y+EXPD
      IF (Y.GT.TU) GO TO 40
      CALL GGUBS (DSEED,1,U)
      IF (U(1).LE.FLMT) GO TO 30
      RLAM = FUNLAM(Y)
      IF (RLAM.LT.0.) GO TO 50
      IF (RLAM.GT.RLAMAX) GO TO 55
      VLMT = RLAM/RLAMAX
      IF (U(1).GT.VLMT) GO TO 35
   30 N = N+1
      R(N) = Y
      IF (N.GE.NUB) GO TO 40
   35 CONTINUE
      GO TO 25
C                                  SELECT APPROPRIATE OUTPUT VALUE FOR
C                                    VECTOR R ACCORDING TO IOPT- THE
C                                    OPTION SWITCH.
   40 IF (IOPT.NE.1) GO TO 9005
      TEMP = R(1)
      R(1) = R(1)-TL
      IF (N.EQ.0.OR.N.EQ.1) GO TO 9005
      DO 45 J=2,N
         R(J) = R(J)-TEMP
         TEMP = TEMP+R(J)
   45 CONTINUE
      GO TO 9005
C                                  FUNLAM IS NEGATIVE
   50 IER = 133
      GO TO 9000
C                                  FUNLAM IS GREATER THAN RLAMAX
   55 IER = 134
      GO TO 9000
C                                  PRINT ERROR MESSAGE
 9000 CONTINUE
      CALL UERTST (IER,'GGNPP ')
 9005 RETURN
      END

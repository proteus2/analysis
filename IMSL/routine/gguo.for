C   IMSL ROUTINE NAME   - GGUO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERATE SET OF ORDER STATISTICS FROM
C                           UNIFORM (0,1) DISTRIBUTION.
C
C   USAGE               - CALL GGUO (DSEED,IFIRST,ILAST,N,R,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT.  DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                IFIRST - INPUT.  FIRST ORDER STATISTIC TO
C                           BE GENERATED.
C                ILAST  - INPUT.  LAST ORDER STATISTIC TO BE
C                           GENERATED.  ILAST MUST BE GREATER THAN
C                           OR EQUAL TO IFIRST.  THE FULL SET OF
C                           ORDER STATISTICS FROM IFIRST TO ILAST
C                           IS GENERATED.  IF ONLY ONE ORDER
C                           STATISTIC IS DESIRED, SET ILAST=IFIRST.
C                N      - INPUT.  HYPOTHETICAL SAMPLE SIZE.
C                R      - OUTPUT.  VECTOR OF LENGTH (ILAST+1-IFIRST)
C                           CONTAINING THE ORDER STATISTICS IN
C                           ASCENDING ORDER.
C                IER    - OUTPUT.  ERROR PARAMETER.
C                         WARNING WITH FIX ERROR
C                           IER = 65  INDICATES IFIRST IS LESS THAN 1
C                             OR ILAST IS GREATER THAN N.  THE COM-
C                             PUTATIONS PROCEED AS IF THE PARAMETER(S)
C                             OUT OF RANGE HAD BEEN AT THE APPROPRIATE
C                             LIMIT(S) OF THE RANGE.
C                         TERMINAL ERROR
C                           IER = 129  INDICATES IFIRST IS GREATER
C                             THAN ILAST.
C
C   REQD. IMSL ROUTINES - GGBTR,GGEXN,GGUBFS,GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGUO   (DSEED,IFIRST,ILAST,N,R,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IFIRST,ILAST,N,IER
      REAL               R(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,J,LM1,IF1,IL1
      REAL               SUM,E(1),B(1),U,UR
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF1 = IFIRST
      IL1 = ILAST
      IF (IFIRST .LE. ILAST) GO TO 5
      IER = 129
      GO TO 9000
    5 IF (ILAST .LE. N) GO TO 10
      IF (IFIRST.GT.N) IF1 = N
      IER = 65
      IL1  = N
   10 IF (IFIRST .GE. 1) GO TO 15
      IF (ILAST.LT.1) IL1 = 1
      IER = 65
      IF1 = 1
   15 L = IL1 - IF1 + 1
      IF (L .LT. N) GO TO 30
C                                  PERFORM RANDOM PARTITIONING
C                                  METHOD FOR L EQUAL TO N
      SUM    = ZERO
      CALL GGEXN (DSEED,ONE,N,R)
      DO 20   J = 1,N
         SUM    = SUM + R(J)
         R(J)   = SUM
   20 CONTINUE
      CALL GGEXN (DSEED,ONE,1,E)
      SUM    = SUM + E(1)
      DO 25   J = 1,N
         R(J)   = R(J)/SUM
   25 CONTINUE
      GO TO 9000
C                                  FOR L LESS THAN N, GENERATE A
C                                  BETA VARIATE AND THEN USE
C                                  DESCENDING SEQUENTIAL METHOD
   30 IF (IF1 .EQ. 1) GO TO 45
      IF (IL1 .EQ. N) GO TO 55
      P = IL1
      Q = N - IL1 + 1
      CALL GGBTR (DSEED,P,Q,1,B)
      R(L) = B(1)
   35 LM1 = L - 1
      IF (LM1.EQ.0) GO TO 9000
      U = IL1 - 1
      DO 40 J=1,LM1
         UR = ONE/U
         R(L-J) = R(L+1-J)*(GGUBFS(DSEED))**UR
         U = U - ONE
   40 CONTINUE
      GO TO 9000
   45 TEMP = ZERO
      DO 50 J=1,L
         UR = ONE/FLOAT(N-J+1)
         R(J) = ONE-(ONE-TEMP)*(GGUBFS(DSEED))**UR
         TEMP = R(J)
   50 CONTINUE
      GO TO 9000
   55 UR = ONE/FLOAT(IL1)
      R(L) = GGUBFS(DSEED)**UR
      GO TO 35
 9000 CONTINUE
      IF (IER.EQ.0) GO TO 9005
      CALL UERTST (IER,'GGUO  ')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - GGNO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERATE SET OF ORDER STATISTICS FROM
C                           NORMAL DISTRIBUTION
C
C   USAGE               - CALL GGNO (DSEED,IFIRST,ILAST,N,R,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
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
C                R      - OUTPUT VECTOR OF LENGTH (ILAST+1-IFIRST)
C                           CONTAINING THE ORDER STATISTICS IN
C                           ASCENDING ORDER.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                         WARNING ERROR (WITH FIX)
C                           IER = 65  INDICATES IFIRST IS LESS THAN 1
C                             OR ILAST IS GREATER THAN N.  THE COM-
C                             PUTATIONS PROCEED AS IF THE PARAMETER(S)
C                             OUT OF RANGE HAD BEEN AT THE APPROPRIATE
C                             LIMIT(S) OF THE RANGE.
C                         TERMINAL ERROR
C                           IER = 129  INDICATES IFIRST IS GREATER
C                             THAN ILAST.
C
C   REQD. IMSL ROUTINES - GGUO,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGNO   (DSEED,IFIRST,ILAST,N,R,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IFIRST,ILAST,N,IER
      REAL               R(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,J,IF1,IL1
      REAL               P,Y
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF1 = IFIRST
      IL1 = ILAST
      IF (IFIRST .GT. ILAST) GO TO 25
      IF (ILAST .LE. N) GO TO 5
      IF (IFIRST.GT.N) IF1 = N
      IER = 65
      IL1 = N
    5 IF (IFIRST .GE. 1) GO TO 10
      IF (ILAST.LT.1) IL1 = 1
      IER = 65
      IF1 = 1
   10 IF (IER .EQ. 0) GO TO 15
      CALL UERTST (IER,6HGGNO  )
   15 L = IL1 - IF1 + 1
      CALL GGUO (DSEED,IF1,IL1,N,R,IER)
      DO 20 J=1,L
         P = R(J)
         CALL MDNRIS(P,Y,IER)
         R(J) = Y
   20 CONTINUE
      GO TO 9005
   25 IER = 129
      CALL UERTST (IER,'GGNO  ')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - GTRTN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - TALLY OF NUMBER OF RUNS UP AND DOWN
C
C   USAGE               - CALL GTRTN (R,N,IOPT,WK,RUNS)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SEQUENCE OF NUMBERS HYPOTHESIZED TO BE
C                           RANDOM.
C                N      - INPUT LENGTH OF VECTOR R
C                IOPT   - INPUT/OUTPUT OPTION PARAMETER. IF THE SEQUENCE
C                           OF RANDOM NUMBERS FITS INTO CORE THEN IOPT
C                           SHOULD BE SET TO 0 ON INPUT AND GTRTN NEED
C                           ONLY BE CALLED ONCE. IF THE SEQUENCE OF
C                           NUMBERS DOES NOT FIT INTO CORE, THE SEQUENCE
C                           SHOULD BE BROKEN INTO PARTS AND GTRTN MAY
C                           BE CALLED AS MANY TIMES AS NEEDED. IN THIS
C                           CASE, IOPT NEED BE SET TO 0 ON THE FIRST
C                           CALL ONLY. IOPT MUST NOT BE CHANGED BY THE
C                           USER BETWEEN SUBSEQUENT CALLS.
C                WK     - WORK VECTOR OF LENGTH 3. WK MUST NOT BE
C                           CHANGED BY THE USER BETWEEN CALLS TO GTRTN.
C                RUNS   - OUTPUT VECTOR OF LENGTH 8 CONTAINING TALLIES.
C                           RUNS(I) CONTAINS THE NUMBER OF RUNS OF
C                           LENGTH I OCCURRING IN THE SEQUENCE. RUNS
C                           MUST NOT BE CHANGED BY THE USER BETWEEN
C                           CALLS TO GTRTN. RUNS OF LENGTH 8 OR MORE
C                           ARE TALLIED IN RUNS(8).
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTRTN (R,N,IOPT,WK,RUNS)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IOPT
      REAL               R(1),RUNS(8),WK(3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IS,J,K,I,LR
      REAL               X,Y
C                                  FIRST EXECUTABLE STATEMENT
      IF(IOPT.NE.0) GO TO 7
      LR = 0
      K = 0
      X = 0.
      IS=1
      J=2
      X=R(1)
      Y=R(2)
         DO 5 I=1,8
    5    RUNS(I)=0.
      K=1
      GO TO 15
    7 LR = WK(1)
      K = WK(2)
      X = WK(3)
      IS = 0
      J = 1
      Y = R(1)
      IF (IOPT.EQ.20) GO TO 20
      GO TO 35
   10 K=0
   15 LR=0
C                                  INCREASING SEQUENCE
      IOPT=20
   20 IF (X .GT. Y) GO TO 25
      IS=IS+1
      LR=LR+1
      J=IS+1
      X=R(IS)
      IF (J .GT. N) GO TO 45
      Y=R(J)
      GO TO 20
   25 IF (LR .NE. 0) GO TO 40
C                                  DECREASING SEQUENCE
   30 IOPT=35
   35 IF (X .LT. Y) GO TO 40
      IS=IS+1
      LR=LR+1
      J=IS+1
      X=R(IS)
      IF (J .GT. N) GO TO 45
      Y=R(J)
      GO TO 35
C                                  TALLY THE RUNS
   40 IF (K .EQ. 1) GO TO 10
      IF (LR .GT. 8) LR=8
      RUNS(LR)=RUNS(LR)+1
      GO TO 15
   45 WK(1) = LR
      WK(2) = K
      WK(3) = X
      RETURN
      END

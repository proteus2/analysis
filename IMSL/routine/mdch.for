C   IMSL ROUTINE NAME   - MDCH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CHI-SQUARED PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDCH (CS,DF,P,IER)
C
C   ARGUMENTS    CS     - INPUT VALUE FOR WHICH THE PROBABILITY IS
C                           COMPUTED. CS MUST BE GREATER THAN OR EQUAL
C                           TO ZERO.
C                DF     - INPUT NUMBER OF DEGREES OF FREEDOM OF THE
C                           CHI-SQUARED DISTRIBUTION. DF MUST BE GREATER
C                           THAN OR EQUAL TO .5 AND LESS THAN OR EQUAL
C                           TO 200,000.
C                P      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE
C                           WHICH FOLLOWS THE CHI-SQUARED DISTRIBUTION
C                           WITH DF DEGREES OF FREEDOM IS LESS THAN OR
C                           EQUAL TO CS.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT CS OR DF WAS
C                             SPECIFIED INCORRECTLY.
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT THE NORMAL PDF
C                             WOULD HAVE PRODUCED AN UNDERFLOW.
C
C   REQD. IMSL ROUTINES - H32/MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,UERTST,
C                           UGETIO
C                       - H36,H48,H60/MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDCH (CS,DF,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               CS,DF,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               PT2
      DOUBLE PRECISION   A,Z,DGAM,EPS,W,W1,B,Z1,HALF,ONE,THRTEN,THRD
      DOUBLE PRECISION   DGAMMA
      REAL               RINFM,X,C
      DATA               EPS/1.0D-6/,HALF/5.D-1/,THRTEN/13.D0/,ONE/1.D0/
      DATA               THRD/.3333333333333333D0/
      DATA               PT2/.2222222E0/
      DATA               RINFM/ZFFFFFFFF/
      FUNC(W,A,Z)=W*DEXP(A*DLOG(Z)-Z)
C                                  FIRST EXECUTABLE STATEMENT
C                                  TEST FOR INVALID INPUT VALUES
      IF (DF .GE. .5 .AND. DF .LE. 2.E5 .AND. CS .GE. 0.0) GO TO 5
      IER=129
      P=RINFM
      GO TO 9000
    5 IER=0
C                                  SET P=0. IF CS IS LESS THAN OR
C                                  EQUAL TO 10.**(-12)
      IF (CS .GT. 1.E-12) GO TO 15
   10 P=0.0
      GO TO 9005
   15 IF(DF.LE.100.) GO TO 20
C                                  USE NORMAL DISTRIBUTION APPROXIMATION
C                                  FOR LARGE DEGREES OF FREEDOM
      IF(CS.LT.2.0) GO TO 10
      X=((CS/DF)**THRD-(ONE-PT2/DF))/SQRT(PT2/DF)
      IF (X .GT. 5.0) GO TO 50
      IF (X .LT. -18.8055) GO TO 55
      CALL MDNOR (X,P)
      GO TO 9005
C                                  INITIALIZATION FOR CALCULATION USING
C                                  INCOMPLETE GAMMA FUNCTION
   20 IF (CS .GT. 200.) GO TO 50
      A=HALF*DF
      Z=HALF*CS
      DGAM = DGAMMA(A)
      W=DMAX1(HALF*A,THRTEN)
      IF (Z .GE. W) GO TO 35
      IF (DF .GT. 25. .AND. CS .LT. 2.) GO TO 10
C                                  CALCULATE USING EQUATION NO. 6.5.29
      W=ONE/(DGAM*A)
      W1=W
         DO 25 I=1,50
         B=I
         W1=W1*Z/(A+B)
         IF (W1 .LE. EPS*W) GO TO 30
         W=W+W1
   25    CONTINUE
   30 P=FUNC(W,A,Z)
      GO TO 9005
C                                  CALCULATE USING EQUATION NO. 6.5.32
   35 Z1=ONE/Z
      B=A-ONE
      W1=B*Z1
      W=ONE+W1
         DO 40 I=2,50
         B=B-ONE
         W1=W1*B*Z1
         IF (W1 .LE. EPS*W) GO TO 45
         W=W+W1
   40    CONTINUE
   45 W=Z1*FUNC(W,A,Z)
      P=ONE-W/DGAM
      GO TO 9005
   50 P=1.0
      GO TO 9005
C                                  WARNING ERROR - UNDERFLOW WOULD HAVE
C                                  OCCURRED
   55 P=0.0
      IER=34
 9000 CONTINUE
      CALL UERTST (IER,6HMDCH  )
 9005 RETURN
      END

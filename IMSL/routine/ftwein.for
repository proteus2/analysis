C   IMSL ROUTINE NAME   - FTWEIN
C
C----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - WIENER FORECAST FOR A STATIONARY STOCHASTIC
C                           PROCESS
C
C   USAGE               - CALL FTWEIN (WHITE,M,N,EPS,X,LP,WA,IER)
C
C   ARGUMENTS    WHITE  - INPUT VARIABLE IN THE RANGE (0,1) WHICH
C                           INDICATES THE EFFECT OF WHITE NOISE ON
C                           E(X**2).  THE PROGRAM ESTIMATES E(X**2) BY
C                           THE INNER PRODUCT OF THE TIME SERIES WITH
C                           ITSELF DIVIDED BY N.  A FINAL ADJUSTMENT
C                           IS MADE BY COMPUTING
C                           E(X**2) = E(X**2)*(1.0+WHITE).
C                M      - INPUT MAXIMUM ALLOWABLE LENGTH OF FORECAST
C                           OPERATOR. M MUST BE GREATER THAN OR EQUAL
C                           TO 1.
C                N      - INPUT LENGTH OF TIME SERIES.
C                EPS    - INPUT VARIABLE WITH A VALUE IN THE RANGE
C                           (0,1).  THE ROUTINE FIRST COMPUTES A
C                           FORECAST OF LENGTH ONE.  IF THE NORMALIZED
C                           MEAN SQUARE ERROR FOR THIS FORECAST OPERATOR
C                           IS LESS THAN EPS, THE ROUTINE TERMINATES.
C                           IF NOT, THE ROUTINE COMPUTES A FORECAST
C                           OPERATOR OF LENGTH 2 AND AGAIN TESTS THE
C                           NORMALIZED MEAN SQUARE ERROR AGAINST EPS,
C                           AND SO FORTH.
C                           THE ROUTINE TERMINATES EITHER WITH
C                           A FORECAST OPERATOR OF LENGTH LP
C                           WHICH SATISFIES THIS CRITERION OR AN
C                           OPERATOR OF MAXIMUM ALLOWABLE LENGTH M
C                           AND AN ERROR SIGNAL.
C                X      - INPUT/OUTPUT VECTOR OF LENGTH N+M.
C                         ON INPUT X(1),...,X(N) CONTAINS THE TIME
C                           SERIES.
C                         ON OUTPUT THE FORECAST OPERATOR WILL OCCUPY
C                           POSITIONS X(N+M),X(N+M-1),...,X(N+M-LP+1)
C                           WHERE LP IS LESS THAN OR EQUAL TO M.
C                           THE USER CAN THEN FORECAST BY THE EQUATION
C                           V(T) = X(N+M)*V(T-1)+X(N+M-1)*V(T-2)+...
C                           X(N+M-LP+1)*V(T-LP), WHERE V(T) IS THE
C                           VALUE OF THE TIME SERIES AT TIME T.
C                LP     - OUTPUT VARIABLE GIVING THE LENGTH OF THE
C                           FORECAST OPERATOR.
C                WA     - WORK VECTOR OF LENGTH 2*M+1
C                IER    - ERROR INDICATOR. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES NO OPERATOR COULD BE
C                             FOUND OF LENGTH LESS THAN OR EQUAL TO M
C                             WHICH PRODUCED A NORMALIZED MEAN SQUARE
C                             ERROR LESS THAN EPS.
C                         TERMINAL ERROR
C                           IER=130 INDICATES THAT M WAS SPECIFIED LESS
C                             THAN 1.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTWEIN (WHITE,M,N,EPS,X,LP,WA,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LP,IER
      REAL               WHITE,EPS,X(1),WA(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,L,LPP1,LPP2,MP1,NPM,NPMMLP,NPMP1,NPMP2
      REAL               C,ZERO
      DOUBLE PRECISION   PIVOT,RHS
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 130
      IF (M .LT. 1) GO TO 9000
C                                  FORM SOME USEFUL INTEGERS
      IER = 0
      MP1 = M+1
      LP = 1
      LPP1 = LP+1
      LPP2 = LPP1+1
      NPM = N+M
      NPMP1 = NPM+1
      NPMP2 = NPMP1+1
C                                  ZERO OUT X(N+1),...,X(N+M)
      DO 1 I=1,M
         X(N+I) = ZERO
    1 CONTINUE
C                                  FORM COVARIANCE
      DO 5 I = 1,MP1
         WA(I) = ZERO
    5 CONTINUE
      DO 15 I = 1,N
         J = I
         DO 10 L = 1,MP1
            WA(L) = WA(L)+X(I)*X(J)
            J = J+1
   10    CONTINUE
   15 CONTINUE
      C = N
      DO 20 I=1,MP1
   20 WA(I) = WA(I)/C
C                                  ADJUST FOR WHITE NOISE
      WA(1) = WA(1)*(1.0+WHITE)
C                                  COMPUTE SOLUTION OF LENGTH 1
      X(NPM) = WA(2)/WA(1)
C                                  CHECK FOR CONVERGENCE
   25 PIVOT = WA(1)
      DO 30 I=2,LPP1
   30 PIVOT = PIVOT-DBLE(X(NPMP2-I))*DBLE(WA(I))
      IF(ABS(SNGL(PIVOT/WA(1))).LE. EPS) GO TO 9005
      IF(LP .GE. M) GO TO 50
C                                  COMPUTE NEW FORECAST
      RHS = WA(LPP2)
      DO 35 I=1,LP
   35 RHS = RHS-DBLE(X(NPMP1-I))*DBLE(WA(LPP2-I))
      PIVOT = RHS/PIVOT
      NPMMLP = NPM-LP
      DO 40 I=1,LP
   40 WA(MP1+I) = X(NPMP1-I) - PIVOT*X(NPMMLP+I)
      DO 45 I=1,LP
   45 X(NPMP1-I) = WA(MP1+I)
      X(NPM-LP) = PIVOT
      LP = LPP1
      LPP1 = LPP2
      LPP2 = LPP2+1
      GO TO 25
   50 IER = 33
 9000 CONTINUE
      CALL UERTST(IER,'FTWEIN')
 9005 RETURN
      END

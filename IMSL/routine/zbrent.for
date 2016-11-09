C   IMSL ROUTINE NAME   - ZBRENT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
C                           GIVEN INTERVAL (BRENT ALGORITHM)
C
C   USAGE               - CALL ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)
C
C   ARGUMENTS    F      - AN EXTERNAL FUNCTION SUBPROGRAM F(X)
C                           PROVIDED BY THE USER WHICH COMPUTES F FOR
C                           ANY X IN THE INTERVAL (A,B). (INPUT)
C                           F MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                           THE CALLING PROGRAM
C                EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT,
C                           B, IS ACCEPTED IF ABS(F(B)) IS LESS THAN OR
C                           EQUAL TO EPS.  EPS MAY BE SET TO ZERO.
C                NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
C                           B, IS ACCEPTED IF THE CURRENT APPROXIMATION
C                           AGREES WITH THE TRUE SOLUTION TO NSIG
C                           SIGNIFICANT DIGITS.
C                A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
C                           AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE
C                           IN SIGN.
C                           ON OUTPUT, BOTH A AND B ARE ALTERED.  B
C                           WILL CONTAIN THE BEST APPROXIMATION TO THE
C                           ROOT OF F. SEE REMARK 1.
C                MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND
C                           ON THE NUMBER OF FUNCTION EVALUATIONS
C                           REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN
C                           WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
C                           EVALUATIONS USED.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THE ALGORITHM FAILED TO
C                             CONVERGE IN MAXFN EVALUATIONS.
C                           IER = 130 INDICATES F(A) AND F(B) HAVE THE
C                             SAME SIGN.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE
C                FOLLOWING,
C                F(A)*F(B) .LE.0,
C                ABS(F(B)) .LE. ABS(F(A)), AND
C                EITHER ABS(F(B)) .LE. EPS OR
C                ABS(A-B) .LE. MAX(ABS(B),0.1)*10.0**(-NSIG).
C                THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES
C                LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE
C                COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED
C                IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL
C                MAGNITUDE.
C            2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN
C                K = (ALOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
C                  D=MIN(OVER X IN (A,B) OF
C                    MAX(ABS(X),0.1)*10.0**(-NSIG)).
C                THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS.
C                RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY
C                ZBRENT EXCEED SQRT(K). D CAN BE COMPUTED AS FOLLOWS,
C                  P = AMIN1(ABS(A),ABS(B))
C                  P = AMAX1(0.1,P)
C                  IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1
C                  D = P*10.0**(-NSIG)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,MAXFN,IER
      REAL               F,EPS,A,B
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IC
      REAL               ZERO,HALF,ONE,THREE,TEN,
     1                   T,FA,FB,C,FC,D,E,TOL,RM,S,P,Q,R,RONE,TEMP
      DATA               ZERO/0.0/,HALF/.5/,ONE/1.0/,THREE/3.0/,
     1                   TEN/10.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      T = TEN**(-NSIG)
      IC = 2
      S = A
      FA = F(S)
      S = B
      FB = F(S)
C                                  TEST FOR SAME SIGN
      IF (FA*FB.GT.ZERO) GO TO 50
    5 C = A
      FC = FA
      D = B-C
      E = D
   10 IF (ABS(FC).GE.ABS(FB)) GO TO 15
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
   15 CONTINUE
      TOL = T*AMAX1(ABS(B),0.1)
      RM = (C-B)*HALF
C                                  TEST FOR FIRST CONVERGENCE CRITERIA
      IF (ABS(FB).LE.EPS) GO TO 40
C                                  TEST FOR SECOND CONVERGENCE CRITERIA
      IF (ABS(C-B).LE.TOL) GO TO 40
C                                  CHECK EVALUATION COUNTER
      IF (IC.GE.MAXFN) GO TO 45
C                                  IS BISECTION FORCED
      IF (ABS(E).LT.TOL) GO TO 30
      IF (ABS(FA).LE.ABS(FB)) GO TO 30
      S = FB/FA
      IF (A.NE.C) GO TO 20
C                                  LINEAR INTERPOLATION
      P = (C-B)*S
      Q = ONE-S
      GO TO 25
C                                  INVERSE QUADRATIC INTERPOLATION
   20 Q = FA/FC
      R = FB/FC
      RONE = R-ONE
      P = S*((C-B)*Q*(Q-R)-(B-A)*RONE)
      Q = (Q-ONE)*RONE*(S-ONE)
   25 IF (P.GT.ZERO) Q = -Q
      IF (P.LT.ZERO) P = -P
      S = E
      E = D
C                                  IF ABS(P/Q).GE.75*ABS(C-B) THEN
C                                     FORCE BISECTION
      IF (P+P.GE.THREE*RM*Q) GO TO 30
C                                  IF ABS(P/Q).GE..5*ABS(S) THEN FORCE
C                                     BISECTION. S = THE VALUE OF P/Q
C                                     ON THE STEP BEFORE THE LAST ONE
      IF (P+P.GE.ABS(S*Q)) GO TO 30
      D = P/Q
      GO TO 35
C                                  BISECTION
   30 E = RM
      D = E
C                                  INCREMENT B
   35 A = B
      FA = FB
      TEMP = D
      IF (ABS(TEMP).LE.HALF*TOL) TEMP = SIGN(HALF*TOL,RM)
      B = B+TEMP
      S = B
      FB = F(S)
      IC = IC+1
      IF (FB*FC.LE.ZERO) GO TO 10
      GO TO 5
C                                  CONVERGENCE OF B
   40 A = C
      MAXFN = IC
      GO TO 9005
C                                  MAXFN EVALUATIONS
   45 IER = 129
      A = C
      MAXFN = IC
      GO TO 9000
C                                  TERMINAL ERROR - F(A) AND F(B) HAVE
C                                  THE SAME SIGN
   50 IER = 130
      MAXFN = IC
 9000 CONTINUE
      CALL UERTST (IER,6HZBRENT)
 9005 RETURN
      END

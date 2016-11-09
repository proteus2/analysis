C   IMSL ROUTINE NAME   - ZXGSN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ONE-DIMENSIONAL UNIMODAL FUNCTION
C                           MINIMIZATION USING THE GOLDEN SECTION
C                           SEARCH METHOD
C
C   USAGE               - CALL ZXGSN (F,A,B,TOL,XMIN,IER)
C
C   ARGUMENTS    F      - A REAL FUNCTION SUBPROGRAM SUPPLIED BY
C                           THE USER. (INPUT)
C                           F MUST BE DECLARED EXTERNAL IN THE CALLING
C                           PROGRAM. F DEFINES THE FUNCTION TO BE
C                           MINIMIZED AND SHOULD BE OF THE FOLLOWING
C                           FORM
C                             F(X)
C                           WHERE X IS THE INDEPENDENT VARIABLE.
C                           F MUST NOT ALTER X.
C                           F IS ASSUMED TO DEFINE A UNIMODAL FUNCTION.
C                           THAT IS, A FUNCTION WITH A UNIQUE MINIMUM
C                           VALUE, XMIN, IN THE INTERVAL DEFINED BY A
C                           AND B. A FUNCTION, F, IS UNIMODAL IF IT
C                           SATISFIES THE FOLLOWING CONDITIONS
C                             FOR ALL X0, X1, AND X2 IN THE INTERVAL
C                             (A,B) INCLUSIVELY, IF X0 IS LESS THAN
C                             X1 AND X1 IS LESS THAN X2, THEN
C                             1. IF F(X0) IS LESS THAN OR EQUAL TO
C                                F(X1), THEN F(X1) IS LESS THAN F(X2).
C                             2. IF F(X1) IS GREATER THAN OR EQUAL
C                                TO F(X2), THEN F(X0) IS GREATER THAN
C                                F(X1).
C                A      - (INPUT/OUTPUT)
C                           ON INPUT, A IS THE LOWER ENDPOINT OF THE
C                           INTERVAL IN WHICH THE MINIMUM OF F IS TO BE
C                           LOCATED.
C                           ON OUTPUT, A IS THE LOWER ENDPOINT OF THE
C                           INTERVAL IN WHICH THE MINIMUM OF F IS
C                           LOCATED.
C                B      - (INPUT/OUTPUT)
C                           ON INPUT, B IS THE UPPER ENDPOINT OF THE
C                           INTERVAL IN WHICH THE MINIMUM OF F IS TO BE
C                           LOCATED.
C                           ON OUTPUT, B IS THE UPPER ENDPOINT OF THE
C                           INTERVAL IN WHICH THE MINIMUM OF F IS
C                           LOCATED.
C                TOL    - THE LENGTH OF THE FINAL SUBINTERVAL
C                           CONTAINING THE MINIMUM. (INPUT)
C                XMIN   - THE APPROXIMATE MINIMUM OF THE FUNCTION F
C                           ON THE ORIGINAL INTERVAL (A,B). (OUTPUT)
C                           ON OUTPUT, WHEN IER=0, THE FOLLOWING
C                           CONDITIONS HOLD
C                           1. (B-A) IS LESS THAN OR EQUAL TO TOL.
C                           2. A IS LESS THAN OR EQUAL TO XMIN AND
C                              XMIN IS LESS THAN OR EQUAL TO B.
C                           3. F(XMIN) IS LESS THAN OR EQUAL TO F(A)
C                              AND F(XMIN) IS LESS THAN OR EQUAL TO
C                              F(B).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 IMPLIES THAT A IS GREATER THAN
C                             OR EQUAL TO B. XMIN IS SET TO A.
C                           IER=130 IMPLIES THAT TOL IS GREATER THAN
C                             OR EQUAL TO THE LENGTH OF THE INPUT
C                             INTERVAL (A,B). XMIN IS SET TO A.
C                           IER=131 IMPLIES THAT THE FUNCTION F IS NOT
C                             UNIMODAL OR APPEARS TO BE NOT UNIMODAL TO
C                             THE ROUTINE DUE TO ROUNDING ERRORS IN THE
C                             EXTERNAL EVALUATION FUNCTION.
C                             WHEN THIS ERROR OCCURS THE FOLLOWING
C                             CONDITIONS HOLD
C                             1. A IS LESS THAN OR EQUAL TO XMIN AND
C                                XMIN IS LESS THAN OR EQUAL TO B.
C                             2. F(XMIN) IS GREATER THAN OR EQUAL TO
C                                F(A) AND F(XMIN) IS GREATER THAN OR
C                                EQUAL TO F(B) (ONLY ONE EQUALITY CAN
C                                HOLD).
C                             FURTHER ANALYSIS OF THE FUNCTION F IS
C                             NECESSARY IN ORDER TO DETERMINE WHETHER
C                             IT IS NOT UNIMODAL IN THE MATHEMATICAL
C                             SENSE OR WHETHER IT APPEARS TO BE NOT
C                             UNIMODAL TO THE ROUTINE DUE TO ROUNDING
C                             ERRORS IN WHICH CASE THE A,B,AND XMIN
C                             RETURNED MAY BE ACCEPTABLE.
C                           IER=132 IMPLIES THAT THE INTERVAL HAS BEEN
C                             REDUCED AS FAR AS NUMERICALLY POSSIBLE.
C                             THIS IS DUE TO TOL BEING TOO SMALL.
C                             WHEN THIS ERROR OCCURS, XMIN GIVES
C                             THE LOCATION OF THE MINIMUM AS
C                             ACCURATELY AS POSSIBLE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IMSL SUBROUTINE ZXGSP IS IDENTICAL IN
C                PURPOSE TO ZXGSN BUT INCLUDES INPUT
C                PARAMETERS P1,P2,IP3,IP4, AND IP5
C                WHICH THE USER MAY REQUIRE TO RELAY DATA
C                FROM THE MAIN PROGRAM TO THE FUNCTION
C                SUBPROGRAM F.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXGSN  (F,A,B,TOL,XMIN,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               F,A,B,TOL,XMIN
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               C,FA,FB,H,V1,V2,FV1,FV2
C                                  FIRST EXECUTABLE STATEMENT
      XMIN = A
      IER = 129
C                                  B MUST BE GREATER THAN A
      IF (B .LE. A) GO TO 9000
      IER = 130
C                                  TOL MUST BE SMALLER THAN THE
C                                    INITIAL INTERVAL
      IF (TOL .GE. (B-A)) GO TO 9000
      IER = 0
C                                  COMPUTE THE FIBONACCI CONSTANT
      C = (3.0-SQRT(5.0))/2.0
C                                  COMPUTE THE INITIAL STEP
      H = C*(B-A)
C                                  COMPUTE THE NEW POINTS
      V1 = A+H
      V2 = B-H
C                                  MAKE THE INITIAL FUNCTION EVALUATIONS
      FA = F(A)
      FB = F(B)
      FV1 = F(V1)
      FV2 = F(V2)
C                                  EACH ITERATION BEGINS HERE
    5 CONTINUE
C                                  HAS THE INTERVAL BECOME TOO SMALL
      IF (A .GE. V1 .OR. V1 .GE. V2 .OR. V2 .GE. B) GO TO 40
C                                  FIND THE CURRENT MINIMUM
      IF (FV1 .GE. FV2) GO TO 10
C                                  V1 IS THE MINIMUM
C                                  CHECK TO SEE IF THE FUNCTION IS
C                                    NOT UNIMODAL
      IF (FV2 .GT. FB) GO TO 25
C                                  UPDATE THE INTERVAL
C                                    V2 BECOMES THE NEW B
      B = V2
C                                  IS THE INTERVAL SUFFICIENTLY SMALL
      IF (TOL .GE. (B-A)) GO TO 15
C                                  REDUCE THE INTERVAL FURTHER
      FB = FV2
      V2 = V1
      FV2 = FV1
      H = C*(B-A)
      V1 = A+H
      FV1 = F(V1)
      GO TO 5
C                                  V2 IS THE MINIMUM
C                                  CHECK TO SEE IF THE FUNCTION IS
C                                    NOT UNIMODAL
   10 IF (FV1 .GT. FA) GO TO 30
C                                  UPDATE THE INTERVAL
C                                    V1 BECOMES THE NEW A
      A = V1
C                                  IS THE INTERVAL SUFFICIENTLY SMALL
      IF (TOL .GE. (B-A)) GO TO 20
C                                  REDUCE THE INTERVAL FURTHER
      FA = FV1
      V1 = V2
      FV1 = FV2
      H = C*(B-A)
      V2 = B-H
      FV2 = F(V2)
      GO TO 5
C                                  CONVERGENCE OBTAINED. V1 OR A
C                                    IS THE MINIMUM
   15 XMIN = V1
      IF (FA .LT. FV1) XMIN = A
      GO TO 9005
C                                  CONVERGENCE OBTAINED. V2 OR B
C                                    IS THE MINIMUM
   20 XMIN = V2
      IF (FB .LT. FV2) XMIN = B
      GO TO 9005
C                                  FUNCTION IS NOT UNIMODAL. RETURN
C                                    THE NECESSARY PARAMETERS
   25 XMIN = V2
      A = V1
      GO TO 35
   30 XMIN = V1
      B = V2
   35 IER = 131
      GO TO 9000
C                                  THE INTERVAL HAS BECOME TOO SMALL
   40 IER = 132
      XMIN = A
      IF (FB .LT. FA) XMIN = B
C
 9000 CONTINUE
      CALL UERTST (IER,6HZXGSN )
 9005 RETURN
      END

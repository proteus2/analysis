C   IMSL ROUTINE NAME   - MERFCI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - INVERSE COMPLEMENTED ERROR FUNCTION
C
C   USAGE               - CALL MERFCI (P,Y,IER)
C
C   ARGUMENTS    P      - INPUT VALUE IN THE EXCLUSIVE RANGE (0.0,2.0)
C                Y      - OUTPUT VALUE OF THE INVERSE COMPLEMENTED
C                           ERROR FUNCTION
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES P LIES OUTSIDE THE LEGAL
C                             RANGE. PLUS OR MINUS MACHINE INFINITY IS
C                             GIVEN AS THE RESULT (SIGN IS THE SIGN OF
C                             THE FUNCTION VALUE OF THE NEAREST LEGAL
C                             ARGUMENT).
C
C   REQD. IMSL ROUTINES - MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MERFCI (P,Y,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               EPS,G0,G1,G2,G3,H0,H1,H2,W,WI,SN,SD
      REAL               SIGMA,X,XINF
      DATA               XINF/Z7FFFFFFF/
      DATA               EPS/Z3C100000/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  ERROR CHECK ON INPUT ARGUMENT
      IF(P .GT. 0.0 .AND. P .LT. 2.0) GO TO 5
      IER = 129
      SIGMA = SIGN(1.0,-P)
      Y = SIGMA * XINF
      GO TO 9000
    5 IF(P.LE.EPS) GO TO 10
      X = 1.0 - P
      CALL MERFI (X,Y,IER)
      GO TO 9005
C                                  P TOO SMALL, COMPUTE Y DIRECTLY
   10 W = SQRT(-ALOG(P+(P-P*P)))
C                                  W GREATER THAN 4., APPROX. F BY A
C                                     RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      Y = W + W*(G0+SN/SD)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMERFCI)
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - RLINCF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - RESPONSE CONTROL USING A FITTED SIMPLE LINEAR
C                           REGRESSION MODEL
C
C   USAGE               - CALL RLINCF (CRIT,IOP,STAT,IER)
C
C   ARGUMENTS    CRIT   - INPUT VECTOR OF LENGTH 9.
C                         CRIT(1) CONTAINS THE INTERCEPT ESTIMATE
C                           (NOT REQUIRED IF IOP(3) IS NOT EQUAL TO 0).
C                         CRIT(2) CONTAINS THE SLOPE ESTIMATE.
C                           CRIT(2) MUST BE NONZERO.
C                         CRIT(3) CONTAINS THE INDEPENDENT
C                           VARIABLE MEAN (NOT REQUIRED IF IOP(3) IS NOT
C                           EQUAL TO 0).
C                         CRIT(4) CONTAINS THE RESPONSE VARIABLE MEAN
C                           (NOT REQUIRED IF IOP(3) IS NOT EQUAL TO 0).
C                         CRIT(5) CONTAINS THE ERROR MEAN SQUARE.
C                         CRIT(6) CONTAINS THE CORRECTED SUM OF SQUARES
C                           FOR THE INDEPENDENT VARIABLE (UNCORRECTED
C                           IF IOP(3) IS NOT EQUAL TO 0).
C                         CRIT(7) CONTAINS THE CONFIDENCE COEFFICIENT
C                           FOR THE INTERVAL ESTIMATE.
C                         CRIT(8) CONTAINS THE THE LIMIT VALUE OF THE
C                           RESPONSE VARIABLE, WHEN ONE-SIDED CONTROL IS
C                           DESIRED. THE LOWER LIMIT OF THE RESPONSE
C                           VARIABLE WHEN TWO-SIDED CONTROL IS DESIRED.
C                         CRIT(9) CONTAINS THE THE UPPER LIMIT OF THE
C                           RESPONSE VARIABLE WHEN TWO-SIDED CONTROL IS
C                           DESIRED. OTHERWISE, CRIT(9) IS UNDEFINED.
C                IOP    - INPUT VECTOR OF LENGTH 5.
C                         IOP(1) SHOULD CONTAIN 0 IF CONTROL OF THE
C                           AVERAGE OF IOP(5) OBSERVATIONS IS DESIRED.
C                           IOP(1) SHOULD CONTAIN ANY OTHER VALUE IF
C                           CONTROL OF THE TRUE MEAN RESPONSE VALUE IS
C                           DESIRED.
C                         IOP(2) SHOULD CONTAIN 0 FOR TWO-SIDED CONTROL.
C                           IOP(2) SHOULD CONTAIN 1 FOR UPPER ONE-SIDED
C                           CONTROL.
C                           IOP(2) SHOULD CONTAIN ANY OTHER VALUE FOR
C                           LOWER ONE-SIDED CONTROL.
C                         IOP(3) SHOULD CONTAIN 0 IF THE FITTED MODEL
C                           HAS BOTH SLOPE AND INTERCEPT PARAMETERS.
C                           IOP(3) SHOULD CONTAIN ANY OTHER VALUE IF
C                           THE FITTED MODEL HAS SLOPE PARAMETER ONLY.
C                         IOP(4) CONTAINS THE NUMBER OF DATA POINTS USED
C                           IN FITTING THE MODEL.
C                         IOP(5) CONTAINS THE NUMBER OF OBSERVATIONS ON
C                           THE RESPONSE WHOSE MEAN IS TO BE CONTROLLED
C                           (REQUIRED ONLY WHEN IOP(1) = 0).
C                STAT   - OUTPUT VECTOR OF LENGTH 2 CONTAINING THE
C                           LIMITS ON THE INDEPENDENT VARIABLE FOR
C                           CONTROLLING THE RESPONSE VARIABLE.
C                         STAT(1) CONTAINS THE LOWER LIMIT (NOT DEFINED
C                           WHEN IOP(2) = 1 AND CRIT(2) IS POSITIVE OR
C                           WHEN IOP(2) IS NOT EQUAL TO 0 OR 1 AND
C                           CRIT(2) IS NEGATIVE).
C                         STAT(2) CONTAINS THE UPPER LIMIT (NOT DEFINED
C                           WHEN IOP(2) = 1 AND CRIT(2) IS NEGATIVE OR
C                           WHEN IOP(2) IS NOT EQUAL TO 0 OR 1 AND
C                           CRIT(2) IS POSITIVE).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 IMPLIES AN ERROR OCCURRED IN MDSTI.
C                           IER=130 IMPLIES THAT THE SLOPE WAS NOT
C                             SIGNIFICANT AT LEVEL (1-CRIT(7)) FOR
C                             TWO-SIDED CONTROL OR (2(1-CRIT(7)))
C                             FOR ONE-SIDED CONTROL.
C                           IER=131 IMPLIES, WHEN IOP(2) =0, THAT
C                             STAT(1) EXCEEDS STAT(2). NO SATISFACTORY
C                             SETTING FOR THE INDEPENDENT VARIABLE
C                             EXISTS TO CONTROL THE RESPONSE
C                             VARIABLE AS DESIRED.
C                           IER=132 IMPLIES THE SLOPE ESTIMATE (CRIT(2))
C                             WAS ZERO.
C
C   REQD. IMSL ROUTINES - MDSTI,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLINCF (CRIT,IOP,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOP(5),IER
      REAL               CRIT(9),STAT(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               QQ,RNDF,X
      REAL               ZERO,HALF,ONE,TWO,RN,RNOB,XL,XP,XQ1,XQ2,
     1                   XR1,XR2,XS1,XS2,Q,CRIT89,CRIT3,CRIT4
      DATA               ZERO/0.0/,HALF/0.5/,ONE/1.0/,TWO/2.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (CRIT(2) - ZERO) 10,5,10
C                                  TERMINAL - SLOPE ESTIMATE = ZERO
    5 IER = 132
      GO TO 9000
   10 IER = 0
      STAT(1) = ZERO
      STAT(2) = ZERO
      CRIT3 = ZERO
      CRIT4 = ZERO
      IF (IOP(3) .NE. 0) GO TO 15
      CRIT3 = CRIT(3)
      CRIT4 = CRIT(4)
   15 RN   = IOP(4)
      IF (IOP(1) .EQ. 0) RNOB = IOP(5)
      RNDF = RN - ONE
      IF (IOP(3) .EQ. 0) RNDF = RNDF - ONE
C                                  GET INVERSE STUDENTS T
      X = ZERO
      IF (CRIT(7) .EQ. HALF .AND. IOP(2) .NE. 0) GO TO 30
      Q = ONE - CRIT(7)
      IF (IOP(2) .EQ. 0) GO TO 20
      Q = TWO*Q
      IF (CRIT(7) .GT. ZERO .AND. CRIT(7) .LT. HALF) Q = -Q + TWO
   20 QQ=Q
      CALL MDSTI(QQ,RNDF,X,IER)
      IF (IER .EQ. 0) GO TO 30
C                                  TERMINAL - ERROR IN MDSTI
      IER = 129
      GO TO 9000
   30 XL = CRIT(2) **2 - X**2*CRIT(5)/CRIT(6)
      IF (XL .GT. ZERO) GO TO 35
C                                  TERMINAL - SLOPE INSIGNIFICANT
      IER = 130
      GO TO 9000
   35 IF (IOP(1) .EQ. 0) GO TO 40
C                                  TRUE MEAN CONTROL OPTION
      XP = ZERO
      IF (IOP(3) .EQ. 0) XP = XL/RN
      GO TO 45
C                                  MEAN OF IOP(5) OBSERVATIONS OPTION
   40 XP = XL/RNOB
      IF (IOP(3) .EQ. 0) XP = XP * (RNOB/RN+ONE)
   45 IF (IOP(2) .EQ. 1) GO TO 55
      XQ1 = CRIT(8) - CRIT4
      XR1 = (XP + XQ1**2/CRIT(6))*CRIT(5)
      XS1 = CRIT3 + CRIT(2)*XQ1/XL
      XR1 = (X/XL)*SQRT(XR1)
      IF (CRIT(2) .GT. ZERO) GO TO 50
      STAT(2) = XS1 - XR1
      GO TO 55
   50 STAT(1) = XS1 + XR1
   55 IF (IOP(2) .NE. 0 .AND. IOP(2) .NE. 1) GO TO 65
      IF (IOP(2) .EQ. 2) CRIT89=CRIT(9)
      IF (IOP(2) .EQ. 1) CRIT89=CRIT(8)
      XQ2 = CRIT89  - CRIT4
      XR2 = (XP + XQ2**2/CRIT(6))*CRIT(5)
      XS2 = CRIT3 + CRIT(2)*XQ2/XL
      XR2 = (X/XL)*SQRT(XR2)
      IF (CRIT(2) .GT. ZERO) GO TO 60
      STAT(1) = XS2 + XR2
      GO TO 65
   60 STAT(2) = XS2 - XR2
C                                  TERMINAL - NO SATISFACTORY SETTING
C                                             FOR INDEPENDENT VARIABLE
   65 IF (IOP(2) .EQ. 0 .AND. STAT(1) .GT. STAT(2)) IER = 131
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST(IER,6HRLINCF)
 9005 RETURN
      END

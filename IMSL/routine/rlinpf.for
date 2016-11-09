
 
C   IMSL ROUTINE NAME   - RLINPF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM/SINGLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INVERSE PREDICTION USING A FITTED SIMPLE
C                           LINEAR REGRESSION MODEL
C
C   USAGE               - CALL RLINPF (CRIT,IOP,STAT,IER)
C
C   ARGUMENTS    CRIT   - INPUT VECTOR OF LENGTH 8.
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
C                         CRIT(8) CONTAINS THE VALUE OF THE RESPONSE
C                           VARIABLE FOR WHICH AN INTERVAL ESTIMATE OF
C                           THE CORRESPONDING INDEPENDENT VARIABLE
C                           VALUE IS DESIRED. WHEN IOP(1) = 0,
C                           CRIT(8) IS THE MEAN OF IOP(6)
C                           OBSERVATIONS.
C                IOP    - INPUT VECTOR OF LENGTH 6.
C                         IOP(1) SHOULD CONTAIN 0 IF AN INTERVAL
C                           ESTIMATE OF THE INDEPENDENT VARIABLE VALUE
C                           FOR THE IOP(6) OBSERVATIONS IS DESIRED.
C                           IOP(1) SHOULD CONTAIN ANY OTHER VALUE IF AN
C                           INTERVAL ESTIMATE OF THE TRUE MEAN
C                           INDEPENDENT VARIABLE VALUE FOR A GIVEN
C                           VALUE (CRIT(8)) OF THE RESPONSE VARIABLE IS
C                           DESIRED.
C                         IOP(2) SHOULD CONTAIN 0 FOR TWO-SIDED INTERVAL
C                           ESTIMATE.
C                           IOP(2) SHOULD CONTAIN 1 FOR UPPER ONE-SIDED
C                           INTERVAL ESTIMATE.
C                           IOP(2) SHOULD CONTAIN ANY OTHER VALUE FOR
C                           LOWER ONE-SIDED INTERVAL ESTIMATE.
C                         IOP(3) SHOULD CONTAIN 0 IF FITTED MODEL HAS
C                           BOTH SLOPE AND INTERCEPT PARAMETERS.
C                           IOP(3) SHOULD CONTAIN ANY OTHER VALUE IF
C                           FITTED MODEL HAS SLOPE PARAMETER ONLY.
C                         IOP(4) CONTAINS THE NUMBER OF DATA POINTS USED
C                           IN FITTING THE MODEL.
C                         IOP(5) CONTAINS THE NUMBER OF DEGREES OF
C                           FREEDOM.
C                         IOP(6) CONTAINS THE NUMBER OF OBSERVATIONS
C                           USED TO OBTAIN CRIT(8) (REQUIRED ONLY WHEN
C                           IOP(1) = 0).
C                STAT   - OUTPUT VECTOR OF LENGTH 3.
C                         STAT(1) CONTAINS THE POINT ESTIMATE OF THE
C                           INDEPENDENT VARIABLE.
C                         STAT(2) CONTAINS THE LOWER LIMIT FOR THE
C                           INDEPENDENT VARIABLE (DEFINED ONLY WHEN
C                           IOP(2) IS NOT EQUAL TO 0 OR 1 AND CRIT(2) IS
C                           POSITIVE, OR WHEN IOP(2) IS 1 AND
C                           CRIT(2) IS NEGATIVE).
C                         STAT(3) CONTAINS THE UPPER LIMIT FOR THE
C                           INDEPENDENT VARIABLE (DEFINED ONLY WHEN
C                           IOP(2) IS NOT EQUAL TO 0 OR 1 AND CRIT(2) IS
C                           NEGATIVE, OR WHEN IOP(2) IS 1 AND
C                           CRIT(2) IS POSITIVE).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 IMPLIES AN ERROR OCCURRED IN MDSTI.
C                             ONLY STAT(1) IS SET.
C                           IER=130 IMPLIES THAT THE SLOPE WAS NOT
C                             SIGNIFICANT AT LEVEL (1-CRIT(7)) FOR A
C                             TWO-SIDED ESTIMATE OR (2(1-CRIT(7)))
C                             FOR A ONE-SIDED ESTIMATE. ONLY
C                             STAT(1) IS SET.
C                           IER=131 IMPLIES THE SLOPE ESTIMATE (CRIT(2))
C                             WAS ZERO.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - MDSTI,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLINPF (CRIT,IOP,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOP(6),IER
      REAL               CRIT(8),STAT(3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               QQ,RNDF,X
      REAL               ONE,TWO,HALF,RN,RNOB,XL,XP,XQ,XR,XS,
     1                   ZERO,Q,CRIT1,CRIT3,CRIT4
      DATA               ZERO/0.0/,HALF/0.5/,ONE/1.0/,TWO/2.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (CRIT(2) - ZERO) 10,5,10
C                                  TERMINAL - SLOPE ESTIMATE = ZERO
   5  IER = 131
      GO TO 9000
  10  STAT(2) = ZERO
      STAT(3) = ZERO
      CRIT1 = ZERO
      CRIT3 = ZERO
      CRIT4 = ZERO
      IF (IOP(3) .NE. 0) GO TO 15
      CRIT1 = CRIT(1)
      CRIT3 = CRIT(3)
      CRIT4 = CRIT(4)
C                                  FLOAT NUMBER OF DATA POINTS
  15  RN   = IOP(4)
C                                  FLOAT DEGREES OF FREEDOM
      RNDF = IOP(5)
C                                  FLOAT NUMBER OF ADDITIONAL OBSER-
C                                    VATIONS
      IF (IOP(1) .EQ. 0) RNOB = IOP(6)
C                                  COMPUTE POINT ESTIMATE
      STAT(1) = (CRIT(8) - CRIT1)/CRIT(2)
C                                  NOW GET INVERSE STUDENTS T
      X = ZERO
      IF (CRIT(7) .EQ. HALF .AND. IOP(2) .NE. 0) GO TO 30
      Q = ONE - CRIT(7)
      IF (IOP(2) .EQ. 0) GO TO 20
      Q = TWO*Q
      IF (CRIT(7) .GT. ZERO .AND. CRIT(7) .LT. HALF) Q = -Q + TWO
   20 QQ = Q
      CALL MDSTI(QQ,RNDF,X,IER)
      IF (IER .EQ. 0) GO TO 30
C                                  TERMINAL - ERROR IN MDSTI
      IER = 129
      GO TO 9000
  30  XL = CRIT(2)**2 - X**2*CRIT(5)/CRIT(6)
      IF (XL .GT. ZERO) GO TO 35
C                                  TERMINAL - SLOPE NOT SIGNIFICANT
      IER = 130
      GO TO 9000
  35  IF (IOP(1) .EQ. 0) GO TO 40
C                                  TRUE MEAN ESTIMATE OPTION
      XP = ZERO
      IF (IOP(3) .EQ. 0) XP = XL/RN
      GO TO 45
C                                  MEAN OF IOP(6) OBSERVATIONS OPTION
  40  XP = XL/RNOB
      IF (IOP(3) .EQ. 0) XP = XP*(RNOB/RN+ONE)
  45  XQ = CRIT(8) - CRIT4
      XR = (XP + XQ**2/CRIT(6))*CRIT(5)
      XS = CRIT3 + CRIT(2) * XQ/XL
      XR = (X/XL)*SQRT(XR)
      IF (IOP(2) .NE. 1) STAT(2) = XS - XR
      IF (IOP(2) .EQ. 0 .OR. IOP(2) .EQ. 1) STAT(3) = XS + XR
 9000 CONTINUE
C                                  IF ERROR, PRINT MESSAGE
      IF (IER .NE. 0) CALL UERTST(IER,6HRLINPF)
 9005 RETURN
      END
 
R; T=0.05/0.44 22:14:29
IF ERROR, PRINT MESSAGE
      IF (IER .NE. 0) CALL UERTST(IER,6HRLINPF)
 9005 RETURN
      END
 
R; T=0.05/0.44 22:14
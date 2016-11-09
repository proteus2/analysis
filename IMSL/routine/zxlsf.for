C   IMSL ROUTINE NAME   - ZXLSF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - ONE-DIMENSIONAL MINIMIZATION OF A SMOOTH
C                           FUNCTION USING SAFEGUARDED QUADRATIC
C                           INTERPOLATION.
C
C   USAGE               - CALL ZXLSF(FUNC,X,STEP,BOUND,XACC,MAXFN,IER)
C
C   ARGUMENTS    FUNC   - A REAL FUNCTION SUBPROGRAM SUPPLIED BY THE
C                           USER.  FUNC MUST BE DECLARED EXTERNAL IN THE
C                           CALLING PROGRAM.  FUNC DEFINES THE FUNCTION
C                           TO BE MINIMIZED AND SHOULD BE OF THE
C                           FOLLOWING FORM
C                             REAL FUNCTION FUNC(X)
C                             REAL X
C                                .
C                                .
C                           WHERE X IS THE INDEPENDENT VARIABLE.  FUNC
C                           MUST NOT ALTER X.
C                X      - THE SCALAR VARIABLE OF THE CALCULATION.  ON
C                           INPUT IT IS THE INITIAL GUESS AT THE
C                           MINIMUM.  ON OUTPUT IT IS THE VALUE OF X
C                           THAT GIVES THE LEAST CALCULATED VALUE OF
C                           FUNC.
C                STEP   - AN ORDER OF MAGNITUDE ESTIMATE OF THE REQUIRED
C                           CHANGE IN X.  (INPUT)
C                BOUND  - A LIMIT, WHICH MUST BE SET TO A POSITIVE
C                           NUMBER, ON THE AMOUNT BY WHICH X MAY BE
C                           CHANGED FROM ITS INITIAL VALUE.  (INPUT)
C                XACC   - THE REQUIRED ABSOLUTE ACCURACY IN THE FINAL
C                           VALUE OF X.  ON A NORMAL RETURN THERE ARE
C                           POINTS ON EITHER SIDE OF X WITHIN A DISTANCE
C                           XACC AT WHICH FUNC IS NO LESS THAN FUNC(X).
C                           (INPUT)
C                MAXFN  - A LIMIT, WHICH MUST BE SET TO A POSITIVE
C                           INTEGER, ON THE NUMBER OF CALLS TO FUNC.
C                           (INPUT)
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                         IER=129 THE VALUE OF BOUND IS SO SMALL THAT NO
C                           CHANGE HAS BEEN MADE TO X.
C                         IER=130 THE FINAL VALUE OF X IS AT A BOUND.
C                           THE MINIMUM IS PROBABLY BEYOND THE BOUND.
C                         IER=131 SUBROUTINE FUNC HAS BEEN CALLED MAXFN
C                           TIMES.
C                         IER=132 COMPUTER ROUNDING ERRORS PREVENT
C                           FURTHER REFINEMENT OF X.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXLSF (FUNC,X,STEP,BOUND,XACC,MAXFN,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MAXFN,IER
      REAL               FUNC,X,STEP,BOUND,XACC
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            INFO,IS,NF
      REAL               ACC,DA,DB,DC,F,FA,FB,FC,FD,H,RHO,SIGMA,ST,TOL,
     *                   X0,XA,XB,XC,XD,XX,TEMP
C                                  FIRST EXECUTABLE STATEMENT
      NF = 0
      GO TO 35
    5 X0 = X
      XB = X
      FB = F
C
C                                  ALTER ANY UNSUITABLE INITIAL
C                                    PARAMETERS, INCLUDING INCREASING
C                                    THE STEP IF NECESSARY IN ORDER THAT
C                                    THE COMPUTER PRECISION GIVES A
C                                    CHANGE IN X.
C
      ACC = AMAX1(0.0,XACC)
      IF (BOUND.LE.ACC) GO TO 135
      ST = STEP
      IF (ABS(ST).GT.BOUND) ST = SIGN(BOUND,ST)
      IF (ST.EQ.0.0) ST = 0.01*BOUND
      IF (ABS(ST).LT.ACC) ST = SIGN(ACC,ST)
      GO TO 15
   10 ST = 5.0*ST
   15 XA = XB+0.5*ABS(ST)
      XC = XB+ABS(ST)
      IF (XA.LE.XB) GO TO 10
      IF (XC.LE.XA) GO TO 10
      IF (BOUND.LT.ABS(ST)) GO TO 135
C
C                                  CALCULATE THE NEXT TRIAL VALUE OF X
C                                  FROM XB AND ST.
C
   20 IS = 1
   25 X = XB+ST
      H = XB+1.5*ST
      IF (ABS(H-X0).LT.BOUND) GO TO 30
      X = X0+SIGN(BOUND,ST)
      IS = 2
C
C                                  CALCULATE THE NEXT VALUE OF F.
C
   30 IF (NF.GE.MAXFN) GO TO 145
   35 NF = NF+1
      F = FUNC(X)
      IF (NF.LE.1) GO TO 5
C
C                                  REVERSE ST IF THE INITIAL STEP SEEMS
C                                  TO BE UPHILL.
C
      IF (NF.GE.3) GO TO 40
      IF (F.LT.FB) GO TO 45
      ST = -ST
      XC = X
      FC = F
      GO TO 20
C
C                                  ENSURE THAT FB IS THE LEAST
C                                  CALCULATED VALUE OF F.
C
   40 XD = XC
      FD = FC
      IF (F.GE.FB) GO TO 60
   45 XC = XB
      FC = FB
      XB = X
      FB = F
C
C                                  CALCULATE AN EXTRA FUNCTION VALUE IF
C                                  X IS AT A BOUND.
C
      IF (IS.GE.4) GO TO 70
      IF (IS.LE.1) GO TO 55
      IF (IS.EQ.3) GO TO 50
      IS = 3
      H = AMAX1(0.9*ACC,0.01*ABS(XB-XC))
      X = XB+SIGN(H,XC-XB)
      IF (ABS(X-XC).LT.ABS(X-XB)) X = 0.5*(XB+XC)
      TEMP = (XB-X)/(XB-XC)
      IF (TEMP.LE.0.0) GO TO 140
      GO TO 30
C
C                                  THIS STAGE IS REACHED WHEN A BRACKET
C                                  IS FOUND NEAR A BOUND.
C
   50 XA = XD
      FA = FD
      IS = 4
      GO TO 85
C
C                                  CALCULATE THE NEXT STEP IN THE SEARCH
C                                  FOR A BRACKET, TRYING TO
C                                    OVERSHOOT THE PREDICTED MINIMUM BY
C                                    THE FACTOR THREE.
C
   55 IF (NF.LE.2) GO TO 25
      RHO = (XB-XC)/(XB-XD)
      SIGMA = (FB-FC)/(FB-FD)
      H = 9.0
      IF (SIGMA.LT.RHO) H = 1.5*(RHO-SIGMA/RHO)/(SIGMA-RHO)
      H = AMIN1(H,9.0)
      H = AMAX1(H,2.0)
      ST = H*ST
      GO TO 25
C
C                                  RETURN IF THE MINIMUM SEEMS TO BE
C                                  BEYOND A BOUND.
C
   60 IF (IS.GE.4) GO TO 65
      IF (IS.EQ.3) GO TO 140
C
C                                  A BRACKET HAS BEEN FOUND SO BRANCH TO
C                                  PREDICT THE MINIMUM.
C
      XA = X
      FA = F
      IS = 4
      GO TO 85
C
C                                  ENSURE THAT XA, XB, XC AND XD ARE
C                                  ORDERED MONOTONICALLY.
C
   65 XC = X
      FC = F
   70 TEMP = (XB-XC)/(XA-XD)
      IF (TEMP.GT.0.0) GO TO 80
   75 H = XA
      XA = XD
      XD = H
      H = FA
      FA = FD
      FD = H
C
C                                  IF THERE ARE THREE CONSECUTIVE EQUAL
C                                  VALUES OF F, ENSURE THAT FB IS
C                                    THE MIDDLE ONE.
C
   80 IF (FA.EQ.FB) GO TO 85
      IF (FD.NE.FB) GO TO 85
      IF (FC.NE.FB) GO TO 85
      XC = XB
      XB = X
      FC = FB
      FB = F
      GO TO 75
C
C                                  USE THE MINIMA OF TWO QUADRATICS TO
C                                  CALCULATE A TOLERANCE ON THE NEXT
C
   85 DA = (FB-FA)/(XB-XA)
      DB = (FC-FB)/(XC-XB)
      IF (IS.GE.5) GO TO 90
      IS = 5
      TOL = 0.01*ABS(XA-XC)
      GO TO 95
   90 DC = (FD-FC)/(XD-XC)
      TOL = 0.0
      IF (DB.EQ.0.0) GO TO 95
      TOL = ABS(XA-XC)
      TEMP = (DC-DB)/(XA-XC)
      IF (TEMP.GE.0.0) GO TO 95
      H = 0.5*ABS(DB*((XD-XB)/(DC-DB)+(XA-XC)/(DB-DA)))
      IF (H.LT.TOL) TOL = H
   95 TOL = AMAX1(TOL,0.9*ACC)
C
C                                  SET X TO THE VALUE THAT MINIMIZES THE
C                                  INTERPOLATING QUADRATIC.
C
      X = XB
      IF (DA.NE.DB) X = 0.5*(DA*(XB+XC)-DB*(XA+XB))/(DA-DB)
C
C                                  ENSURE THAT ABS(XA-XB).GE.ABS(XB-XC).
C
      TEMP = (XA-XB)/(XB-XC)
      IF (TEMP.GE.1.0) GO TO 100
      H = XA
      XA = XC
      XC = H
      H = FA
      FA = FC
      FC = H
C
C                                  TEST FOR CONVERGENCE.
C
  100 IF (ABS(XA-XB).LE.ACC) GO TO 150
C
C                                  IF ABS(XA-XB).LE.2.9*ACC, CHOOSE THE
C                                  NEXT X TO AVOID AN UNNECESSARY
C                                    FUNCTION EVALUATION.
C
      TEMP = (XA-XB)/(XB-XC)
      IF (TEMP.GT.10.0) GO TO 110
      IF (ABS(XA-XB).GT.2.9*ACC) GO TO 105
      X = 0.5*(XA+XB)
      IF (ABS(X-XB).LE.ACC) GO TO 115
      X = 0.67*XB+0.33*XA
      GO TO 115
C
C                                  IF (XA-XB)/(XB-XC).LE.10, ENSURE THAT
C                                  THE DISTANCE FROM X TO XB IS
C                                    NORMALLY AT LEAST TOL.
C
  105 IF (ABS(X-XB).GE.TOL) GO TO 115
      X = XB+SIGN(TOL,X-XB)
      IF (ABS(X-XC).LT.ABS(X-XB)) X = XB+SIGN(TOL,XA-XB)
      IF (ABS(X-XA).LT.ABS(X-XB)) X = 0.5*(XA+XB)
      GO TO 115
C
C                                  WHEN (XA-XB)/(XB-XC).GT.10, ENSURE
C                                  THAT X IS IN THE LONGER INTERVAL,
C                                    AND TRY TO OVERSHOOT THE MINIMUM.
C
  110 H = (X-XB)*SIGN(3.0,XA-XB)
      H = AMAX1(H,TOL,ABS(XB-XC))
      H = AMIN1(H,0.1*ABS(XA-XB))
      X = XB+SIGN(H,XA-XB)
C
C                                  TEST WHETHER ROUNDING ERRORS MAKE THE
C                                  NEW X THE SAME AS A PREVIOUS X.
C
  115 IF (X.EQ.XB) GO TO 120
      IF ((XA-X)/(XA-XB).LE.0.0) GO TO 120
      IF ((X-XC)/(XA-XB).GT.0.0) GO TO 30
C
C                                  SEEK A NEW VALUE OF X BETWEEN XA AND
C                                  XB, BUT ROUNDING ERRORS MAY NOT
C                                    ALLOW ONE.
C
  120 X = XA
  125 XX = 0.5*(X+XB)
      TEMP = (XX-XB)/(XA-XB)
      IF (TEMP.LE.0.0) GO TO 130
      TEMP = (X-XX)/(XA-XB)
      IF (TEMP.LE.0.0) GO TO 130
      X = XX
      GO TO 125
  130 IF (X.NE.XA) GO TO 30
C
C                                  RETURN FROM THE SUBROUTINE.
C
      INFO = 5
      GO TO 155
  135 INFO = 2
      GO TO 155
  140 INFO = 3
      GO TO 155
  145 INFO = 4
      GO TO 155
  150 INFO = 1
  155 X = XB
      IER = 0
      IF (INFO.EQ.1) RETURN
      IER = 127+INFO
      CALL UERTST(IER,6HZXLSF )
      RETURN
      END

C   IMSL ROUTINE NAME   - MDCHI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INVERSE CHI-SQUARED PROBABILITY DISTRIBUTION
C                           FUNCTION
C
C   USAGE               - CALL MDCHI (P,DF,X,IER)
C
C   ARGUMENTS    P      - INPUT PROBABILITY IN THE EXCLUSIVE RANGE
C                           (0,1)
C                DF     - INPUT NUMBER OF DEGREES OF FREEDOM. DF MUST BE
C                           IN THE EXCLUSIVE RANGE (.5,200000.).
C                X      - OUTPUT CHI-SQUARED VALUE, SUCH THAT A RANDOM
C                           VARIABLE, DISTRIBUTED AS CHI-SQUARED WITH
C                           DF DEGREES OF FREEDOM, WILL BE LESS THAN OR
C                           EQUAL TO X WITH PROBABILITY P.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE BOUNDS WHICH
C                             ENCLOSED P COULD NOT BE FOUND WITHIN 20
C                             (NCT) ITERATIONS
C                           IER = 130 INDICATES THAT AN ERROR OCCURRED
C                             IN IMSL ROUTINE MDNRIS (P IS NOT A VALID
C                             VALUE)
C                           IER = 131 INDICATES THAT AN ERROR OCCURRED
C                             IN IMSL ROUTINE MDCH
C                           IER = 132 INDICATES THAT THE VALUE X COULD
C                             NOT BE FOUND WITHIN 50 (ITMAX) ITERATIONS,
C                             SO THAT THE ABSOLUTE VALUE OF P1-P WAS
C                             LESS THAN OR EQUAL TO EPS. (P1 IS THE
C                             CALCULATED PROBABILITY AT X, EPS = .0001)
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MDNRIS,MERFI,MERRC=ERFC,
C                           MGAMAD=DGAMMA,UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MDNRIS,MERFI,
C                           MERRC=ERFC,MGAMA=GAMMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDCHI   (P,DF,X,IER)
C
      DATA               EPS/.0001/,ITMAX/50/,NCT/20/,NSIG/5/
C                                  FIRST EXECUTABLE STATEMENT
      IC = 0
      IER = 0
C                                  ESTIMATE STARTING X
      CALL MDNRIS (P,XP,IER)
      IF (IER .NE. 0) GO TO 80
      DFF = .2222222
      DFF = DFF/DF
      D = SQRT(DFF)
      X  = DF*(1. -DFF + XP*D)**3
C                                  IS THE CASE ASYMPTOTICALLY NORMAL
      IF (DF .GE. 40.) GO TO 9005
C                                  FIND BOUNDS (IN X) WHICH ENCLOSE P
      NCNT = 0
      IST = 0
      ISW = 0
      DX = X * .125
    5 IF (X) 10,15,20
   10 X = 0.0
      DX = -DX
      GO TO 20
   15 DX = .1
   20 CALL MDCH (X,DF,P1,IER)
      DX = DX + DX
      IF (IER .NE. 0) GO TO 85
      CSS = X
      NCNT = NCNT + 1
      IF (NCNT .GT. NCT) GO TO 90
      IF (P1 - P) 25,9005,30
   25 X = X + DX
      ISW = 1
      IF (IST .EQ. 0) GO TO 5
      GO TO 35
   30 X = X - DX
      IST = 1
      IF (ISW .EQ. 0) GO TO 5
      XR = CSS
      XL = X
      GO TO 40
   35 XL = CSS
      XR = X
C                                  PREPARE FOR ITERATION TO FIND X
   40 EPSP = 10.**(-NSIG)
      IF (XL .LT. 0.) XL = 0.0
      CALL MDCH (XL,DF,P1,IER)
      IF (IER .NE. 0) GO TO 85
      FXL = P1 - P
      CALL MDCH (XR,DF,P1,IER)
      IF (IER .NE. 0) GO TO 85
      FXR = P1-P
      IF (FXL*FXR .NE. 0.0) GO TO 45
      X = XR
      IF (FXL .EQ. 0.0) X = XL
      GO TO 9005
   45 IF (DF .LE. 2. .OR. P .GT. .98) GO TO 50
C                                  REGULA FALSI METHOD
      X = XL + FXL*(XR-XL)/(FXL-FXR)
      GO TO 55
C                                  BISECTION METHOD
   50 X = (XL+XR) * .5
   55 CALL MDCH (X,DF,P1,IER)
      IF (IER .NE. 0) GO TO 85
      FCS = P1-P
      IF (ABS(FCS) .GT. EPS) GO TO 60
      GO TO 9005
   60 IF (FCS * FXL .GT. 0.0) GO TO 65
      XR = X
      FXR = FCS
      GO TO 70
   65 XL = X
      FXL = FCS
   70 IF (XR-XL .GT. EPSP*ABS(XR)) GO TO 75
      GO TO 9005
   75 IC = IC+1
      IF (IC .LE. ITMAX) GO TO 45
      IER = 132
      GO TO 9000
C                                  ERROR RETURNED FROM MDNRIS
   80 IER = 130
      GO TO 9000
C                                  ERROR RETURNED FROM MDCH
   85 IER = 131
      GO TO 9000
   90 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMDCHI )
 9005 RETURN
      END

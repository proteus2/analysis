C   IMSL ROUTINE NAME   - ZFALSE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZERO OF A FUNCTION GIVEN AN INTERVAL
C                           CONTAINING THE ZERO
C
C   USAGE               - CALL ZFALSE (F,EPS,NSIG,XL,XR,XAPP,ITMAX,IER)
C
C   ARGUMENTS    F      - A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM
C                           SUPPLIED BY THE USER. (INPUT)
C                           F MUST BE DECLARED EXTERNAL IN THE CALLING
C                           PROGRAM. F DEFINES THE FUNCTION FOR WHICH
C                           THE ROOT IS TO BE FOUND.
C                EPS    - FIRST CONVERGENCE CRITERION. (INPUT)
C                           A ROOT, XAPP, IS ACCEPTED IF
C                           ABS(F(XAPP)) .LE. EPS.
C                NSIG   - SECOND CONVERGENCE CRITERION. (INPUT)
C                           A ROOT, XAPP, IS ACCEPTED IF TWO
C                           SUCCESSIVE ITERATES AGREE TO NSIG
C                           SIGNIFICANT DIGITS.
C                XL     - A VALUE KNOWN TO BE TO THE LEFT OF THE
C                           ROOT. (INPUT)
C                XR     - A VALUE KNOWN TO BE TO THE RIGHT OF THE
C                           ROOT. (INPUT)
C                         NOTE - THE USER SHOULD BE CERTAIN THAN F(XL)
C                           AND F(XR) DO NOT AGREE IN SIGN.
C                XAPP   - THE APPROXIMATE VALUE OF THE TRUE ROOT WHICH
C                           LIES BETWEEN XL AND XR. (OUTPUT)
C                ITMAX  - ITERATION INDICATOR. (INPUT AND OUTPUT)
C                           ON INPUT, ITMAX IS THE MAXIMUM NUMBER OF
C                           ITERATIONS TO BE TAKEN TO FIND THE ROOT.
C                           ON OUTPUT, ITMAX IS THE NUMBER OF ITERATIONS
C                           THE ROUTINE USED TO FIND THE ROOT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 35, XL GREATER THAN XR. COMPUTATIONS
C                             PERFORMED TREATING XR AS BEING LEFT OF
C                             THE ROOT AND XL AS BEING RIGHT OF THE
C                             ROOT.
C                         TERMINAL ERROR
C                           IER = 129, INDICATES F(XL) AND F(XR) HAVE
C                             THE SAME SIGN
C                           IER = 130, NUMBER OF ITERATIONS HAS REACHED
C                             THE MAXIMUM ALLOWABLE NUMBER.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZFALSE (F,EPS,NSIG,XL,XR,XAPP,ITMAX,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,ITMAX,IER
      REAL               F,EPS,XL,XR,XAPP
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IC
      REAL               EPSP,FXAPP,FPREV,FXL,FXR,XLL,XRR,ZERO,TEN,HALF
      DATA               ZERO/0.0/,TEN/10.0/,HALF/0.5/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IC = 0
      XLL = AMIN1(XL,XR)
      XRR = AMAX1(XL,XR)
      IF(XL .NE. XLL)IER = 35
      EPSP = TEN**(-NSIG)
      FXL = F(XLL)
      FPREV = FXL
      FXR = F(XRR)
      IF (FXL .EQ. ZERO .OR. FXR .EQ. ZERO) GO TO 10
      IF (SIGN(1.0,FXL) .NE. SIGN(1.0,FXR)) GO TO 15
C                                  TERMINAL ERROR
      IER = 129
      GO TO 40
   10 CONTINUE
C                                  FXL OR FXR = 0
      XAPP = XRR
      IF (FXL .EQ. ZERO) XAPP = XLL
      GO TO 40
   15 CONTINUE
C                                  COMPUTE APPROXIMATE ROOT
      XAPP = XLL+FXL*(XRR-XLL)/(FXL-FXR)
      FXAPP = F(XAPP)
      IF (ABS(FXAPP) .GT. EPS) GO TO 20
      GO TO 40
C                                  DETERMINE WHETHER APPROXIMATE ROOT
C                                  LIES BETWEEN XAPP AND XLL OR XAPP
C                                  AND XRR
   20 IF (FXAPP*FXL .GT. ZERO) GO TO 25
      XRR = XAPP
      FXR = FXAPP
      IF (FPREV*FXR .GT. ZERO) FXL = HALF*FXL
      FPREV = FXR
      GO TO 30
   25 XLL = XAPP
      FXL = FXAPP
      IF (FPREV*FXL .GT. ZERO) FXR = HALF*FXR
      FPREV = FXL
C                                  DETERMINE IF (XRR-XLL .GT. EPSP)
   30 IF (XRR-XLL .GT. EPSP*ABS(XRR)) GO TO 35
      GO TO 40
C                                  CONTINUE FOR ITMAX ITERATIONS
   35 IC = IC+1
      IF (IC .LE. ITMAX) GO TO 15
      IER = 130
   40 ITMAX = IC
      IF (IER .NE. 0) GO TO 9000
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HZFALSE)
 9005 RETURN
      END

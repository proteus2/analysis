C   IMSL ROUTINE NAME   - MDBETA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - BETA PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDBETA (X,A,B,P,IER)
C
C   ARGUMENTS    X      - VALUE TO WHICH FUNCTION IS TO BE INTEGRATED.
C                           X MUST BE IN RANGE (0,1) INCLUSIVE.(INPUT)
C                A      - INPUT (1ST) PARAMETER (MUST BE GREATER THAN 0)
C                B      - INPUT (2ND) PARAMETER (MUST BE GREATER THAN 0)
C                P      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE
C                           FROM A BETA DISTRIBUTION HAVING PARAMETERS
C                           A AND B WILL BE LESS THAN OR EQUAL TO X.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES X IS NOT IN RANGE (0,1)
C                             INCLUSIVE
C                           IER = 130 INDICATES A AND/OR B LESS THAN OR
C                             EQUAL TO 0
C
C   REQD. IMSL ROUTINES - H32/MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDBETA  (X,A,B,P,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               X,A,B,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               AA,BB,TEMP
      DOUBLE PRECISION   PS,PX,Y,P1,DA,XINT,CNT,WH,XB,DB,C,EPS,EPS1,
     *                   ALEPS,TOT,PQ,D4,ONE,ZERO
      DOUBLE PRECISION   DLGAMA
      INTEGER            INT1,IB
C                                  MACHINE PRECISION
      DATA               EPS/1.D-6/
C                                  SMALLEST POSITIVE NUMBER
C                                  SAFELY REPRESENTABLE
      DATA               EPS1/1.D-78/
C                                  NATURAL LOG OF EPS1
      DATA               ALEPS/-179.6016D0/
      DATA               ONE/1.0D0/,ZERO/0.0D0/
C                                  CHECK RANGES OF THE ARGUMENTS
C                                  FIRST EXECUTABLE STATEMENT
      Y = X
      IF ((X .LE. 1.0) .AND. (X .GE. 0.0)) GO TO 5
      IER = 129
      GO TO 9000
    5 IF ((A .GT. 0.0) .AND. (B .GT. 0.0)) GO TO 10
      IER = 130
      GO TO 9000
   10 IER = 0
      AA = A
      BB = B
      IF (X .GT. 0.5) GO TO 15
      INT1 = 0
      GO TO 20
C                                  SWITCH ARGUMENTS FOR MORE EFFICIENT
C                                  USE OF THE POWER SERIES
   15 INT1 = 1
      TEMP = AA
      AA = BB
      BB = TEMP
      Y = ONE-Y
   20 IF (X .NE. 0. .AND. X .NE. 1.) GO TO 25
C                                  SPECIAL CASE - X IS 0. OR 1.
      P = 0.E0
      GO TO 60
   25 TEMP = AINT(BB)
      PS = BB-TEMP
      IF (BB .EQ. TEMP) PS = ONE
      DA = AA
      DB = BB
      PX = DA*DLOG(Y)
      D4 = DLOG(DA)
      PQ = DLGAMA(DA+DB)
      P1 = DLGAMA(DA)
      C = DLGAMA(DB)
      XB = PX+DLGAMA(PS+DA)-DLGAMA(PS)-D4-P1
C                                  SCALING
      IB = XB/ALEPS
      XINT = ZERO
C                                  FIRST TERM OF A DECREASING SERIES
C                                  WILL UNDERFLOW
      IF (IB .NE. 0) GO TO 35
      XINT = DEXP(XB) * 1.D10
      CNT = XINT*DA
C                                  CNT WILL EQUAL DEXP(TEMP)*(1.D0-PS)I*
C                                  P*Y**I/FACTORIAL(I)
      WH = ZERO
   30 WH = WH+ONE
      CNT = CNT*(WH-PS)*Y/WH
      XB = CNT/(DA+WH)
      XINT = XINT+XB
      IF (XB/EPS .GT. XINT) GO TO 30
      XINT = XINT * 1.D-10
   35 TOT = ZERO
      IF (DB .LE. ONE) GO TO 55
      XB = PX+DB*DLOG(ONE-Y)+PQ-P1-DLOG(DB)-C
C                                  SCALING
      TEMP = XB/ALEPS
      TEMP = AINT(TEMP)
      IF(TEMP.LT.0.0E0) TEMP = 0.0E0
      C = ONE/(ONE-Y)
      PS = TEMP
      CNT = DEXP(XB-PS*ALEPS)
      PS = DB
      WH = DB
   40 WH = WH-ONE
      IF (WH .LE. ZERO) GO TO 55
      PX = (PS*C)/(DA+WH)
      IF (PX .GT. ONE) GO TO 45
      IF (CNT/EPS .LE. TOT .OR. CNT .LE. EPS1/PX) GO TO 55
   45 CNT = CNT*PX
      IF (CNT .LE. ONE) GO TO 50
C                                  RESCALE
      TEMP = TEMP - 1.0E0
      CNT = CNT*EPS1
   50 PS = WH
      IF (TEMP .EQ. 0.0E0) TOT = TOT+CNT
      GO TO 40
   55 P = TOT+XINT
   60 IF (INT1 .NE. 0) P = 1.-P
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMDBETA)
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - MDBNOR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - BIVARIATE NORMAL PROBABILITY DISTRIBUTION
C                           FUNCTION
C
C   USAGE               - CALL MDBNOR (X,Y,RHO,P,IER)
C
C   ARGUMENTS    X      - INPUT UPPER LIMIT OF INTEGRATION FOR THE
C                           FIRST VARIABLE
C                Y      - INPUT UPPER LIMIT OF INTEGRATION FOR THE
C                           SECOND VARIABLE
C                RHO    - INPUT CORRELATION COEFFICIENT
C                P      - OUTPUT PROBABILITY THAT THE FIRST VARIABLE
C                           IS LESS THAN OR EQUAL TO X AND THAT THE
C                           SECOND VARIABLE IS LESS THAN OR EQUAL TO Y
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                         IER = 129 INDICATES THE ABSOLUTE VALUE OF RHO
C                             IS GREATER THAN OR EQUAL TO ONE
C
C   REQD. IMSL ROUTINES - MDNOR,MDTNF,MERRC=ERFC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDBNOR (X,Y,RHO,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               X,Y,RHO,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IAX,IAY,IND
      REAL               C1,EPS,F1,XY,AX,AY,TY,TX,QX,QY
      DATA               C1/1.0E0/
      DATA               EPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (ABS(RHO) .LT. C1) GO TO 5
C                                  TERMINAL - RHO OUT OF RANGE
      IER = 129
      GO TO 9000
    5 F1 = 1.0/SQRT(1.0 - RHO**2)
      XY = X*Y
      IAX = 0
      IAY = 0
      IND = 0
      IF (XY .EQ. 0.) GO TO 10
      AX = F1*(Y/X - RHO)
      AY = F1*(X/Y - RHO)
      GO TO 25
   10 IF (X .NE. 0.) GO TO 15
      IF (Y .NE. 0.) GO TO 20
C                                                                2 1/2
C                                  FOR X=Y=0 AX=AY=(1-RHO)/(1-RHO )
      AX = F1*(1.0 - RHO)
      AY = AX
      GO TO 25
C                                  FOR Y=0,X LESS THAN 0     TY = -1/4
C                                  FOR Y=0,X GREATER THAN 0  TY =  1/4
   15 TY = 0.25
      IF (X .LT. 0.0) TY = -TY
      AX = -F1*RHO
      IND = 1
      GO TO 25
C                                  FOR X=0,Y LESS THAN 0     TX = -1/4
C                                  FOR X=0,Y GREATER THAN 0  TX =  1/4
   20 TX = 0.25
      IF (Y .LT. 0.0) TX = -TX
      AY = -F1*RHO
      GO TO 35
   25 IF (AX .GE. 0.0) GO TO 30
      IAX = 1
      AX = -AX
   30 CALL MDTNF(X,AX,EPS,TX)
      IF (IAX .NE. 0) TX = -TX
      IF (IND .NE. 0) GO TO 45
   35 IF (AY .GE. 0.0) GO TO 40
      IAY = 1
      AY = -AY
   40 CALL MDTNF(Y,AY,EPS,TY)
      IF (IAY .NE. 0) TY = -TY
   45 IF (X .GT. 0.0) GO TO 50
      CALL MDNOR(X,QX)
      GO TO 55
   50 CALL MDNOR(-X,QX)
      QX = 1.- QX
   55 IF (Y .GT. 0.0) GO TO 60
      CALL MDNOR(Y,QY)
      GO TO 65
   60 CALL MDNOR(-Y,QY)
      QY = 1.- QY
C                                  NOW EVALUATE P
   65 P = 0.5*(QX + QY) - TX - TY
      IF (XY .LE. 0.0 .AND.(XY .NE. 0.0 .OR. X+Y .LT. 0.0)) P = P - 0.5
      P = AMIN1(AMAX1(0.0,P),1.0)
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST(IER,6HMDBNOR)
 9005 RETURN
      END

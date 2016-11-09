C   IMSL ROUTINE NAME   - GGEXT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - RANDOM DEVIATE GENERATOR FOR A MIXTURE
C                           OF TWO EXPONENTIALS.
C
C   USAGE               - CALL GGEXT (DSEED,P,XM1,XM2,NR,R,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                P      - INPUT. MIXING PARAMETER P IN THE DENSITY
C                           FUNCTION.  P RELATES TO XM1, 1-P TO XM2.
C                XM1    - INPUT. MEAN OF 1ST EXPONENTIAL DISTRIBUTION.
C                           XM1 MUST BE POSITIVE AND .GE. XM2.
C                XM2    - INPUT. MEAN OF 2ND EXPONENTIAL DISTRIBUTION.
C                           XM2 MUST BE POSITIVE AND .LE. XM1.
C                NR     - INPUT. NUMBER OF DEVIATES TO BE GENERATED.
C                R      - OUTPUT.  VECTOR OF LENGTH NR CONTAINING THE
C                           MIXED EXPONENTIAL DEVIATES.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                           IER=130  IMPLIES THAT P IS NEGATIVE.
C                           IER=131  IMPLIES THAT P IS GREATER THAN
C                             XM1/(XM1-XM2).
C                           IER=132  IMPLIES THAT EITHER XM1 OR XM2 IS
C                             NEGATIVE.
C                           IER=133  IMPLIES XM2 IS GREATER THAN XM1.
C
C   REQD. IMSL ROUTINES - GGUBS,GGUBFS,GGEXN
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGEXT (DSEED,P,XM1,XM2,NR,R,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P,R(NR),XM1,XM2
      INTEGER            N,IER
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               A,B,C,PMIN,PMAX,U,U1,U2,U3,XM12,XM,ONE,ZERO
      INTEGER            I
C                                  DATA INITIALIZATION
      DATA               ONE/1.0/, ZERO/0.0/,
     1                   PMIN/0.015625/, PMAX/0.984375/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  ERROR CHECKING
      IF (XM2.GT.XM1) IER = 133
      IF (XM1.LT.ZERO.OR.XM2.LT.ZERO) IER = 132
      IF (IER.NE.0) GO TO 9000
      IF (XM1.EQ.XM2.OR.P.EQ.ONE) GO TO 5
      XM12 = XM1/(XM1-XM2)
      IF (P.GT.XM12.AND.P.GT.ONE) IER = 131
      IF (P.LT.ZERO) IER = 130
      IF (IER.NE.0) GO TO 9000
C                                  FILL R WITH UNIFORM RANDOM DEVIATES
      CALL GGUBS (DSEED,NR,R)
C                                  CONTROL LOGIC
      IF (P.GT.ZERO.AND.P.LT.ONE) GO TO 10
      IF (P.LE.XM12.AND.P.GT.ONE) GO TO 60
    5 XM = XM1
      IF (P.EQ.ZERO) XM = XM2
C                                  X IS EXPONENTIAL WITH MEAN XM
      CALL GGEXN (DSEED,XM,NR,R)
      RETURN
C                                  CASE I
   10 CONTINUE
      IF (P.LT.PMIN.OR.P.GT.PMAX) GO TO 35
C                                  CASE I.A
      B = -ALOG(P)
      C = -ALOG(ONE-P)
      DO 30 I=1,NR
         U = R(I)
         IF (P-U) 20,15,15
   15    R(I) = -XM1*(ALOG(U)+B)
         GO TO 25
   20    R(I) = -XM2*(ALOG(U-P)+C)
   25    CONTINUE
   30 CONTINUE
      RETURN
C                                  CASE I.B
   35 CONTINUE
      DO 55 I=1,NR
         U1 = R(I)
         U2 = GGUBFS(DSEED)
         IF (P-U1) 45,40,40
   40    R(I) = -XM1*ALOG(U2)
         GO TO 50
   45    R(I) = -XM2*ALOG(U2)
   50    CONTINUE
   55 CONTINUE
      RETURN
C                                  CASE II
   60 CONTINUE
      A = P - (P-ONE)*XM1/XM2
      A1 = ONE-A
      C = -ALOG(A1)
      IF (A.LT.PMIN.OR.A.GT.PMAX) GO TO 85
C                                  CASE II.A
      DO 80 I=1,NR
         U1 = R(I)
         U2 = GGUBFS(DSEED)
         IF (A1-U1) 70,65,65
   65    R(I) = -XM1*ALOG(U2)-XM2*(ALOG(U1)+C)
         GO TO 75
   70    R(I) = -XM1*ALOG(U2)
   75    CONTINUE
   80 CONTINUE
      RETURN
C                                  CASE II.B
   85 CONTINUE
      DO 105 I=1,NR
         U1 = R(I)
         U2 = GGUBFS(DSEED)
         U3 = GGUBFS(DSEED)
         IF (A1-U1) 95,90,90
   90    R(I) = -XM1*ALOG(U2)-XM2*ALOG(U3)
         GO TO 100
   95    R(I) = -XM1*ALOG(U2)
  100    CONTINUE
  105 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GGEXT ')
 9005 RETURN
      END

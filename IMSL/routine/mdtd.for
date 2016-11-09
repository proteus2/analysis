C   IMSL ROUTINE NAME   - MDTD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - STUDENTS T PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDTD (TVAL,DF,Q,IER)
C
C   ARGUMENTS    TVAL   - INPUT NON-NEGATIVE CONSTANT. MDTD COMPUTES
C                           THE PROBABILITY, Q, THAT TVAL WILL BE
C                           EXCEEDED IN ABSOLUTE VALUE.
C                DF     - INPUT DEGREES OF FREEDOM. DF MUST BE GREATER
C                           THAN OR EQUAL TO 1.
C                Q      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE
C                           FOLLOWING STUDENTS T DISTRIBUTION WILL
C                           EXCEED TVAL IN ABSOLUTE VALUE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT DF (DEGREES OF
C                             FREEDOM) IS LESS THAN 1.
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MERRC=ERFC,MLGAMD=DLGAMA,UERTST,
C                           UGETIO
C                       - H36,H48,H60/MDBETA,MERRC=ERFC,MLGAMA=ALGAMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDTD (TVAL,DF,Q,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               TVAL,DF,Q
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N
      REAL               AN,A,B,CON1,T,W,XJ,Y,Z
      DATA               CON1/.6366198/
C                                  FIRST EXECUTABLE STATEMENT
      IF (DF.LT.1.0) GO TO 50
      IER = 0
      T = TVAL*TVAL
      IF (DF.LE.T) GO TO 45
      T = TVAL
      AN = DF
      N = AN
      T = T*T
      Y = T/AN
      B = 1.0+Y
      IF (AN.NE.N .OR. (AN.GE.20.0 .AND. T.LT.AN) .OR. AN.GT.200.0) GO
     *TO 20
      IF (AN.LT.20.0 .AND. T.LT.4.0) GO TO 30
C                                  *TAIL* SERIES EXPANSION FOR LARGE
C                                      T-VALUES
      A = 1.0
      Y = AN
      XJ = 0.0
      Z = 0.0
    5 IF (A.EQ.Z) GO TO 10
      XJ = XJ+2.0
      Z = A
      Y = Y*(XJ-1.0)/(B*XJ)
      A = A+Y/(AN+XJ)
      GO TO 5
   10 IF (AN.LE.1.0 .OR. A.LT.1.0E-60) GO TO 15
      A = (AN-1.0)/(B*AN)*A
      AN = AN-2.0
      GO TO 10
   15 IF (AN.NE.0.0) A = SQRT(B)*CON1*A/B
      Q = A
      GO TO 9005
C                                  ASYMPTOTIC SERIES FOR LARGE AN
   20 W = B-1.0
      IF (W.NE.0.0) Y = Y*(ALOG(B)/W)
      A = AN-.5
      B = 48.0*A*A
      Y = Y*A
      Y = (((((-0.4*Y-3.3)*Y-24.0)*Y-85.5)/(0.8*(Y*Y)+100.0+B)+Y+3.0)/
     *B+1.0)*SQRT(Y)
      IF (Y.LT.18.8125) GO TO 25
      Q = 0.0
      GO TO 9005
   25 Q = ERFC(Y*.7071068)
      GO TO 9005
C                                  NESTED SUMMATION OF *COSINE* SERIES
   30 Y = SQRT(Y)
      A = Y
      IF (AN.EQ.1.0) A = 0.0
   35 AN = AN-2.0
      IF (AN.LE.1.0) GO TO 40
      A = (AN-1.0)/(B*AN)*A+Y
      GO TO 35
   40 IF (AN.EQ.0.0) A = A/SQRT(B)
      IF (AN.NE.0.0) A = (ATAN(Y)+A/B)*CON1
      Q = 1.0-A
      GO TO 9005
   45 T = DF/(DF+T)
      A = .5*DF
      B = .5
      CALL MDBETA (T,A,B,Q,IER)
      GO TO 9005
   50 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMDTD  )
 9005 RETURN
      END

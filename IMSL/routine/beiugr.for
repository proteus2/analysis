C   IMSL ROUTINE NAME   - BEIUGR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ESTIMATION OF  BASIC STATISTICAL PARAMETERS
C                           USING UNGROUPED DATA
C
C   USAGE               - CALL BEIUGR (Y,N,IOPT,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           RESPONSES. IF IOPT(4)=1, THE ELEMENTS OF Y
C                           MUST BE IN INCREASING ORDER OF MAGNITUDE.
C                N      - NUMBER OF RESPONSES IN THE SAMPLE. (INPUT)
C                IOPT   - INPUT VECTOR OF LENGTH FIVE INDICATING WHICH
C                           ESTIMATES ARE TO BE RETURNED.
C                           IF IOPT(1) = 1 THE ARITHMETIC MEAN ESTIMATE
C                             IS RETURNED.
C                           IF IOPT(2) = 1 THE GEOMETRIC MEAN ESTIMATE
C                             IS RETURNED.
C                           IF IOPT(3) = 1 THE HARMONIC MEAN ESTIMATE
C                             IS RETURNED.
C                           IF IOPT(4) = 1 THE MEDIAN ESTIMATE
C                             IS RETURNED.
C                           IF IOPT(5) = 1 THE VARIANCE ESTIMATE
C                             IS RETURNED.
C                STAT   - OUTPUT VECTOR OF LENGTH FIVE CONTAINING THE
C                           ESTIMATES INDICATED BY IOPT. IOPT(I) = 0
C                           IMPLIES STAT(I) IS SET EQUAL TO ZERO, FOR
C                           I = 1,2,...,5.
C                           STAT(1) = ARITHMETIC MEAN ESTIMATE
C                           STAT(2) = GEOMETRIC MEAN ESTIMATE
C                           STAT(3) = HARMONIC MEAN ESTIMATE
C                           STAT(4) = MEDIAN ESTIMATE
C                           STAT(5) = VARIANCE ESTIMATE
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 IMPLIES THE NUMBER OF INPUT
C                             OBSERVATIONS IS LESS THAN 2.
C                           IER=130 IMPLIES AT LEAST ONE OBSERVATION IS
C                             LESS THAN OR EQUAL TO ZERO WHEN IOPT(2)=1.
C                           IER=131 IMPLIES AT LEAST ONE OBSERVATION IS
C                             LESS THAN OR EQUAL TO ZERO WHEN IOPT(3)=1.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEIUGR (Y,N,IOPT,STAT,IER)
C
      INTEGER            IOPT(1)
      DOUBLE PRECISION   SUM
      REAL               Y(1),STAT(1),AY,AN,ZERO,ONE
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER=0
      AN = ONE/N
      DO 5 I=1,5
         STAT(I) = ZERO
    5 CONTINUE
C                                  CHECK FOR ERRORS
      DO 10 I=1,5
         IF (IOPT(I).EQ.1.AND.N.LT.2) GO TO 15
   10 CONTINUE
      GO TO 20
   15 IER=129
      GO TO 9000
   20 DO 30 I=1,N
         IF (Y(I) .GT. ZERO) GO TO 30
         IF (Y(I) .LE. ZERO) GO TO 25
         IF (IOPT(2).EQ.1) GO TO 40
         GO TO 30
   25    IF (IOPT(2).EQ.1) GO TO 40
         IF (IOPT(3).EQ.1) GO TO 35
   30 CONTINUE
      GO TO 45
   35 IER=131
      GO TO 9000
   40 IER=130
      GO TO 9000
   45 IF (IOPT(1).NE.1.AND.IOPT(5).NE.1) GO TO 55
C                                  COMPUTE ARITHMETIC MEAN
      SUM = 0.D0
      DO 50 I=1,N
         SUM=SUM+Y(I)
   50 CONTINUE
      STAT(1)=SUM*AN
   55 IF (IOPT(3).NE.1) GO TO 65
C                                  COMPUTE HARMONIC MEAN
      SUM = 0.D0
      DO 60 I=1,N
         AY = ONE/Y(I)
         SUM=SUM+AY
   60 CONTINUE
      STAT(3)=SUM*AN
      STAT(3) = ONE/STAT(3)
   65 IF (IOPT(4).NE.1) GO TO 75
C                                  COMPUTE MEDIAN
      N2=N/2
      IF (N.NE.(N2+N2)) GO TO 70
C                                  N IS EVEN
      STAT(4)=(Y(N2)+Y(N2+1))*.5
      GO TO 75
C                                  N IS ODD
   70 NN = (N+1)/2
      STAT(4) = Y(NN)
   75 IF (IOPT(5).NE.1) GO TO 85
C                                  COMPUTE VARIANCE
      SUM = 0.D0
      DO 80  I=1,N
         AY=(Y(I)-STAT(1))
         SUM=SUM+AY*AY
   80 CONTINUE
      STAT(5)=SUM/(N-1)
   85 IF (IOPT(2).NE.1) GO TO 9005
C                                  COMPUTE GEOMETRIC MEAN
      SUM = 0.D0
      DO 95 I=1,N
         SUM=SUM+ALOG(Y(I))
   95 CONTINUE
      STAT(2)=SUM*AN
      STAT(2)=EXP(STAT(2))
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'BEIUGR')
 9005 RETURN
      END

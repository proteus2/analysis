C   IMSL ROUTINE NAME   - BEIGRP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ESTIMATION OF BASIC STATISTICAL PARAMETERS
C                           USING GROUPED DATA
C
C   USAGE               - CALL BEIGRP (F,Y,K,YLM,WID,IOPT,STAT,IER)
C
C   ARGUMENTS    F      - INPUT VECTOR OF LENGTH K CONTAINING THE
C                           FREQUENCY FOR EACH CLASS INTERVAL.
C                Y      - VECTOR OF LENGTH K FOR WORK STORAGE.
C                         ON OUTPUT, Y CONTAINS THE CLASS MARKS.
C                K      - NUMBER OF GROUPS INTO WHICH RESPONSES
C                           ARE CLASSIFIED. (INPUT)
C                YLM    - LOWER LIMIT OF THE LOWEST CLASS INTERVAL.
C                           (INPUT)
C                WID    - WIDTH OF EACH CLASS INTERVAL. (INPUT)
C                IOPT   - INPUT VECTOR OF LENGTH EIGHT INDICATING WHICH
C                           ESTIMATES ARE TO BE RETURNED.
C                           IF IOPT(1) = 1 THE ARITHMETIC MEAN ESTIMATE
C                             IS RETURNED.
C                           IF IOPT(2) = 1 THE GEOMETRIC MEAN ESTIMATE
C                             IS RETURNED.
C                           IF IOPT(3) = 1 THE HARMONIC MEAN ESTIMATE IS
C                             RETURNED.
C                           IF IOPT(4) = 1 THE MEDIAN ESTIMATE IS
C                             RETURNED.
C                           IF IOPT(5) = 1 THE MODE ESTIMATE IS
C                             RETURNED.
C                           IF IOPT(6) = 1 THE VARIANCE ESTIMATE IS
C                             RETURNED.
C                           IF IOPT(7) = 1 THE THIRD CENTRAL MOMENT
C                             ESTIMATE IS RETURNED.
C                           IF IOPT(8) = 1 THE FOURTH CENTRAL MOMENT
C                             ESTIMATE IS RETURNED.
C                STAT   - OUTPUT VECTOR OF LENGTH EIGHT CONTAINING THE
C                           ESTIMATES INDICATED BY IOPT. IOPT(I)
C                           NOT EQUAL TO ONE IMPLIES STAT(I) IS SET
C                           EQUAL TO ZERO, I=1,2,...,8.
C                           STAT(1) = ARITHMETIC MEAN ESTIMATE
C                           STAT(2) = GEOMETRIC MEAN ESTIMATE
C                           STAT(3) = HARMONIC MEAN ESTIMATE
C                           STAT(4) = MEDIAN ESTIMATE
C                           STAT(5) = MODE ESTIMATE
C                           STAT(6) = VARIANCE ESTIMATE
C                           STAT(7) = THIRD CENTRAL MOMENT ESTIMATE
C                           STAT(8) = FOURTH CENTRAL MOMENT ESTIMATE
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 IMPLIES THE NUMBER OF INPUT GROUPS,
C                             K, IS LESS THAN 2.
C                           IER=130 IMPLIES YLM IS LESS THAN ZERO WHEN
C                             IOPT(2) = 1.
C                           IER=131 IMPLIES YLM IS LESS THAN ZERO WHEN
C                             IOPT(3) = 1.
C                         WARNING ERROR
C                           IER=36 IMPLIES THAT THE HIGHEST CLASS
C                             FREQUENCY OCCURS FOR TWO OR MORE CLASSES
C                             WHEN IOPT(5)=1.  THUS THE MODE ESTIMATE
C                             IS NOT RETURNED.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEIGRP (F,Y,K,YLM,WID,IOPT,STAT,IER)
C
      INTEGER            IOPT(1)
      REAL               F(1),STAT(1),Y(1)
      DOUBLE PRECISION   SUM
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER AND STAT VECTOR
      IER = 0
      DO 5  I=1,8
         STAT(I) = 0.0
    5 CONTINUE
C                                  CHECK FOR ERRORS
      IF (K .GE. 2) GO TO 10
      IER = 129
      GO TO 9000
   10 IF (IOPT(2) .EQ. 1 .AND. YLM .LT. 0.0) GO TO 15
      IF (IOPT(3) .EQ. 1 .AND. YLM .LT. 0.0) GO TO 30
      GO TO 20
   15 IER = 130
      GO TO 9000
C                                  EVALUATE THE CLASS MARKS ARRAY
   20 Y(1) = YLM + WID * .5
      KM1 = K-1
      DO 25  J=1,KM1
         Y(J+1) = Y(J) + WID
   25 CONTINUE
      GO TO 35
   30 IER = 131
      GO TO 9000
C                                  COMPUTE THE SUM OF ALL F
   35 AN = 0.0
      DO 40 I=1,K
         AN = F(I) + AN
   40 CONTINUE
      AND2 = AN*.5
      AN=1.0/AN
      IF (IOPT(1) .NE. 1 .AND. IOPT(6) .NE. 1 .AND. IOPT(7) .NE. 1 .AND.
     1    IOPT(8) .NE. 1) GO TO 50
C                                  COMPUTE ARITHMETIC MEAN
      SUM = 0.0
      DO 45 I=1,K
         SUM = SUM + F(I) * Y(I)
   45 CONTINUE
      STAT(1) = SUM*AN
   50 IF (IOPT(2) .NE. 1) GO TO 60
C                                  COMPUTE GEOMETRIC MEAN
      SUM = 0.0
      DO 55 I=1,K
         SUM = SUM + ALOG(Y(I)) * F(I)
   55 CONTINUE
      STAT(2) = SUM*AN
      STAT(2) = EXP(STAT(2))
   60 IF (IOPT(3) .NE. 1) GO TO 70
C                                  COMPUTE HARMONIC MEAN
      SUM = 0.0
      DO 65 I=1,K
         SUM = SUM + F(I)/Y(I)
   65 CONTINUE
      STAT(3) = SUM*AN
      STAT(3)  = 1.0 / STAT(3)
   70 IF (IOPT(4) .NE. 1) GO TO 85
C                                  COMPUTE THE MEDIAN
      FJ = 0.0
      DO 75 I=1,K
         FJ = FJ + F(I)
         IF (FJ .LT. AND2) GO TO 75
         J = I-1
         FJ = FJ - F(I)
         GO TO 80
   75 CONTINUE
   80 FJ1 = YLM + J * WID
      STAT(4) = FJ1 + WID * ((AND2 - FJ)/F(J+1))
   85 IF (IOPT(5) .NE. 1) GO TO 95
C                                  FIND THE MAXIMUM OF ALL F, F(J),
C                                  AND CHECK IF F(J) IS UNIQUE
      FMAX = 0.0
      J=1
      IEND=1+(K-1)
      DO 86 I=1,IEND
         T=ABS(F(I))
         IF (T .LE. FMAX) GO TO 86
         FMAX=T
         J=I
   86 CONTINUE
      J=(J-1)+1
      ICOUNT = 0
      DO 90 I=1,K
         IF (F(I) - FMAX) 90,87,90
   87    ICOUNT = ICOUNT + 1
   90 CONTINUE
      IF (ICOUNT .GT. 1) IER = 36
      STAT(5) = Y(J)
   95 IF (IOPT(6) .NE. 1) GO TO 105
C                                  COMPUTE VARIANCE
      SUM = 0.0
      DO 100 I=1,K
         SUM = SUM + F(I)*(Y(I)-STAT(1))*(Y(I)-STAT(1))
  100 CONTINUE
      STAT(6) = SUM*AN
  105 IF (IOPT(7) .NE. 1) GO TO 115
C                                  COMPUTE THIRD CENTRAL MOMENT
      SUM = 0.0
      DO 110 I=1,K
         SUM = SUM + F(I) * (Y(I) - STAT(1)) **3
  110 CONTINUE
      STAT(7) = SUM*AN
  115 IF (IOPT(8) .NE. 1) GO TO 125
C                                  COMPUTE FOURTH CENTRAL MOMENT
      SUM = 0.0
      DO 120 I=1,K
      SUM = SUM + F(I) * (Y(I)-STAT(1))**4
  120 CONTINUE
      STAT(8) = SUM*AN
  125 IF (IER .EQ. 0) GO TO 9005
      IF (IER.EQ.36) STAT(5)=0.0
 9000 CONTINUE
      CALL UERTST (IER,'BEIGRP')
 9005 RETURN
      END

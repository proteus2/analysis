C   IMSL ROUTINE NAME   - NHINC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INCLUDANCE TEST
C
C   USAGE               - CALL NHINC (X,N,Y,M,I1,I2,STAT,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING
C                           SAMPLE 1. ON OUTPUT, X IS
C                           ORDERED MONOTONICALLY INCREASING.
C                N      - INPUT LENGTH OF X, SAMPLE 1.
C                Y      - INPUT VECTOR OF LENGTH M CONTAINING
C                           SAMPLE 2.
C                M      - INPUT LENGTH OF Y, SAMPLE 2.
C                I1     - INPUT INDEX OF THE ELEMENT OF THE ORDERED X
C                           VECTOR WHICH IS TO BE THE LOW ENDPOINT OF
C                           THE RANGE CONSIDERED. I1 MUST BE GREATER
C                           THAN OR EQUAL TO 1 AND LESS THAN I2.
C                I2     - INPUT INDEX OF THE ELEMENT OF THE ORDERED X
C                           VECTOR WHICH IS TO BE THE HIGH ENDPOINT OF
C                           THE RANGE CONSIDERED. I2 MUST BE GREATER
C                           THAN I1 AND LESS THAN OR EQUAL TO N.
C                STAT   - INPUT/OUTPUT VECTOR OF LENGTH 4.
C                         ON INPUT, STAT(1) MUST CONTAIN EPS, WHICH
C                           IS USED TO JUDGE TIES. IF A Y ELEMENT IS
C                           WITHIN EPS OF X(I1) OR X(I2), A TIE WILL
C                           BE COUNTED. SEE REMARKS.
C                         ON OUTPUT, STAT CONTAINS THE RESULTANT
C                           STATISTICS FROM THE TEST.
C                         STAT(1) CONTAINS THE NUMBER OF TIES DETECTED.
C                         STAT(2) CONTAINS THE NUMBER OF (UNTIED)
C                           ELEMENTS OF Y LYING OUTSIDE
C                           THE RANGE (X(I1),X(I2)), INCLUSIVELY.
C                         STAT(3) CONTAINS THE PROBABILITY OF STAT(2) OR
C                           MORE ELEMENTS OF Y LYING OUTSIDE THE RANGE
C                           (X(I1),X(I2)) IF THE NULL HYPOTHESIS IS TRUE
C                           (COUNTING TIES AS BEING IN THE RANGE).
C                         STAT(4) CONTAINS THE PROBABILITY OF
C                           STAT(1)+STAT(2) OR MORE ELEMENTS OF Y LYING
C                           OUTSIDE THE RANGE (X(I1),X(I2)) IF THE NULL
C                           HYPOTHESIS IS TRUE (COUNTING TIES AS BEING
C                           OUTSIDE THE RANGE).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS THAT I1 IS NOT LESS THAN I2
C                           IER=130 MEANS I1 IS LESS THAN 1
C                           IER=131 MEANS THAT I2 IS GREATER THAN N
C                           IER=132 MEANS THAT AN ERROR OCCURRED
C                             IN ROUTINE MDHYP
C                         WARNING ERROR
C                           IER=37 MEANS THAT AT LEAST 1 TIE WAS
C                             DETECTED
C
C   REQD. IMSL ROUTINES - MDHYP,UERTST,UGETIO,VSRTA
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  ON INPUT TO NHINC, STAT(1) MUST CONTAIN A VALUE,
C                CALL IT EPS, WHICH IS USED TO DETERMINE THE
C                OCCURRENCE OF TIES. A TIE IS COUNTED WHEN EITHER
C                OF THE FOLLOWING CONDITIONS IS TRUE
C                (A)  Y(I) IS IN THE INTERVAL (X(I1)-EPS,X(I1)+EPS),
C                     INCLUSIVELY, FOR I=1,2,...,M.
C                (B)  Y(I) IS IN THE INTERVAL (X(I2)-EPS,X(I2)+EPS),
C                     INCLUSIVELY, FOR I=1,2,...,M.
C            2.  IF I1=1 AND I2=N, NHINC TESTS THE HYPOTHESIS THAT
C                THE Y POPULATION LIES, IN EQUAL PROPORTION TO THE
C                X POPULATION, BETWEEN THE ENDPOINTS OF THE X SAMPLE.
C            3.  IF I1=(N+1)/4 AND I2=3*(N+1)/4, THE FIRST AND THIRD
C                QUARTILE ESTIMATES OF THE X POPULATION ARE BEING
C                CONSIDERED. THE ACTUAL NULL HYPOTHESIS IS THAT X
C                AND Y ARE INDEPENDENT SAMPLES DRAWN FROM THE SAME
C                POPULATION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NHINC  (X,N,Y,M,I1,I2,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,I1,I2,IER
      REAL               X(1),Y(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IB,IER2,K,L,ND,NOUT,NTY
      REAL               AHIGH,ALOW,BHIGH,BLOW,EPS
      DOUBLE PRECISION   PEQK,PLEK
C                                  FIRST EXECUTABLE STATEMENT
      IF (I1 .LT. I2) GO TO 5
      IER = 129
      GO TO 9000
    5 IF (I1 .GE. 1) GO TO 10
      IER = 130
      GO TO 9000
   10 IF (I2 .LE. N) GO TO 15
      IER = 131
      GO TO 9000
   15 IER = 0
C                                  SORT X VECTOR
      CALL VSRTA(X,N)
      NTY = 0
      NOUT = 0
C                                  SPECIFY BOUNDS FOR TIES
      EPS = STAT(1)
      ALOW = X(I1) - EPS
      AHIGH = X(I1) + EPS
      BLOW = X(I2) - EPS
      BHIGH = X(I2) + EPS
C                                  COUNT THE ELEMENTS OF Y OUTSIDE OF
C                                    RANGE AND TIES
      DO 30 I=1,M
         IF ((Y(I) .LE. BHIGH) .AND. (Y(I) .GE. ALOW)) GO TO 20
         NOUT = NOUT + 1
         GO TO 30
   20    IF ((Y(I) .GE. BLOW) .AND. (Y(I) .LE. BHIGH)) GO TO 25
         IF ((Y(I) .LT. ALOW) .OR. (Y(I) .GT. AHIGH)) GO TO 30
   25    NTY = NTY + 1
   30 CONTINUE
      STAT(1) = NTY
      STAT(2) = NOUT
      IF (NTY .GT. 0) IER = 37
C                                  INVOKE MDHYP  FOR STAT(3) AND (4)
      L = N + M
C                                  CONSIDER TIES AS OUTSIDE RANGE
      I = 3
      IB = NOUT
   35 ND = N - I2 + I1 + IB
      K = IB - 1
      IF (K .GT. -1) GO TO 40
      STAT(I) = 1.0
      GO TO 50
   40 CALL MDHYP(K,M,L,ND,PEQK,PLEK,IER2)
C                                  CHECK FOR ERROR IN MDHYP
      IF (IER2 .EQ. 0) GO TO 45
      IER = 132
      GO TO 9000
   45 STAT(I) = 1.D0 - PLEK
   50 IF (I .EQ. 4) GO TO 55
C                                  CONSIDER TIES AS INSIDE RANGE
      I = 4
      IB = NOUT + NTY
      GO TO 35
   55 IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HNHINC )
 9005 RETURN
      END

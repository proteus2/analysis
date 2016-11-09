C   IMSL ROUTINE NAME   - NBSDL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE      NBSDL  - COX AND STUART SIGN TEST FOR TRENDS IN
C                           DISPERSION AND LOCATION
C
C   USAGE               - CALL NBSDL (IOPT,X,N,K,IDS,EPS,NSTAT,PSTAT,
C                           IER)
C
C   ARGUMENTS    IOPT   - INPUT OPTION PARAMETER.
C                         IF IOPT IS EQUAL TO ZERO, THE COX AND STUART
C                           TEST FOR TRENDS IN DISPERSION
C                           (THE S2 STATISTIC) IS INVOKED.
C                         IF IOPT IS NONZERO, THE COX AND STUART TEST
C                           FOR TRENDS IN LOCATION
C                           (THE S3 STATISTIC) IS INVOKED.
C                X      - INPUT/OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           OBSERVATIONS, IN TIME ORDER.
C                         IF THE ROUTINE IS CALLED WITH IOPT EQUAL TO
C                           ZERO, ON OUTPUT X WILL CONTAIN THE N/K
C                           DISPERSION STATISTICS FOR EACH OF THE K-SIZE
C                           SETS OF THE ORIGINAL X VECTOR.
C                         IF THE ROUTINE IS CALLED WITH IOPT NONZERO,
C                           THE X VECTOR IS NOT ALTERED UPON OUTPUT.
C                N      - INPUT NUMBER OF OBSERVATIONS.
C                K      - INPUT NUMBER OF CONSECUTIVE X ELEMENTS TO BE
C                           USED TO MEASURE DISPERSION.
C                           REQUIRED ONLY WHEN IOPT IS EQUAL TO ZERO.
C                IDS    - INPUT INDICATOR OF DISPERSION MEASURE
C                           TECHNIQUE.
C                           REQUIRED ONLY WHEN IOPT IS EQUAL TO ZERO.
C                         IF IDS IS EQUAL TO ZERO, THE RANGE WILL BE
C                           USED AS THE MEASURE OF DISPERSION.
C                         IF IDS IS NONZERO, THE SUMS OF SQUARES WILL
C                           BE USED AS THE MEASURE OF DISPERSION.
C                EPS    - INPUT VALUE TO BE USED TO DETERMINE WHEN
C                           VALUES IN X ARE EQUAL. IF THE ABSOLUTE VALUE
C                           OF THE DIFFERENCE BETWEEN TWO ELEMENTS OF X
C                           IS LESS THAN OR EQUAL TO EPS, A TIE IS
C                           COUNTED.
C                NSTAT  - OUTPUT VECTOR OF LENGTH 8.
C                           SEE THE MANUAL DOCUMENT FOR USAGE OF NSTAT.
C                         NSTAT(1) CONTAINS THE NUMBER OF NEGATIVE
C                           DIFFERENCES USED IN CALCULATING THE S2
C                           STATISTIC.
C                         NSTAT(2) CONTAINS THE NUMBER OF POSITIVE
C                           DIFFERENCES USED IN CALCULATING THE S2
C                           STATISTIC.
C                         NSTAT(3) CONTAINS THE NUMBER OF ZERO
C                           DIFFERENCES USED IN CALCULATING THE S2
C                           STATISTIC.
C                         NSTAT(4) CONTAINS THE NUMBER OF DIFFERENCES
C                           USED TO CALCULATE PSTAT(1) THROUGH PSTAT(4)
C                           (SEE THE DESCRIPTION OF PSTAT BELOW)
C                         NSTAT(5) CONTAINS THE NUMBER OF NEGATIVE
C                           DIFFERENCES USED IN CALCULATING THE S3
C                           STATISTIC.
C                         NSTAT(6) CONTAINS THE NUMBER OF POSITIVE
C                           DIFFERENCES USED IN CALCULATING THE S3
C                           STATISTIC.
C                         NSTAT(7) CONTAINS THE NUMBER OF ZERO
C                           DIFFERENCES USED IN CALCULATING THE S3
C                           STATISTIC.
C                         NSTAT(8) CONTAINS THE NUMBER OF DIFFERENCES
C                           USED TO CALCULATE PSTAT(5) THROUGH PSTAT(8)
C                           (SEE THE DESCRIPTION OF PSTAT BELOW)
C                PSTAT  - OUTPUT VECTOR OF LENGTH 8.
C                           SEE THE MANUAL DOCUMENT FOR USAGE OF PSTAT.
C                         PSTAT(1) THROUGH PSTAT(4) ARE ASSOCIATED
C                           WITH THE S2 STATISTIC.
C                         PSTAT(1) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(1)+NSTAT(3) OR MORE NEGATIVE SIGNS
C                           (TIES CONSIDERED NEGATIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         PSTAT(2) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(2) OR MORE POSITIVE SIGNS
C                           (TIES CONSIDERED NEGATIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         PSTAT(3) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(1)+NSTAT(3) OR MORE NEGATIVE SIGNS
C                           (TIES CONSIDERED POSITIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         PSTAT(4) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(2) OR MORE POSITIVE SIGNS
C                           (TIES CONSIDERED POSITIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         PSTAT(5) THROUGH PSTAT(8) ARE ASSOCIATED WITH
C                           THE S3 STATISTIC.
C                         PSTAT(5) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(5)+NSTAT(7) OR MORE NEGATIVE SIGNS
C                           (TIES CONSIDERED NEGATIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         PSTAT(6) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(6) OR MORE POSITIVE SIGNS
C                           (TIES CONSIDERED NEGATIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         PSTAT(7) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(5)+NSTAT(7) OR MORE NEGATIVE SIGNS
C                           (TIES CONSIDERED POSITIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         PSTAT(8) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(6) OR MORE POSITIVE SIGNS
C                           (TIES CONSIDERED POSITIVE)
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES NSTAT(4) IS LESS THAN 1
C                           IER=130 INDICATES NSTAT(8) IS LESS THAN 1
C                           IER=131 INDICATES AN ERROR OCCURRED
C                             IN IMSL ROUTINE MDBIN
C                         WARNING ERROR
C                           IER=36 INDICATES K IS GREATER THAN 5
C                           IER=37 INDICATES NSTAT(4) IS LESS THAN 8
C                           IER=38 INDICATES NSTAT(8) IS LESS THAN 8
C                           IER=39 INDICATES AT LEAST 1 TIE WAS DETECTED
C
C   REQD. IMSL ROUTINES - MDBIN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      SINCE IT IS POSSIBLE THAT MORE THAN ONE OF THE
C                CONDITIONS WHICH SET THE WARNING ERROR (IER WITH A
C                VALUE OF 36, 37, 38, OR 39), THE USER SHOULD BE
C                AWARE THAT IER WILL CONTAIN THE VALUE INDICATIVE
C                OF THE CONDITION DETECTED LAST DURING EXECUTION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NBSDL  (IOPT,X,N,K,IDS,EPS,NSTAT,PSTAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,N,K,IDS,NSTAT(1),IER
      REAL               X(1),EPS,PSTAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IND,J,JHIGH,JLOW,K1,L,NEG,NNEG,NNPOS,NPOS,
     1                   NTY,NUM,N1
      REAL               DIF,PEQK,PLEK,SUMX,SUMX2,XK
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IOPT .NE. 0) GO TO 34
C                                  TRENDS IN DISPERSION
C                                  CHECK K
      IF (K .GT. 5) IER=36
      N1 = N/K
C                                  CHECK TYPE OF DISPERSION MEASURE
      IF (IDS .EQ. 0) GO TO 15
C                                  SUMS OF SQUARES STATISTIC USED
      XK = K
      DO 10 I=1,N1
         JLOW = (I-1) * K + 1
         JHIGH = JLOW + K - 1
C                                  INITIALIZE SUMS FOR EACH SERIES
         SUMX = 0.0
         SUMX2 = 0.0
         DO 5 J=JLOW,JHIGH
            SUMX = SUMX + X(J)
            SUMX2 = SUMX2 + X(J)*X(J)
    5    CONTINUE
C                                  TRANSFER NEW STATISTIC TO X-VECTOR
         X(I) = SUMX2 - SUMX*SUMX/XK
   10 CONTINUE
      GO TO 35
C                                  RANGE STATISTIC USED
   15 DO 30 I=1,N1
         JLOW = (I-1) * K + 2
         JHIGH = JLOW + K - 2
C                                  INITIALIZE EXTREMES
         SUMX = X(JLOW-1)
         SUMX2 = X(JLOW-1)
         DO 25 J=JLOW,JHIGH
            IF (X(J) .GT. SUMX) GO TO 20
            IF (X(J) .GE. SUMX2) GO TO 25
            SUMX2 = X(J)
            GO TO 25
   20       SUMX = X(J)
   25    CONTINUE
C                                  TRANSFER RANGE STATISTIC TO X-VECTOR
         X(I) = SUMX - SUMX2
   30 CONTINUE
      GO TO 35
C
   34 CONTINUE
C                                  TRENDS IN LOCATION
      N1 = N
C                                  COMPUTE NSTAT(4). BEGIN S2 CALC.
   35 K1 = N1/2
      NSTAT(4) = K1
C                                  CHECK NSTAT(4)
      IF (NSTAT(4) .GE. 1) GO TO 40
      IER = 129
      GO TO 9000
C                                  CHECK SIZE OF N1
   40 IF (N1 .LT. 16) IER=37
C                                  CHECK FOR N1 ODD OR EVEN
      L = K1
      IF ((N1/2)*2 .NE. N1) L=K1+1
C                                  COUNT TIES AND NEGATIVE DIFFERENCES
      NTY = 0
      NEG = 0
      DO 50 I=1,K1
         DIF = X(I) - X(I+L)
         IF (ABS(DIF) .GT. EPS) GO TO 45
C                                  INCREMENT TIE COUNT
         NTY = NTY + 1
         GO TO 50
   45    IF (DIF .GE. 0.0) GO TO 50
C                                  INCREMENT NEGATIVE COUNT
         NEG = NEG + 1
   50 CONTINUE
      NUM = K1
      IND = 0
   55 NPOS = NUM - NTY - NEG
      NSTAT(IND+1) = NEG
      NSTAT(IND+2) = NPOS
      NSTAT(IND+3) = NTY
      XK = 0.5
C                                  COUNT TIES AS NEGATIVE DIFFERENCES
      NNEG = NUM - NPOS
      NNPOS = NPOS
      JLOW = IND + 1
      JHIGH = JLOW + 2
C                                  LOOP TO COMPUTE TWO PAIRS OF PROBS.
C                                  I IS INITIATED TWICE FOR 8 PROBS. TOT
      DO 80 I=JLOW,JHIGH,2
C                                  FIND MOST EFFICIENT BINOM. ARGUMENTS
         IF (NNEG .GT. NNPOS) GO TO 65
C                                  USE BINOM. DISTRIBUTION OF NEGATIVES
         CALL MDBIN(NNEG,NUM,XK,PLEK,PEQK,J)
         IF (J .LT. 128) GO TO 60
C                                  ERROR IN MDBIN
         IER = 131
         GO TO 9000
C                                  PROB OF THIS NUM. OR MORE OF +
   60    PSTAT(I+1) = PLEK
C                                  PROB OF THIS NUM. OR MORE OF -
         PSTAT(I) = 1.0 - PLEK + PEQK
         GO TO 75
C                                  USE BINOM. DISTRIBUTION OF POSITIVES
   65    CALL MDBIN(NNPOS,NUM,XK,PLEK,PEQK,J)
         IF (J .LT. 128) GO TO 70
C                                  ERROR IN MDBIN
         IER = 131
         GO TO 9000
C                                  PROB OF THIS NUM. OR MORE OF +
   70    PSTAT(I+1) = 1.0 - PLEK + PEQK
C                                  PROB OF THIS NUM. OR MORE OF -
         PSTAT(I) = PLEK
C                                  BYPASS SWITCHES FOR I=IHIGH
C                                  DO NOT RECYCLE WHEN NUM. OF TIES = 0
   75    IF ((I .EQ. JHIGH) .OR. (NTY .EQ. 0)) GO TO 85
C                                  SWITCH THE COUNT TO POSITIVE TALLY
         NNEG = NEG
         NNPOS = NUM - NEG
   80 CONTINUE
   85 IF (NTY .GT. 0) GO TO 90
C                                  PSTAT(3) AND PSTAT(4) OR PSTAT(7) AND
C                                  PSTAT(8) WHEN NUMBER OF TIES = 0
      PSTAT(JHIGH) = PSTAT(JLOW)
      PSTAT(JHIGH+1) = PSTAT(JLOW+1)
      GO TO 95
C                                  WARNING WHEN AT LEAST 1 TIE DETECTED
   90 IER = 39
   95 IF (IND .EQ. 4) GO TO 115
C                                  COMPUTE NSTAT(8). BEGIN S3 CALC.
      K1 = (2*N1)/3
      L = (N1 + 2)/3
      NSTAT(8) = L
C                                  CHECK NSTAT(8)
      IF (L .GE. 1) GO TO 100
      IER = 130
      GO TO 9000
C                                  CHECK FOR NSTAT(8) LESS THAN 8
  100 IF (L .LT. 8) IER=38
C                                  COUNT TIES AND NEGATIVE DIFFERENCES
      NTY = 0
      NEG = 0
      DO 110 I=1,L
         SUMX = X(I) - X(I+K1)
         IF (ABS(SUMX) .GT. EPS) GO TO 105
C                                  INCREMENT TIE COUNT
         NTY = NTY + 1
         GO TO 110
  105    IF (SUMX .GE. 0.0) GO TO 110
C                                  INCREMENT NEGATIVE COUNT
         NEG = NEG + 1
  110 CONTINUE
      NUM = L
      IND = 4
C                                  TRANSFER TO COMPUTE PROBABILITIES
      GO TO 55
  115 IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HNBSDL )
 9005 RETURN
      END

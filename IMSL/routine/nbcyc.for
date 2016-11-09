C   IMSL ROUTINE NAME   - NBCYC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NOETHERS TEST FOR CYCLICAL TREND
C
C   USAGE               - CALL NBCYC (X,N,EPS,NSTAT,P,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           OBSERVATIONS IN TIME ORDER
C                N      - INPUT NUMBER OF OBSERVATIONS. N MUST BE
C                           GREATER THAN OR EQUAL TO 3.
C                EPS    - INPUT VALUE TO BE USED IN DETERMINING WHEN
C                           VALUES IN X ARE EQUAL. IF THE ABSOLUTE VALUE
C                           OF THE DIFFERENCE BETWEEN TWO DIFFERENT
C                           VALUES OF X IS LESS THAN OR EQUAL TO EPS,
C                           THEN A TIE BETWEEN THOSE TWO ELEMENTS
C                           EXISTS.
C                NSTAT  - OUTPUT VECTOR OF LENGTH 6.
C                         NSTAT(1) CONTAINS THE NUMBER OF SEQUENCES
C                           USED TO DETERMINE THE PRESENCE OF CYCLICAL
C                           TREND, CONSIDERING THAT THE TYING
C                           OBSERVATION HAS BEEN IGNORED AND THE
C                           SEQUENCE SHORTENED.
C                         NSTAT(2) CONTAINS THE NUMBER OF MONOTONIC
C                           SEQUENCES IN THE SET DEFINED BY NSTAT(1).
C                         NSTAT(3) CONTAINS THE NUMBER OF MONOTONIC
C                           SEQUENCES WHERE TIED THREESOMES WERE COUNTED
C                           AS NONMONOTONIC.
C                         NSTAT(4) CONTAINS THE NUMBER OF MONOTONIC
C                           SEQUENCES WHERE TIED THREESOMES WERE COUNTED
C                           AS MONOTONIC.
C                         NSTAT(5) CONTAINS THE NUMBER OF TIES DETECTED
C                           IN FORMING THE NSTAT(1) SEQUENCES.
C                         NSTAT(6) CONTAINS THE NUMBER OF TIES DETECTED
C                           IN FORMING THE NSTAT(3) (AND NSTAT(4))
C                           SEQUENCES.
C                P      - OUTPUT VECTOR OF LENGTH 3.
C                         P(1) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(2) OR MORE MONOTONIC SEQUENCES
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                           IF NSTAT(1) IS LESS THAN 1, P(1) IS SET
C                           TO MACHINE INFINITY.
C                         P(2) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(3) OR MORE MONOTONIC SEQUENCES
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                         P(3) CONTAINS THE PROBABILITY OF OBTAINING
C                           NSTAT(4) OR MORE MONOTONIC SEQUENCES
C                           IF THE NULL HYPOTHESIS IS TRUE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES N IS LESS THAN 3
C                         WARNING (WITH FIX) ERROR
C                           IER=66 INDICATES NSTAT(1) IS LESS THAN 1.
C                             P(1) WILL BE SET TO MACHINE INFINITY.
C                         WARNING ERROR
C                           IER=35 INDICATES NSTAT(1) OR THE TOTAL
C                             NUMBER OF SEQUENCES USED FOR DETERMINING
C                             NSTAT(3) AND NSTAT(4) IS LESS THAN 8.
C                           IER=36 INDICATES AT LEAST 1 TIE WAS DETECTED
C
C   REQD. IMSL ROUTINES - MDBIN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF N IS GREATER THAN OR EQUAL TO 3 BUT NSTAT(1) IS
C                LESS THAN 1, P(1) WILL BE SET TO MACHINE INFINITY AND
C                IMSL ROUTINE MDBIN WILL NOT BE INVOKED. HOWEVER, THE
C                REMAINING STATISTICS AND ASSOCIATED PROBABILITIES WILL
C                BE DETERMINED AND RETURNED AS DESCRIBED.
C            2.  THE USER IS WARNED WHEN NSTAT(1) OR THE TOTAL NUMBER
C                OF SEQUENCES USED FOR DETERMINING NSTAT(3) AND NSTAT(4)
C                IS LESS THAN 8, AND WHEN AT LEAST ONE TIE IS DETECTED.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NBCYC  (X,N,EPS,NSTAT,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSTAT(1),IER
      REAL               X(1),EPS,P(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IER1,IER3,J,K,NMON,NTIE,NTIE1,NTOT
      REAL               RINFP,DIF1,DIF2,PR
      DATA               RINFP/Z7FFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IER1 = 0
      IER3 = 0
      IF (N .GE. 3) GO TO 5
      IER = 129
      GO TO 9000
C                                  INITIALIZE VARIABLES
    5 PR = 0.3333333
      NTOT = 0
      NMON = 0
      NTIE = 0
      I = 1
      J = 2
      K = 3
C                                  CHECK MID-VALUE FOR TIE
   10 DIF1 = ABS(X(I) - X(J)) - EPS
      DIF2 = ABS(X(K) - X(J)) - EPS
      IF ((DIF1 .LE. 0.0) .OR. (DIF2 .LE. 0.0)) GO TO 25
C                                  CHECK FOR MONOTONIC SEQUENCE
      IF ((X(I) .GT. X(J)) .AND. (X(J) .GT. X(K))) GO TO 30
      IF ((X(I) .LT. X(J)) .AND. (X(J) .LT. X(K))) GO TO 30
C                                  INCREMENT TOTAL SET NO. AND INDICES
   15 NTOT = NTOT + 1
      I = K+1
      J = I+1
   20 K = J+1
      IF (K .GT. N) GO TO 35
      GO TO 10
C                                  TIE DETECTED
   25 NTIE = NTIE + 1
      J = K
      GO TO 20
C                                  MONOTONIC SEQUENCE DETECTED
   30 NMON = NMON + 1
      GO TO 15
C                                  ALL SEQUENCES CHECKED VIA METHOD 1
   35 NSTAT(1) = NTOT
      NSTAT(2) = NMON
      IF (NTIE .GT. 0) GO TO 40
C                                  NO TIES, THUS METHOD 2 WOULD HAVE
C                                  SAME RESULTS AS METHOD 1
      NSTAT(3) = NMON
      NSTAT(4) = NMON
      NSTAT(5) = 0
      NSTAT(6) = 0
      IF (NTOT .LT. 8) IER = 35
      CALL MDBIN(NMON,NTOT,PR,DIF1,DIF2,I)
      P(1) = 1.0 - DIF1 + DIF2
      P(2) = P(1)
      P(3) = P(1)
      IF (IER .EQ. 0) GO TO 9005
      GO TO 9000
C                                  METHOD 2 - USED WHEN TIES DETECTED
   40 IER3 = 36
      IF (NTOT .LT. 1) IER1=66
      NTOT = 0
      NMON = 0
      NTIE1 = 0
      NTOT = N/3
      J = 3 * NTOT
      DO 55 I=1,J,3
         DIF1 = ABS(X(I) - X(I+1)) - EPS
         DIF2 = ABS(X(I+2) - X(I+1)) - EPS
         IF ((DIF1 .LE. 0.0) .OR. (DIF2 .LE. 0.0)) GO TO 45
         IF ((X(I) .GT. X(I+1)) .AND. (X(I+1) .GT. X(I+2))) GO TO 50
         IF ((X(I) .LT. X(I+1)) .AND. (X(I+1) .LT. X(I+2))) GO TO 50
         GO TO 55
   45    NTIE1 = NTIE1 + 1
         GO TO 55
   50    NMON = NMON + 1
   55 CONTINUE
      NSTAT(3) = NTOT
      NSTAT(4) = NMON
      NSTAT(5) = NTOT
      NSTAT(6) = NMON + NTIE1
      J = 1
      IF (IER1 .EQ. 0) GO TO 60
C                                  NSTAT(1) LESS THAN 1
      P(1) = RINFP
      J = 2
C                                  LOOP FOR OBTAINING PROBABILITIES
   60 DO 65 I=J,3
         CALL MDBIN(NSTAT(I+I),NSTAT(I+I-1),PR,DIF1,DIF2,K)
         P(I) = 1.0 - DIF1 + DIF2
   65 CONTINUE
      NSTAT(3) = NMON
      NSTAT(4) = NSTAT(6)
      NSTAT(5) = NTIE
      NSTAT(6) = NTIE1
      IF (NTOT .LT. 8) IER = 35
 9000 CONTINUE
      IF (IER1 .NE. 0) CALL UERTST (IER1,6HNBCYC )
      IF (IER .NE. 0)  CALL UERTST (IER,6HNBCYC )
      IF (IER3 .NE. 0) CALL UERTST (IER3,6HNBCYC )
      IER = MAX0(IER1,IER,IER3)
 9005 RETURN
      END

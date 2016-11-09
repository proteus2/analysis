C   IMSL ROUTINE NAME   - NBSIGN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SIGN TEST (FOR PERCENTILES)
C
C   USAGE               - CALL NBSIGN (X,N,Q,P,NSIGN,PROB,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SAMPLE OF INDEPENDENT OBSERVATIONS.
C                N      - INPUT SAMPLE SIZE. N MUST BE GREATER THAN OR
C                           EQUAL TO 2.
C                Q      - INPUT HYPOTHESIZED PERCENTILE BOUNDARY OF
C                           THE POPULATION FROM WHICH X WAS DRAWN.
C                P      - INPUT. Q IS THE 100*P PERCENTILE OF THE
C                           POPULATION. P SHOULD BE IN THE RANGE (0,1).
C                NSIGN  - OUTPUT VECTOR OF LENGTH 3.
C                         NSIGN(1) CONTAINS THE NUMBER OF NEGATIVE
C                           DIFFERENCES IN X(I)-Q, FOR I=1,2,...,N.
C                         NSIGN(2) CONTAINS THE NUMBER OF POSITIVE
C                           DIFFERENCES IN X(I)-Q, FOR I=1,2,...,N.
C                         NSIGN(3) CONTAINS THE NUMBER OF ZERO
C                           DIFFERENCES (TIES) IN X(I)-Q, FOR
C                           I=1,2,...,N.
C                PROB   - OUTPUT VECTOR OF LENGTH 4.
C                         PROB(1) CONTAINS THE PROBABILITY OF
C                           NSIGN(1)+NSIGN(3) NEGATIVE DIFFERENCES OR
C                           MORE (CONSIDERING TIES NEGATIVE).
C                         PROB(2) CONTAINS THE PROBABILITY OF NSIGN(2)
C                           POSITIVE DIFFERENCES OR MORE (CONSIDERING
C                           TIES NEGATIVE).
C                         PROB(3) CONTAINS THE PROBABILITY OF NSIGN(1)
C                           NEGATIVE DIFFERENCES OR MORE (CONSIDERING
C                           TIES POSITIVE).
C                         PROB(4) CONTAINS THE PROBABILITY OF
C                           NSIGN(2)+NSIGN(3) POSITIVE DIFFERENCES OR
C                           MORE (CONSIDERING TIES POSITIVE).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT THE SIZE OF THE
C                             SAMPLE, N, WAS LESS THAN 2.
C                           IER=130 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL ROUTINE MDBIN.
C                         WARNING ERROR
C                           IER=35 INDICATES AT LEAST ONE TIE WAS
C                             DETECTED.
C
C   REQD. IMSL ROUTINES - MDBIN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE USER IS WARNED IF ONE OR MORE TIES ARE DETECTED.
C            2.  IN JUDGING THE HYPOTHESES, ONE SHOULD BE CONCERNED
C                WITH THE NUMBER OF TIES (NBSIGN(3)), AS WELL AS THE
C                REJECTION PROBABILITIES PROVIDED IN PROB.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NBSIGN (X,N,Q,P,NSIGN,PROB,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSIGN(1),IER
      REAL               X(1),Q,P,PROB(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IER2,NEG,NNEG,NNPOS,NPOS,NTY
      REAL               PCOMP,PEQK,PLEK
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .GE. 2) GO TO 5
C                                  SAMPLE SIZE TOO SMALL
      IER = 129
      GO TO 9000
    5 IER = 0
      NEG = 0
      NTY = 0
C                                  ACCUMULATE NEG DIFFERENCES AND TIES
      DO 20  I=1,N
         IF(X(I)-Q)10,15,20
   10    NEG = NEG + 1
         GO TO 20
   15    NTY = NTY + 1
   20 CONTINUE
      NPOS = N - NEG - NTY
      NSIGN(1) = NEG
      NSIGN(2) = NPOS
      NSIGN(3) = NTY
C                                  COUNT TIES AS NEGATIVE DIFFERENCES
      NNEG = N - NPOS
      NNPOS = NPOS
C                                  LOOP TO COMPUTE 2 PAIRS OF PROBS.
      DO 45 I=1,3,2
C                                  FIND MOST EFFICIENT BINOM. ARGUMENTS
         IF(NNEG .GT. NNPOS) GO TO 30
C                                  USE BINOM. DISTRIBUTION OF NEGATIVES
         CALL MDBIN(NNEG,N,P,PLEK,PEQK,IER2)
         IF (IER2. LT. 128) GO TO 25
C                                  ERROR IN MDBIN
         IER = 130
         GO TO 9000
C                                  PROB OF THIS NUM. OF +DIFFS OR MORE
   25    PROB(I+1) = PLEK
C                                  PROB OF THIS NUM. OF -DIFFS OR MORE
         PROB(I) = 1.0 - PLEK + PEQK
         GO TO 40
C                                  USE BINOM. DISTRIBUTION OF POSITIVES
   30    PCOMP = 1.0 - P
         CALL MDBIN(NNPOS,N,PCOMP,PLEK,PEQK,IER2)
         IF (IER2 .LT. 128) GO TO 35
C                                  ERROR IN MDBIN
         IER = 130
         GO TO 9000
C                                  PROB OF THIS NUM. OF +DIFFS OR MORE
   35    PROB(I+1) = 1.0 - PLEK + PEQK
C                                  PROB OF THIS NUM. OF -DIFFS OR MORE
         PROB(I) = PLEK
C                                  BYPASS SWITCHES FOR I=3
C                                  DO NOT RECYCLE WHEN NUM. OF TIES = 0
   40    IF ((I .EQ. 3) .OR. (NTY .EQ. 0)) GO TO 50
C                                  SWITCH TIE COUNT TO POSITIVE TALLY
         NNEG = NEG
         NNPOS = N - NEG
   45 CONTINUE
   50 IF (NTY .GT. 0) GO TO 55
C                                  PROB(3)AND(4) WHEN NUM. OF TIES = 0
      PROB(3) = PROB(1)
      PROB(4) = PROB(2)
      GO TO 9005
   55 IER = 35
 9000 CONTINUE
      CALL UERTST(IER,6HNBSIGN)
 9005 RETURN
      END

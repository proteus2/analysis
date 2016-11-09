C   IMSL ROUTINE NAME   - NKS1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - KOLMOGOROV-SMIRNOV ONE-SAMPLE TEST
C
C   USAGE               - CALL NKS1 (PDF,X,N,PDIF,IER)
C
C   ARGUMENTS    PDF    - USER SUPPLIED SUBROUTINE WHICH MUST DETERMINE
C                           THE ORDINATE OF THE THEORETICAL DISTRIBUTION
C                           FUNCTION AT AN OBSERVATION (THE AREA OF THE
C                           THEORETICAL DENSITY FROM MINUS INFINITY TO
C                           THE OBSERVATION). PDF MUST BE DECLARED
C                           EXTERNAL IN THE CALLING PROGRAM AND MUST
C                           HAVE TWO ARGUMENTS, Y AND F, WHERE Y IS THE
C                           INPUT VALUE AND F IS THE OUTPUT ORDINATE.
C                           Y MUST NOT BE ALTERED IN THE SUBROUTINE.
C                X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           INDEPENDENT OBSERVATIONS (THE SAMPLE).
C                           PRIOR TO CALLING THE ROUTINE,
C                           X MUST BE SORTED INTO ASCENDING ORDER.
C                N      - INPUT NUMBER OF OBSERVATIONS (SAMPLE SIZE).
C                PDIF   - OUTPUT VECTOR OF LENGTH 6.
C                         PDIF(1) CONTAINS D, THE MAXIMUM OF D+,D-
C                           (SEE DESCRIPTION BELOW)
C                         PDIF(2) CONTAINS D+, THE SUPREMUM OF THE
C                           DIFFERENCES BETWEEN THE SAMPLE DISTRIBUTION
C                           FUNCTION AND THE THEORETICAL DISTRIBUTION
C                           FUNCTION.
C                         PDIF(3) CONTAINS D-, THE SUPREMUM OF THE
C                           DIFFERENCES BETWEEN THE THEORETICAL
C                           DISTRIBUTION FUNCTION AND THE SAMPLE
C                           DISTRIBUTION FUNCTION.
C                         PDIF(4) CONTAINS Z, THE STATISTIC USED TO
C                           OBTAIN THE PROBABILITIES.
C                         PDIF(5) CONTAINS THE PROBABILITY OF THE
C                           STATISTIC EXCEEDING Z IF THE HYPOTHESIS OF
C                           EQUALITY IS TRUE AND THE ALTERNATIVE
C                           IS ONE SIDED.
C                         PDIF(6) CONTAINS THE PROBABILITY OF EXCEEDING
C                           Z WHEN THE HYPOTHESIS OF EQUALITY IS TRUE
C                           AND THE ALTERNATIVE IS TWO SIDED.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 IMPLIES THAT THE SAMPLE SIZE IS LESS
C                             THAN 80.
C   REQD. IMSL ROUTINES - MDSMR,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE USER MAY WISH TO USE ROUTINES IN CHAPTER V OF THE
C                IMSL LIBRARY TO SORT THE OBSERVATIONS IN ASCENDING
C                ORDER.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NKS1   (PDF,X,N,PDIF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      REAL               X(N),PDIF(6)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,J,N1
      REAL               E,EE,F,FN,S,SS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF(N.GE.80) GO TO 5
C                                  SAMPLE SIZE WARNING
      IER = 33
    5 PDIF(2) = 0.
      PDIF(3) = 0.
      FN = N
      FN = 1./FN
      CALL PDF(X(1),F)
      S = 0.
      N1 = N - 1
      II = 1
C                                  CHECK FOR TIES
   10 DO 15 I=II,N1
         J = I
         IF(X(J).LT.X(J+1))GO TO 25
   15 CONTINUE
   20 J = N
   25 II = J + 1
C                                  RETAIN SAMPLE PDF VALUES
C                                  EVALUATE SAMPLE PDF
      SS = J*FN
C                                  EVALUATE THEORETICAL PDF
      E = F-S
      EE = SS-F
      IF (J.EQ.N) GO TO 30
      CALL PDF(X(J+1),F)
      S=SS
C                                  RETAIN MAX, MIN VALUES
   30 PDIF(2) = AMAX1(PDIF(2),EE)
      PDIF(3) = AMAX1(PDIF(3),E)
C                                  CHECK FOR LAST SAMPLE POINT
      IF (II-N)10,20,35
   35 PDIF(1) = AMAX1(PDIF(2),PDIF(3))
      PDIF(4) = PDIF(1)/SQRT(FN)
      CALL MDSMR(PDIF(4),E,S)
      PDIF(5) = 1. - E
      PDIF(6) = 1. - S
      IF(IER.EQ.0)GO TO 9005
C                                  CHECK ERROR INDICATOR
 9000 CONTINUE
      CALL UERTST(IER,6HNKS1  )
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - NBQT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COCHRAN Q TEST
C
C   USAGE               - CALL NBQT (X,N,M,IA,Q,PQ,IER)
C
C   ARGUMENTS    X      - INPUT N BY M MATRIX OF DICHOTOMIZED DATA,
C                           CONTAINING N READINGS OF ZERO OR ONE ON M
C                           VARIABLES (SAMPLES).
C                N      - INPUT NUMBER OF OBSERVATION SETS.
C                M      - INPUT NUMBER OF VARIABLES (SAMPLES). M IS THE
C                           NUMBER OF READINGS OF ZERO OR ONE IN EACH
C                           OBSERVATION SET. M MUST BE GREATER THAN
C                           OR EQUAL TO 3.
C                IA     - INPUT ROW DIMENSION OF THE MATRIX X EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                Q      - OUTPUT Q STATISTIC
C                PQ     - OUTPUT PROBABILITY OF EXCEEDING Q IF THE
C                           HYPOTHESIS OF EQUALITY OF UNDERLYING
C                           POPULATIONS IS TRUE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL ROUTINE MDCH
C                           IER=130 INDICATES THAT THE NUMBER OF
C                             SAMPLES, M, IS LESS THAN 3.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,
C                           UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE INPUT DATA MUST CONSIST OF ZEROS AND ONES ONLY.
C                FOR EXAMPLE, THE DATA MAY BE PASS-FAIL INFORMATION
C                ON M QUESTIONS ASKED OF N PEOPLE (1 INDICATING PASS,
C                0 INDICATING FAIL), OR TEST RESPONSES OF N INDIVI-
C                DUALS TO M DIFFERENT CONDITIONS. THE RESULTING
C                STATISTIC IS DISTRIBUTED APPROXIMATELY AS CHI-SQUARED
C                WITH M-1 DEGREES OF FREEDOM IF N IS NOT TOO SMALL.
C                N GREATER THAN OR EQUAL TO 5*M IS A CONSERVATIVE
C                RECOMMENDATION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NBQT   (X,N,M,IA,Q,PQ,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IA,IER
      REAL               X(IA,1),Q,PQ
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J
      REAL               FM1,QP,SUMC,SUMC1,S1,S2,S3
C                                  FIRST EXECUTABLE STATEMENT
      IER = 130
      IF(M.LT.3) GO TO 9000
      IER = 0
C                                  COMPUTE COLUMN TOTALS
      SUMC = 0.0
      DO 10 J=1,M
         SUMC1 = 0.0
         DO 5 I=1,N
    5    SUMC1 = SUMC1+X(I,J)
   10 SUMC = SUMC + SUMC1*SUMC1
C                                  COMPUTE ROW TOTALS
      S1 = 0.0
      S2 = 0.0
      DO 20 I=1,N
         S3 = 0.0
         DO 15 J=1,M
   15    S3 = S3 + X(I,J)
         S2 = S2 + S3
   20 S1 = S1+S3*S3
C                                  COMPUTE COCHRAN Q VALUE + P
      FM1 = M-1.
      Q = FM1*(M*SUMC-S2*S2)/(M*S2-S1)
      CALL MDCH(Q,FM1,QP,IER)
      IF (IER.GT.127) GO TO 9000
      PQ = 1.-QP
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HNBQT  )
 9005 RETURN
      END

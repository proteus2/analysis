C   IMSL ROUTINE NAME   - GTPST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - PAIRS TEST OR GOODS SERIAL TEST
C
C   USAGE               - CALL GTPST (A,IA,K,CS,Q,STD,IER)
C
C   ARGUMENTS    A      - INPUT K BY K TALLY MATRIX. A MAY BE OBTAINED
C                           AS OUTPUT FROM IMSL ROUTINE GTPR.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                K      - INPUT NUMBER OF SUBDIVISIONS OF THE UNIT LINE
C                           USED FOR TALLYING.
C                CS     - OUTPUT CHI-SQUARE STATISTIC
C                Q      - OUTPUT PROBABILITY OF EXCEEDING CHI-SQUARE
C                           IF THE NULL HYPOTHESIS OF UNIFORMITY IS
C                           TRUE
C                STD    - OUTPUT STANDARDIZED CHI-SQUARE STATISTIC
C                           STD=(CS-F)/SQRT(2*F) WHERE F=K*K-1 = THE
C                           DEGREES OF FREEDOM FOR THE CHI-SQUARE TEST.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT AN ERROR OCCURRED
C                             IN IMSL SUBROUTINE MDCH
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MDNRIS,MERRC=ERFC,
C                           MGAMAD=DGAMMA,UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MDNRIS,MERRC=ERFC,
C                           MGAMA=GAMMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  NORMALLY, TO USE THE PAIRS TEST WITH AN N LENGTH
C                SEQUENCE R IN MIND, ONE WOULD CALL GTCN TO DETERMINE
C                THE SQUARE OF THE APPROXIMATE K. THEN CALL GTPR TO
C                TALLY THE PAIRS (R(I),R(I+1)) WHERE I=1,3,...,N-1
C                AND CALL GTPST TO OBTAIN THE STATISTICS FOR CONSIDER-
C                ING THE HYPOTHESIS.
C            2.  IF THE USER IS FORCED BY CORE LIMITATION TO ENTER
C                GTPR MULTIPLE TIMES, THEN THE USER SHOULD NOTE THAT
C                THE LAST PAIR TALLIED IS ALWAYS (R(N-L),R(N)). THUS
C                EACH TIME GTPR IS ENTERED L-1 ELEMENTS OF THE SEQUENCE
C                ARE LEFT OUT OF THE TALLY. THIS, HOWEVER, IS OF LITTLE
C                CONCERN TO THE USER BECAUSE GTPST ACTUALLY SUMS THE
C                TALLY MATRIX TO OBTAIN THE CORRECT SAMPLE SIZE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTPST  (A,IA,K,CS,Q,STD,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,K,IER
      REAL               A(IA,1),CS,Q,STD
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J
      REAL               DF,S,SS,DQ
C                                  FIRST EXECUTABLE STATEMENT
      DF=K*K
C                                  SUM TALLY MATRIX
      S=0.0
         DO 5 I=1,K
            DO 5 J=1,K
    5       S=S+A(I,J)
      S=S/DF
C                                  CALCULATE CHI-SQUARE STATISTIC
      CS=0.0
         DO 10 I=1,K
            DO 10 J=1,K
            SS=A(I,J)-S
   10       CS=CS+SS*SS
      CS=CS/S
C                                  CALCULATE STANDARDIZED CHI-SQUARE
C                                  STATISTIC
      DF=DF-1
      STD=(CS-DF)/SQRT(DF+DF)
C                                  CALCULATE PROBABILITY
      CALL MDCH (CS,DF,DQ,IER)
      Q=1.0-DQ
      IF (IER .LT. 128) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GTPST ')
 9005 RETURN
      END

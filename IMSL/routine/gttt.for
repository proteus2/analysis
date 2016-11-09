C   IMSL ROUTINE NAME   - GTTT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TRIPLETS TEST
C
C   USAGE               - CALL GTTT (A,IA1,IA2,K,CS,Q,STD,IER)
C
C   ARGUMENTS    A      - INPUT THREE-DIMENSIONAL TALLY MATRIX. A IS
C                           DIMENSIONED K BY K BY K. A MAY BE OBTAINED
C                           FROM IMSL ROUTINE GTTRT.
C                IA1    - INPUT FIRST DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IA2    - INPUT SECOND DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                K      - INPUT NUMBER OF SUBDIVISIONS OF THE UNIT LINE
C                           USED FOR TALLYING. K SHOULD BE GREATER THAN
C                           OR EQUAL TO 2 AND LESS THAN OR EQUAL TO 58
C                CS     - OUTPUT CHI-SQUARE STATISTIC
C                Q      - OUTPUT PROBABILITY OF EXCEEDING CS IF THE
C                           NULL HYPOTHESIS OF UNIFORMITY IS TRUE
C                STD    - OUTPUT STANDARDIZED CHI-SQUARE STATISTIC
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES AN ERROR OCCURRED IN
C                             IMSL ROUTINE MDCH.
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
C   REMARKS      NORMALLY, TO PERFORM THE TRIPLETS TEST FOR A SEQUENCE
C                R OF LENGTH N, ONE WOULD CALL GTCN TO DETERMINE THE
C                CUBE OF AN APPROXIMATE K, CALL GTTRT TO TALLY THE
C                TRIPLETS IN A AND THEN CALL GTTT TO OBTAIN THE
C                STATISTICS FOR EXAMINING THE HYPOTHESIS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTTT  (A,IA1,IA2,K,CS,Q,STD,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA1,IA2,K,IER
      REAL               A(IA1,IA2,1),CS,Q,STD
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,L
      REAL               DF,S,SS,DQ
C                                  FIRST EXECUTABLE STATEMENT
      DF=K*K*K
C                                  SUM THE TALLY MATRIX
      S=0.0
         DO 5 I=1,K
            DO 5 J=1,K
               DO 5 L=1,K
    5          S=S+A(I,J,L)
      S=S/DF
C                                  CALCULATE CHI-SQUARE STATISTIC
      CS=0.0
         DO 10 I=1,K
            DO 10 J=1,K
               DO 10 L=1,K
               SS=A(I,J,L)-S
   10          CS=CS+SS*SS
      CS=CS/S
C                                  CALCULATE STANDARDIZED CHI-SQUARE
C                                  STATISTIC
      DF=DF-1.
      STD=(CS-DF)/SQRT(DF+DF)
C                                  CALCULATE PROBABILITY
      CALL MDCH (CS,DF,DQ,IER)
      Q=1.0-DQ
      IF (IER .LT. 128) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GTTT  ')
 9005 RETURN
      END

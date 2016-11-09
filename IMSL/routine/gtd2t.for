C   IMSL ROUTINE NAME   - GTD2T
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - THE D-SQUARE TEST
C
C   USAGE               - CALL GTD2T (COUNT,K,E,CS,STD,Q,IER)
C
C   ARGUMENTS    COUNT  - INPUT VECTOR WHICH CONTAINS THE TALLY OF
C                           DISTANCES (BETWEEN PAIRS OF POINTS IN THE
C                           UNIT SQUARE) WHICH OCCURRED IN RESPECTIVE
C                           K CATEGORIES. COUNT IS USUALLY THE OUTPUT
C                           VECTOR FROM IMSL ROUTINE GTDDU AND IS OF
C                           LENGTH K.
C                K      - INPUT NUMBER OF CELLS INTO WHICH THE SQUARED
C                           DISTANCE INTERVAL (0,2) HAS BEEN DIVIDED.
C                           DEGREES OF FREEDOM FOR CHI-SQUARE TEST ARE
C                           K-1. K MUST BE GREATER THAN ONE.
C                E      - OUTPUT EXPECTED VALUE OF EACH CELL OF COUNT,
C                           GIVEN THE NULL HYPOTHESIS OF UNIFORMITY IS
C                           TRUE.
C                CS     - OUTPUT CHI-SQUARE STATISTIC
C                STD    - OUTPUT CHI-SQUARE STATISTIC, STANDARDIZED
C                Q      - OUTPUT PROBABILITY OF OBTAINING A STATISTIC
C                           AT LEAST AS GREAT AS CS, GIVEN THE NULL
C                           HYPOTHESIS OF UNIFORMITY IS TRUE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES K IS LESS THAN 2
C                             (INSUFFICIENT NUMBER OF SUBDIVISIONS OF
C                             DISTANCE)
C                         WARNING ERROR WITH FIX
C                           IER = 66 INDICATES THAT E IS LESS THAN 1.
C                             Q WILL BE SET TO MACHINE INFINITY AND
C                             IMSL ROUTINE MDCH WILL NOT BE INVOKED.
C                         WARNING ERROR
C                           IER = 35 INDICATES E LESS THAN 5
C                             (SMALLER THAN ADVISABLE)
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
C-----------------------------------------------------------------------
C
      SUBROUTINE GTD2T  (COUNT,K,E,CS,STD,Q,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,IER
      REAL               COUNT(1),E,CS,STD,Q
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IER2
      REAL               RINFP,SUM,DF,QUE
C                                  RINFP=LGST POS REAL -SINGLE PRECISION
      DATA               RINFP/Z7FFFFFFF/
C                                  CHECK K=NUMBER OF CELLS
C                                  FIRST EXECUTABLE STATEMENT
      IF (K .GE. 2) GO TO 5
      IER = 129
      GO TO 9000
    5 IER = 0
C                                  ACCUMULATE TALLY TOTAL
      SUM = 0.0
      DO 10 I = 1,K
         SUM = SUM + COUNT(I)
   10 CONTINUE
C                                  COMPUTE EXPECTED VALUE OF COUNT
      E = SUM/K
C                                  CHECK FOR E LESS THAN 5
      IF (E .LT. 5.0) IER=35
C                                  COMPUTE CHI-SQ STATISTIC
      SUM = 0.0
      DO 15 I=1,K
         SUM = SUM + (COUNT(I)-E)**2
   15 CONTINUE
      CS = SUM/E
C                                  COMPUTE STANDARDIZED CHI-SQ
      DF = K-1
      STD = (CS-DF)/SQRT(DF+DF)
C                                  BYPASS CALL FOR E LESS THAN 1
      IF (E .LT. 1.0) GO TO 20
      CALL MDCH(CS,DF,QUE,IER2)
      Q = 1.-QUE
      IF (IER .EQ. 0) GO TO 9005
      GO TO 9000
   20 IER = 66
      Q = RINFP
 9000 CONTINUE
      CALL UERTST(IER,'GTD2T ')
 9005 RETURN
      END

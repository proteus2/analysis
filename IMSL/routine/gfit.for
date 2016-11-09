C   IMSL ROUTINE NAME   - GFIT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CHI-SQUARED GOODNESS OF FIT TEST
C
C   USAGE               - CALL GFIT(CDF,K,OBS,N,CELLS,COMP,CS,IDF,Q,IER)
C
C   ARGUMENTS    CDF    - THEORETICAL CUMLATIVE DISTRIBUTION FUNCTION
C                           AGAINST WHICH SAMPLE IS TESTED.  CDF IS A
C                           USER-SUPPLIED EXTERNAL SUBROUTINE WHICH
C                           PROVIDES THE ORDINATE OF THE HYPOTHESIZED
C                           CUMULATIVE DISTRIBUTION FUNCTION AT POINTS
C                           X. GFIT CALLS CDF AS FOLLOWS,
C                                    CALL CDF (X,P)
C                           P IS THE RESULTANT ORDINATE AT X. P MUST NOT
C                           BE OUTSIDE THE RANGE OF THE CLOSED INTERVAL
C                           (0,1). CDF MUST APPEAR IN AN EXTERNAL STATE-
C                           MENT IN THE CALLING PROGRAM.
C                K      - INPUT NUMBER OF EQUIPROBABLE CELLS INTO WHICH
C                           THE OBSERVATIONS ARE TO BE TALLIED.
C                OBS    - INPUT VECTOR OF LENGTH N, CONTAINING THE
C                           SAMPLE.
C                N      - INPUT LENGTH OF THE VECTOR OBS.
C                CELLS  - OUTPUT VECTOR OF LENGTH K CONTAINING COUNTS
C                           OF OBSERVATIONS WHICH FALL INTO THE
C                           EQUIPROBABLE CATEGORIES.
C                COMP   - OUTPUT VECTOR OF LENGTH K CONTAINING THE
C                           COMPONENTS OF THE CHI-SQUARED STATISTIC,
C                           K*(CELLS-N/K)**2/N.
C                CS     - OUTPUT CHI-SQUARED STATISTIC.
C                IDF    - INPUT IDF CONTAINS THE NUMBER OF HYPOTHESIZED
C                           DISTRIBUTION FUNCTION PARAMETERS ESTIMATED
C                           FROM THE DATA, OBS.  IDF SHOULD BE ZERO
C                           IF NONE WERE ESTIMATED.  ON OUTPUT, IDF
C                           CONTAINS THE DEGREES OF FREEDOM OF THE
C                           STATISTIC CS.
C                Q      - OUTPUT PROBABILITY OF THE CHI-SQUARED
C                           STATISTIC EXCEEDING CS IF THE NULL
C                           HYPOTHESIS IS TRUE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 MEANS AN ERROR OCCURRED IN MDCH.
C                           IER = 130 MEANS K WAS LESS THAN 2.
C                         WARNING ERROR
C                           IER = 35 MEANS AN EXPECTED CELL VALUE IS
C                             LESS THAN 5.
C                           IER = 36 MEANS AN EXPECTED CELL VALUE IS
C                             LESS THAN 1. IN THIS CASE, Q IS SET TO
C                             NEGATIVE MACHINE INFINITY.
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
      SUBROUTINE GFIT (CDF,K,OBS,N,CELLS,COMP,CS,IDF,Q,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,N,IDF,IER
      REAL               OBS(1),CELLS(1),COMP(1),CS,Q
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,IET
      REAL               X,P,EN,EK,CELL1,E,CS1,DF,Q1
      DATA               RINFM/ZFFFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF(K.LT.2) GO TO 30
C                                  INITIALIZE COUNTS
      DO 5 I=1,K
    5 CELLS(I) = 0.0
C                                  GENERATE COUNTS
      DO 10 I=1,N
         X = OBS(I)
         CALL CDF(X,P)
         J = K*P + 1
         IF(P.EQ.1.0) J = K
   10 CELLS(J) = CELLS(J) + 1.0
C                                  DETERMINE EXPECTED VALUES OF CELLS
      EN = N
      EK = K
      E = EN/EK
      DO 15 I=1,K
         CELL1 = CELLS(I)-E
C                                  FIND COMPONENTS OF STATISTICS
   15 COMP(I) = EK*(CELL1*CELL1)/EN
      CS1 = 0.0
      DO 20 I=1,K
   20 CS1 = CS1 + COMP(I)
      CS = CS1
      IDF = K-1-IDF
      IF(E.LT.5.0) IER = 35
      IF(E.LT.1.0) GO TO 25
C                                  FIND PROBABILITY
      DF=IDF
      CALL MDCH (CS1,DF,Q1,IET)
      IF(IET.GT.127) GO TO 35
      Q = 1.0 - Q1
      IF(IER - 32) 9005,9000,9000
   25 IER = 36
      Q = RINFM
      GO TO 9000
   30 IER = 130
      GO TO 9000
   35 IER = IET
 9000 CONTINUE
      CALL UERTST(IER,'GFIT  ')
 9005 RETURN
      END

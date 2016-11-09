C   IMSL ROUTINE NAME   - GTPOK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - THE POKER TEST
C
C   USAGE               - CALL GTPOK (HANDS,HE,K1,CSOBS,K2,STAT,IER)
C
C   ARGUMENTS    HANDS  - INPUT VECTOR OF LENGTH K1 CONTAINING COUNTS
C                           OF TYPES OF POKER HANDS. (OUTPUT VECTOR
C                           HSAVE FROM IMSL ROUTINE GTPL MAY BE USED
C                           AS INPUT VECTOR HANDS.)
C                HE     - INPUT/OUTPUT VECTOR OF LENGTH K1. ON INPUT,
C                           HE CONTAINS EXPECTED VALUES OF HANDS.(OUT-
C                           PUT VECTOR HE FROM IMSL ROUTINE GTPL MAY
C                           BE USED AS INPUT.) ON OUTPUT, HE CONTAINS
C                           EXPECTED VALUES OF THE CELLS IN THE VECTOR
C                           HANDS.
C                K1     - INPUT LENGTH OF VECTORS HANDS AND HE.
C                CSOBS  - INPUT TALLY VECTOR OF LENGTH K2 CONTAINING
C                           CHI-SQUARE STATISTICS. (OUTPUT VECTOR CSOBS
C                           FROM IMSL ROUTINE GTPL.)
C                K2     - INPUT LENGTH OF VECTOR CSOBS.
C                STAT   - OUTPUT VECTOR OF LENGTH 8 CONTAINING TEST
C                           STATISTICS, WHERE
C                           STAT(1) = CHI-SQUARE STATISTIC BASED ON
C                             THE TYPE OF POKER HANDS.
C                           STAT(2) = STANDARDIZED CHI-SQUARE STATISTIC
C                             BASED ON TYPE OF POKER HANDS.
C                           STAT(3) = PROBABILITY THAT THE CHI-SQUARE
C                             STATISTIC ABOVE WOULD BE EXCEEDED, IF THE
C                             UNIFORM HYPOTHESIS IS TRUE.
C                           STAT(4) = DEGREES OF FREEDOM FOR STAT(1).
C                           STAT(5) = CHI-SQUARE STATISTIC FROM
C                             GOODNESS OF FIT TEST ON CSOBS.
C                           STAT(6) = STANDARDIZED CHI-SQUARE STATISTIC
C                             FROM GOODNESS OF FIT TEST ON CSOBS.
C                           STAT(7) = PROBABILITY THAT THE CHI-SQUARE
C                             STATISTIC OF STAT(5) WOULD BE EXCEEDED IF
C                             THE UNIFORM HYPOTHESIS IS TRUE.
C                           STAT(8) = DEGREES OF FREEDOM FOR STAT(5).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES AN ERROR OCCURED IN MDCH
C                         WARNING ERROR
C                           IER = 34 INDICATES AN EXPECTED VALUE LESS
C                             THAN FIVE WAS CALCULATED INTO HE.
C                           IER = 35 INDICATES AN EXPECTED VALUE LESS
C                             THAN ONE WAS CALCULATED INTO HE. STAT(3)
C                             WAS SET TO MACHINE INFINITY.
C                           IER = 36 INDICATES THAT THE EXPECTED VALUE
C                             ON WHICH THE GOODNESS OF FIT TEST WAS
C                             BASED WAS LESS THAN 5.
C                           IER = 37 INDICATES THAT THE EXPECTED VALUE
C                             ON WHICH THE GOODNESS OF FIT TEST WAS
C                             BASED WAS LESS THAN 1.  STAT(7) WAS SET
C                             TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,
C                           UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,
C                           MGAMA=GAMMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTPOK(HANDS,HE,K1,CSOBS,K2,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K1,K2,IER
      REAL               HANDS(1),HE(1),CSOBS(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,IET
      REAL               RINFP,G1M1,S,HTEMP2,HTEMP1,HCS,HCD,Q,CSE
      DATA               RINFP/Z7FFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      G1M1 = K1-1
      S = 0.0
C                                  COUNT TOTAL NUMBER OF HANDS
      DO 5 I=1,K1
    5 S = S + HANDS(I)
C                                  FIND EXPECTED NUMBER OF HANDS
      DO 10 J=1,K1
         HE(J) = HE(J) * S
         IF(HE(J).LT.5.0) IER = 34
         IF(HE(J).LT.1.0) IER = 35
   10 CONTINUE
      HTEMP2 = 0.0
C                                  CALCULATED TOTAL HAND CHI-SQUARED
      DO 15 I=1,K1
         HTEMP1 = HANDS(I)-HE(I)
   15 HTEMP2 = HTEMP2+(HTEMP1*HTEMP1)/HE(I)
C                                  CHI-SQUARED STANDARDIZED
      HCS = (HTEMP2-G1M1)/SQRT(G1M1+G1M1)
      IF(IER.EQ.35) GO TO 20
      HCD = HTEMP2
C                                  OBTAIN PROBABILITY WITH CHI-SQUARED
C                                     STATISTIC
      CALL MDCH  (HCD,G1M1,Q,IET)
      IF(IET.GT.127) GO TO 50
      GO TO 25
   20 Q = RINFP
   25 STAT(1) = HTEMP2
      STAT(2) = HCS
      STAT(3) = 1.0-Q
      STAT(4) = G1M1
      HTEMP2 = 0.0
      G1M1 = K2-1
C                                  EXPECTED NUMBER FOR GOODNESS OF FIT
C                                     TEST
      DO 30 I=1,K2
   30 HTEMP2 = HTEMP2+CSOBS(I)
      CSE = HTEMP2/K2
      IF(CSE.LT.5.0) IER = 36
      IF(CSE.LT.1.0) IER = 37
      HTEMP2 = 0.0
C                                  CHI-SQUARED FOR GOODNESS OF FIT TEST
      DO 35 I=1,K2
         HTEMP1 = CSOBS(I)-CSE
   35 HTEMP2 = HTEMP2+HTEMP1*HTEMP1
      STAT(5) = HTEMP2/CSE
      IF(IER.EQ.37) GO TO 40
      HCD = STAT(5)
C                                  PROBABILITY FOR GOODNESS OF FIT TEST
      CALL MDCH  (HCD,G1M1,Q,IET)
      IF(IET.GT.127) GO TO 50
      GO TO 45
   40 Q = RINFP
C                                  STANDARDIZED CHI-SQUARED STATISTIC
   45 STAT(6) = (STAT(5)-G1M1)/SQRT(G1M1+G1M1)
      STAT(7) = 1.0-Q
      STAT(8) = G1M1
      IF(IER) 9000,9005,9000
   50 IER = IET
 9000 CONTINUE
      CALL UERTST(IER,'GTPOK ')
 9005 RETURN
      END

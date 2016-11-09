C   IMSL ROUTINE NAME   - BESTA2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - COMPUTATIONS OF CONFIDENCE INTERVALS AND
C                           OTHER BASIC STATISTICS USING OUTPUT FROM
C                           IMSL ROUTINE BESTAT
C
C   USAGE               - CALL BESTA2 (NVAR,IOPT,CLM,CLV,STATS,IS,
C                           RNG,CV,CIM,CIV,IER)
C
C   ARGUMENTS    NVAR   - INPUT NUMBER OF VARIABLES.
C                IOPT   - INPUT OPTIONS VECTOR OF LENGTH 4, CONTAINING
C                           A 1 IN EACH POSITION CORRESPONDING TO A
C                           STATISTIC TO BE COMPUTED FOR EACH VARIABLE
C                           IOPT(1) = 1 INDICATES RANGES ARE TO BE
C                               COMPUTED.
C                           IOPT(2) = 1 INDICATES COEFFICIENTS  OF
C                               VARIATION ARE TO BE COMPUTED.
C                           IOPT(3) = 1 INDICATES CONFIDENCE INTERVALS
C                               FOR THE MEANS ARE TO BE COMPUTED.
C                           IOPT(4) = 1 INDICATES CONFIDENCE INTERVALS
C                               FOR THE VARIANCES ARE TO BE COMPUTED.
C                CLM    - INPUT CONFIDENCE LEVEL FOR THE INTERVAL
C                           ESTIMATE OF THE MEAN (ASSUMING NORMALITY).
C                           CLM IS A PERCENTAGE IN THE EXCLUSIVE RANGE
C                           0 TO 100.
C                CLV    - INPUT CONFIDENCE LEVEL FOR THE INTERVAL
C                           ESTIMATE OF THE VARIANCE (ASSUMING
C                           NORMALITY). CLV IS A PERCENTAGE IN THE
C                           EXCLUSIVE RANGE 0 TO 100.
C                STATS  - INPUT 8 BY NVAR MATRIX CONTAINING IN EACH
C                           ROW STATISTICS OF EACH OF THE VARIABLES.
C                             STATS(1,*) CONTAINS MEANS,
C                             STATS(2,*) CONTAINS VARIANCES,
C                             STATS(3,*) CONTAINS STANDARD DEVIATIONS,
C                             STATS(4,*) IS NOT USED BY BESTA2,
C                             STATS(5,*) IS NOT USED BY BESTA2,
C                             STATS(6,*) CONTAINS MINIMA,
C                             STATS(7,*) CONTAINS MAXIMA, AND
C                             STATS(8,*) CONTAINS NUMBERS (COUNTS).
C                           THESE STATISTICS ARE IN THE SAME ORDER AS
C                           PRODUCED BY THE IMSL ROUTINE BESTAT.
C                IS     - INPUT ROW DIMENSION OF STATS EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                RNG    - OUTPUT VECTOR OF LENGTH NVAR CONTAINING THE
C                           RANGES.
C                CV     - OUTPUT VECTOR OF LENGTH NVAR CONTAINING THE
C                           COEFFICIENTS OF VARIATION WHEN THEY ARE
C                           DEFINED. IF THE COEFFICIENT IS NOT DEFINED
C                           FOR A GIVEN VARIABLE, CV CONTAINS A
C                           ZERO IN THE CORRESPONDING POSITION.
C                CIM    - OUTPUT VECTOR OF LENGTH 2*NVAR CONTAINING
C                           THE CONFIDENCE LIMITS FOR THE MEANS
C                           (ASSUMING NORMALITY). FOR THE (I)TH
C                           VARIABLE, THE LOWER LIMIT IS IN THE
C                           (2*I-1)TH POSITION OF CIM AND THE UPPER
C                           LIMIT IS IN THE (2*I)TH  POSITION.
C                CIV    - OUTPUT VECTOR OF LENGTH 2*NVAR CONTAINING
C                           THE CONFIDENCE LIMITS FOR THE VARIANCES
C                           (ASSUMING NORMALITY). FOR THE (I)TH
C                           VARIABLE, THE LOWER LIMIT IS IN THE
C                           (2*I-1)TH POSITION OF CIM AND THE UPPER
C                           LIMIT IS IN THE (2*I)TH POSITION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=34 INDICATES THAT SOME VARIABLE(S)
C                             HAD FEWER THAN TWO OBSERVATIONS.
C                             THE PERTINENT STATISTICS ARE
C                             SET TO NEGATIVE MACHINE INFINITY.
C                           IER=35 INDICATES A NEGATIVE RANGE WAS
C                             COMPUTED BECAUSE OF INCORRECT INPUT IN
C                             STATS.THE RANGE IS SET TO NEGATIVE
C                             MACHINE INFINITY.
C                           IER=36 INDICATES THAT A MEAN WAS SO CLOSE
C                             TO ZERO THAT A COEFFICIENT OF VARIATION
C                             COULD NOT BE COMPUTED. THE CORRESPONDING
C                             COEFFICIENT OF VARIATION IS SET TO ZERO.
C                             THIS ERROR IS ALSO GENERATED IF A
C                             STANDARD DEVIATION IS NONPOSITIVE.
C                           IER=37 INDICATES THAT TWO OR MORE OF THE
C                             ABOVE ERROR SITUATIONS OCCURRED.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT CLM IS OUT OF RANGE.
C                             CIM IS NOT COMPUTED.
C                           IER=130 INDICATES THAT CLV IS OUT OF RANGE.
C                             CIV IS NOT COMPUTED.
C                           IER=131 INDICATES THAT AN ERROR OCCURRED IN
C                             COMPUTING A CONFIDENCE LIMIT FOR THE MEAN.
C                             THE CORRESPONDING LIMITS ARE SET TO
C                             NEGATIVE MACHINE INFINITY.
C                           IER=132 INDICATES THAT AN ERROR OCCURRED IN
C                             COMPUTING A CONFIDENCE LIMIT FOR THE
C                             VARIANCE. THE CORRESPONDING LIMITS ARE
C                             SET TO NEGATIVE MACHINE INFINITY.
C                           IER=133 INDICATES THAT ERRORS OCCURRED IN
C                             COMPUTING A CONFIDENCE LIMIT FOR THE MEAN
C                             AND THE VARIANCE. THE CORRESPONDING LIMITS
C                             ARE SET TO NEGATIVE MACHINE INFINITY.
C
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDCHI,MDNOR,MDNRIS,MDSTI,MERFI,
C                           MERRC=ERFC,MGAMAD=DGAMMA,UERSET,UERTST,
C                           UGETIO
C                       - H36,H48,H60/MDCH,MDCHI,MDNOR,MDNRIS,MDSTI,
C                           MERFI,MERRC=ERFC,MGAMA=GAMMA,UERSET,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BESTA2 (NVAR,IOPT,CLM,CLV,STATS,IS,RNG,CV,CIM,CIV,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NVAR,IS,IER,IOPT(1)
      REAL               CLM,CLV,STATS(IS,1),RNG(NVAR),CV(NVAR),CIM(1),
     *                   CIV(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JER,KER,LER,LEVOLD
      REAL               ABSM,HALFL,HUND,ONE,PT5,XETA,XINF,XINFM,XMEAN,
     *                   XSD,ZERO
      REAL               Q,Q2,DF,TVAL,X1,X2
      DATA               ZERO,ONE,HUND,PT5/0.0,1.0,100.,0.5/
      DATA               XINF/Z7FFFFFFF/
      DATA               XINFM/ZFFFFFFFF/
      DATA               XETA/Z00100000/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
      KER = 0
      CALL UERSET (0,LEVOLD)
      IF (IOPT(1).NE.1) GO TO 10
C                                  COMPUTE RANGES
      DO 5 I = 1,NVAR
         RNG(I) = STATS(7,I) - STATS(6,I)
         IF (RNG(I).GE.ZERO) GO TO 5
         RNG(I) = XINFM
         IER=35
    5 CONTINUE
   10 IF (IOPT(2).NE.1) GO TO 35
C                                  COMPUTE COEFFICIENTS OF VARIATION
      JER = 0
      DO 30 I = 1,NVAR
         XMEAN=STATS(1,I)
         ABSM=ABS(XMEAN)
         XSD=STATS(3,I)
         IF (XSD.LE.ZERO) GO TO 20
         IF(ABSM.GT.ONE) GO TO 15
         Q = XINF*ABSM
         IF(XSD.LT.Q) GO TO 25
         GO TO 20
   15    Q = XETA*ABSM
         IF(XSD.GT.Q) GO TO 25
   20    JER=36
         IF(IER.GT.0) JER = 37
         CV(I)=ZERO
         GOTO 30
   25    CV(I) = XSD / XMEAN
   30 CONTINUE
      IF (JER.NE.0) IER = JER
   35 IF (IOPT(3).NE.1) GO TO 60
C                                  COMPUTE CONFIDENCE INTERVALS FOR MEAN
      Q=ONE-CLM/HUND
      IF(Q.GT.ZERO.AND.Q.LT.ONE) GOTO 40
      IER=129
      GOTO 60
   40 JER = 0
      DO 55 I = 1,NVAR
         DF=STATS(8,I)-ONE
         IF(DF.LT.PT5) GOTO 45
         IF (STATS(3,I).LT.ZERO) GO TO 50
         CALL MDSTI (Q,DF,TVAL,KER)
         IF(KER.NE.0) GOTO 50
         HALFL = TVAL*STATS(3,I)/SQRT(STATS(8,I))
         CIM(2*I-1) = STATS(1,I) - HALFL
         CIM(2*I) = STATS(1,I) + HALFL
         GOTO 55
   45    JER=34
         IF(IER.GT.0) JER = 37
         CIM(2*I-1) = XINFM
         CIM(2*I) = XINFM
         GOTO 55
   50    CIM(2*I-1) = XINFM
         CIM(2*I) = XINFM
         IER = 131
   55 CONTINUE
      IF (JER.NE.0) IER = JER
   60 IF (IOPT(4).NE.1) GO TO 85
C                                  COMPUTE CONFIDENCE INTERVALS
C                                  FOR VARIANCES
      JER=0
      Q=PT5*(ONE-CLV/HUND)
      Q2=ONE-Q
      IF(Q.GT.ZERO.AND.Q.LT.PT5) GOTO 65
      IER=130
      GOTO 85
   65 DO 80 I = 1,NVAR
         DF=STATS(8,I)-ONE
         IF(DF.LT.PT5) GOTO 70
         IF(STATS(2,I).LT.ZERO) GO TO 75
         CALL MDCHI (Q,DF,X1,LER)
         CALL MDCHI (Q2,DF,X2,LER)
         IF(LER.NE.0) GO TO 75
         CIV(2*I-1) = STATS(2,I)*DF/X2
         CIV(2*I) = STATS(2,I)*DF/X1
         GOTO 80
   70    JER=34
         IF(IER.GT.34) JER = 37
         CIV(2*I-1) = XINFM
         CIV(2*I) = XINFM
         GOTO 80
   75    KER=132
         CIV(2*I-1)=XINFM
         CIV(2*I)=XINFM
   80 CONTINUE
      IF (JER.NE.0) IER = JER
      IF(IER.NE.131.AND.KER.NE.0)IER=132
      IF(IER.NE.132.AND.KER.NE.0)IER=133
   85 CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BESTA2')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - BENSON
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - VARIANCE INFERENCES USING A SAMPLE FROM A
C                           NORMAL POPULATION WITH KNOWN MEAN
C
C   USAGE               - CALL BENSON (Y,N,IOP,CRIT,M,VAR,STAT,NDF,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SAMPLE.
C                N      - SAMPLE SIZE. (INPUT)
C                IOP    - INPUT VECTOR OF LENGTH 2 INDICATING, FOR THE
C                           VARIANCE, HYPOTHESIS TEST AND INTERVAL
C                           ESTIMATE TO BE COMPUTED.
C                           IOP(1) = I IMPLIES THE FOLLOWING TEST WHEN
C                             I=1, TWO-TAILED.
C                             I=2, UPPER ONE-TAILED.
C                             I=3, LOWER ONE-TAILED.
C                           I NOT EQUAL TO 1,2, OR 3 IMPLIES TEST IS NOT
C                           DESIRED. SEE REMARKS.
C                           IOP(2) = I IMPLIES THE FOLLOWING INTERVAL
C                           ESTIMATE WHEN
C                             I=1, TWO-SIDED.
C                             I=2, UPPER ONE-SIDED.
C                             I=3, LOWER ONE-SIDED.
C                           I NOT EQUAL TO 1,2, OR 3 IMPLIES ESTIMATE IS
C                           NOT DESIRED.
C                CRIT   - INPUT VECTOR OF LENGTH 3 CONTAINING CONSTANTS
C                           REQUIRED FOR INTERVAL ESTIMATE AND HYPOTHE-
C                           SIS TEST. THE I-TH ELEMENT OF CRIT CONTAINS,
C                           WHEN
C                           I = 1, CONFIDENCE COEFFICIENT, BETWEEN
C                             0.0 AND 1.0 EXCLUSIVELY, FOR INTERVAL
C                             ESTIMATE OF THE VARIANCE (REQUIRED ONLY
C                             WHEN IOP(2) = 1, 2, OR 3).
C                           I = 2, HYPOTHESIZED VALUE OF VARIANCE
C                             (REQUIRED ONLY WHEN IOP(1) = 1,2, OR 3).
C                           I = 3, KNOWN MEAN OF POPULATION.
C                M      - NUMBER OF RESPONSES IN EARLIER SAMPLE
C                           YIELDING ADDITIONAL ESTIMATE OF VARIANCE. M
C                           LESS THAN 2 IMPLIES NO ADDITIONAL ESTIMATE
C                           IS AVAILABLE. (INPUT)
C                VAR    - IF M GREATER THAN OR EQUAL TO 2, VAR CONTAINS
C                           ADDITIONAL ESTIMATE OF VARIANCE ON INPUT
C                           AND COMBINED ESTIMATE ON OUTPUT. IF M IS
C                           LESS THAN 2, VAR IS UNDEFINED ON INPUT AND
C                           CONTAINS ESTIMATE OF VARIANCE ON OUTPUT.
C                STAT   - OUTPUT VECTOR OF LENGTH 3 CONTAINING RESULTS
C                           OF HYPOTHESIS TEST AND INTERVAL ESTIMATE
C                           FOR VARIANCE. THE I-TH ELEMENT OF STAT
C                           CONTAINS, WHEN
C                           I = 1, PROBABILITY OF COMPUTED OR MORE
C                             EXTREME VALUE OF TEST STATISTIC FOR
C                             HYPOTHESIS TEST ON VARIANCE. (DEFINED ONLY
C                             WHEN IOP(1) = 1,2, OR 3).
C                           I = 2, LOWER LIMIT FOR VARIANCE (DEFINED
C                             ONLY WHEN IOP(2) = 1 OR 3).
C                           I = 3, UPPER LIMIT FOR VARIANCE (DEFINED
C                             ONLY WHEN IOP(2) = 1 OR 2).
C                NDF    - NUMBER OF DEGREES OF FREEDOM IN VARIANCE
C                           ESTIMATE. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS NUMBER OF RESPONSES SPECIFIED
C                             IS LESS THAN 1.
C                           IER=130 MEANS AN ERROR OCCURRED IN MDCH.
C                           IER=131 MEANS AN ERROR OCCURRED IN MDCHI.
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/MDCH,MDCHI,MDNOR,MDNRIS,MERFI,
C                           MERRC=ERFC,MGAMAD=DGAMMA,UERTST,UGETIO
C                       - SINGLE(H36,H48,H60)/MDCH,MDCHI,MDNOR,MDNRIS,
C                           MERFI,MERRC=ERFC,MGAMA=GAMMA,UERTST,UGETIO
C                       - DOUBLE/MDCH,MDCHI,MDNOR,MDNRIS,MERFI,
C                           MERRC=ERFC,MGAMAD=DGAMMA,UERTST,UGETIO,
C                           VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      LETTING (M,S) BE THE TRUE MEAN AND KNOWN VARIANCE OF
C                THE POPULATION, RESPECTIVELY, AND DEFINING
C                           S0 = CRIT(2)
C                THE TABLE BELOW CLARIFIES HYPOTHESIS TESTING OPTIONS.
C
C                      IOP(1)            NULL            ALTERNATIVE
C                                     HYPOTHESIS         HYPOTHESIS
C                      ------         ----------         -----------
C                        1             S.EQ.S0             S.NE.S0
C                        2             S.LE.S0             S.GT.S0
C                        3             S.GE.S0             S.LT.S0
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BENSON (Y,N,IOP,CRIT,M,VAR,STAT,NDF,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,NDF,IER,IOP(2)
      REAL               Y(N),CRIT(3),VAR,STAT(3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IOP1,IOP2,I
      REAL               HALF,ONE,SSCP,TWO,XM,XNM,XN,ZERO
      DOUBLE PRECISION   CHISQ,SS
      DATA               ZERO,HALF,ONE,TWO/0.0,0.5,1.0,2.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.GE.1) GO TO 5
C                                  TERMINAL ERROR NUMBER OF RESPONSES
C                                  (N) SPECIFIED .LT. 1
      IER = 129
      GO TO 9000
    5 IER = 0
      STAT(1) = ZERO
      STAT(2) = ZERO
      STAT(3) = ZERO
      IOP1 = IOP(1)
      IOP2 = IOP(2)
      IF (M.LE.1) M = 0
      XN = N
      XM = M
      NDF = N+M
      XNM = XN+XM
      TEMP1 = CRIT(1)
      TEMP2 = XNM
C                                  COMPUTE VARIANCE
      SS = 0.D0
      DO 10 I=1,N
         CHISQ = Y(I)-CRIT(3)
         SS = SS+CHISQ*CHISQ
   10 CONTINUE
      SSCP = SS/XN
      IF (M.LT.2) VAR = SSCP
      IF (M.GE.2) VAR = (XN*SSCP+XM*VAR)/XNM
      IF (IOP1.LE.0.OR.IOP1.GT.3) GO TO 15
      CHISQ = (XNM*VAR)/CRIT(2)
      CHI = CHISQ
      CALL MDCH (CHI,TEMP2,P,IER)
      IF (IER.NE.0) GO TO 30
      STAT(1) = P
      IF (IOP1.EQ.2) STAT(1) = ONE-P
      IF (IOP1.NE.1) GO TO 15
      STAT(1) = P+P
      IF (P.GE.HALF) STAT(1) = TWO-STAT(1)
   15 IF (IOP2.LE.0.OR.IOP2.GT.3) GO TO 9005
      TEMP2 = CRIT(1)
      TEMP1 = ONE-CRIT(1)
      IF (IOP2.NE.1) GO TO 20
      TEMP1 = TEMP1*HALF
      TEMP2 = (ONE+CRIT(1))*HALF
   20 IF (IOP2.EQ.2) GO TO 25
      P = XNM
      CALL MDCHI (TEMP2,P,X2,IER)
      IF (IER.NE.0) GO TO 35
      STAT(2) = (XNM*VAR)/X2
   25 IF (IOP2.EQ.3) GO TO 9005
      P = XNM
      CALL MDCHI (TEMP1,P,X3,IER)
      IF (IER.NE.0) GO TO 35
      STAT(3) = (XNM*VAR)/X3
      GO TO 9005
   30 IER = 130
      GO TO 9000
   35 IER = 131
 9000 CONTINUE
      CALL UERTST (IER,'BENSON')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - BEMSON
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - MEAN AND VARIANCE INFERENCES USING A SAMPLE
C                           FROM A NORMAL POPULATION
C
C   USAGE               - CALL BEMSON (Y,N,IOP,CRIT,M,PAR,STAT,NDF,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SAMPLE.
C                N      - SAMPLE SIZE. (INPUT)
C                IOP    - INPUT VECTOR OF LENGTH 4 INDICATING HYPOTHESIS
C                           TESTS AND INTERVAL ESTIMATES TO BE COMPUTED.
C                           IOP(1) = I IMPLIES THE FOLLOWING TEST OF THE
C                           MEAN, WHEN
C                             I=1, TWO-TAILED.
C                             I=2, UPPER ONE-TAILED.
C                             I=3, LOWER ONE-TAILED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES TEST IS
C                             NOT DESIRED. SEE REMARKS.
C                           IOP(2) = I IMPLIES THE FOLLOWING INTERVAL
C                           ESTIMATE OF THE MEAN, WHEN
C                             I=1, TWO-SIDED.
C                             I=2, UPPER ONE-SIDED.
C                             I=3, LOWER ONE-SIDED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES ESTIMATE
C                             IS NOT DESIRED.
C                           IOP(3) = I IMPLIES THE FOLLOWING TEST OF THE
C                           VARIANCE, WHEN
C                             I=1, TWO-TAILED.
C                             I=2, UPPER ONE-TAILED.
C                             I=3, LOWER ONE-TAILED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES TEST IS
C                             NOT DESIRED. SEE REMARKS.
C                           IOP(4) = I IMPLIES THE FOLLOWING INTERVAL
C                           ESTIMATE OF THE VARIANCE, WHEN
C                             I=1, TWO-SIDED.
C                             I=2, UPPER ONE-SIDED.
C                             I=3, LOWER ONE-SIDED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES ESTIMATE
C                             IS NOT DESIRED.
C                CRIT   - INPUT VECTOR OF LENGTH 4 CONTAINING CONSTANTS
C                           REQUIRED FOR INTERVAL ESTIMATES AND HYPO-
C                           THESIS TESTS. THE I-TH ELEMENT OF CRIT CON-
C                           TAINS, WHEN
C                           I=1, CONFIDENCE COEFFICIENT, BETWEEN
C                             0.0 AND 1.0 EXCLUSIVELY, FOR INTERVAL
C                             ESTIMATE OF THE MEAN (REQUIRED ONLY WHEN
C                             IOP(2)=1,2, OR 3).
C                           I=2, CONFIDENCE COEFFICIENT, BETWEEN
C                             0.0 AND 1.0 EXCLUSIVELY, FOR INTERVAL
C                             ESTIMATE OF THE VARIANCE (REQUIRED ONLY
C                             WHEN IOP(4)=1, 2, OR 3).
C                           I=3, HYPOTHESIZED VALUE OF MEAN (REQUIRED
C                             ONLY WHEN IOP(1)=1,2, OR 3).
C                           I=4, HYPOTHESIZED VALUE OF VARIANCE (RE-
C                             QUIRED ONLY WHEN IOP(3)=1,2, OR 3).
C                M      - NUMBER OF RESPONSES IN EARLIER SAMPLE YIELDING
C                           ADDITIONAL PARAMETER ESTIMATES. M LESS THAN
C                           2 IMPLIES NO ADDITIONAL INFORMATION
C                           IS AVAILABLE. (INPUT)
C                PAR    - VECTOR OF LENGTH 2 CONTAINING ADDITIONAL PARA-
C                           METER ESTIMATES ON INPUT AND COMBINED ESTI-
C                           MATES ON OUTPUT. THE I-TH ELEMENT OF PAR
C                           CONTAINS, WHEN
C                           I=1, IF M IS GREATER THAN OR EQUAL TO
C                             2, ADDITIONAL ESTIMATE OF MEAN ON INPUT
C                             AND COMBINED ESTIMATE ON OUTPUT. IF
C                             M IS LESS THAN 2, UNDEFINED ON INPUT
C                             AND ESTIMATE OF MEAN ON OUTPUT.
C                           I=2, IF M IS GREATER THAN OR EQUAL TO
C                             2, ADDITIONAL ESTIMATE OF VARIANCE ON
C                             INPUT AND COMBINED ESTIMATE ON OUTPUT. IF
C                             M IS LESS THAN 2, UNDEFINED ON INPUT
C                             AND ESTIMATE OF VARIANCE ON OUTPUT.
C                STAT   - OUTPUT VECTOR OF LENGTH 6 CONTAINING PROBABIL-
C                           ITIES (ASSUMING NULL HYPOTHESES ARE TRUE) OF
C                           COMPUTED OR MORE EXTREME VALUES OF TEST STA-
C                           TISTICS AND CONTAINING INTERVAL ESTIMATES.
C                           THE I-TH ELEMENT OF STAT CONTAINS, WHEN
C                           I=1, PROBABILITY FOR HYPOTHESIS TEST ON MEAN
C                             (DEFINED ONLY WHEN IOP(1)=1,2, OR 3).
C                           I=2, PROBABILITY FOR HYPOTHESIS TEST ON
C                             VARIANCE (DEFINED ONLY WHEN IOP(3)=1,2,
C                             OR 3).
C                           I=3, LOWER LIMIT FOR MEAN (DEFINED ONLY WHEN
C                             IOP(2) = 1 OR 3).
C                           I=4, UPPER LIMIT FOR MEAN (DEFINED ONLY WHEN
C                             IOP(2)=1 OR 2).
C                           I=5, LOWER LIMIT FOR VARIANCE (DEFINED ONLY
C                             WHEN IOP(4) = 1 OR 3).
C                           I=6, UPPER LIMIT FOR VARIANCE (DEFINED ONLY
C                             WHEN IOP(4)=1 OR 2).
C                NDF    - NUMBER OF DEGREES OF FREEDOM IN VARIANCE
C                           ESTIMATE. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS NUMBER OF RESPONSES SPECIFIED
C                             IS LESS THAN 2.
C                           IER=130 MEANS EITHER CRIT(1) OR CRIT(2)
C                             WAS NOT IN THE INTERVAL (0,1), BUT THE
C                             CORRESPONDING CONFIDENCE INTERVAL(S) WAS
C                             (WERE) REQUESTED,  OR THAT CRIT(4) WAS
C                             NEGATIVE AND A TEST THAT THE VARIANCE
C                             IS EQUAL TO CRIT(4) WAS REQUESTED.
C                           IER=131 MEANS PAR(2) WAS NEGATIVE YET
C                             PAR(2) WAS A PRIOR ESTIMATE OF THE
C                             VARIANCE.
C                           IER=132 MEANS CONFIDENCE LIMIT(S) ON MEAN
C                             WAS (WERE) SO EXTREME THAT IF THE LOWER
C                             LIMIT (STAT(3)) WAS REQUESTED, IT WAS SET
C                             TO NEGATIVE MACHINE INFINITY AND IF THE
C                             UPPER LIMIT (STAT(4)) WAS REQUESTED, IT
C                             WAS SET TO MACHINE INFINITY.
C                           IER=133 MEANS ROUNDING ERRORS AND/OR
C                             EXTREME VALUES OF THE CONFIDENCE
C                             COEFFEICIENTS PREVENTED COMPUTATION OF
C                             STAT(5) AND/OR STAT(6).
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/MDBETA,MDCH,MDCHI,MDNOR,MDNRIS,
C                           MDSTI,MDTD,MERFI,MERRC=ERFC,MGAMAD=DGAMMA,
C                           MLGAMD=DLGAMA,UGETIO,UERTST,UGETIO
C                       - SINGLE(H36,H48,H60)/MDBETA,MDCH,MDCHI,MDNOR,
C                           MDNRIS,MDSTI,MDTD,MERFI,MERRC=ERFC,
C                           MGAMA=GAMMA,MLGAMA=ALGAMA,UERSET,UERTST,
C                           UGETIO
C                       - DOUBLE/MDBETA,MDCH,MDCHI,MDNOR,MDNRIS,MDSTI,
C                           MDTD,MERFI,MERRC=ERFC,MGAMAD=DGAMMA,
C                           MLGAMD=DLGAMA,UERSET,UERTST,UGETIO,VXADD,
C                           VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      LETTING (M,S) BE THE TRUE MEAN AND VARIANCE FOR THE
C                POPULATION, RESPECTIVELY, AND DEFINING
C                            M0 = CRIT(3)
C                            S0 = CRIT(4)
C                THE TABLE BELOW CLARIFIES HYPOTHESIS TESTING OPTIONS.
C
C                      IOP(1)            NULL            ALTERNATIVE
C                                     HYPOTHESIS         HYPOTHESIS
C                      ------         ----------         -----------
C                        1              M.EQ.M0            M.NE.M0
C                        2              M.LE.M0            M.GT.M0
C                        3              M.GE.M0            M.LT.M0
C
C
C                      IOP(3)            NULL            ALTERNATIVE
C                                     HYPOTHESIS         HYPOTHESIS
C                      ------         ----------         -----------
C                        1              S.EQ.S0            S.NE.S0
C                        2              S.LE.S0            S.GT.S0
C                        3              S.GE.S0            S.LT.S0
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEMSON (Y,N,IOP,CRIT,M,PAR,STAT,NDF,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IOP(4),M,NDF,IER
      REAL               Y(N),CRIT(4),PAR(2),STAT(6)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IOP1,IOP2,IOP3,IOP4,I,KER,LEVEL,LEVOLD
      REAL               HALF,ONE,S2SS,SDN,TC,XINF,XINFM,
     1                   TWO,VAR,XMN,XM,XNDF,XN,YMN,ZERO
      REAL               CHI
      DOUBLE PRECISION   CHISQ,SSCP,SS,YCOM,YC
      DATA               ZERO,HALF,ONE,TWO/0.0,.5,1.0,2.0/
      DATA               XINF/Z7FFFFFFF/
      DATA               XINFM/ZFFFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      LEVEL = 0
      CALL UERSET(LEVEL,LEVOLD)
      IOP1 = IOP(1)
      IOP2 = IOP(2)
      IOP3 = IOP(3)
      IOP4 = IOP(4)
      IF (M.GE.2) GO TO 10
      PAR(1) = 0.0
      PAR(2) = 0.0
      GO TO 20
   10 IF (PAR(2) .GE. 0.0) GO TO 20
      IER = 131
      GO TO 9000
   20 TEMP1 = CRIT(1)
      TEMP2 = CRIT(2)
      IF (IOP2.LT.1 .OR. IOP2.GT.3) GO TO 30
      IF (TEMP1.GT.0.0 .AND. TEMP1.LT.1.0) GO TO 30
      IER = 130
      GO TO 9000
   30 IF (IOP4.LT.1 .OR. IOP4.GT.3) GO TO 40
      IF (TEMP2.GT.0.0 .AND. TEMP2.LT.1.0) GO TO 40
      IER = 130
      GO TO 9000
   40 IF (IOP3.LT.1 .OR. IOP3.GT.3) GO TO 50
      IF (CRIT(4).GE.0.0) GO TO 50
      IER = 130
      GO TO 9000
   50 IF (N.GE.2) GO TO 60
      IER = 129
      GO TO 9000
   60 CONTINUE
      DO 70 I=1,6
         STAT(I) = ZERO
   70 CONTINUE
      NDF = N-1
      YMN = PAR(1)
      VAR = PAR(2)
      XNDF = NDF
      XMN = M+N
      XMN1 = XMN-ONE
      IF (M.LT.2) XMN1 = XNDF
      XM = M
      XN = N
C                                  COMPUTE MEAN FROM SAMPLE
      YC = 0.D0
      DO 80 I=1,N
         YC = YC+Y(I)
   80 CONTINUE
      YCOM = YC/XN
C                                  COMPUTE VARIANCE FROM SAMPLE
      SS = 0.D0
      DO 90 I=1,N
         CHISQ = Y(I)-YCOM
         SS = SS+CHISQ*CHISQ
   90 CONTINUE
      SSCP = SS/XNDF
      SDN = DSQRT(SSCP/XN)
      PAR(1) = YCOM
      PAR(2) = SSCP
      IF (M.LE.1) GO TO 100
C                                  COMPUTE MEAN AND VARIANCE USING
C                                  USER*S INPUT DATA
      S2SS = ((XN*XM)/XMN)*(YCOM-YMN)**2
      PAR(1) = (XN*YCOM+XM*YMN)/XMN
      NDF = NDF+M
      XNDF = NDF
      PAR(2) = ((XN-ONE)*SSCP+(XM-ONE)*VAR+S2SS)/XNDF
      XN = XMN
      SDN = SQRT(PAR(2)/XN)
  100 IF (IOP1.LE.0.OR.IOP1.GT.3) GO TO 110
C                                  COMPUTE T PROBABILITY
      TC = (PAR(1)-CRIT(3))/SDN
      TEMP3 = XNDF
      TEMP4 = TC
      CALL MDTD (TEMP4,TEMP3,R,KER)
      STAT(1) = R
      IF (IOP1.EQ.2) TC = -TC
      IF (IOP1.EQ.1) GO TO 110
      STAT(1) = HALF*STAT(1)
      IF (TC.GE.ZERO) STAT(1) = ONE-STAT(1)
C                                  COMPUTE INTERVAL FOR MEAN
  110 IF (IOP2.GT.3.OR.IOP2.LE.0) GO TO 150
      X = 0.
      TD = ONE-CRIT(1)
      IF (IOP2.EQ.1) GO TO 120
      IF (CRIT(1).EQ.HALF) GO TO 140
      TD = TD+TD
      IF (CRIT(1).LT.HALF) TD = CRIT(1)+CRIT(1)
  120 CALL MDSTI (TD,XMN1,X,KER)
      IF (KER.EQ.0) GO TO 130
      IER = 132
      IF (IOP2.NE.2) STAT(3) = XINFM
      IF (IOP2.NE.3) STAT(4) = XINF
      GO TO 150
  130 IF (IOP2.NE.1.AND.CRIT(1).LT.HALF) X = -X
  140 IF (IOP2.NE.2) STAT(3) = PAR(1)-X*SDN
      IF (IOP2.NE.3) STAT(4) = PAR(1)+X*SDN
  150 IF (IOP3.LE.0.OR.IOP3.GT.3) GO TO 160
      CHISQ = (XN-ONE)*PAR(2)/CRIT(4)
      CHI = CHISQ
      TEMP3 = XNDF
      CALL MDCH (CHI,TEMP3,P,KER)
      STAT(2) = P
      IF (IOP3.EQ.2) STAT(2) = ONE-P
      IF (IOP3.NE.1) GO TO 160
      STAT(2) = P+P
      IF (P.GE.HALF) STAT(2) = TWO-STAT(2)
  160 IF (IOP4.GT.3.OR.IOP4.LE.0) GO TO 9005
      TEMP1 = ONE-CRIT(2)
      IF (IOP4.EQ.1) TEMP1 = TEMP1*HALF
      IF (IOP4.GT.2) GO TO 170
      CALL MDCHI (TEMP1,XMN1,X6,KER)
      IF (KER.NE.0) GO TO 180
      STAT(6) = ((XN-ONE)*PAR(2))/X6
      IF (IOP4.EQ.1) TEMP2 = (TEMP2+ONE)*HALF
  170 IF (IOP4.EQ.2) GO TO 9005
      CALL MDCHI (TEMP2,XMN1,X5,KER)
      IF (KER.NE.0) GO TO 180
      STAT(5) = ((XN-ONE)*PAR(2))/X5
      GO TO 9005
C                                  TERMINAL ERROR ON RETURN FROM MDCHI
  180 IER = 133
 9000 CONTINUE
      CALL UERTST (IER,'BEMSON')
 9005 CALL UERSET(LEVOLD,LEVEL)
      RETURN
      END

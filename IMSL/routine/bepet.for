C   IMSL ROUTINE NAME   - BEPET
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - MEAN AND VARIANCE INFERENCES USING SAMPLES
C                           FROM EACH OF TWO NORMAL POPULATIONS WITH
C                           EQUAL VARIANCES
C
C   USAGE               - CALL BEPET  (Y,N,IOP,CRIT,STAT,NDF,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N(1)+N(2) CONTAINING
C                           FIRST SAMPLE FOLLOWED BY SECOND SAMPLE.
C                N      - INPUT VECTOR OF LENGTH 2 CONTAINING SAMPLE
C                           SIZES. THE I-TH ELEMENT OF N CONTAINS, WHEN
C                           I = 1, SIZE OF SAMPLE ONE.
C                           I = 2, SIZE OF SAMPLE TWO.
C                IOP    - INPUT VECTOR OF LENGTH 4 INDICATING HYPOTHESIS
C                           TESTS AND INTERVAL ESTIMATES TO BE COMPUTED.
C                           IOP(1) = I IMPLIES THE FOLLOWING TEST OF THE
C                           DIFFERENCE BETWEEN THE MEANS, WHEN
C                             I=1, TWO-TAILED.
C                             I=2, UPPER ONE-TAILED.
C                             I=3, LOWER ONE-TAILED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES TEST IS
C                             NOT DESIRED. SEE REMARKS.
C                           IOP(2) = I IMPLIES THE FOLLOWING INTERVAL
C                           ESTIMATE OF THE DIFFERENCE BETWEEN THE
C                           MEANS, WHEN
C                             I=1, TWO-SIDED.
C                             I=2, UPPER ONE-SIDED.
C                             I=3, LOWER ONE-SIDED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES ESTIMATE
C                             IS NOT DESIRED. SEE REMARKS.
C                           IOP(3) = I IMPLIES THE FOLLOWING TEST OF
C                           THE VARIANCE, WHEN
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
C                CRIT   - INPUT VECTOR OF LENGTH 5 CONTAINING CONSTANTS
C                           REQUIRED FOR INTERVAL ESTIMATES AND HYPOTHE-
C                           SIS TESTS. THE I-TH ELEMENT OF CRIT
C                           CONTAINS, WHEN
C                           I = 1, CONFIDENCE COEFFICIENT, BETWEEN
C                             0.0 AND 1.0 EXCLUSIVELY, FOR INTERVAL
C                             ESTIMATE INVOLVING THE MEANS. REQUIRED
C                             ONLY WHEN IOP(2) = 1,2, OR 3.
C                           I = 2, CONFIDENCE COEFFICIENT, BETWEEN
C                             0.0 AND 1.0 EXCLUSIVELY, FOR INTERVAL
C                             ESTIMATE OF THE VARIANCE. REQUIRED ONLY
C                             WHEN IOP(4) = 1, 2, OR 3.
C                           I = 3, HYPOTHESIZED VALUE OF VARIANCE. RE-
C                             QUIRED ONLY WHEN IOP(3) = 1,2, OR 3.
C                           I = 4, MULTIPLICATIVE CONSTANT. REQUIRED
C                             ONLY WHEN IOP(1) OR IOP(2) EQUALS 1,2,
C                             OR 3. NORMALLY, CRIT(4) = 1.0 IS
C                             APPROPRIATE. SEE REMARKS.
C                           I = 5, ADDITIVE CONSTANT. REQUIRED ONLY
C                             WHEN IOP(1) OR IOP(2) EQUALS 1,2, OR 3.
C                             NORMALLY, CRIT(5) = 0.0 IS APPROPRIATE.
C                             SEE REMARKS.
C                STAT   - OUTPUT VECTOR OF LENGTH 9 CONTAINING RESULTS
C                           OF INFERENCES ABOUT MEANS AND VARIANCE.
C                           LET H0 AND H1 BE NULL AND ALTERNATIVE
C                           HYPOTHESES, RESPECTIVELY. THE
C                           I-TH ELEMENT OF STAT CONTAINS, WHEN
C                           I = 1, ESTIMATE OF FIRST MEAN.
C                           I = 2, ESTIMATE OF SECOND MEAN.
C                           I = 3, ESTIMATE OF VARIANCE.
C                           I = 4, FOR HYPOTHESIS TEST ON THE MEANS, THE
C                             PROBABILITY ASSUMING H0 IS TRUE, OF THE
C                             COMPUTED OR MORE EXTREME VALUE OF THE
C                             TEST STATISTIC. DEFINED ONLY WHEN IOP(1)
C                             = 1,2,OR 3.
C                           I = 5, FOR HYPOTHESIS TEST ON THE VARIANCE,
C                             THE PROBABILITY, ASSUMING H0 IS TRUE, OF
C                             THE COMPUTED OR MORE EXTREME VALUE OF THE
C                             TEST STATISTIC. DEFINED ONLY WHEN IOP(3)
C                             = 1,2, OR 3.
C                           I = 6, LOWER LIMIT FOR INTERVAL ESTIMATE
C                             INVOLVING THE MEANS. DEFINED ONLY WHEN
C                             IOP(2) = 1 OR 3.
C                           I = 7, UPPER LIMIT FOR INTERVAL ESTIMATE
C                             INVOLVING THE MEANS. DEFINED ONLY WHEN
C                             IOP(2) = 1 OR 2.
C                           I = 8, LOWER LIMIT FOR INTERVAL ESTIMATE OF
C                             VARIANCE. DEFINED ONLY WHEN IOP(4) = 1 OR
C                             3.
C                           I = 9, UPPER LIMIT FOR INTERVAL ESTIMATE OF
C                             VARIANCE. DEFINED ONLY WHEN IOP(4) = 1 OR
C                             2.
C                NDF    - NUMBER OF DEGREES OF FREEDOM IN VARIANCE
C                           ESTIMATE. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS TOTAL NUMBER OF RESPONSES
C                             SPECIFIED IS LESS THAN 3.
C                           IER=130 MEANS AN ERROR OCCURRED IN MDTD.
C                           IER=131 MEANS AN ERROR OCCURRED IN MDCH.
C                           IER=132 MEANS AN ERROR OCCURRED IN MDSTI.
C                           IER=133 MEANS AN ERROR OCCURRED IN MDCHI.
C
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/MDBETA,MDCH,MDCHI,MDNOR,MDNRIS,
C                           MDSTI,MDTD,MERFI,MERRC=ERFC,MGAMAD=DGAMMA,
C                           MLGAMD=DLGAMA,UERTST,UGETIO
C                       - SINGLE(H36,H48,H60)/MDBETA,MDCH,MDCHI,MDNOR,
C                           MDNRIS,MDSTI,MDTD,MERFI,MERRC=ERFC,
C                           MGAMA=GAMMA,MLGAMA=ALGAMA,UERTST,UGETIO
C                       - DOUBLE/MDBETA,MDCH,MDCHI,MDNOR,MDNRIS,MDSTI,
C                           MDTD,MERFI,MERRC=ERFC,MGAMAD=DGAMMA,
C                           MLGAMD=DLGAMA,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  LETTING (M1,S) AND (M2,S) BE THE TRUE MEANS AND
C                VARIANCES FOR POPULATIONS ONE AND TWO, RESPECTIVELY,
C                AND DEFINING
C                            M3 = CRIT(4)*M2+CRIT(5)
C                            S0 = CRIT(3)
C                THE TABLE BELOW CLARIFIES HYPOTHESIS TESTING OPTIONS.
C
C                      IOP(1)            NULL            ALTERNATIVE
C                                     HYPOTHESIS         HYPOTHESIS
C                      ------         ----------         -----------
C                        1             M1.EQ.M3           M1.NE.M3
C                        2             M1.LE.M3           M1.GT.M3
C                        3             M1.GE.M3           M1.LT.M3
C
C
C                      IOP(3)            NULL            ALTERNATIVE
C                                     HYPOTHESIS         HYPOTHESIS
C                      ------         ----------         -----------
C                        1             S .EQ.S0           S .NE.S0
C                        2             S .LE.S0           S .GT.S0
C                        3             S .GE.S0           S .LT.S0
C
C            2.  THE MEANS INTERVAL ESTIMATE IS FOR M1-M3.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEPET  (Y,N,IOP,CRIT,STAT,NDF,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOP(4),N(2),NDF,IER
      REAL               Y(1),CRIT(5),STAT(9)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   CHISQ,SS,YC,S,TL
      REAL               ZERO,HALF,ONE,TWO,XN1
      REAL               Q,XN2,XNDF,VAR1,VAR2,TC
      DATA               ZERO,HALF,ONE,TWO/0.0,0.5,1.0,2.0/
C                                  FIRST EXECUTABLE STATEMENT
      ISW = 0
    5 IF (N(1)+N(2).LT.3) GO TO 75
      DO 10 I=1,9
         STAT(I) = ZERO
   10 CONTINUE
      IER = 0
      XN1 = N(1)
      XN2 = N(2)
      IOP1 = IOP(1)
      IOP2 = IOP(2)
      IOP3 = IOP(3)
      IOP4 = IOP(4)
      N1 = N(1)
      N2 = N(2)
      NDF = N1+N2-2
      XNDF = NDF
      TEMP3 = XNDF
C                                  COMPUTE FIRST MEAN
      YC = 0.D0
      DO 15 I=1,N1
         YC = YC+Y(I)
   15 CONTINUE
      STAT(1) = YC/XN1
C                                  COMPUTE SECOND MEAN
      YC = 0.D0
      DO 20 I=1,N2
         YC = YC+Y(I+N1)
   20 CONTINUE
      STAT(2) = YC/XN2
C                                  COMPUTE FIRST VARIANCE
      SS = 0.D0
      DO 25 I=1,N1
         CHISQ = Y(I)-STAT(1)
         SS = SS+CHISQ*CHISQ
   25 CONTINUE
      VAR1 = SS
C                                  COMPUTE SECOND VARIANCE
      SS = 0.D0
      DO 30 I=1,N2
         CHISQ = Y(I+N1)-STAT(2)
         SS = SS+CHISQ*CHISQ
   30 CONTINUE
      VAR2 = SS
      STAT(3) = (VAR1+VAR2)/XNDF
      TL = STAT(1)-CRIT(4)*STAT(2)-CRIT(5)
      S = SQRT(STAT(3)*(ONE/XN1+(CRIT(4)*CRIT(4))/XN2))
   35 IF (IOP1.LE.0.OR.IOP1.GT.3) GO TO 40
C                                  COMPUTE T PROBABILITY
      TC = TL/S
      TEMP4 = TC
      CALL MDTD (TEMP4,TEMP3,R,IER)
      IF (IER.NE.0) GO TO 80
      STAT(4) = R
      IF (IOP1.EQ.1) GO TO 40
      IF (IOP1.EQ.2) TC = -TC
      STAT(4) = HALF*STAT(4)
      IF (TC.GE.ZERO) STAT(4) = ONE-STAT(4)
   40 IF (IOP2.LE.0.OR.IOP2.GT.3) GO TO 55
C                                  COMPUTE INTERVALS FOR MEAN
      X = ZERO
      TD = ONE-CRIT(1)
      IF (IOP2.EQ.1) GO TO 45
      IF (CRIT(1).EQ.HALF) GO TO 50
      TD = TD+TD
      IF (CRIT(1).LT.HALF) TD = CRIT(1)+CRIT(1)
   45 CALL MDSTI (TD,TEMP3,X,IER)
      IF (IER.NE.0) GO TO 90
      IF (IOP2.EQ.2.AND.CRIT(1).LT.HALF) X = -X
   50 IF (IOP2.NE.2) STAT(6) = TL-X*S
      IF (IOP2.NE.3) STAT(7) = TL+X*S
   55 IF (IOP3.LE.0.OR.IOP3.GT.3) GO TO 60
      CHI = ((XNDF)*STAT(3))/CRIT(3)
      CALL MDCH (CHI,TEMP3,P,IER)
      IF (IER.NE.0) GO TO 85
      STAT(5) = P
      IF (IOP3.EQ.2) STAT(5) = ONE-P
      IF (IOP3.NE.1) GO TO 60
      STAT(5) = P+P
      IF (P.GE.HALF) STAT(5) = TWO-STAT(5)
   60 IF (IOP4.GT.3.OR.IOP4.LE.0) GO TO 100
      TEMP1 = ONE-CRIT(2)
      TEMP2 = CRIT(2)
      IF (IOP4.NE.1) GO TO 65
      TEMP1 = (ONE-CRIT(2))*HALF
      TEMP2 = (ONE+CRIT(2))*HALF
   65 IF (IOP4.EQ.3) GO TO 70
      CALL MDCHI (TEMP1,TEMP3,X9,IER)
      IF (IER.NE.0) GO TO 95
      STAT(9) = (XNDF*STAT(3))/X9
   70 IF (IOP4.EQ.2) GO TO 100
      CALL MDCHI (TEMP2,TEMP3,X8,IER)
      IF (IER.NE.0) GO TO 95
      STAT(8) = (XNDF*STAT(3))/X8
      GO TO 100
C                                  TERMINAL ERROR N1 + N2 .LT. 3
   75 IER = 129
      GO TO 9000
C                                  TERMINAL ERROR ON RETURN FROM MDTD
   80 IER = 130
      GO TO 9000
C                                  TERMINAL ERROR ON RETURN FROM MDCH
   85 IER = 131
      GO TO 9000
C                                  TERMINAL ERROR ON RETURN FROM MDSTI
   90 IER = 132
      GO TO 9000
C                                  TERMINAL ERROR ON RETURN FROM MDCHI
   95 IER = 133
      GO TO 9000
  100 IF (ISW.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BEPET ')
 9005 RETURN
      END

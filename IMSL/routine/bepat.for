C   IMSL ROUTINE NAME   - BEPAT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - MEAN AND VARIANCE INFERENCES USING SAMPLES
C                           FROM EACH OF TWO NORMAL POPULATIONS WITH
C                           UNEQUAL VARIANCES
C
C   USAGE               - CALL BEPAT  (Y,N,IOP,CRIT,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N(1)+N(2) CONTAINING
C                           FIRST SAMPLE FOLLOWED BY SECOND SAMPLE.
C                N      - INPUT VECTOR OF LENGTH 2 CONTAINING SAMPLE
C                           SIZES. THE I-TH ELEMENT OF N CONTAINS, WHEN
C                           I = 1, SIZE OF SAMPLE ONE.
C                           I = 2, SIZE OF SAMPLE TWO.
C                IOP    - INPUT VECTOR OF LENGTH 5 INDICATING HYPOTHE-
C                           SIS TESTS AND INTERVAL ESTIMATES TO BE
C                           COMPUTED. IOP(1)=I IMPLIES THE FOLLOWING
C                           TEST OF THE DIFFERENCE BETWEEN THE MEANS.
C                           WHEN
C                             I=1, TWO-TAILED.
C                             I=2, UPPER ONE-TAILED.
C                             I=3, LOWER ONE-TAILED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES TEST IS
C                             NOT DESIRED. SEE REMARKS.
C                           IOP(2) = I IMPLIES THE FOLLOWING INTERVAL
C                           ESTIMATE OF THE DIFFERENCE BETWEEN THE
C                           MEANS. WHEN
C                             I=1, TWO-SIDED.
C                             I=2, UPPER ONE-SIDED.
C                             I=3, LOWER ONE-SIDED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES ESTIMATE
C                             IS NOT DESIRED. SEE REMARKS.
C                           IOP(3) = I IMPLIES THE FOLLOWING TEST OF
C                           THE DIFFERENCE BETWEEN THE VARIANCES. WHEN
C                             I=1, TWO-TAILED.
C                             I=2, UPPER ONE-TAILED.
C                             I=3, LOWER ONE-TAILED.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES TEST IS
C                             NOT DESIRED. SEE REMARKS.
C                           IOP(4) = I IMPLIES AN INTERVAL ESTIMATE
C                           OF THE FOLLOWING RATIO OF VARIANCES. WHEN
C                             I=1, FIRST OVER SECOND.
C                             I=2, SECOND OVER FIRST.
C                             I NOT EQUAL TO 1,2, OR 3 IMPLIES ESTIMATE
C                             IS NOT DESIRED. SEE REMARKS.
C                           IOP(5) = I IMPLIES (REQUIRED ONLY IF IOP(4)
C                           EQUALS 1 OR 2) THE FOLLOWING TYPE OF
C                           ESTIMATE FOR THE RATIO SPECIFIED BY IOP(4).
C                           WHEN
C                             I=1, TWO-SIDED.
C                             I=2, UPPER ONE-SIDED.
C                             I NOT EQUAL TO 1 OR 2, LOWER ONE-SIDED.
C                             SEE REMARKS.
C                CRIT   - INPUT VECTOR OF LENGTH 5 CONTAINING CONSTANTS
C                             REQUIRED FOR INTERVAL ESTIMATES AND HYPO-
C                             THESIS TESTS. THE I-TH ELEMENT OF CRIT
C                             CONTAINS WHEN
C                           I = 1, CONFIDENCE COEFFICIENT, BETWEEN
C                             0.0 AND 1.0 EXCLUSIVELY, FOR INTERVAL
C                             ESTIMATE INVOLVING THE MEANS. REQUIRED
C                             ONLY WHEN IOP(2) = 1,2, OR 3.
C                           I = 2, CONFIDENCE COEFFICIENT, BETWEEN
C                             0.0 AND 1.0 EXCLUSIVELY, FOR INTERVAL
C                             ESTIMATE OF RATIO OF VARIANCES.
C                             REQUIRED ONLY WHEN IOP(4) = 1 OR 2.
C                           I = 3, MULTIPLICATIVE CONSTANT. REQUIRED
C                             ONLY WHEN IOP(1) OR IOP(2) EQUALS
C                             1,2, OR 3. NORMALLY, CRIT(3)=1.0
C                             IS APPROPRIATE. SEE REMARKS.
C                           I = 4, ADDITIVE CONSTANT. REQUIRED ONLY
C                             WHEN IOP(1) OR IOP(2) EQUALS 1,2, OR 3.
C                             NORMALLY, CRIT(4)=0.0 IS APPROPRIATE. SEE
C                             REMARKS.
C                           I = 5, ADDITIONAL MULTIPLICATIVE CONSTANT.
C                             REQUIRED ONLY WHEN IOP(3)=1,2, OR 3
C                             OR WHEN IOP(4)=1 OR 2. NORMALLY,
C                             CRIT(5)=1.0 IS APPROPRIATE. SEE REMARKS.
C                STAT   - OUTPUT VECTOR OF LENGTH 11 CONTAINING RESULTS
C                           OF INFERENCES ABOUT MEANS AND VARIANCES.
C                           LET H0 AND H1 BE NULL AND ALTERNATIVE
C                           HYPOTHESES, RESPECTIVELY. THE
C                           I-TH ELEMENT OF STAT CONTAINS WHEN
C                           I = 1, ESTIMATE OF FIRST MEAN.
C                           I = 2, ESTIMATE OF SECOND  MEAN.
C                           I = 3, ESTIMATE OF FIRST VARIANCE.
C                           I = 4, ESTIMATE OF SECOND VARIANCE.
C                           I = 5, FOR HYPOTHESIS TEST ON THE MEANS,
C                             THE APPROXIMATE PROBABILITY, ASSUMING H0
C                             IS TRUE, OF THE COMPUTED OR MORE EXTREME
C                             VALUE OF THE TEST STATISTIC. DEFINED ONLY
C                             WHEN IOP(1)=1,2, OR 3.
C                           I = 6, FOR HYPOTHESIS TEST ON THE VARIANCES,
C                             THE PROBABILITY, ASSUMING H0 IS TRUE, OF
C                             THE COMPUTED OR MORE EXTREME VALUE OF THE
C                             TEST STATISTIC. DEFINED ONLY WHEN IOP(3)=
C                             1,2, OR 3.
C                           I = 7, LOWER LIMIT FOR INTERVAL ESTIMATE
C                             INVOLVING THE MEANS. DEFINED ONLY WHEN
C                             IOP(2) = 1 OR 3.
C                           I = 8, UPPER LIMIT FOR INTERVAL ESTIMATE
C                             INVOLVING THE MEANS. DEFINED ONLY WHEN
C                             IOP(2) = 1 OR 2.
C                           I = 9, LOWER LIMIT FOR INTERVAL ESTIMATE OF
C                             RATIO OF VARIANCES. DEFINED ONLY WHEN BOTH
C                             IOP(4) = 1 OR 2 AND IOP(5) NOT EQUAL TO 2.
C                           I = 10, UPPER LIMIT FOR INTERVAL ESTIMATE
C                             OF RATIO OF VARIANCES. DEFINED ONLY WHEN
C                             BOTH IOP(4) = 1 OR 2 AND IOP(5) = 1 OR 2.
C                           I = 11, NUMBER OF DEGREES OF FREEDOM
C                             ASSOCIATED WITH APPROXIMATE TEST AND
C                             INTERVAL ESTIMATE INVOLVING THE MEANS.
C                             IN GENERAL, STAT(11) WILL NOT BE AN
C                             INTEGER.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS NUMBER OF RESPONSES SPECIFIED
C                             IN FIRST SAMPLE IS LESS THAN 2.
C                           IER=130 MEANS NUMBER OF RESPONSES SPECIFIED
C                             IN SECOND SAMPLE IS LESS THAN 2
C                           IER=131 MEANS AN ERROR OCCURRED IN MDBETA.
C                           IER=132 MEANS AN ERROR OCCURRED IN MDFI.
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/MDBETA,MDBETI,MDFI,MLGAMD=DLGAMA,
C                           UERTST,UGETIO
C                       - SINGLE(H36,H48,H60)/MDBETA,MDBETI,MDFI,
C                           MLGAMA=ALGAMA,UERTST,UGETIO
C                       - DOUBLE/MDBETA,MDBETI,MDFI,MLGAMD=DLGAMA,
C                           UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  LETTING (M1,S1) AND (M2,S2) BE THE TRUE MEANS AND
C                VARIANCES FOR POPULATIONS ONE AND TWO, RESPECTIVELY,
C                AND DEFINING
C                            M3 = CRIT(3)*M2+CRIT(4)
C                            S3 = CRIT(5)*S2
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
C                        1             S1.EQ.S3           S1.NE.S3
C                        2             S1.LE.S3           S1.GT.S3
C                        3             S1.GE.S3           S1.LT.S3
C
C            2.  THE MEANS INTERVAL ESTIMATE IS FOR M1-M3.
C            3.  THE VARIANCES INTERVAL ESTIMATE IS FOR S1/S3 OR
C                S3/S1 AS IOP(4) IS 1 OR 2, RESPECTIVELY.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEPAT  (Y,N,IOP,CRIT,STAT,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOP(1),N(1),IER
      REAL               Y(1),CRIT(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   CHISQ,SD,VDV,VDV1,TL,F,TC
      REAL               XN1,XN2,SS,YC,ZERO,HALF
      REAL               ONE,TWO,CRIT3,CRIT4,CRIT5
      DATA               ZERO,HALF,ONE,TWO/0.0,0.5,1.0,2.0/
      DATA               SHALF,ONF/0.5,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      ISW = 0
      CRITF2 = CRIT(2)
      CRIT3 = CRIT(3)
      CRIT4 = CRIT(4)
      CRIT5 = CRIT(5)
      IF (N(1).LT.2) GO TO 140
C                                  TERMINAL ERROR N2 .LT. 2
      IF (N(2).GT.1) GO TO 5
      IER = 130
      GO TO 9000
    5 IER = 0
      DO 10 I=1,11
   10 STAT(I) = ZERO
      XN1 = N(1)
      XN2 = N(2)
      IOP1 = IOP(1)
      IOP2 = IOP(2)
      IOP3 = IOP(3)
      IOP4 = IOP(4)
      IOP5 = IOP(5)
      N1 = N(1)
      N2 = N(2)
C                                  COMPUTE FIRST SAMPLE MEAN
      YC = 0.D0
      DO 15 I=1,N1
   15 YC = YC+Y(I)
      STAT(1) = YC/XN1
C                                  COMPUTE SECOND SAMPLE MEAN
      YC = 0.D0
      DO 20 I=1,N2
   20 YC = YC+Y(I+N1)
      STAT(2) = YC/XN2
C                                  COMPUTE FIRST SAMPLE VARIANCE
      SS = 0.D0
      DO 25 I=1,N1
         CHISQ = Y(I)-STAT(1)
   25 SS = SS+CHISQ*CHISQ
      STAT(3) = SS/(XN1-ONE)
C                                  COMPUTE SECOND SAMPLE VARIANCE
      SS = 0.D0
      DO 30 I=1,N2
         CHISQ = Y(I+N1)-STAT(2)
   30 SS = SS+CHISQ*CHISQ
      STAT(4) = SS/(XN2-ONE)
      TL = STAT(1)-CRIT(3)*STAT(2)-CRIT(4)
      SD = STAT(3)/XN1+(CRIT(3)*CRIT(3))*STAT(4)/XN2
      F = SD**2/(((STAT(3)/XN1)**2/(XN1-ONE))+((CRIT(3)**2*STAT(4)/XN2)
     1**2)/(XN2-ONE))
      STAT(11) = F
      SD = DSQRT(SD)
      IF (IOP1.LE.0.OR.IOP1.GT.3) GO TO 45
C                                  COMPUTE T PROBABILITY
      TC = TL/SD
      TT = F/(F+TC**2)
      FF = F*0.5D0
      CALL MDBETA (TT,FF,SHALF,R,IEM)
      IF (IEM.NE.0) GO TO 130
      STAT(5) = R
      GO TO (45,35,40), IOP1
   35 TC = -TC
   40 STAT(5) = HALF*STAT(5)
      IF (TC.LT.ZERO) GO TO 45
      STAT(5) = ONE-STAT(5)
   45 IF (IOP2.GT.3.OR.IOP2.LE.0) GO TO 65
C                                  COMPUTE INTERVALS VALUES FOR MEANS
      IF (IOP2.NE.1.AND.CRIT(1).LE.HALF) GO TO 50
      TT = CRIT(1)
      IF (IOP2.NE.1) TT = CRIT(1)+CRIT(1)-ONE
      STATF = STAT(11)
      CALL MDFI (TT,ONF,STATF,X,JER)
      IF (JER.NE.0) GO TO 135
      X = SQRT(X)
      GO TO 60
   50 IF (CRIT(1).LT.HALF) GO TO 55
      X = ZERO
      GO TO 60
   55 PF = ONE-TWO*CRIT(1)
      CALL MDFI (PF,ONF,STATF,X,JER)
      IF (JER.NE.0) GO TO 135
      X = SQRT(X)
      IF (IOP2.EQ.2) X = -X
   60 IF (IOP2.NE.2) STAT(7) = TL-X*SD
      IF (IOP2.NE.3) STAT(8) = TL+X*SD
   65 N1 = N1-1
      N2 = N2-1
      VDV = STAT(3)/(CRIT(5)*STAT(4))
      VDV1 = (STAT(4)*CRIT(5))/STAT(3)
      IF (IOP3.LE.0.OR.IOP3.GT.3) GO TO 85
      SVDV = N2/(N2+N1*VDV)
      CALL MDBETA (SVDV,N2*SHALF,N1*SHALF,Q,JER)
      IF (JER.NE.0) GO TO 130
      GO TO (70,75,80), IOP3
   70 STAT(6) = TWO*AMIN1(Q,1.0-Q)
      GO TO 85
   75 STAT(6) = Q
      GO TO 85
   80 STAT(6) = ONE-Q
   85 IF (IOP4.LE.0.OR.IOP4.GT.2) GO TO 145
      FN1 = N1
      FN2 = N2
      IF (IOP5.NE.1.AND.IOP5.NE.2) GO TO 90
      PP = (CRIT(2)+ONE)*HALF
      CALL MDFI (PP,FN1,FN2,X2,JER)
      IF (JER.NE.0) GO TO 135
      CALL MDFI (PP,FN2,FN1,X3,JER)
      IF (JER.NE.0) GO TO 135
   90 GO TO (95,110), IOP4
   95 IF (IOP5.NE.1) GO TO 100
      STAT(9) = VDV/X2
      STAT(10) = VDV*X3
      GO TO 145
  100 IF (IOP5.NE.2) CALL MDFI (CRITF2,FN1,FN2,X4,JER)
      IF (IOP5.EQ.2) CALL MDFI (CRITF2,FN2,FN1,X5,JER)
      IF (JER.NE.0) GO TO 135
      IF (IOP5.NE.2) GO TO 105
      STAT(10) = VDV*X5
      GO TO 145
  105 STAT(9) = VDV/X4
      GO TO 145
  110 IF (IOP5.EQ.2) CALL MDFI (CRITF2,FN1,FN2,X4,JER)
      IF (IOP5.NE.2) CALL MDFI (CRITF2,FN2,FN1,X5,JER)
      IF (IOP5.NE.1.AND.IOP5.NE.2) GO TO 125
      GO TO (115,120), IOP5
  115 STAT(9) = VDV1/X3
      STAT(10) = VDV1*X2
      GO TO 145
  120 STAT(10) = VDV1*X4
      GO TO 145
  125 STAT(9) = VDV1/X5
      GO TO 145
C                                  TERMINAL ERROR ON RETURN FROM MDBETA
  130 IER = 131
      CRIT(3) = CRIT3
      CRIT(4) = CRIT4
      CRIT(5) = CRIT5
      GO TO 9000
C                                  TERMINAL ERROR ON RETURN FROM MDFI
  135 IER = 132
      CRIT(3) = CRIT3
      CRIT(4) = CRIT4
      CRIT(5) = CRIT5
      GO TO 9000
C                                  TERMINAL ERROR N1 .LT. 2
  140 IER = 129
      GO TO 9000
  145 IF (ISW.EQ.0) GO TO 9005
      CRIT(3) = CRIT3
      CRIT(4) = CRIT4
      CRIT(5) = CRIT5
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BEPAT ')
 9005 RETURN
      END

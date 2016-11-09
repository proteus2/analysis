C   IMSL ROUTINE NAME   - RLMUL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MULTIPLE LINEAR REGRESSION ANALYSIS
C
C   USAGE               - CALL RLMUL (A,XYBAR,N,M,ALFA,ANOVA,B,IB,VARB,
C                           IER)
C
C   ARGUMENTS    A      - ON INPUT, A IS AN (M+1) BY (M+1) SYMMETRIC
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE. IT
C                           CONTAINS THE CORRECTED SUMS OF SQUARES AND
C                           CROSS PRODUCTS OF THE INDEPENDENT AND
C                           DEPENDENT VARIABLES. A VECTOR OF LENGTH
C                           (M+1)*(M+2)/2 IS REQUIRED. IMSL SUBROUTINE
C                           BECOVM MAY BE USED TO GENERATE THE INPUT
C                           MATRIX A.
C                         ON OUTPUT, A IS DESTROYED.
C                XYBAR  - INPUT VECTOR OF LENGTH M+1 CONTAINING THE
C                           VARIABLE MEANS WITH THE RESPONSE
C                           (DEPENDENT) VARIABLE MEAN LAST.
C                N      - INPUT NUMBER OF DATA POINTS. N MUST BE
C                           GREATER THAN M+1.
C                M      - INPUT NUMBER OF INDEPENDENT VARIABLES. M MUST
C                           BE POSITIVE AND LESS THAN N-1.
C                ALFA   - INPUT RISK LEVEL IN THE INTERVAL (0,1) FOR
C                           THE REGRESSION COEFFICIENT INTERVAL
C                           ESTIMATES. THE CHOICE 0.05 IS A COMMON ONE.
C                ANOVA  - OUTPUT VECTOR OF LENGTH 14 CONTAINING ANALYSIS
C                           OF VARIANCE TABLE ENTRIES AND OTHER
C                           INFORMATION.
C                         ANOVA(1), ANOVA(2), AND ANOVA(3) CONTAIN
C                           THE DEGREES OF FREEDOM FOR THE
C                           REGRESSION, RESIDUAL, AND CORRECTED TOTAL
C                           SOURCES OF VARIATION, RESPECTIVELY.
C                         ANOVA(4), ANOVA(5), AND ANOVA(6) CONTAIN
C                           THE SUMS OF SQUARES CORRESPONDING
C                           TO THE DEGREES OF FREEDOM IN LOCATIONS 1,
C                           2, AND 3.
C                         ANOVA(7) AND ANOVA(8) CONTAIN THE MEAN
C                           SQUARES FOR THE REGRESSION
C                           AND RESIDUAL SOURCES OF VARIATION,
C                           RESPECTIVELY.
C                         ANOVA(9) AND ANOVA(10) CONTAIN THE COMPUTED
C                           F-VALUE AND CORRESPONDING TAIL AREA OF THE
C                           F DISTRIBUTION, RESPECTIVELY, FOR TESTING
C                           THE NULL HYPOTHESIS THAT ALL REGRESSION
C                           COEFFICIENT PARAMETERS ARE ZERO.
C                         ANOVA(11) CONTAINS THE PERCENTAGE OF
C                           VARIATION EXPLAINED BY THE ESTIMATED MODEL.
C                         ANOVA(12) CONTAINS THE STANDARD DEVIATION OF
C                           RESIDUALS
C                         ANOVA(13) CONTAINS THE RESIDUAL STANDARD
C                           DEVIATION EXPRESSED AS A PERCENTAGE OF THE
C                           RESPONSE MEAN.
C                         ANOVA(14) CONTAINS THE NUMBER OF DECIMAL
C                           DIGITS OF ACCURACY IN THE REGRESSION
C                           COEFFICIENT ESTIMATES, DEFINED WHEN WARNING
C                           ERROR IER=33 OCCURS. OTHERWISE ANOVA(14) IS
C                           SET TO 4.0.
C                B      - OUTPUT (M+1) BY 7 MATRIX CONTAINING INFERENCE
C                           RESULTS FOR EACH REGRESSION COEFFICIENT AND
C                           THE INTERCEPT. EACH ROW CORRESPONDS TO A
C                           COEFFICIENT WITH INTERCEPT RESULTS APPEARING
C                           IN THE LAST ROW.
C                         B(I,1) FOR I=1,...,M+1 CONTAIN THE PARAMETER
C                           ESTIMATES.
C                         B(I,2) FOR I=1,...,M+1 CONTAIN THE LOWER
C                           BOUNDS FOR INTERVAL ESTIMATES
C                         B(I,3) FOR I=1,...,M+1 CONTAIN THE UPPER
C                           BOUNDS FOR INTERVAL ESTIMATES
C                         B(I,4) FOR I=1,...,M+1 CONTAIN THE STANDARD
C                           ERRORS OF PARAMETER ESTIMATES.
C                         B(I,5) FOR I=1,...,M CONTAIN THE ADJUSTED
C                           SUMS OF SQUARES
C                         B(I,6) FOR I=1,...,M CONTAIN THE PARTIAL
C                           F-VALUES.
C                         B(I,7) FOR I=1,...,M CONTAIN THE F
C                           DISTRIBUTION TAIL AREAS CORRESPONDING TO
C                           THE PARTIAL F-VALUES.
C                         THE LAST THREE ELEMENTS OF ROW M+1 ARE WORK
C                           SPACE.
C                IB     - INPUT ROW DIMENSION OF THE MATRIX B EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                VARB   - OUTPUT VECTOR OF LENGTH M*(M+1)/2 CONTAINING
C                           THE M BY M SYMMETRIC INVERSE OF THE
C                           INFORMATION MATRIX (PART OF INPUT A).
C                           VARB*ANOVA(8) IS THE ESTIMATED VARIANCE-
C                           COVARIANCE MATRIX FOR THE REGRESSION
C                           COEFFICIENT ESTIMATES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT FEWER THAN THREE
C                             SIGNIFICANT DIGITS OF ACCURACY WERE
C                             RETAINABLE IN COMPUTING THE COEFFICIENT
C                             ESTIMATES. THE ACTUAL NUMBER OF DIGITS IS
C                             RETURNED IN ANOVA(14).
C                           IER=34 INDICATES THAT THE LAST ELEMENT OF
C                             INPUT A WAS ZERO (I.E. THE RESPONSE
C                             VARIABLE WAS CONSTANT) OR A PERFECT FIT
C                             TO THE DATA WAS OBTAINED. IN SUCH A CASE,
C                             OUTPUT ELEMENTS 7 THROUGH 13 OF ANOVA AND
C                             ALL COLUMNS OF B EXCEPTING 1 AND 5 HAVE
C                             BEEN SET TO ZERO. OTHER OUTPUT INFORMATION
C                             IS CORRECT.
C                         TERMINAL ERROR
C                           IER=131 INDICATES AT LEAST ONE OF THE
C                             PARAMETERS N,M, OR ALFA WAS SPECIFIED
C                             INCORRECTLY.
C                           IER=132 INDICATES THE PROBLEM WAS TOO ILL-
C                             CONDITIONED TO BE SOLVED CORRECTLY BY
C                             RLMUL. THE USER MAY WISH TO CONSIDER
C                             EXECUTION OF RLMUL IN HIGHER PRECISION,
C                             REVIEW OF THE MODEL FOR DELETION AND/OR
C                             AUGMENTATION IN THE INDEPENDENT VARIABLE
C                             SET, OR SCALING OF VARIABLES TO SIMILAR
C                             RANGES.
C                           IER=133 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL SUBROUTINE MDBETA OR MDBETI. USAGE OF
C                             RLMUL IS NOT APPROPRIATE IN THIS CASE.
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/LUELMP,MDBETA,MDBETI,
C                           MLGAMD=DLGAMA,UERTST,UGETIO,VMULFS
C                       - SINGLE(H36,H48,H60)/LUELMP,MDBETA,MDBETI,
C                           MLGAMA=ALGAMA,UERTST,UGETIO,VMULFS
C                       - DOUBLE/LUELMP,MDBETA,MDBETI,MLGAMD=DLGAMA,
C                           UERTST,UGETIO,VMULFS,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IN A USER APPLICATION PROGRAM CALLING RLMUL, IT
C                MAY BE DESIRABLE, PRIOR TO CALLING RLMUL, TO UTILIZE
C                IMSL ROUTINES RLFITI (LACK OF FIT), BDTRGI
C                (TRANSFORMATION), BECORI (CORRELATION), BECOVM
C                (CROSS PRODUCTS (INFORMATION) MATRIX), AND RLSUM
C                (SELECTION OF A SUBSET PROBLEM FOR ANALYSIS).
C                SUBSEQUENT TO THE RLMUL CALL, IMSL ROUTINES
C                RLRES (RESIDUAL ANALYSIS), RLPRDI (PREDICTION
C                (EVALUATION) ANALYSIS), AND USPLT (PLOTTING) (USPLTD
C                FOR DOUBLE PRECISION) MAY BE UTILIZED.
C            2.  ALTHOUGH M MAY BE AS LARGE AS N-2, TYPICAL USAGE
C                WOULD REQUIRE N GREATER THAN M+10.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLMUL (A,XYBAR,N,M,ALFA,ANOVA,B,IB,VARB,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IB,IER
      REAL               A(1),XYBAR(1),ALFA,ANOVA(14),B(IB,7),VARB(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IP1,IP,IQ,IR,I,JER,J,KK,K,LM,L,M1,M2,MEND,MM1,
     1                   MND,MROW,NKK1,NM1
      REAL               ALFA1,DF1,DF2,F,P
      REAL               ABIG,ABSX,EEPS,EPS4,EPS,FOUR,HUND,ONE,
     1                   P5,REPS,RN,RPEXE,SIXTN,TEMP,X,ZERO
      DOUBLE PRECISION   SUM
      DATA               REPS/Z3C100000/
      DATA               RPEXE/174.6731/
      DATA               EPS4/1.E+3/,ZERO/0./,HUND/100.0/,P5/.5/
      DATA               ONE/1.0/,SIXTN/16.0/,FOUR/4.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      NM1 = N-1
      IF (M.GE.NM1 .OR. N.LT.2 .OR. ALFA.LE.ZERO .OR. ALFA.GE.ONE) GO
     *TO 100
      M1 = M+1
      M2 = M+2
      MND = (M*M1)/2
      MROW = MND+1
      RN = ONE/(M1*SIXTN)
      MM1 = M-1
      MEND = (M1*M2)/2
      ANOVA(14) = FOUR
      EEPS = REPS*EPS4
      EPS = REPS
      IP = 1
C                                  DECOMPOSITION OF MATRIX A INTO L*
C                                  L-TRANSPOSE
      DO 35 I=1,M
         IQ = IP
         IR = 1
         DO 30 J=1,I
            X = A(IP)
            IF (J.EQ.1) GO TO 10
            DO 5 K=IQ,IP1
               X = X-A(K)*A(IR)
               IR = IR+1
    5       CONTINUE
   10       IF (I.NE.J) GO TO 20
            TEMP = X*RN
            IF (A(IP)+TEMP.LE.A(IP)) GO TO 95
            TEMP = EEPS*A(IP)
            IF (X.GT.TEMP) GO TO 15
C                                  WARNING - FEWER THAN 3 SIGNIFICANT
C                                  DIGITS OF ACCURACY WERE RETAINABLE
            IER = 33
            ANOVA(14) = AMIN1(ANOVA(14),-ALOG10(REPS*A(IP)/X))
   15       A(IP) = ONE/SQRT(X)
            GO TO 25
   20       A(IP) = X*A(IR)
   25       IP1 = IP
            IP = IP+1
            IR = IR+1
   30    CONTINUE
   35 CONTINUE
C                                  WARNING - PERFECT FIT
      IF (A(MEND).LE.ZERO) IER = 34
C                                  COMPUTE PARAMETER ESTIMATES
      CALL LUELMP(A,A(MROW),M,B(1,1))
      SUM = XYBAR(M1)
      KK = M
      DO 40 I=1,M
         SUM = SUM-DBLE(XYBAR(I))*DBLE(B(I,1))
   40 CONTINUE
      B(M1,1) = SUM
C                                  STORE DEGREES OF FREEDOM AND SUMS OF
C                                  SQUARES
      ANOVA(1) = KK
      NKK1 = NM1-KK
      ANOVA(2) = NKK1
      ANOVA(3) = NM1
      SUM = 0.D0
      DO 45 I=1,M
         SUM = SUM+DBLE(B(I,1))*DBLE(A(MND+I))
   45 CONTINUE
      ANOVA(4) = SUM
      ANOVA(6) = A(MEND)
      IF (ANOVA(6).LT.ZERO) ANOVA(6) = ZERO
      ANOVA(5) = ANOVA(6)-ANOVA(4)
      IF (ANOVA(5).LE.EPS*HUND*A(MEND)) IER = 34
      IF (KK.NE.0 .AND. IER.NE.34) GO TO 55
      ANOVA(5) = ZERO
      DO 50 I=7,13
         ANOVA(I) = ZERO
   50 CONTINUE
      GO TO 60
C                                  COMPUTE MEAN SQUARES
   55 ANOVA(7) = ANOVA(4)/ANOVA(1)
      ANOVA(8) = ANOVA(5)/ANOVA(2)
C                                  COMPUTE F VALUES
      ANOVA(9) = ANOVA(7)/ANOVA(8)
      F = ONE-NKK1/(NKK1+KK*ANOVA(9))
      DF1 = KK*P5
      DF2 = NKK1*P5
      CALL MDBETA(F,DF1,DF2,P,JER)
      IF (JER.NE.0) GO TO 105
      ANOVA(10) = ONE-P
      ANOVA(11) = HUND*ANOVA(4)/ANOVA(6)
      ANOVA(12) = SQRT(ANOVA(8))
      ANOVA(13) = ZERO
      ABSX = ABS(XYBAR(M1))
      IF (ABSX.EQ.ZERO) GO TO 60
      ABSX = ALOG(ANOVA(12))-ALOG(ABSX)
      IF (ABSX.GT.RPEXE) GO TO 60
      ANOVA(13) = HUND*ANOVA(12)/XYBAR(M1)
   60 K = 0
      DF1 = P5
      DF2 = ANOVA(2)*P5
      ALFA1 = ONE-ALFA
      CALL MDBETI(ALFA1,DF1,DF2,F,JER)
      IF (JER.NE.0) GO TO 105
      X = NKK1*F/(ONE-F)
      X = SQRT(X)
      DO 70 I=1,M1
         DO 65 J=2,7
            B(I,J) = ZERO
   65    CONTINUE
   70 CONTINUE
C                                  COMPUTE THE INVERSE OF THE
C                                  INFORMATION MATRIX
      L = 1
      LM = M
      IF (IER.NE.34) TEMP = ONE/ANOVA(8)
      DF2 = NKK1*P5
      DO 80 I=1,M
         DO 75 J=L,LM
            VARB(J) = ZERO
   75    CONTINUE
         VARB(L+I-1) = ONE
         CALL LUELMP(A,VARB(L),M,VARB(L))
         L = L+I
         LM = L+MM1
C                                  COMPUTE STANDARD ERRORS OF PARAMETER
C                                  ESTIMATES
         K = K+I
         SUM = DBLE(B(I,1))*DBLE(B(I,1))
         B(I,5) = SUM/VARB(K)
         IF (IER.EQ.34) GO TO 80
C                                  COMPUTE UPPER AND LOWER BOUNDS FOR
C                                  INTERVAL ESTIMATES
         B(I,4) = SQRT(ANOVA(8)*VARB(K))
         ABIG = X*B(I,4)
         B(I,2) = B(I,1)-ABIG
         B(I,3) = B(I,1)+ABIG
         B(I,6) = SUM*TEMP/VARB(K)
         F = ONE-NKK1/(NKK1+B(I,6))
         CALL MDBETA(F,DF1,DF2,P,JER)
         IF (JER.NE.0) GO TO 105
         B(I,7) = ONE-P
   80 CONTINUE
      IF (IER.EQ.34) GO TO 90
      CALL VMULFS(XYBAR,VARB,1,M,1,A,1)
      SUM = 0.D0
      DO 85 I=1,M
         SUM = SUM+DBLE(A(I))*DBLE(XYBAR(I))
   85 CONTINUE
      ABIG = SUM
      B(M1,4) = SQRT(ANOVA(8)/N+ANOVA(8)*ABIG)
      ABIG = B(M1,4)*X
      B(M1,2) = B(M1,1)-ABIG
      B(M1,3) = B(M1,1)+ABIG
   90 IER = MAX0(IER,JER)
      IF (IER.EQ.0) GO TO 9005
      GO TO 9000
C                                  PROBLEM IS TOO ILL-CONDITIONED
   95 IER = 132
      GO TO 9000
C                                  INPUT PARAMETERS SPECIFIED
C                                  INCORRECTLY
  100 IER = 131
      GO TO 9000
C                                  AN ERROR OCCURRED IN IMSL ROUTINE
C                                  MDFD OR MDBETA
  105 IER = 133
 9000 CONTINUE
      CALL UERTST(IER,6HRLMUL )
 9005 RETURN
      END

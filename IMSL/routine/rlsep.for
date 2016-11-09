C   IMSL ROUTINE NAME   - RLSEP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SELECTION OF A REGRESSION MODEL USING A
C                           FORWARD STEPWISE ALGORITHM, AND COMPUTATION
C                           OF THE USUAL ANALYSIS OF VARIANCE TABLE
C                           ENTRIES - EASY TO USE VERSION
C
C   USAGE               - CALL RLSEP (XY,N,M,IX,ALFA,IJOB,IND,ANOVA,
C                           XYB,IB,VARB,IER)
C
C   ARGUMENTS    XY     - INPUT N BY M+1 DATA MATRIX CONTAINING THE
C                           SETTINGS FOR THE INDEPENDENT VARIABLE IN
C                           THE FIRST M COLUMNS, AND THE CORRESPONDING
C                           RESPONSES IN COLUMN M+1.
C                         ON OUTPUT, XY IS DESTROYED.
C                N      - INPUT NUMBER OF DATA POINTS.
C                M      - INPUT NUMBER OF INDEPENDENT VARIABLES UNDER
C                           CONSIDERATION.
C                IX     - INPUT ROW DIMENSION OF THE MATRIX XY EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                ALFA   - INPUT VECTOR OF LENGTH 2 CONTAINING
C                           SIGNIFICANCE LEVELS FOR ENTERING AND
C                           DELETING VARIABLES.
C                         ALFA(1) CONTAINS THE LEVEL FOR ENTRANCE, IN
C                           THE INCLUSIVE RANGE (0,1). THE CHOICE 0.05
C                           IS A COMMON ONE.
C                         ALFA(2) CONTAINS THE LEVEL FOR DELETION, IN
C                           THE INCLUSIVE RANGE (0,1). ALFA(2) MUST NOT
C                           BE LESS THAN ALFA(1).
C                IJOB   - INPUT OPTION VECTOR OF LENGTH 2.
C                         IJOB(1) CONTAINS THE LACK OF FIT OPTION.
C                           WHEN IJOB(1) IS ZERO, A LACK OF FIT TEST
C                           BASED ON REPLICATE OBSERVATIONS IS NOT
C                           PERFORMED.
C                           WHEN IJOB(1) IS NONZERO, A LACK OF FIT
C                           TEST WILL BE PERFORMED USING ONLY THOSE
C                           REPLICATE OBSERVATIONS THAT APPEAR
C                           CONTIGUOUSLY IN MATRIX XY.
C                         IJOB(2) CONTAINS THE PARTIAL F-TEST OPTION.
C                           WHEN IJOB(2) IS ZERO, ONLY THE OVERALL
C                           F-TEST FOR TERMS IN THE MODEL IS PERFORMED.
C                           WHEN IJOB(2) IS NON-ZERO, A PARTIAL
C                           F-TEST FOR EACH TERM IN THE MODEL IS
C                           ALSO PERFORMED.
C                IND    - INPUT/OUTPUT VECTOR OF LENGTH 2*M+1.
C                         ON INPUT, THE FIRST M LOCATIONS MUST BE ZERO
C                           EXCEPT THAT IND(I) MUST BE NONZERO IF THE
C                           I-TH INDEPENDENT VARIABLE IS TO BE FORCED
C                           INTO THE MODEL, I=1,2,...,M. THE REMAINING
C                           M+1 LOCATIONS ARE WORK SPACE.
C                         ON OUTPUT, THE SECOND M LOCATIONS INDICATE
C                           WHICH TERMS ARE IN THE MODEL.  THE I-TH
C                           INDEPENDENT VARIABLE IS OR IS NOT IN THE
C                           MODEL AS IND(M+I) HAS THE VALUE 1 OR 0,
C                           RESPECTIVELY, I=1,2,...,M. THE LAST
C                           LOCATION IS USED AS WORK SPACE.
C                ANOVA  - OUTPUT VECTOR OF LENGTH THE MAXIMUM OF
C                           (16,M+1). ONLY THE FIRST 6 ELEMENTS OF ANOVA
C                           ARE DEFINED IF NO INDEPENDENT VARIABLES
C                           ENTER THE MODEL. FOR THE SELECTED MODEL,
C                         ANOVA(1), ANOVA(2), AND ANOVA(3) CONTAIN THE
C                           DEGREES OF FREEDOM FOR THE REGRESSION,
C                           RESIDUAL, AND CORRECTED TOTAL SOURCES OF
C                           VARIATION, RESPECTIVELY.
C                         ANOVA(4), ANOVA(5), AND ANOVA(6) CONTAIN THE
C                           THE SUMS OF SQUARES CORRESPONDING
C                           TO THE DEGREES OF FREEDOM IN LOCATIONS
C                           1, 2, AND 3.
C                         ANOVA(7) AND ANOVA(8) CONTAIN THE MEAN SQUARES
C                           FOR THE REGRESSION AND RESIDUAL SOURCES OF
C                           VARIATION, RESPECTIVELY.
C                         ANOVA(9) AND ANOVA(10) CONTAIN THE COMPUTED F
C                           VALUE AND CORRESPONDING TAIL AREA OF THE F
C                           DISTRIBUTION, RESPECTIVELY, FOR TESTING
C                           THE NULL HYPOTHESIS THAT ALL REGRESSION
C                           COEFFICIENT PARAMETERS ARE ZERO.
C                         ANOVA(11) CONTAINS THE PERCENTAGE OF VARIATION
C                           EXPLAINED.
C                         ANOVA(12) CONTAINS THE STANDARD DEVIATION OF
C                           THE RESIDUALS.
C                         ANOVA(13) CONTAINS THE RESIDUAL STANDARD
C                           DEVIATION EXPRESSED AS A PERCENTAGE OF THE
C                           RESPONSE MEAN.
C                         ANOVA(14), ANOVA(15), AND ANOVA(16) CONTAIN
C                           THE SUM OF SQUARES, COMPUTED
C                           F VALUE, AND CORRESPONDING TAIL AREA OF
C                           THE F DISTRIBUTION, RESPECTIVELY, FOR
C                           LACK OF FIT. THESE 3 LOCATIONS ARE
C                           DEFINED ONLY IF IJOB(1) IS NONZERO AND
C                           IF THE WARNING ERROR(IER = 33) DOES NOT
C                           OCCUR.
C                         NOTE THAT IF M IS GREATER THAN 15, THE
C                           REMAINING LOCATIONS OF ANOVA ARE USED AS
C                           WORK STORAGE.
C                XYB    - OUTPUT M+1 BY 5 MATRIX CONTAINING SUMMARY
C                           STATISTICS FOR THE INPUT DATA AND THE
C                           SELECTED MODEL.
C                         COLUMN 1 CONTAINS THE INDEPENDENT VARIABLE
C                           MEANS FOLLOWED BY THE RESPONSE MEAN.
C                         COLUMN 2 CONTAINS THE REGRESSION COEFFI-
C                           CIENTS. FOR TERMS NOT IN THE MODEL,
C                           CORRESPONDING COEFFICIENTS ARE SET TO
C                           ZERO. THE LAST ELEMENT IS THE
C                           INTERCEPT ESTIMATE.
C                         COLUMNS 3,4,5 ARE DEFINED ONLY IF IJOB(2)
C                           IS NONZERO, AND THEN CONTAIN THE
C                           ADJUSTED SUM OF SQUARES, F-VALUE, AND
C                           CORRESPONDING TAIL AREA OF THE F-
C                           DISTRIBUTION, RESPECTIVELY, FOR EACH
C                           TERM IN THE MODEL. ELEMENTS
C                           CORRESPONDING TO TERMS NOT IN THE MODEL
C                           ARE SET TO ZERO. THE LAST ELEMENT OF
C                           EACH COLUMN IS UNDEFINED.
C                IB     - INPUT ROW DIMENSION OF THE MATRIX XYB EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                VARB   - OUTPUT VECTOR OF LENGTH (M+1)*(M+2)/2
C                           CONTAINING THE INVERSE OF THE CORRECTED
C                           SUMS OF SQUARES AND CROSS PRODUCTS
C                           MATRIX CORRESPONDING TO VARIABLES IN THE
C                           MODEL. ASSUMING Q (LESS THAN OR EQUAL TO
C                           M) VARIABLES ARE IN THE MODEL, THE Q BY Q
C                           SYMMETRIC INVERSE WILL APPEAR IN SYMMETRIC
C                           STORAGE MODE, IN THE FIRST Q*(Q+1)/2
C                           LOCATIONS OF VARB. THE REMAINING LOCATIONS
C                           ARE WORK STORAGE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES A WARNING ERROR OCCURRED IN
C                             IMSL ROUTINE RLFITI WHEN THE LACK OF FIT
C                             TEST WAS REQUESTED (IJOB(1)).
C                           IER=36 INDICATES THAT THE LAST COLUMN OF
C                             INPUT XY WAS CONSTANT (I.E. THE RESPONSE
C                             VARIABLE WAS CONSTANT) OR A PERFECT FIT
C                             TO THE DATA WAS OBTAINED, POSSIBLY
C                             OCCURRING PRIOR TO NORMAL (AS EXPECTED BY
C                             THE USER) TERMINATION OF RLSEP. OUTPUT
C                             ELEMENTS 7 THROUGH 16 OF ANOVA AND COLUMNS
C                             4 AND 5 OF XYB HAVE BEEN SET TO ZERO.
C                         TERMINAL ERROR
C                           IER=130 INDICATES A TERMINAL ERROR OCCURRED
C                             IN IMSL ROUTINE BECOVM.
C                           IER=131 INDICATES A TERMINAL ERROR OCCURRED
C                             IN IMSL ROUTINE RLSTP.
C
C   REQD. IMSL ROUTINES - SINGLE/BECOVM,MDFD,MERRC=ERFC,RLFITI,RLSTP,
C                           RLSUBM,UERTST,UGETIO
C                       - DOUBLE/BECOVM,MDFD,MERRC=ERFC,RLFITI,RLSTP,
C                           RLSUBM,UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  PRIOR TO CALLING RLSEP, THE USER CAN TRANSFORM
C                THE DATA BY USING IMSL ROUTINE BDTRGI, COMPUTE A
C                CORRELATION MATRIX BY CALLING IMSL ROUTINE BECORI,
C                AND/OR GENERATE THE VARIABLE SETTINGS FOR A FULL
C                QUADRATIC MODEL IN THE INDEPENDENT VARIABLES
C                BY CALLING IMSL ROUTINE RLGQMI.
C            2.  AFTER CALLING RLSEP, THE USER CAN COMPUTE INTERVAL
C                ESTIMATES FOR THE RESPONSE BY CALLING IMSL ROUTINES
C                BECVLI AND RLPRDI, AND/OR IF RLGQMI WAS USED PRIOR
C                TO CALLING RLSEP, DECODE THE FITTED MODEL BY
C                CALLING IMSL ROUTINE RLDCQM.
C            3.  THE PARTIAL F-TEST CALCULATIONS CONTROLLED BY
C                IJOB(2), TEST THE NULL HYPOTHESIS THAT EACH
C                REGRESSION COEFFICIENT PARAMETER IS ZERO, GIVEN
C                THAT ALL OTHER SELECTED TERMS ARE IN THE MODEL.
C            4.  THE USER MAY USE THIS ROUTINE FOR AN ORDINARY
C                MULTIPLE REGRESSION ANALYSIS BY FORCING ALL M
C                INDEPENDENT VARIABLES INTO THE MODEL.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSEP  (XY,N,M,IX,ALFA,IJOB,IND,ANOVA,XYB,IB,VARB,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IX,IJOB(2),IND(1),IB,IER
      REAL               XY(IX,1),ALFA(2),ANOVA(1),XYB(IB,5),VARB(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NBR(6),I,IOPT,J,JER,K,KK,L,MEND,M1,NDF,NKK1,
     1                   NKK1DF,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,
     2                   I13,I14,I15,I16
      REAL               F,P
      REAL               SS,ONE,HUND,ZERO,TEMP,SST,SS1(1)
      DOUBLE PRECISION   SUM
      DATA               ONE/1.0/,HUND/100./,ZERO/0./
      DATA               NBR(4)/1/,NBR(5)/1/,NBR(6)/1/,I1/1/,I2/2/,
     1                   I3/3/,I4/4/,I5/5/,I6/6/,I7/7/,I8/8/,I9/9/,
     2                   I10/10/,I11/11/,I12/12/,I13/13/,I14/14/,
     3                   I15/15/,I16/16/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IOPT = 0
      JER = 0
      M1 = M+1
      MEND = (M1*(M+2))/2
      NBR(1) = M1
      NBR(2) = N
      NBR(3) = N
      IF (IJOB(1) .NE. 0) CALL RLFITI(XY,N,M,IX,XY(1,M1),1,IX,SS1,NDF,
     *   JER)
      IF (IJOB(1) .NE. 0) SS = SS1(1)
      CALL BECOVM(XY,IX,NBR,ANOVA,XYB(1,1),VARB,IER)
      IF (IER .NE. 0) GO TO 80
      K = 0
      DO 5 I = 1,M
         IND(M1+I) = 0
         IF (IND(I) .EQ. 0) GO TO 5
         IND(I) = 1
         K = K+1
    5 CONTINUE
      IND(M1) = 0
      SST = VARB(MEND)
      CALL RLSTP (VARB,M,N,ALFA(1),ALFA(2),IND,IND(M1),XYB(1,2),IOPT,
     1   IER)
   15 IF (IER .GT. 128) GO TO 85
      SUM = XYB(M1,1)
      KK = 0
      DO 25 I = 1,M
         IF (IND(M+I) .EQ. -1) GO TO 20
         SUM = SUM-DBLE(XYB(I,1))*DBLE(XYB(I,2))
         KK = KK+1
         GO TO 25
   20    IND(M+I) = 0
   25 CONTINUE
      XYB(M1,2) = SUM
      NKK1 = N-KK-1
      ANOVA(I1) = KK
      ANOVA(I2) = NKK1
      ANOVA(I3) = N-1
      ANOVA(I5) = VARB(MEND)
      IF (ANOVA(I5) .LT. ZERO) ANOVA(I5) = ZERO
      ANOVA(I6) = SST
      IF (ANOVA(I6) .LT. ZERO) ANOVA(I6) = ZERO
      ANOVA(I4) = ANOVA(I6)-ANOVA(I5)
      IF (KK .NE. 0 .AND. IER .NE. 37) GO TO 35
      DO 30 I = 7,13
         ANOVA(I) = ZERO
   30 CONTINUE
      GO TO 40
   35 ANOVA(I7) = ANOVA(I4)/ANOVA(I1)
      ANOVA(I8) = ANOVA(I5)/ANOVA(I2)
      ANOVA(I9) = ANOVA(I7)/ANOVA(I8)
      F = ANOVA(I9)
      CALL MDFD(F,KK,NKK1,P,IER)
      ANOVA(I10) = ONE-P
      ANOVA(I12) =  SQRT(ANOVA(I8))
      ANOVA(I13) = HUND*ANOVA(I12)/XYB(M1,1)
      ANOVA(I11) = HUND*ANOVA(I4)/SST
      IF (IJOB(1) .NE. 0 .AND. JER.EQ. 0) GO TO 45
   40 ANOVA(I14) = ZERO
      ANOVA(I15) = ZERO
      ANOVA(I16) = ZERO
      GO TO 50
   45 ANOVA(I14) = ANOVA(I5)-SS
      NKK1DF = NKK1-NDF
      ANOVA(I15) = ANOVA(I14)*NDF/(SS*NKK1DF)
      F = ANOVA(I15)
      CALL MDFD(F,NKK1DF,NDF,P,IER)
      ANOVA(I16) = ONE-P
   50 DO 60 I = 1,M1
         DO 55 J = 3,5
            XYB(I,J) = ZERO
   55    CONTINUE
   60 CONTINUE
      IF (KK .EQ. 0) GO TO 75
      IF (IJOB(2) .EQ. 0) GO TO 70
      L = 0
      IF (IER .NE. 37) TEMP = ONE/ANOVA(I8)
      DO 65 I = 1,M
         L = L+I
         IF (IND(M+I) .EQ. 0) GO TO 65
         SUM = DBLE(XYB(I,2))*DBLE(XYB(I,2))/VARB(L)
         XYB(I,3) = SUM
         IF (IER .EQ. 37) GO TO 65
         XYB(I,4) = SUM*TEMP
         F = XYB(I,4)
         CALL MDFD(F,1,NKK1,P,IER)
         XYB(I,5) = ONE-P
   65 CONTINUE
   70 CALL RLSUBM(VARB,M1,IND(M1),VARB,KK)
   75 IF (JER .NE. 0) IER = 33
      IF (IER .EQ. 37) IER = 36
      IF (IER .EQ. 0) GO TO 9005
      GO TO 9000
   80 IER = 130
      GO TO 9000
   85 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,6HRLSEP )
 9005 RETURN
      END

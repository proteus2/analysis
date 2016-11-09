C   IMSL ROUTINE NAME   - SSRAND
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SIMPLE RANDOM SAMPLING WITH CONTINUOUS
C                           DATA - INFERENCES REGARDING THE POPULATION
C                           MEAN AND TOTAL USING RATIO OR REGRESSION
C                           ESTIMATION
C
C   USAGE               - CALL SSRAND (Y,IY,IOPT,NBR,ALPHA,TEMP,XBAR,B,
C                           STAT,IER)
C
C   ARGUMENTS    Y      - INPUT NBR(2) BY 2 SUBMATRIX OF THE MATRIX
C                           (CALL IT YY) CONTAINING THE ENTIRE RANDOM
C                           SAMPLE.  THE AUXILIARY VARIABLE SETTINGS ARE
C                           IN COLUMN ONE WITH CORRESPONDING VARIABLE
C                           OF INTEREST SETTINGS IN COLUMN TWO.  THE
C                           LAST SUBMATRIX IN YY MAY HAVE FEWER THAN
C                           NBR(2) ROWS.  SEE EXAMPLE.
C                IY     - INPUT ROW DIMENSION OF MATRIX Y EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IOPT   - INPUT ESTIMATION OPTION.
C                         IF IOPT IS NEGATIVE, ONLY INFERENCE ABOUT THE
C                           POPULATION RATIO IS DESIRED.
C                         IF IOPT = 0, RATIO ESTIMATION IS TO BE USED
C                           FOR INFERENCE ABOUT THE POPULATION MEAN,
C                           TOTAL, AND RATIO.
C                         IF IOPT IS POSITIVE, REGRESSION ESTIMATION IS
C                           TO BE USED FOR INFERENCE ABOUT THE POPULA-
C                           TION MEAN AND TOTAL.
C                NBR    - INPUT VECTOR OF LENGTH 7.  NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF PAIRS OF OBSERVATIONS IN YY.
C                           I=2, NUMBER OF PAIRS OF OBSERVATIONS IN
C                             EACH SUBMATRIX Y, NOT INCLUDING THE LAST
C                             SUBMATRIX WHERE THE NUMBER MAY BE LESS
C                             THAN OR EQUAL TO NBR(2).  HOWEVER NBR(2)
C                             SHOULD BE THE SAME FOR ALL CALLS.
C                           I=3, THE NUMBER OF THE SUBMATRIX STORED IN
C                             Y. SEE REMARKS.
C                           I=4, THE TEMPORARY MEAN INDICATOR.  IF
C                             NBR(4) = 0, THE USER SUPPLIES TEMPORARY
C                             MEANS IN TEMP.  OTHERWISE, THE FIRST ROW
C                             OF YY (OR FIRST ROW OF Y WHEN NBR(3) = 1)
C                             IS UTILIZED.
C                           I=5, NUMBER OF PAIRS OF ELEMENTS IN THE
C                             SAMPLED POPULATION.
C                           I=6, OPTION FOR POPULATION MEAN OF THE AUX-
C                             ILIARY VARIATE, DEFINED ONLY WHEN IOPT IS
C                             NEGATIVE.  IF NBR(6) = 0, THIS MEAN IS
C                             KNOWN AND IS INPUT VIA XBAR.  FOR NON-ZERO
C                             NBR(6) THE MEAN IS UNKNOWN.
C                           I=7, OPTION FOR THE REGRESSION COEFFICIENT,
C                             DEFINED ONLY WHEN IOPT IS POSITIVE.  IF
C                             NBR(7) = 0, THE REGRESSION COEFFICIENT IS
C                             PREASSIGNED VIA B.  FOR NON-ZERO NBR(7)
C                             THIS PARAMETER IS ESTIMATED FROM THE DATA.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT CON-
C                           FIDENCE INTERVALS FOR THE MEAN, TOTAL, AND
C                           RATIO PARAMETERS.  THE VALUE 0.05 IS A
C                           COMMON CHOICE.
C                TEMP   - INPUT VECTOR OF LENGTH 2. IF NBR(4) = 0, TEMP
C                           MUST CONTAIN TEMPORARY MEANS FOR THE TWO
C                           COLUMNS OF YY, RESPECTIVELY. OTHERWISE TEMP
C                           IS UNDEFINED.
C                XBAR   - INPUT POPULATION MEAN OF THE AUXILIARY
C                           VARIATE. NOT REQUIRED WHEN IOPT IS NEGATIVE
C                           AND NBR(6) IS NON-ZERO.
C                B      - INPUT PREASSIGNED REGRESSION COEFFICIENT.
C                           REQUIRED ONLY WHEN IOPT IS POSITIVE AND
C                           NBR(7) = 0.
C                STAT   - OUTPUT VECTOR OF LENGTH 18. NOTE THAT STAT(1)
C                         THROUGH STAT(8) ARE NOT DEFINED WHEN IOPT IS
C                         NEGATIVE, AND THAT STAT(9) THROUGH STAT(12)
C                         ARE NOT DEFINED WHEN IOPT IS POSITIVE.
C                         STAT(I) CONTAINS, FOR
C                           I=1, ESTIMATE OF THE MEAN.
C                           I=2, ESTIMATE OF THE TOTAL.
C                           I=3, VARIANCE ESTIMATE OF THE MEAN ESTIMATE.
C                           I=4, VARIANCE ESTIMATE OF THE TOTAL ESTIMATE
C                           I=5, LOWER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=6, UPPER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=7, LOWER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=8, UPPER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=9, ESTIMATE OF THE RATIO.
C                           I=10,VARIANCE ESTIMATE OF THE RATIO ESTIMATE
C                           I=11,LOWER CONFIDENCE LIMIT FOR THE RATIO.
C                           I=12,UPPER CONFIDENCE LIMIT FOR THE RATIO.
C                           I=13, ESTIMATE (EXPRESSED AS A PERCENTAGE)
C                             OF THE COEFFICIENT OF VARIATION OF THE
C                             MEAN, TOTAL, AND RATIO ESTIMATES THAT
C                             ARE DEFINED, AS CONTROLLED BY IOPT.
C                           I=14,ESTIMATE (EXPRESSED AS A PERCENTAGE)
C                             OF THE COEFFICIENT OF VARIATION OF THE
C                             MEAN OF THE AUXILIARY VARIATE (COLUMN
C                             ONE OF YY).
C                           I=15,ESTIMATE (EXPRESSED AS A PERCENTAGE)
C                             OF THE COEFFICIENT OF VARIATION OF THE
C                             MEAN OF THE VARIABLE OF INTEREST (COLUMN
C                             TWO OF YY).
C                           I=16,AVERAGE OF COLUMN ONE OF YY.
C                           I=17,AVERAGE OF COLUMN TWO OF YY.
C                           I=18,ESTIMATE OF THE REGRESSION COEFFICIENT.
C                             DEFINED ONLY WHEN IOPT IS POSITIVE AND
C                             NBR(7) IS NON-ZERO.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33 INDICATES THAT EITHER THE SAMPLE
C                             SIZE (NBR(1)) DOES NOT EXCEED 30, OR THAT
C                             ONE OR BOTH OF THE COEFFICIENTS OF
C                             VARIATION (STAT(14) AND STAT(15)) FOR THE
C                             TWO COLUMNS OF YY EXCEED 10 PERCENT. THE
C                             RESULTS DEPEND ON THE ASSUMPTION OF
C                             APPROXIMATE NORMALITY OF THE PARAMETER
C                             ESTIMATES AND THAT ASSUMPTION MAY BE
C                             INVALID IN THIS CASE.
C                         TERMINAL ERROR
C                           IER = 130 INDICATES NBR(1) IS LESS THAN 3 OR
C                             THAT NBR(1) EXCEEDS NBR(5).
C                           IER = 131 INDICATES NBR(3) IS LESS THAN ONE
C                             OR THAT NBR(2)*(NBR(3)-1) EXCEEDS NBR(1)
C                             (I.E., THE NUMBER OF THE SUBMATRIX
C                             (NBR(3)) IMPLIES THE NUMBER OF OBSERVA-
C                             TIONS IN YY (NBR(1)) IS NOT SUFFICIENT).
C                           IER = 132 INDICATES THAT NBR(2) IS LESS THAN
C                             1 OR THAT ALPHA IS NOT IN THE EXCLUSIVE
C                             INTERVAL (0,1).
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      BETWEEN THE FIRST CALL AND THE LAST CALL (M-TH CALL)
C                TO SSRAND ONLY NBR(3) MAY BE MODIFIED AND IT SHOULD
C                FOLLOW THE PATTERN 1,2,...,M. THOUGH THIS PATTERN IS
C                THE OBVIOUS ONE TO FOLLOW, IT IS NOT NECESSARY IN ITS
C                ENTIRETY. FOR CALLS 2,3,...,M-1, NBR(3) MAY TAKE ANY
C                VALUE IN THE SET (2,3,...,M-1). ON THE FIRST CALL
C                NBR(3) MUST EQUAL 1, AND ON THE M-TH CALL NBR(3) MUST
C                EQUAL M.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSRAND (Y,IY,IOPT,NBR,ALPHA,TEMP,XBAR,B,STAT,IER)
C
      DIMENSION          Y(IY,1),NBR(1),TEMP(1),STAT(1),TA(10)
      DOUBLE PRECISION   SUM1,SUM2,SUM3,SUM4,SUM5
      EQUIVALENCE        (TA(1),SUM1),(TA(3),SUM2),(TA(5),SUM3),
     1                   (TA(7),SUM4),(TA(9),SUM5)
      DATA               ZERO,ONE,TWO,HUND/0.,1.,2.,100./
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  BEGIN FIRST CALL
      IF(NBR(3) - 1) 40,5,35
C                                  WARNING ERROR - NUMBER PAIRS OBSERVA-
C                                                  TIONS IS TOO SMALL
    5 IF(NBR(1) .LE. 30) IER = 33
C                                  TERMINAL ERROR - NUMBER PAIRS OBSER-
C                                                   VATIONS NOT IN RANGE
      IF(NBR(1) .GE. 3 .AND. NBR(1) .LE. NBR(5)) GO TO 10
      IER = 130
      GO TO 9000
C                                  TERMINAL ERROR - NUMBER PAIRS OBSERVA
C                                                   TIONS IN SUBMATRIX
C                                                   NEGATIVE OR ZERO
   10 IF(NBR(2) .LT. 1) GO TO 15
C                                  TERMINAL ERROR - ALPHA NOT IN RANGE
      IF(ALPHA .GT. 0.0 .AND. ALPHA .LT. 1.0) GO TO 20
   15 IER = 132
      GO TO 9000
C                                  INITIALIZE CONSTANTS
   20 SN = NBR(1)
      CN = NBR(5)
      F = ONE-SN/CN
      DO 25 I=1,18
         STAT(I) = ZERO
   25 CONTINUE
C                                  DEFINE TEMP(I),MSTR,NSTR
      IF (NBR(4) .EQ. 0) GO TO 30
      TEMP(1) = Y(1,1)
      TEMP(2) = Y(1,2)
C                                  BEGIN INTERMEDIATE AND LAST CALL
   30 STAT(11) = FLOAT((NBR(1)+NBR(2)-1)/NBR(2))
      STAT(12) = FLOAT(NBR(2))
   35 MSTR = STAT(11)
      NSTR = STAT(12)
      IF (NBR(3) - MSTR) 50,45,40
C                                  TERMINAL ERROR - DATA SET EXCEEDS
C                                                   SPECIFIED RANGE
   40 IER = 131
      GO TO 9000
   45 NSTR = NBR(1) - (MSTR-1)*NSTR
C                                  COMPUTE MEANS AND TEMPORARY CORRECTED
C                                  SUMS OF SQUARES AND CROSS PRODUCTS.
   50 DO 55 I = 1,10
         TA(I) = STAT(I)
   55 CONTINUE
      DO 60 I = 1,NSTR
         T1 = Y(I,1)
         T2 = Y(I,2)
         SUM4 = SUM4 + DBLE(T1)
         SUM5 = SUM5 + DBLE(T2)
         T1 = T1 - TEMP(1)
         T2 = T2 - TEMP(2)
         SUM1 = SUM1 + DBLE(T1)*DBLE(T1)
         SUM2 = SUM2 + DBLE(T2)*DBLE(T2)
         SUM3 = SUM3 + DBLE(T1)*DBLE(T2)
   60 CONTINUE
      IF (NBR(3) .EQ. MSTR) GO TO 70
      DO 65 I = 1,10
         STAT(I) = TA(I)
   65 CONTINUE
      GO TO 9000
C                                  END INTERMEDIATE CALL; FINISH LAST
C                                  CALL BY BEGINNING PRELIMINARY CALCU-
C                                  LATIONS OF STAT(I)
   70 STAT(3) = SUM1
      STAT(4) = SUM2
      STAT(5) = SUM3
      STAT(16) = SUM4/SN
      STAT(17) = SUM5/SN
      SY = STAT(17) - TEMP(2)
      SX = STAT(16) - TEMP(1)
      UX = SN*SX
      SXY = STAT(5) - UX*SY
      SX = STAT(3) - UX*SX
      SY = STAT(4) - SN*SY*SY
      UX = SUM5 - SN*TEMP(2)
      UY = STAT(4) + TEMP(2)*(SUM5 + UX)
      UXY = STAT(5) + TEMP(2)*SUM4 + TEMP(1)*UX
      UX = STAT(3) + TEMP(1)*(SUM4 + SUM4 - SN*TEMP(1))
      STAT(15) = ONE/(SN*(SN-ONE))
      RHAT = SUM5/SUM4
      ST = UY + RHAT*(RHAT*UX-UXY-UXY)
      SG = SY - SXY*SXY/SX
      P = ONE - ALPHA/TWO
      CALL MDNRIS(P,TA2,IER)
      IF (IOPT) 80,75,90
   75 STAT(1) = RHAT*XBAR
      STAT(2) = STAT(1)*CN
C                                  CALCULATE COMMON STAT(I) IF IOPT IS
C                                  NEGATIVE OR ZERO
C                                  DEFINE XBAR
   80 IF (NBR(6) .NE. 0 .AND. IOPT .LT. 0) XBAR = STAT(16)
      STAT(3) = (ST*F)*STAT(15)
      STAT(9) = RHAT
      STAT(10) = STAT(3)/(XBAR*XBAR)
      STAT(12) = SQRT(STAT(10))
      STAT(13) = HUND*STAT(12)/STAT(9)
      D = TA2*STAT(12)
      STAT(11) = STAT(9) - D
      STAT(12) = STAT(9) + D
      IF (IOPT .GE. 0) GO TO 100
      DO 85 I = 1,8
         STAT(I) = ZERO
   85 CONTINUE
      GO TO 110
C                                  CALCULATE COMMON STAT(I) IF IOPT IS
C                                  ZERO OR POSITIVE
   90 D = ZERO
      IF (NBR(7) .EQ.0) SGZ = SY + B * (B*SX-SXY-SXY)
C                                  DEFINE B
      IF (NBR(7) .EQ. 0) GO TO 95
      B = SXY/SX
      SGZ = SG
      D = ONE
   95 D = D + ONE
      STAT(1) = STAT(17) + B*(XBAR-STAT(16))
      STAT(2) = CN*STAT(1)
      STAT(3) = (SGZ*F)/(SN*(SN-D))
      STAT(13) = SQRT(STAT(3))*HUND/STAT(1)
      STAT(18) = B
  100 STAT(4) = CN*CN*STAT(3)
      STAT(6) = SQRT(STAT(3))*TA2
      STAT(5) = STAT(1) - STAT(6)
      STAT(6) = STAT(1) + STAT(6)
      STAT(8) = TA2*SQRT(STAT(4))
      STAT(7) = STAT(2) - STAT(8)
      STAT(8) = STAT(2) + STAT(8)
      IF(IOPT .LE. 0) GO TO 110
      DO 105 I = 9,12
         STAT(I) = ZERO
  105 CONTINUE
  110 STAT(14) = HUND * (SQRT(SX*STAT(15)))/XBAR
      STAT(15) = HUND * (SQRT(SY*STAT(15)))/STAT(17)
C                                  WARNING ERROR - ONE OR BOTH COEF-
C                                                  FICIENTS OF VARIATION
C                                                  EXCEED 10 PERCENT
      IF (STAT(14) .GT. 10.0 .OR. STAT(15) .GT. 10.0) IER = 33
 9000 CONTINUE
      IF(IER .EQ. 0) GO TO 9005
      CALL UERTST (IER,6HSSRAND)
 9005 RETURN
      END

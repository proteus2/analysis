C   IMSL ROUTINE NAME   - BESTAT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - COMPUTATIONS OF BASIC UNIVARIATE STATISTICS
C                           FROM DATA POSSIBLY CONTAINING MISSING
C                           VALUES, WITH WEIGHTING ON OPTION.
C
C   USAGE               - CALL BESTAT (X,IX,WT,NBR,XMISS,STATS,IS,WK,
C                           IER)
C
C   ARGUMENTS    X      - INPUT. X IS AN NBR(3) BY NBR(1) SUBMATRIX
C                           OF THE MATRIX (CALL IT XX) OF DATA FOR
C                           WHICH BASIC UNIVARIATE STATISTICS ARE
C                           DESIRED. THE LAST SUBMATRIX IN XX MAY
C                           HAVE FEWER THAN NBR(3) ROWS. IF ALL THE
C                           DATA ARE AVAILABLE AT ONCE, X CAN BE THE
C                           FULL DATA SET.
C                IX     - INPUT, ROW DIMENSION OF X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                WT     - INPUT VECTOR OF LENGTH NBR(3) IF WEIGHTING
C                           IS DESIRED. IN THIS CASE WT CONTAINS THE
C                           NONNEGATIVE WEIGHTS CORRESPONDING TO THE
C                           ROWS OF X. IF EQUAL WEIGHTING IS DESIRED,
C                           WT IS OF LENGTH 1 AND WT(1) MAY BE SET
C                           TO ANY NEGATIVE QUANTITY ON THE FIRST
C                           CALL TO BESTAT FOR A GIVEN PROBLEM, AND THEN
C                           OTHER ELEMENTS OF WT ARE NOT USED. IF WT(1)
C                           IS NONNEGATIVE, ALL THE ELEMENTS OF WT MUST
C                           BE NONNEGATIVE.
C                NBR    - INPUT VECTOR OF LENGTH 5. NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF VARIABLES.
C                           I=2, NUMBER OF OBSERVATIONS PER VARIABLE
C                             IN XX.
C                           I=3, NUMBER OF OBSERVATIONS PER VARIABLE IN
C                             EACH SUBMATRIX X, NOT INCLUDING THE LAST
C                             SUBMATRIX WHERE THE NUMBER MAY BE LESS
C                             THAN OR EQUAL TO NBR(3). HOWEVER, NBR(3)
C                             SHOULD BE THE SAME FOR ALL CALLS.
C                           I=4, THE NUMBER OF THE SUBMATRIX STORED IN
C                             X. SEE EXAMPLE.
C                             NBR(4) CAN ALSO BE USED AS A NO-DATA
C                             INDICATOR. IF NBR(4) IS NEGATIVE, NO DATA
C                             FROM X ARE USED AND THE ONLY COMPUTATIONS
C                             THAT ARE PERFORMED ARE TO COMPUTE THE
C                             STANDARD DEVIATIONS AND COEFFICIENTS OF
C                             SKEWNESS AND KURTOSIS FROM THE QUANTITIES
C                             IN STATS, WHICH HAVE BEEN ACCUMULATED ON
C                             PREVIOUS CALLS TO BESTAT. THIS METHOD IS
C                             GENERALLY NOT NEEDED, BUT IT IS USEFUL
C                             WHEN AT THE TIME OF A CALL TO BESTAT, IT
C                             WAS NOT KNOWN THAT X WAS THE LAST SUB-
C                             MATRIX OF XX. IN THIS CASE, NBR(2) WAS
C                             INITIALLY SET TO A LARGER VALUE THAN THE
C                             ACTUAL NUMBER OF OBSERVATIONS AND THEN A
C                             FINAL CALL TO BESTAT FOR POST-PROCESSING
C                             IS MADE WITH NBR(2) SET CORRECTLY AND
C                             NBR(4) SET TO A NEGATIVE VALUE.
C                           I=5, THE MISSING VALUE OPTION.
C                             IF NBR(5) = 0, ALL DATA IN X ARE ASSUMED
C                               TO BE VALID.
C                             OTHERWISE, THE VALUES IN XMISS ARE
C                               INTERPRETED AS MISSING VALUE CODES AND
C                               ANY X VALUE EQUAL TO THE CORRESPONDING
C                               VALUE IN XMISS IS EXCLUDED FROM THE
C                               COMPUTATIONS.
C                               IF NBR(5) IS POSITIVE, THE EXCLUSION IS
C                                 ELEMENTWISE.
C                               IF NBR(5) IS NEGATIVE, THE EXCLUSION IS
C                                 LISTWISE. (THE ENTIRE ROW OF X IS
C                                 EXCLUDED IF ANY OF THE VALUES OF THE
C                                 ROW IS EQUAL TO THE CORRESPONDING
C                                 VALUE IN XMISS.)
C                XMISS  - INPUT VECTOR OF LENGTH NBR(1) REQUIRED IF
C                           NBR(5) .NE. 0. IN THIS CASE, ANY VALUE IN
C                           THE ITH COLUMN OF THE DATA MATRIX X WHICH
C                           IS EQUAL TO THE CORRESPONDING VALUE IN THE
C                           ITH POSITION OF VECTOR XMISS IS ASSUMED
C                           TO BE INVALID OR MISSING. THE INVALID OR
C                           MISSING DATA ARE HANDLED BY EITHER ELEMENT-
C                           WISE OR LISTWISE DELETION AS SPECIFIED IN
C                           NBR(5).
C                STATS  - INPUT/OUTPUT 8 BY NBR(1) MATRIX CONTAINING IN
C                           EACH ROW STATISTICS ON ALL OF THE VARIABLES.
C                             STATS(1,*) CONTAINS MEANS,
C                             STATS(2,*) CONTAINS VARIANCES,
C                             STATS(3,*) CONTAINS STANDARD DEVIATIONS,
C                             STATS(4,*) CONTAINS COEFFICIENTS OF
C                                        SKEWNESS,
C                             STATS(5,*) CONTAINS COEFFICIENTS OF
C                                        EXCESS (KURTOSIS),
C                             STATS(6,*) CONTAINS MINIMA,
C                             STATS(7,*) CONTAINS MAXIMA,
C                             STATS(8,*) CONTAINS NUMBERS (COUNTS) OF
C                                        NONMISSING OBSERVATIONS.
C                IS     - INPUT, ROW DIMENSION OF STATS EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                WK     - WORK VECTOR OF LENGTH 3*(NBR(1)+1). WK SHOULD
C                           NOT BE CHANGED BETWEEN CALLS TO BESTAT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=34 INDICATES THAT FEWER THAN TWO VALID
C                             OBSERVATIONS WERE PRESENT FOR SOME
C                             VARIABLE. THE PERTINENT STATISTICS ARE
C                             SET TO NEGATIVE MACHINE INFINITY.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NBR(3)*(NBR(4)-1)
C                             EXCEEDS NBR(2) (WHEN NBR(4) IS POSITIVE).
C                           IER=130 INDICATES THAT NBR(1) IS LESS THAN
C                             1 OR NBR(2) IS LESS THAN 2 OR NBR(3)
C                             EXCEEDS NBR(2) OR NBR(3) IS LESS THAN 1.
C                           IER=131 INDICATES THAT NBR(4)=1 AND WT(1) IS
C                             NONNEGATIVE (INPLYING WEIGHTS ARE TO BE
C                             USED) AND THAT SOME WEIGHTS ARE NEGATIVE.
C
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BESTAT (X,IX,WT,NBR,XMISS,STATS,IS,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,IS,IER,NBR(1)
      REAL               X(IX,1),WT(1),XMISS(1),STATS(IS,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,MM,N1,N2,N3,NN
      REAL               CNT,CNT2,CNT3,CNTJ,CNTJ2,CNTJ3,CNTJM1,CNTM1,
     *                   FOUR,OMRAT,OMRAT2,OMRAT3,ONE,RAT,RAT2,RAT3,
     *                   SIX,STAT2J,STAT3J,SUMWT,THREE,TJ,
     *                   TJ2,TJ3,TJ4,WI,XIJ,XINFM,ZERO
      DATA               ZERO,ONE/0.0,1.0/
      DATA               THREE,FOUR,SIX/3.0,4.0,6.0/
      DATA               XINFM/ZFFFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
      N1 = NBR(1)
      N2 = 2*N1
      N3 = 3*N1
      IF (NBR(1) .GE. 1 .AND. NBR(2) .GE. 2 .AND. NBR(3) .LE. NBR(2)
     1   .AND. NBR(3) .GE. 1) GO TO 5
      IER = 130
      GO TO 9000
    5 IF (NBR(4) .LT. 0) GO TO 175
      IF (NBR(4) .EQ. 0) GO TO 10
      IF (NBR(3)*(NBR(4)-1) .LE. NBR(2)) GO TO 15
C                                  TERMINAL ERROR - INVALID NBR VALUES
   10 IER = 129
      GO TO 9000
   15 MM = (NBR(2)+NBR(3)-1)/NBR(3)
      IF (NBR(4) .NE. 1) GO TO 30
C                                  INITIALIZE WK,STATS
      WK(N3+1) = ZERO
      WK(N3+2) = ZERO
      WK(N3+3) = WT(1)
      DO 25 J=1,N1
         WK(J) = X(1,J)
         IF(NBR(5).NE.0) WK(J)=XMISS(J)
         WK(N1+J) = ZERO
         WK(N2+J) = ZERO
         DO 20 I=1,5
            STATS(I,J) = ZERO
   20    CONTINUE
         STATS(6,J) = X(1,J)
         STATS(7,J) = X(1,J)
         IF(NBR(5).NE.0) STATS(6,J)=XMISS(J)
         IF(NBR(5).NE.0) STATS(7,J)=XMISS(J)
         STATS(8,J) = ZERO
   25 CONTINUE
C                                  ALL ENTRIES EXCEPT FIRST AND SPECIAL
C                                    NO OPTS BEGIN HERE
   30 NN = NBR(3)
      CNT = WK(N3+1)
      SUMWT = WK(N3+2)
      IF (NBR(4) .EQ. MM) NN = NBR(2)-(MM-1)*NBR(3)
      IF (WK(N3+3).GE.ZERO) GO TO 95
C
C                                  ACCUMULATE MEANS AND ADJUSTED SUMS
C                                  OF SQUARES AND TEST FOR CONSTANT
C                                  OBSERVATION SET.
C
C
C                                  NO WEIGHTING
C
      IF(NBR(5).GT.0) GOTO 70
C                                  NO MISSING VALUES OR ELSE LISTWISE
C                                  DELETION OPTION
      DO 60 I=1,NN
         IF(NBR(5).EQ.0) GOTO 40
C                                  CHECK TO SEE IF ANY MISSING VALUES
         DO 35 J=1,N1
            IF(X(I,J).EQ.XMISS(J)) GOTO 60
   35    CONTINUE
   40    CNTM1 = CNT
         CNT = CNTM1 + ONE
         CNT2 = CNT*CNT
         CNT3 = CNT*CNT2
         DO 55 J=1,N1
            XIJ=X(I,J)
            IF(NBR(5).EQ.0) GO TO 45
            IF(WK(J).NE.XMISS(J)) GO TO 45
            STATS(6,J)=XIJ
            STATS(7,J)=XIJ
            WK(J)=XIJ
            GO TO 50
   45       IF(XIJ.LT.STATS(6,J))STATS(6,J)=XIJ
            IF(XIJ.GT.STATS(7,J))STATS(7,J)=XIJ
            IF (XIJ.NE.WK(J)) WK(N1+J)=ONE
   50       TJ = XIJ-STATS(1,J)
            TJ2 = TJ*TJ
            STATS(1,J) = STATS(1,J) + (TJ/CNT)
            STATS(2,J) = STATS(2,J) + TJ2*(CNTM1)/CNT
            STATS(4,J) = STATS(4,J) +
     *                TJ*(TJ2*(CNT2-ONE)/CNT2-THREE*STATS(2,J)/CNT)
            STATS(5,J) = STATS(5,J) +
     *                TJ*(TJ*(TJ2*(CNT3-ONE)/CNT3-SIX*STATS(2,J)/CNT2)-
     *                      FOUR*STATS(4,J)/CNT)
   55    CONTINUE
   60 CONTINUE
      DO 65 J=1,N1
         STATS(8,J) = CNT
   65 CONTINUE
      GO TO 170
   70 CONTINUE
C                                  ELEMENTWISE DELETION OPTION
      DO 90 I=1,NN
         DO 85 J=1,N1
            XIJ=X(I,J)
            IF(XIJ.EQ.XMISS(J)) GOTO 85
            IF(WK(J).NE.XMISS(J)) GOTO 75
            STATS(6,J)=XIJ
            STATS(7,J)=XIJ
            WK(J)=XIJ
            GO TO 80
   75       IF(XIJ.LT.STATS(6,J))STATS(6,J)=XIJ
            IF(XIJ.GT.STATS(7,J))STATS(7,J)=XIJ
            IF (XIJ.NE.WK(J)) WK(N1+J)=ONE
   80       CNTJM1 = STATS(8,J)
            CNTJ = CNTJM1 + ONE
            CNTJ2 = CNTJ*CNTJ
            CNTJ3 = CNTJ*CNTJ2
            STATS(8,J) = CNTJ
            TJ = XIJ-STATS(1,J)
            TJ2=TJ*TJ
            STATS(1,J) = STATS(1,J)+ TJ/CNTJ
            STATS(2,J) = STATS(2,J) + TJ2*CNTJM1/CNTJ
            STATS(4,J) = STATS(4,J) +
     *                   TJ*(TJ2*(CNTJ2-ONE)/CNTJ2-THREE*STATS(2,J)
     *                      /CNTJ)
            STATS(5,J) = STATS(5,J) +
     *                TJ*(TJ*(TJ2*(CNTJ3-ONE)/CNTJ3-SIX*STATS(2,J)
     *                      /CNTJ2)- FOUR*STATS(4,J)/CNTJ)
   85    CONTINUE
   90 CONTINUE
      GOTO 170
C
C                                  WEIGHTING
C
   95 CONTINUE
      IF(NBR(5).GT.0) GOTO 140
C                                  NO MISSING VALUES OR ELSE LISTWISE
C                                  DELETION OPTION
      DO 130 I=1,NN
         IF(NBR(5).EQ.0) GOTO 105
         DO 100 J=1,N1
            IF(X(I,J).EQ.XMISS(J)) GOTO 130
  100    CONTINUE
  105    WI = WT(I)
         IF (WI.GE.ZERO) GO TO 110
         IER = 131
         GO TO 9000
  110    IF (WI.EQ.ZERO) GO TO 130
         SUMWT = SUMWT+WI
         RAT = WI/SUMWT
         CNT = CNT+ONE
         DO 125 J=1,N1
            XIJ=X(I,J)
            IF(NBR(5).EQ.0) GO TO 115
            IF(WK(J).NE.XMISS(J))GO TO 115
            STATS(6,J)=XIJ
            STATS(7,J)=XIJ
            WK(J)=XIJ
            GO TO 120
  115       IF(XIJ.LT.STATS(6,J))STATS(6,J)=XIJ
            IF(XIJ.GT.STATS(7,J))STATS(7,J)=XIJ
            IF (XIJ.NE.WK(J)) WK(N1+J)=ONE
  120       TJ = XIJ-STATS(1,J)
            TJ2=TJ*TJ
            TJ3=TJ2*TJ
            TJ4=TJ3*TJ
            RAT2=RAT*RAT
            RAT3=RAT2*RAT
            OMRAT=ONE-RAT
            OMRAT2=ONE-RAT2
            OMRAT3=ONE-RAT3
            STATS(1,J) = STATS(1,J)+(TJ*RAT)
            STATS(2,J) = STATS(2,J) + TJ2*OMRAT*WI
            STATS(4,J) = STATS(4,J)-THREE*STATS(2,J)*RAT*TJ+
     *                   WI*OMRAT2*TJ3
            STATS(5,J) = STATS(5,J)-FOUR*STATS(4,J)*RAT*TJ-SIX*
     *                   STATS(2,J)*RAT2*TJ2+WI*OMRAT3*TJ4
  125    CONTINUE
  130 CONTINUE
      DO 135 J=1,N1
         STATS(8,J) = CNT
  135 CONTINUE
      GOTO 170
  140 CONTINUE
C                                  ELEMENTWISE DELETION OPTION
      DO 165 I=1,NN
         WI = WT(I)
         IF (WI.GE.ZERO) GO TO 145
         IER = 131
         GO TO 9000
  145    IF (WI.EQ.ZERO) GO TO 165
         DO 160 J=1,N1
            XIJ=X(I,J)
            IF(XIJ.EQ.XMISS(J)) GOTO 160
            IF(WK(J).NE.XMISS(J)) GOTO 150
            STATS(6,J)=XIJ
            STATS(7,J)=XIJ
            WK(J)=XIJ
            GO TO 155
  150       IF(XIJ.LT.STATS(6,J))STATS(6,J)=XIJ
            IF(XIJ.GT.STATS(7,J))STATS(7,J)=XIJ
            IF (XIJ.NE.WK(J)) WK(N1+J)=ONE
  155       WK(N2+J) = WK(N2+J)+WI
            STATS(8,J)=STATS(8,J) + ONE
            TJ = X(I,J)-STATS(1,J)
            TJ2=TJ*TJ
            TJ3=TJ2*TJ
            TJ4=TJ3*TJ
            RAT=WI/WK(N2+J)
            RAT2=RAT*RAT
            RAT3=RAT2*RAT
            OMRAT=ONE-RAT
            OMRAT2=ONE-RAT2
            OMRAT3=ONE-RAT3
            STATS(1,J) = STATS(1,J)+(TJ*RAT)
            STATS(2,J) = STATS(2,J) + TJ2*OMRAT*WI
            STATS(4,J) = STATS(4,J)-THREE*STATS(2,J)*RAT*TJ+
     *                   WI*OMRAT2*TJ3
            STATS(5,J) = STATS(5,J)-FOUR*STATS(4,J)*RAT*TJ-SIX*
     *                   STATS(2,J)*RAT2*TJ2+WI*OMRAT3*TJ4
  160    CONTINUE
  165 CONTINUE
  170 CONTINUE
      IF (NBR(4).NE.MM) GO TO 9005
  175 CONTINUE
C                                  POSTPROCESSING BEGINS HERE
C                                  COMPUTE THE VARIANCES AND THE
C                                  STANDARD DEVIATIONS
      DO 190 J=1,N1
         IF (STATS(8,J).GT.ONE) GO TO 180
         IER = 34
         STATS(2,J) = XINFM
         STATS(3,J) = XINFM
         STATS(4,J) = XINFM
         STATS(5,J) = XINFM
         GO TO 190
  180    IF (WK(N1+J).EQ.ZERO.OR.STATS(2,J).LE.ZERO) GO TO 185
         STAT2J = STATS(2,J)/STATS(8,J)
         STAT3J = SQRT(STAT2J)
         STATS(2,J) = STATS(2,J)/(STATS(8,J)-ONE)
         STATS(3,J) = SQRT(STATS(2,J))
         STATS(4,J) = STATS(4,J)/(STATS(8,J)*STAT3J**3)
         STATS(5,J) = STATS(5,J)/(STATS(8,J)*STAT2J**2) - THREE
         GO TO 190
  185     STATS(2,J) = ZERO
         STATS(3,J) = ZERO
         STATS(4,J) = ZERO
         STATS(5,J) = ZERO
  190 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BESTAT')
 9005 WK(N3+1) = CNT
      WK(N3+2) = SUMWT
      RETURN
      END

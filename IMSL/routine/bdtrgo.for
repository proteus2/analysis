C   IMSL ROUTINE NAME   - BDTRGO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TRANSGENERATION OF THE COLUMNS OF A MATRIX
C                           (OUT-OF-CORE VERSION)
C
C   USAGE               - CALL BDTRGO (X,NT,ITRG,IT,C,IER)
C
C   ARGUMENTS    X      - ON INPUT, X IS A ROW OF THE MATRIX TO BE
C                           TRANSGENERATED. X HAS LENGTH M.
C                         ON OUTPUT, X IS THE TRANSGENERATED ROW OF THE
C                           MATRIX.
C                NT     - NUMBER OF TRANSGENERATIONS TO BE PERFORMED.
C                           (INPUT)
C                ITRG   - NT BY 4 INPUT TRANSGENERATION MATRIX.
C                IT     - ROW DIMENSION OF MATRIX ITRG EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - INPUT VECTOR OF CONSTANTS REQUIRED BY THE
C                           TRANSGENERATIONS. C HAS LENGTH NT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=128 + I INDICATES AN ERROR OCCURRED IN
C                             THE I-TH TRANSGENERATION.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  X MUST BE A SINGLE DIMENSION ARRAY LONG ENOUGH TO HOLD
C                ADDITIONAL COMPONENTS RESULTING FROM TRANSGENERATIONS.
C            2.  TO TRANSGENERATE THE COLUMNS OF AN N BY M MATRIX X,
C                ONE MUST READ EACH ROW OF THE MATRIX INTO A SINGLE
C                DIMENSION ARRAY Y AND CALL BDTRGO N TIMES USING THE
C                SAME ITRG MATRIX FOR EACH CALL.
C            3.  TRANSGENERATIONS AVAILABLE VIA BDTRGO FOLLOW. THE
C                LETTER C REFERS TO ELEMENTS OF VECTOR C IN THE
C                ARGUMENT LIST.
C
C                CODE  TRANSGENERATION               RESTRICTION
C                ----  ---------------               -----------
C                1     SQRT(X(I)) INTO X(K)          X(I) NON-NEGATIVE
C                2     LOG(10) X(I) INTO X(K)        X(I) POSITIVE
C                3     LN X(I) INTO X(K)             X(I) POSITIVE
C                4     E**X(I) INTO X(K)
C                5     SIN X(I) INTO X(K)            X(I) RADIANS
C                6     COS X(I) INTO X(K)            X(I) RADIANS
C                7     SIN X(I) INTO X(K)            X(I) DEGREES
C                8     COS X(I) INTO X(K)            X(I) DEGREES
C                9     ARCSIN(X(I)) INTO X(K)        X(I) IN INCLUSIVE
C                      (RADIANS)                     INTERVAL (-1,1)
C                10    ARCTAN(X(I)) INTO X(K)
C                      (RADIANS)
C                11    ARCSINH(X(I)) INTO X(K)
C                12    X(I) INTO X(K)
C                13    X(I) + C INTO X(K)
C                14    X(I) * C INTO X(K)
C                15    C/X(I) INTO X(K)              X(I).NE.0
C                16    X(I)**C INTO X(K)             X(I) NON-NEGATIVE
C                      (0**0 IS DEFINED AS ZERO)
C                17    C**X(I) INTO X(K)             C NON-NEGATIVE
C                      (0**0 IS DEFINED AS ZERO)
C                18    1 INTO X(K), IF X(I).EQ.C
C                      0 INTO X(K), IF X(I).NE.C
C                19    1 INTO X(K), IF X(I).GE.C
C                      0 INTO X(K), IF X(I).LT.C
C                20    1 INTO X(K), IF X(I).GT.C
C                      0 INTO X(K), IF X(I).LE.C
C                21    C INTO X(K)
C                22    X(I)+X(J) INTO X(K)
C                23    X(I)-X(J) INTO X(K)
C                24    X(I)*X(J) INTO X(K)
C                25    X(I)/X(J) INTO X(K)           X(J).NE.0
C                26    X(I)**X(J) INTO X(K)          X(I) NON-NEGATIVE
C                      (0**0 IS DEFINED AS ZERO)
C                27    C INTO X(K), IF X(I).EQ.X(J)
C                      0 INTO X(K), IF X(I).NE.X(J)
C                28    C INTO X(K), IF X(I).GE.X(J)
C                      0 INTO X(K), IF X(I).LT.X(J)
C                29    C INTO X(K), IF X(I).GT.X(J)
C                      0 INTO X(K), IF X(I).LE.X(J)
C                30    X(J) INTO X(K), IF X(I).EQ.C
C                      X(K) UNCHANGED, IF X(I).NE.C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDTRGO (X,NT,ITRG,IT,C,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               X(1),C(1)
      INTEGER            NT,IT,ITRG(IT,4),IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TD,T
      REAL               CON,ZERO,TEMP,EPS
      INTEGER            I,LL,II,JJ,KK
      DATA               ZERO/0.0/
      DATA               EPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DO 190 I=1,NT
         LL = ITRG(I,4)
         IF (LL.LE.0.OR.LL.GT.30) GO TO 195
         II = ITRG(I,1)
         JJ = ITRG(I,2)
         KK = ITRG(I,3)
         IF (LL.GE.13.AND.LL.LE.21.OR.LL.GE.27.AND.LL.LE.30) CON = C(I)
         GO TO (5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,85,100,110,
     *   120,125,130,135,140,145,150,160,170,175,180,185), LL
    5    IF (X(II).LT.ZERO) GO TO 195
         X(KK) = SQRT(X(II))
         GO TO 190
   10    IF (X(II).LE.ZERO) GO TO 195
         X(KK) = ALOG10(X(II))
         GO TO 190
   15    IF (X(II).LE.ZERO) GO TO 195
         X(KK) = ALOG(X(II))
         GO TO 190
   20    X(KK) = EXP(X(II))
         GO TO 190
   25    X(KK) = SIN(X(II))
         GO TO 190
   30    X(KK) = COS(X(II))
         GO TO 190
   35    X(KK) = SIN(X(II)*.0174533)
         GO TO 190
   40    X(KK) = COS(X(II)*.0174533)
         GO TO 190
   45    TEMP = X(II)
         IF (TEMP.LT.-1..OR.TEMP.GT.1.) GO TO 195
C                                  CALCULATE THE ARCSIN
         X(KK) = ATAN2(TEMP,SQRT(1.-TEMP*TEMP))
         GO TO 190
   50    X(KK) = ATAN(X(II))
         GO TO 190
   55    TD = X(II)
         X(KK) = DLOG(TD+DSQRT(TD*TD+1.D0))
         GO TO 190
   60    X(KK) = X(II)
         GO TO 190
   65    X(KK) = X(II)+CON
         GO TO 190
   70    X(KK) = CON*X(II)
         GO TO 190
   75    IF (X(II)) 80,195,80
   80    X(KK) = CON/X(II)
         GO TO 190
   85    TEMP = X(II)
         IF (TEMP) 195,95,90
   90    X(KK) = TEMP**CON
         GO TO 190
   95    X(KK) = ZERO
         GO TO 190
  100    IF (CON) 195,95,105
  105    X(KK) = CON**X(II)
         GO TO 190
  110    IF (ABS(X(II)-CON).GT.EPS*ABS(CON)) GO TO 95
  115    X(KK) = 1.
         GO TO 190
  120    T = X(II)
         IF (T.GT.CON) GO TO 115
         IF (DABS(T-CON).LE.EPS*ABS(CON)) GO TO 115
         GO TO 95
  125    IF (X(II).GT.CON) GO TO 115
         GO TO 95
  130    X(KK) = CON
         GO TO 190
  135    X(KK) = X(II)+X(JJ)
         GO TO 190
  140    X(KK) = X(II)-X(JJ)
         GO TO 190
  145    X(KK) = X(II)*X(JJ)
         GO TO 190
  150    IF (X(JJ)) 155,195,155
  155    X(KK) = X(II)/X(JJ)
         GO TO 190
  160    IF (X(II)) 195,95,165
  165    X(KK) = X(II)**X(JJ)
         GO TO 190
  170    IF (ABS(X(II)-X(JJ)).GT.EPS*ABS(X(JJ))) GO TO 95
         GO TO 130
  175    TEMP = X(II)
         T = X(JJ)
         IF (TEMP.GT.T) GO TO 130
         IF (DABS(TEMP-T).GT.EPS*DABS(T)) GO TO 95
         GO TO 130
  180    IF (X(II).GT.X(JJ)) GO TO 130
         GO TO 95
  185    IF (ABS(X(II)-CON).GT.EPS*ABS(CON)) GO TO 190
         X(KK) = X(JJ)
  190 CONTINUE
      GO TO 9005
  195 IER = 128+I
      CALL UERTST (IER,'BDTRGO')
 9005 RETURN
      END

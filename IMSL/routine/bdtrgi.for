C   IMSL ROUTINE NAME   - BDTRGI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - TRANSGENERATION OF THE COLUMNS OF A MATRIX
C                           (IN-CORE VERSION)
C
C   USAGE               - CALL BDTRGI (X,N,IX,NT,ITRG,IT,C,INFER,IER)
C
C   ARGUMENTS    X      - ON INPUT, X IS THE N BY M MATRIX TO BE
C                           TRANSGENERATED.
C                         ON OUTPUT, X IS THE TRANSGENERATED MATRIX.
C                           THE NUMBER OF COLUMNS MAY OR MAY NOT HAVE
C                           BEEN CHANGED.
C                N      - NUMBER OF ROWS IN THE INPUT MATRIX X. (INPUT)
C                IX     - ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NT     - NUMBER OF TRANSGENERATIONS TO BE PERFORMED.
C                           (INPUT)
C                ITRG   - NT BY 4 INPUT TRANSGENERATION MATRIX. SEE
C                           ALGORITHM.
C                IT     - ROW DIMENSION OF MATRIX ITRG EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - INPUT VECTOR OF CONSTANTS REQUIRED BY THE
C                           TRANSGENERATIONS. C HAS LENGTH NT.
C                INFER  - OUTPUT VECTOR OF LENGTH 2 CONTAINING ERROR
C                           INFORMATION.
C                         IF IER IS NOT EQUAL TO 0,
C                           INFER(1) CONTAINS THE NUMBER OF THE
C                             TRANSGENERATION IN WHICH THE ERROR
C                             OCCURRED.
C                           INFER(2) CONTAINS THE NUMBER OF THE ROW IN
C                             THE X MATRIX IN WHICH THE ERROR OCCURRED.
C                           OTHERWISE INFER(1)=INFER(2)=0.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT A NON-EXISTENT
C                             TRANSGENERATION CODE APPEARED.
C                           IER=130 INDICATES THAT A COMPONENT OF X WAS
C                             NOT IN THE ALLOWABLE RANGE.
C                           IER=131 INDICATES THAT A CONSTANT WAS NOT
C                             IN THE ALLOWABLE TRANSGENERATION RANGE.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE COLUMN DIMENSION OF X MUST BE SUFFICIENTLY LARGE
C                TO ACCOMMODATE ADDITIONAL COLUMNS RESULTING FROM
C                TRANSGENERATIONS.
C            2.  TRANSGENERATIONS AVAILABLE VIA BDTRGI FOLLOW. THE
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
C                31    X(I)-XX(I) INTO X(K)
C                      XX(I) = SUM FROM J=1 TO N OF
C                      X(J,I)/N
C                32    X(I)/S(I) INTO X(K)           S(I).NE.0
C                      S(I)**2 = SUM FROM J=1 TO N
C                      OF (X(J,I)-XX(I))**2/(N-1)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDTRGI (X,N,IX,NT,ITRG,IT,C,INFER,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IX,NT,IT,ITRG(IT,1),INFER(2),IER
      REAL               X(IX,1),C(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TD,T,Z
      REAL               EPS,TEMP,CON,ZERO,ONEM,ONE
      DATA               ZERO/0.0/,ONEM/-1.0/,ONE/1.0/
      DATA               EPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      INFER(1) = 0
      INFER(2) = 0
      IER = 0
      DO 380 I=1,NT
         LL = ITRG(I,4)
         IF (LL.GT.0.AND.LL.LE.32) GO TO 5
         J = 0
         IER = 129
         GO TO 395
    5    II = ITRG(I,1)
         JJ = ITRG(I,2)
         KK = ITRG(I,3)
         IF (LL.GE.13.AND.LL.LE.21.OR.LL.GE.27.AND.LL.LE.30) CON = C(I)
         GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,1
     1   75,185,200,215,230,240,250,260,270,280,295,310,325,335,345,360)
     2   , LL
   10    DO 15 J=1,N
            IF (X(J,II).LT.ZERO) GO TO 385
   15    X(J,KK) = SQRT(X(J,II))
         GO TO 380
   20    DO 25 J=1,N
            IF (X(J,II).LE.ZERO) GO TO 385
   25    X(J,KK) = ALOG10(X(J,II))
         GO TO 380
   30    DO 35 J=1,N
            IF (X(J,II).LE.ZERO) GO TO 385
   35    X(J,KK) = ALOG(X(J,II))
         GO TO 380
   40    DO 45 J=1,N
   45    X(J,KK) = EXP(X(J,II))
         GO TO 380
   50    DO 55 J=1,N
   55    X(J,KK) = SIN(X(J,II))
         GO TO 380
   60    DO 65 J=1,N
   65    X(J,KK) = COS(X(J,II))
         GO TO 380
   70    DO 75 J=1,N
   75    X(J,KK) = SIN(X(J,II)*.0174533)
         GO TO 380
   80    DO 85 J=1,N
   85    X(J,KK) = COS(X(J,II)*.0174533)
         GO TO 380
   90    DO 95 J=1,N
            TEMP = X(J,II)
            IF (TEMP.LT.ONEM.OR.TEMP.GT.ONE) GO TO 385
C                                  CALCULATE THE ARCSIN
   95    X(J,KK) = ATAN2(TEMP,SQRT(1.-TEMP*TEMP))
         GO TO 380
  100    DO 105 J=1,N
  105    X(J,KK) = ATAN(X(J,II))
         GO TO 380
  110    DO 115 J=1,N
            TD = X(J,II)
  115    X(J,KK) = DLOG(TD+DSQRT(TD*TD+1.D0))
         GO TO 380
  120    DO 125 J=1,N
  125    X(J,KK) = X(J,II)
         GO TO 380
  130    DO 135 J=1,N
  135    X(J,KK) = X(J,II)+CON
         GO TO 380
  140    DO 145 J=1,N
  145    X(J,KK) = CON*X(J,II)
         GO TO 380
  150    DO 155 J=1,N
            IF (X(J,II).EQ.ZERO) GO TO 385
  155    X(J,KK) = CON/X(J,II)
         GO TO 380
  160    DO 170 J=1,N
            TEMP = X(J,II)
            IF (TEMP.LT.ZERO) GO TO 385
            IF (TEMP.EQ.ZERO) GO TO 165
            X(J,KK) = TEMP**CON
            GO TO 170
  165       X(J,KK) = ZERO
  170    CONTINUE
         GO TO 380
  175    J = 0
         IF (CON.LT.ZERO) GO TO 390
         IF (CON.EQ.ZERO) GO TO 230
         DO 180 J=1,N
  180    X(J,KK) = CON**X(J,II)
         GO TO 380
  185    TEMP = ABS(CON)
         DO 195 J=1,N
            IF (ABS(X(J,II)-CON).LE.EPS*TEMP) GO TO 190
            X(J,KK) = ZERO
            GO TO 195
  190       X(J,KK) = ONE
  195    CONTINUE
         GO TO 380
  200    TEMP = ABS(CON)
         DO 210 J=1,N
            T = X(J,II)
            IF (T.GT.CON) GO TO 205
            IF (DABS(T-CON).LE.EPS*TEMP) GO TO 205
            X(J,KK) = ZERO
            GO TO 210
  205       X(J,KK) = ONE
  210    CONTINUE
         GO TO 380
  215    DO 225 J=1,N
            IF (X(J,II).GT.CON) GO TO 220
            X(J,KK) = ZERO
            GO TO 225
  220       X(J,KK) = ONE
  225    CONTINUE
         GO TO 380
  230    DO 235 J=1,N
  235    X(J,KK) = CON
         GO TO 380
  240    DO 245 J=1,N
  245    X(J,KK) = X(J,II)+X(J,JJ)
         GO TO 380
  250    DO 255 J=1,N
  255    X(J,KK) = X(J,II)-X(J,JJ)
         GO TO 380
  260    DO 265 J=1,N
  265    X(J,KK) = X(J,II)*X(J,JJ)
         GO TO 380
  270    DO 275 J=1,N
            TEMP = X(J,JJ)
            IF (TEMP.EQ.ZERO) GO TO 385
  275    X(J,KK) = X(J,II)/X(J,JJ)
         GO TO 380
  280    DO 290 J=1,N
            IF (X(J,II).LT.ZERO) GO TO 385
            IF (X(J,II).EQ.ZERO) GO TO 285
            X(J,KK) = X(J,II)**X(J,JJ)
            GO TO 290
  285       X(J,KK) = ZERO
  290    CONTINUE
         GO TO 380
  295    DO 305 J=1,N
            IF (ABS(X(J,II)-X(J,JJ)).LE.EPS*ABS(X(J,JJ))) GO TO 300
            X(J,KK) = ZERO
            GO TO 305
  300       X(J,KK) = CON
  305    CONTINUE
         GO TO 380
  310    DO 320 J=1,N
            TEMP = X(J,II)
            T = X(J,JJ)
            IF (TEMP.GT.T) GO TO 315
            IF (DABS(TEMP-T).LE.EPS*DABS(T)) GO TO 315
            X(J,KK) = ZERO
            GO TO 320
  315       X(J,KK) = CON
  320    CONTINUE
         GO TO 380
  325    DO 330 J=1,N
            TEMP = ZERO
            IF (X(J,II).GT.X(J,JJ)) TEMP = CON
  330    X(J,KK) = TEMP
         GO TO 380
  335    TEMP = ABS(CON)
         DO 340 J=1,N
            IF (ABS(X(J,II)-CON).GT.EPS*TEMP) GO TO 340
            X(J,KK) = X(J,JJ)
  340    CONTINUE
         GO TO 380
  345    TD = 0.D0
         DO 350 J=1,N
  350    TD = TD+X(J,II)
         TD = TD/N
         DO 355 J=1,N
  355    X(J,KK) = X(J,II)-TD
         GO TO 380
  360    TD = 0.D0
         DO 365 J=1,N
  365    TD = TD+X(J,II)
         TD = TD/N
         T = 0.D0
         DO 370 J=1,N
            Z = X(J,II)-TD
  370    T = T+Z*Z
         T = T/(N-1)
         T = DSQRT(T)
         J = 0
         IF (T.EQ.0.D0) GO TO 385
         DO 375 J=1,N
  375    X(J,KK) = X(J,II)/T
  380 CONTINUE
      GO TO 9005
  385 IER = 130
      GO TO 395
  390 IER = 131
  395 INFER(1) = I
      INFER(2) = J
 9000 CONTINUE
      CALL UERTST (IER,'BDTRGI')
 9005 RETURN
      END

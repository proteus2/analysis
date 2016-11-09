C   IMSL ROUTINE NAME   - FTRDIF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TRANSFORMATIONS, DIFFERENCES AND SEASONAL
C                           DIFFERENCES OF A TIME SERIES FOR MODEL
C                           IDENTIFICATION
C
C   USAGE               - CALL FTRDIF (ID1,ID2,IP,IS,LZ,Z,SHIFT,LW,IER)
C
C   ARGUMENTS    ID1    - INPUT ORDER OF NON-SEASONAL DIFFERENCE.
C                           ID1 MUST BE GREATER THAN OR EQUAL TO ZERO.
C                ID2    - INPUT ORDER OF SEASONAL DIFFERENCE.
C                           ID2 MUST BE GREATER THAN OR EQUAL TO ZERO.
C                IP     - INPUT TRANSFORMATION EXPONENT.
C                         IF IP IS EQUAL TO ZERO, A LOGARITHMIC
C                           TRANSFORMATION OF THE FORM
C                           Z(I) = LOG(Z(I)+SHIFT) IS PERFORMED ON THE
C                           INPUT SERIES, WHERE LOG IS THE BASE E
C                           LOGARITHM FUNCTION. SEE DESCRIPTION OF SHIFT
C                           BELOW.
C                         IF IP IS NOT EQUAL TO ZERO, AN EXPONENTIAL
C                           TRANSFORMATION OF THE FORM
C                           Z(I) = Z(I)**IP IS PERFORMED ON THE INPUT
C                           SERIES.
C                IS     - INPUT LENGTH OF SEASONAL PERIOD.
C                           IS MUST BE GREATER THAN OR EQUAL TO ZERO.
C                LZ     - INPUT LENGTH OF THE TIME SERIES.
C                           LZ MUST BE GREATER THAN OR EQUAL TO TWO.
C                Z      - INPUT AND OUTPUT VECTOR OF LENGTH LZ.
C                         ON INPUT, Z CONTAINS THE TIME SERIES.
C                         ON OUTPUT, Z CONTAINS THE TRANSFORMED AND
C                           DIFFERENCED TIME SERIES.
C                         GENERALLY THE OUTPUT TIME SERIES WILL BE
C                            SHORTER THAN THE INPUT TIME SERIES.
C                SHIFT  - OUTPUT CONSTANT.
C                         IF IP IS EQUAL TO ZERO, SHIFT WILL BE SET TO
C                           1.0-ZMIN WHERE ZMIN IS THE MINIMUM OF Z(I),
C                           I=1,...,LZ.
C                         IF IP IS NOT EQUAL TO ZERO, SHIFT IS NOT
C                           USED AND WILL BE SET TO ZERO.
C                LW     - OUTPUT TIME SERIES LENGTH.
C                           LW = LZ-ID1-IS*ID2. THE PARAMETERS LZ,ID1,
C                           ID2, AND IS SHOULD BE SPECIFIED SO THAT
C                           LZ IS GREATER THAN OR EQUAL TO ONE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT ID1 OR ID2 IS LESS
C                             THAN ZERO.
C                           IER=130 INDICATES THAT LZ IS LESS THAN TWO.
C                           IER=131 INDICATES THAT LW IS LESS THAN ONE.
C                           IER=132 INDICATES THAT IS IS LESS THAN ZERO
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      FTRDIF SHOULD BE USED WHEN THE GIVEN TIME SERIES IS
C                EXPONENTIAL, LOGARITHMIC, OR CONTAINS RAMPS OR
C                PLATEAUS. THE OUTPUT TIME SERIES WILL NOT EXHIBIT
C                THESE CHARACTERISTICS AND HENCE IS BETTER FOR
C                STOCHASTIC MODELING TECHNIQUES IN THE SENSE THAT IT
C                CAN POSSIBLY BE FITTED BY AUTOREGRESSIVE AND MOVING
C                AVERAGE MODELS.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTRDIF (ID1,ID2,IP,IS,LZ,Z,SHIFT,LW,IER)
C             SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ID1,ID2,IP,IS,LZ,LW,IER
      REAL               Z(LZ),SHIFT
C             SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L
      REAL               ZERO,ONE
      DATA               ZERO/0.0/,ONE/1.0/
C             FIRST EXECUTABLE STATEMENT
      IER = 129
      IF(ID1 .LT. 0 .OR. ID2 .LT. 0) GO TO 9000
      IER = 130
      IF (LZ .LT. 2) GO TO 9000
      IER = 131
      K = LZ-ID1-(IS*ID2)
      IF (K .LT. 1) GO TO 9000
      IER = 132
      IF (IS .LT. 0) GO TO 9000
      IER = 0
      K = IS+1
      SHIFT = ZERO
      IF (IP .NE. 0) GO TO 15
C           FIND MINIMUM VALUE IN Z
      SHIFT = Z(1)
      DO 5  J=2,LZ
         IF (SHIFT .GT. Z(J)) SHIFT = Z(J)
    5 CONTINUE
      SHIFT = -SHIFT+ONE
      DO 10  J=1,LZ
   10 Z(J) = ALOG(Z(J)+SHIFT)
      GO TO 25
   15 IF (IP .EQ. 1) GO TO 25
      DO 20 I=1,LZ
   20 Z(I) = Z(I)**IP
   25 LW = LZ
      IF (ID2 .EQ. 0) GO TO 40
C                TRANSFORM TO TIME SERIES
      DO 35  L=1,ID2
         DO 30  I=K,LW
   30    Z(I-IS) = Z(I)-Z(I-IS)
         LW = LW - IS
   35 CONTINUE
   40 IF (ID1 .EQ. 0) GO TO 9005
      DO 50  L=1,ID1
         DO 45  I=2,LW
   45    Z(I-1) = Z(I)-Z(I-1)
         LW = LW-1
   50 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'FTRDIF')
 9005 RETURN
      END

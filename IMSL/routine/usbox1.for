C   IMSL ROUTINE NAME   - USBOX1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY ROUTINE USBOX
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USBOX1 (X,N,FIVNO,SC,ISTRT,MAXL)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,ISTRT,MAXL
      REAL               X(1),FIVNO(5),SC
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,IL1,IL2,IS,IS1,ITEMP,I1,I2,I3,I4,I5,IJ,J
      REAL               TEMP(3,129),SD(2),CSTAR(3),K(10),BLANK,DASH,
     *                   ONEH,PLUS,SD1,STEP,X1
      DATA               K /1H ,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1H+/,
     *                   PLUS /1H+/,CSTAR /1H*,1HO,1HX/,DASH /1H-/,
     *                   BLANK /1H /,ONEH /1.5/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZATIONS
C                                  BLANK OUT MATRIX
      CALL UGETIO(1,NIN,NOUT)
      DO 10 I=1,3
         DO 5 J=1,MAXL
    5    TEMP(I,J) = BLANK
   10 CONTINUE
      X1 = X(1)
C                                  IF ALL DATA EQUAL,PLOT ONE POINT
      ITEMP = (FIVNO(5)-FIVNO(1))/SC
      IF (ITEMP.LE.0) GO TO 15
C                                  IF N.LE.5 PLOT ALL POINTS
      IF (N.GT.5) GO TO 25
   15 J = 0
      IS1 = 0
      DO 20 I=1,N
         IS = ISTRT+(X(I)-X1)/SC
         J = J+1
         IF (IS.NE.IS1) J = 1
         IS1 = IS
         IF (J.GT.10) J = 10
         TEMP(2,IS) = CSTAR(3)
         TEMP(3,IS) = K(J)
   20 CONTINUE
      GO TO 90
   25 STEP = ONEH*(FIVNO(4)-FIVNO(2))
C                                  BEGIN FILLING IN MIDDLE LINE
C                                  START TWO STEPS BELOW, THEN ONE
C                                  STEP BELOW.
      IS1 = 0
      I = 1
      J = 0
      SD(2) = FIVNO(2)-STEP
      SD(1) = SD(2)-STEP
      DO 35 II=1,2
   30    IF (X(I).GT.SD(II)) GO TO 35
         IS = ISTRT+(X(I)-X1)/SC
         J = J+1
         IF (IS.NE.IS1) J = 1
         IS1 = IS
         IF (J.GT.10) J = 10
         TEMP(2,IS) = CSTAR(II)
         TEMP(3,IS) = K(J)
         I = I+1
         GO TO 30
   35 CONTINUE
C                                  DETERMINE ENDS OF WHISKERS, ENDS
C                                  OF BOX AND MEDIAN
      I1 = ISTRT+(X(I)-X1)/SC
      I2 = ISTRT+(FIVNO(2)-X1)/SC
      I3 = ISTRT+(FIVNO(3)-X1)/SC
      I4 = ISTRT+(FIVNO(4)-X1)/SC
      I = (3*N)/4-2
      SD1 = FIVNO(4)+STEP
   40 I = I+1
      IF (I.GT.N) GO TO 45
      IF (X(I).LT.SD1) GO TO 40
   45 I5 = ISTRT+(X(I-1)-X1)/SC
      IF (I.GT.N) GO TO 60
C                                  IDENTIFY POINTS BEYOND FENCE AND
C                                  BEYOND TWO STEPS
      SD1 = SD1+STEP
   50 IF (X(I).GT.SD1) GO TO 55
      IS = ISTRT+(X(I)-X1)/SC
      J = J+1
      IF (IS.NE.IS1) J = 1
      IS1 = IS
      IF (J.GT.10) J = 10
      TEMP(2,IS) = CSTAR(2)
      TEMP(3,IS) = K(J)
      IF (I.GE.N) GO TO 60
      I = I+1
      GO TO 50
   55 IS = ISTRT+(X(I)-X1)/SC
      J = J+1
      IF (IS.NE.IS1) J = 1
      IS1 = IS
      IF (J.GT.10) J = 10
      TEMP(2,IS) = CSTAR(3)
      TEMP(3,IS) = K(J)
      IF (I.GE.N) GO TO 60
      I = I+1
      GO TO 55
C                                  FILL IN BOX AND MEDIAN
   60 TEMP(2,I2) = PLUS
      TEMP(2,I4) = PLUS
      TEMP(2,I1) = PLUS
      TEMP(2,I5) = PLUS
      TEMP(1,I2) = PLUS
      TEMP(1,I4) = PLUS
      TEMP(3,I2) = PLUS
      TEMP(3,I4) = PLUS
      TEMP(2,I3) = CSTAR(1)
      IF (I4.LE.I2) GO TO 70
      IL1 = I2+1
      IL2 = I4-1
      DO 65 II=IL1,IL2
         TEMP(1,II) = DASH
         TEMP(3,II) = DASH
   65 CONTINUE
C                                  FILL IN WHISKERS
   70 IL1 = I1+1
      IF (I2.LE.IL1) GO TO 80
      IL2 = I2-1
      DO 75 II=IL1,IL2
   75 TEMP(2,II) = DASH
   80 IL1 = I4+1
      IF (I5.LE.IL1) GO TO 90
      IL2 = I5-1
      DO 85 II=IL1,IL2
   85 TEMP(2,II) = DASH
C
C                                  NOW PRINT THE PLOT
   90 CONTINUE
      IF (MAXL.EQ.80) WRITE (NOUT,95) ((TEMP(I,J),J=1,79),I=1,3)
      IF (MAXL.EQ.129) WRITE (NOUT,100) ((TEMP(I,J),J=1,128),I=1,3)
   95 FORMAT (1X, 79A1)
  100 FORMAT (1X, 128A1)
      RETURN
      END

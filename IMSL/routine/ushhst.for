C   IMSL ROUTINE NAME   - USHHST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - PRINT A HORIZONTAL HISTOGRAM
C
C   USAGE               - CALL USHHST (T,N,IOPT,IER)
C
C   ARGUMENTS    T      - VECTOR OF LENGTH N CONTAINING THE
C                           FREQUENCIES (COUNTS). (INPUT)
C                           ELEMENTS MUST BE NON-NEGATIVE.
C                N      - LENGTH OF T VECTOR, (NUMBER OF BARS TO PRINT).
C                           N MUST BE A POSITIVE INTEGER. (INPUT)
C                IOPT   - OPTION VECTOR OF LENGTH 5. (INPUT)
C                           IOPT(I) DENOTES,
C                           I=1, SPACING BETWEEN HORIZONTAL HISTOGRAM
C                             LINES.  0, 1, OR 2 SPACES ARE ALLOWED.
C                           I=2, ZERO WILL CAUSE A FULL (HORIZONTAL)
C                             PAGE HISTOGRAM. IOPT(2)=1 WILL LIMIT THE
C                             WIDTH TO EIGHT AND ONE HALF INCHES.
C                           I=3, THE UPPER LIMIT OF THE NUMBER OF
C                             LINES TO PRINT WITHIN THE HISTOGRAM
C                             PER PAGE.  AFTER THAT NUMBER OF LINES
C                             IS PRINTED, THE ROUTINE WILL SKIP TO A
C                             NEW PAGE TO CONTINUE PRINTING.
C                             IF IOPT(3)=0, THEN THE MAXIMUM NUMBER OF
C                             LINES COINCIDES WITH THE STANDARD PRINTER
C                             PAGE.
C                           I=4, IF ZERO, THEN, IF MULTIPLE PAGES ARE
C                             REQUIRED, THE FREQUENCY LINE (BOTTOM)
C                             AND THE CLASS LINE (TOP) ARE REPEATED FOR
C                             EACH PAGE.  IF NONZERO, THE CLASS AND
C                             FREQUENCY WILL BE PRINTED ON THE FIRST AND
C                             LAST PAGE OF THE HISTOGRAM, RESPECTIVELY.
C                           I=5, IF ZERO, SKIP TO NEW PAGE BEFORE
C                             PRINTING FIRST LINE.  IF NONZERO, TWO
C                             SPACES ARE SKIPPED AND PRINTING BEGINS ON
C                             THE SAME PAGE.  IOPT(5) SHOULD BE NONZERO
C                             IF THE USER WISHES TO PRINT A TITLE ABOVE
C                             THE HISTOGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING WITH FIX ERROR
C                           IER=65 MEANS THAT IOPT(1) IS NOT 0, 1, OR 2.
C                             THE ZERO OPTION IS USED FOR IOPT(1).
C                           IER=66 MEANS THAT IOPT(2) IS NOT 0 OR 1.
C                             THE ZERO OPTION IS USED FOR IOPT(2).
C                         TERMINAL ERROR
C                           IER=132 MEANS THAT THE LENGTH OF T IS NOT
C                             POSITIVE.
C                           IER=133 MEANS THAT THE MAXIMUM ELEMENT OF T
C                             IS LESS THAN ONE.  THE BODY OF THE
C                             HISTOGRAM IS BLANK.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO.  SEE THE UGETIO DOCUMENT FOR DETAILS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USHHST (T,N,IOPT,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,IOPT(5),N
      REAL               T(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IAK,IIK,IIL,IJ,IK,ILL,IL,IM,IN,IOP,IOPT1,IP,
     *                   I,J1,J2,JAK,JI,JJ,JK,J,K,NIN,NI,NL,NM,NN,
     *                   NOUT,NP1
      REAL               AA,AK,CHAR,HYP,TMAX,WK(161),WW
      DATA               AA/1H /,HYP/1H*/,WW/4H----/,CHAR/1HI/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO (1,NIN,NOUT)
      IER = 0
      IOPT1 = IOPT(1)
      IF (N.LE.0) GO TO 75
      IF (IOPT1.GE.0 .AND. IOPT1.LE.2) GO TO 5
      IER = 65
      IOPT1 = 0
    5 NM = N
      NP1 = N+1
      NN = NM*(IOPT1+1)-IOPT1
      TMAX = 0.0
C                                  FIND WIDTH OF HISTOGRAM
      DO 10 I=1,N
         IF (TMAX.LT.T(I)) TMAX = T(I)
   10 CONTINUE
      AK = AINT(.99168+.00833*TMAX)
      IF (IOPT(2).EQ.1) AK = AINT(.93334+.01667*TMAX)
      IF (AK.LT.1.) AK = 1.0
      IF (IOPT(2).LT.0 .OR. IOPT(2).GT.1) IER = 66
      IAK = AK
      IF (TMAX.GE.1.0) GO TO 15
      IER = 133
      GO TO 9000
   15 IL = TMAX/AK
      ILL = IL+3
      IIL = IL+2
      J = IL/4+1
C                                  INITIALIZE A TO BLANKS
      DO 20 I=2,IIL
         WK(I) = AA
   20 CONTINUE
      J1 = ILL+1
      J2 = ILL+J
      DO 25 I=J1,J2
         WK(I) = WW
   25 CONTINUE
C                                  OUTPUT PRELIMINARY LINE
      WK(1) = HYP
      WK(ILL) = HYP
      NL = 54
      IF (IOPT(3).NE.0) NL = IOPT(3)
      IF (IOPT(5).NE.0) WRITE (NOUT,80)
      IF (IOPT(5).EQ.0) WRITE (NOUT,110)
      WRITE (NOUT,85) (WK(I),I=J1,J2)
      IJ = N+1
      IP = IOPT1+1
      IOP = IOPT1
      DO 65 IK=1,NN,NL
         IIK = NL
         DO 50 IM=1,IIK,IP
            IJ = IJ-1
            JAK = T(IJ)/AK
            IF (JAK.LT.1) GO TO 35
            DO 30 IN=1,JAK
               WK(IN+1) = CHAR
   30       CONTINUE
   35       JI = 1
            WRITE (NOUT,90) IJ, (WK(K),K=1,ILL)
            DO 40 JJ=2,IIL
               WK(JJ) = AA
   40       CONTINUE
            IF (IJ.EQ.1) GO TO 70
            IF (IOPT1.EQ.0) GO TO 50
            DO 45 JK=1,IOP
               WRITE (NOUT,95) (WK(NI),NI=1,ILL)
   45       CONTINUE
   50    CONTINUE
         IF (IOPT(4).NE.0) GO TO 60
         WRITE (NOUT,55) WW, WW, (WK(I),I=J1,J2)
         WRITE (NOUT,100) (I,I=5,IL,5)
         WRITE (NOUT,105) IAK
         WRITE (NOUT,85) (WK(I),I=J1,J2)
   55    FORMAT (1H , 32A4, A3)
   60    IF (IK.LT.NN .AND. IOPT(4).NE.0) WRITE (NOUT,110)
   65 CONTINUE
   70 WRITE (NOUT,55) WW, WW, (WK(I),I=J1,J2)
      WRITE (NOUT,100) (I,I=5,IL,5)
      WRITE (NOUT,105) IAK
      IF (IER) 9005, 9005, 9000
   75 IER = 132
 9000 CONTINUE
      CALL UERTST (IER,6HUSHHST)
 9005 RETURN
   80 FORMAT (1H0)
   85 FORMAT (1H , 8HCLASS --, 30A4, A3)
   90 FORMAT (1H , I4, 2X, 126A1)
   95 FORMAT (1H , 6X, 126A1)
  100 FORMAT (1H , 9HFREQUENCY, 24(I3, 2X))
  105 FORMAT (1H , 11X, 30HONE FREQUENCY UNIT IS EQUAL TO, I5, 6H COUNT,
     *8H UNIT(S))
  110 FORMAT (1H1)
      END

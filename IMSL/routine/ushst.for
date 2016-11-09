C   IMSL ROUTINE NAME   - USHST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - PRINT A VERTICAL HISTOGRAM
C
C   USAGE               - CALL USHST (T,N,ISP,IER)
C
C   ARGUMENTS    T      - REAL VECTOR OF LENGTH N CONTAINING THE
C                           FREQUENCIES (COUNTS). ELEMENTS MUST BE
C                           NON-NEGATIVE. (INPUT)
C                N      - LENGTH OF T, (NUMBER OF BARS TO PRINT).
C                           IF N EXCEEDS 100/(ISP+1), N=100/(ISP+1) IS
C                           USED.  N MUST BE A POSITIVE INTEGER. (INPUT)
C                ISP    - SPACING BETWEEN HISTOGRAM BARS. ISP MAY EQUAL
C                           0,1, OR 4. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=34 MEANS THAT N*(ISP+1) IS LESS THAN 1
C                             OR GREATER THAN 100.  THE WIDTH OF THE
C                             HISTOGRAM IS SET TO 100, AND 100/(ISP+1)
C                             BARS ARE PRINTED. THE NUMBER OF CLASS
C                             INTERVALS WILL BE PRINTED COMPLETELY IF
C                             ISP.NE.0 AND WILL ALWAYS BE PRINTED UP TO
C                             AND INCLUDING 100/(ISP+1) EVEN THOUGH THE
C                             HISTOGRAM BODY IS ONLY 100 SPACES WIDE.
C                             IF THE CONDITION OF IER=35 ALSO OCCURS,
C                             IER IS SET TO 35.
C                           IER=35 MEANS THAT ISP IS OUT OF ITS RANGE.
C                             ISP=0 IS USED.
C                         TERMINAL ERROR
C                           IER=132 MEANS THAT THE LENGTH OF T IS NOT
C                             POSITIVE.
C                           IER=133 MEANS THAT THE MAXIMUM ELEMENT OF T
C                             IS LESS THAN 1. THE BODY OF THE HISTOGRAM
C                             IS BLANK.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF N*(ISP+1) .LE. 100, THE HISTOGRAM WILL BE PRINTED
C                CORRECTLY.
C            2.  IF THE MAXIMUM FREQUENCY IS .GT. 9999, THE FREQUENCY
C                COLUMN WILL CONTAIN **** ON SOME LINES.
C            3.  OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO. SEE THE UGETIO DOCUMENT FOR DETAILS.
C            4.  EACH HISTOGRAM IS PRINTED ON A NEW PAGE. A USER MIGHT
C                SUPPLY A TITLE FOR A HISTOGRAM BY PRINTING LINES OF
C                TEXT AFTER CALLING USHST . THE LENGTH OF THE TEXT
C                SHOULD BE MIN(N(ISP+1)+12,112).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USHST  (T,N,ISP,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,ISP,IER
      REAL               T(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IDEL,ISP2,IWK,J,J1,J2,L,M,NIN,NM,NN,NOUT
      REAL               AA,AK,CHAR,FREQ(3),HYP,TMAX,WK(129),WW
      DATA               AA/1H /,HYP/1H*/,WW/4H----/,FREQ(1)/4HFREQ/
      DATA               CHAR/1HI/,FREQ(2)/4HUENC/,FREQ(3)/4HY---/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      IER = 0
      ISP2 = ISP
      IF (N.LE.0) GO TO 70
C                                  CHECK HISTOGRAM WIDTH
      M = ISP2+1
      NM = N*M
      L = NM + 4
      IF (L.LE.104.AND.L.GT.0) GO TO 5
      L = 104
      IER = 34
C                                  CHECK BAR SPACING
    5 IF (ISP2.EQ.4) GO TO 15
      IF (ISP2.EQ.1) GO TO 20
      IF (ISP2.EQ.0) GO TO 10
      IER = 35
      M = 1
      L = N+4
      ISP2 = 0
   10 IDEL = 5
      GO TO 25
   15 IDEL = 1
      GO TO 25
   20 IDEL = 2
C                                  FIND MAXIMUM FREQUENCY
   25 TMAX = 0.0
      NN = N
      IF (NM.GT.100) NN = 100/M
      DO 30 I=1,NN
         IF (TMAX.LT.T(I)) TMAX = T(I)
   30 CONTINUE
      J = L/4-1
      IF (TMAX.GE.1.0) GO TO 35
      IER = 133
      GO TO 75
C                                  INITIALIZE A TO BLANKS,
C                                     W TO DASHES (--)
   35 DO 40 I=2,L
   40 WK(I) = AA
      WK(1) = HYP
      J1=L+1
      J2=L+J
      DO 45 I=J1,J2
   45 WK(I) = WW
C                                  OUTPUT PRELIMINARY LINE
      WRITE (NOUT,100) FREQ,(WK(I),I=J1,J2)
      AK = AINT(.98001+.02*TMAX)
      IF (AK.EQ.0.) AK = 1.0
      IWK = TMAX
      IAK = AK
      IF (MOD(IWK,IAK).NE.0) TMAX = TMAX + AMOD(TMAX,AK)
   50 IWK = TMAX
      DO 55 I=1,NN
         IF (T(I) .GE. FLOAT(IWK)) WK(I*M+1) = CHAR
   55 CONTINUE
      WK(NN*M+4) = HYP
      WRITE (NOUT,80) IWK,(WK(I),I=1,L)
      TMAX = TMAX - AK
      IF (TMAX.GE.1.0) GO TO 50
   60 WRITE (NOUT,90) WW,WW,WW,(WK(I),I=J1,J2)
      IF(ISP2.EQ.1) GO TO 65
      WRITE (NOUT,95) (I,I=IDEL,NN,IDEL)
      GO TO 75
   65 WRITE (NOUT,85) (I,I=IDEL,NN,IDEL)
      GO TO 75
   70 IER=132
   75 IF (IER.EQ.0) GO TO 9005
   80 FORMAT (1H ,I4,2X,120A1)
   85 FORMAT (6H CLASS,2X,25I4)
   90 FORMAT (1H ,32A4)
   95 FORMAT (6H CLASS,2X,20I5)
  100 FORMAT (1H1,32A4)
 9000 CONTINUE
      CALL UERTST(IER,6HUSHST )
 9005 RETURN
      END

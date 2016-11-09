C   IMSL ROUTINE NAME   - USHST2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - PRINT A VERTICAL HISTOGRAM, PLOTTING TWO
C                           FREQUENCIES WITH ONE BAR OF THE HISTOGRAM
C
C   USAGE               - CALL USHST2 (T,U,N,ISP,IER)
C
C   ARGUMENTS    T      - REAL VECTOR OF LENGTH N CONTAINING THE
C                           FREQUENCIES (COUNTS). ELEMENTS MUST BE
C                           NON-NEGATIVE. (INPUT)
C                U      - REAL VECTOR OF LENGTH N CONTAINING THE
C                           FREQUENCIES (NEW COUNTS). ELEMENTS MUST BE
C                           NON-NEGATIVE. (INPUT)
C                N      - LENGTH OF T AND U, (NUMBER OF BARS TO PRINT).
C                           IF N EXCEEDS 100/(ISP+1), N=100/(ISP+1) IS
C                           USED.  N MUST BE A POSITIVE INTEGER. (INPUT)
C                ISP    - SPACING BETWEEN HISTOGRAM BARS. (INPUT)
C                           ISP=0,1, OR 4 IS ALLOWED.
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
C                             THE ZERO OPTION IS USED FOR ISP.
C                         TERMINAL ERROR
C                           IER=132 MEANS THAT N IS NOT POSITIVE.
C                           IER=133 MEANS THAT THE MAXIMUM VALUE IN T
C                             AND U IS LESS THAN ONE.  THE BODY OF THE
C                             HISTOGRAM IS BLANK.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF N*(ISP+1) IS .LE. 100, THE HISTOGRAM WILL BE
C                PRINTED CORRECTLY.
C            2.  IF THE MAXIMUM FREQUENCY IS .GT. 9999, THE FREQUENCY
C                COLUMN WILL CONTAIN **** ON SOME LINES.
C            3.  OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO.  SEE THE UGETIO DOCUMENT FOR DETAILS.
C            4.  EACH HISTOGRAM IS PRINTED ON A NEW PAGE. A USER MIGHT
C                SUPPLY A TITLE FOR A HISTOGRAM BY PRINTING LINES OF
C                TEXT AFTER CALLING USHST2. THE LENGTH OF THE TEXT
C                SHOULD BE MIN(N(ISP+1)+12,112).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USHST2 (T,U,N,ISP,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,ISP,IER
      REAL               T(N),U(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IDEL,IM1,ISP2,IWK,J,J1,J2,L,M,NIN,NM,NN,NOUT
      REAL               AA,AK,CHAR,CHM,CHP,FIW,FREQ(3),HYP,TMAX,
     *                   WK(129),WW
      DATA               AA/1H /,HYP/1H*/,WW/4H----/,FREQ(1)/4HFREQ/
      DATA               CHAR/1HI/,FREQ(2)/4HUENC/,FREQ(3)/4HY---/
      DATA               CHP/1H+/,CHM/1H-/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      IER = 0
      IF (N.LE.0) GO TO 75
C                                  ROUND DATA TO NEAREST INTEGER
      DO 5 I=1,N
         T(I) = AINT(T(I)+.5)
         U(I) = AINT(U(I)+.5)
    5 CONTINUE
C                                  INITIALIZE ISP2
      ISP2 = ISP
C                                  CHECK HISTOGRAM WIDTH
      M = ISP2+1
      NM = N*M
      L = NM+4
      IF (L.LE.104 .AND. L.GT.0) GO TO 10
      L = 104
      IER = 34
C                                  CHECK BAR SPACING
   10 IF (ISP2.EQ.4) GO TO 20
      IF (ISP2.EQ.1) GO TO 25
      IF (ISP2.EQ.0) GO TO 15
      IER = 35
      M = 1
      L = N+4
      ISP2 = 0
   15 IDEL = 5
      GO TO 30
   20 IDEL = 1
      GO TO 30
   25 IDEL = 2
C                                  FIND MAXIMUM FREQUENCY
   30 TMAX = 0.0
      NN = N
      IF (NM.GT.100) NN = 100/M
      DO 35 I=1,NN
         IF (TMAX.LT.T(I)) TMAX = T(I)
         IF (TMAX.LT.U(I)) TMAX = U(I)
   35 CONTINUE
      J = L/4-1
      IF (TMAX.GE.1.0) GO TO 40
      IER = 133
      GO TO 65
C                                  INITIALIZE A TO BLANKS,
C                                     W TO DASHES (--)
   40 DO 45 I=2,L
   45 WK(I) = AA
      WK(1) = HYP
      J1 = L+1
      J2 = L+J
      DO 50 I=J1,J2
   50 WK(I) = WW
C                                  OUTPUT PRELIMINARY LINE
      WRITE (NOUT,105) FREQ, (WK(I),I=J1,J2)
      AK = AINT(.98001+.02*TMAX)
      IWK = TMAX
      IAK = AK
      IF (MOD(IWK,IAK).NE.0) TMAX = TMAX + AMOD(TMAX,AK)
   55 IWK = TMAX
      FIW = IWK
      IM1 = 1
      DO 60 I=1,NN
         IM1 = IM1+M
         IF ((T(I).LT.FIW) .AND. (U(I).LT.FIW)) GO TO 60
         WK(IM1) = CHP
         IF (T(I).GT.U(I)) WK(IM1) = CHM
         IF (TMAX.LE.AMIN1(T(I),U(I))) WK(IM1) = CHAR
   60 CONTINUE
      WK(NN*M+4) = HYP
      WRITE (NOUT,85) IWK, (WK(I),I=1,L)
      TMAX = TMAX-AK
      IF (TMAX.GT.0.0) GO TO 55
   65 WRITE (NOUT,95) WW, WW, WW, (WK(I),I=J1,J2)
      IF (ISP2.EQ.1) GO TO 70
      WRITE (NOUT,100) (I,I=IDEL,NN,IDEL)
      GO TO 80
   70 WRITE (NOUT,90) (I,I=IDEL,NN,IDEL)
      GO TO 80
   75 IER = 132
   80 IF (IER.EQ.0) GO TO 9005
   85 FORMAT (1H , I4, 2X, 120A1)
   90 FORMAT (6H CLASS, 2X, 25I4)
   95 FORMAT (1H , 32A4)
  100 FORMAT (6H CLASS, 2X, 20I5)
  105 FORMAT (1H1, 32A4)
 9000 CONTINUE
      CALL UERTST(IER,6HUSHST2)
 9005 RETURN
      END

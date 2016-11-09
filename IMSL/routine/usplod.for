C   IMSL ROUTINE NAME   - USPLOD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINTER PLOT OF UP TO TEN FUNCTIONS
C                           (DOUBLE PRECISION VERSION)
C
C   USAGE               - CALL USPLOD (X,Y,IY,N,M,INC,ITITLE,NTITLE,
C                           IXLABL,NXLABL,IYLABL,NYLABL,RANGE,ICHAR,
C                           IOPT,IER)
C
C   ARGUMENTS    X      - DOUBLE PRECISION VECTOR OF LENGTH N CONTAINING
C                           THE INDEPENDENT VARIABLE VALUES. (INPUT)
C                Y      - DOUBLE PRECISION MATRIX OF DIMENSION N BY M
C                           CONTAINING THE M SETS OF FUNCTION VALUES.
C                           (INPUT)
C                IY     - ROW DIMENSION OF MATRIX Y EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - LENGTH OF THE X VECTOR. (INPUT)
C                M      - NUMBER OF FUNCTIONS TO BE PLOTTED. (INPUT)
C                INC    - DISPLACEMENT BETWEEN ELEMENTS OF THE VECTOR X
C                           TO BE USED. USPLOD PLOTS X(1+(I-1)*INC) FOR
C                           I=1,...N. (INPUT)
C                ITITLE - CHARACTER STRING USED AS THE PLOT TITLE.
C                           (INPUT) THE LENGTH OF ITITLE MUST NOT
C                           EXCEED 72.
C                NTITLE - LENGTH OF ITITLE. (INPUT)
C                           IF NTITLE IS 0 THE TITLE IS LEFT BLANK.
C                IXLABL - CHARACTER STRING USED TO LABEL THE X AXIS.
C                           (INPUT) THE LENGTH OF IXLABL MUST NOT
C                           EXCEED 36.
C                NXLABL - LENGTH OF IXLABL. (INPUT) IF NXLABL IS 0
C                           THE X-AXIS LABEL IS LEFT BLANK.
C                IYLABL - CHARACTER STRING USED TO LABEL THE Y AXIS.
C                           (INPUT) THE LENGTH OF IYLABL MUST NOT
C                           EXCEED 36.
C                NYLABL - LENGTH OF IYLABL. (INPUT) IF NYLABL IS 0
C                           THE Y-AXIS LABEL IS LEFT BLANK.
C                RANGE  - DOUBLE PRECISION VECTOR OF LENGTH 4 SPECIFYING
C                           MIN AND MAX RANGES FOR X,Y AXES. (INPUT)
C                           (MIN X, MAX X, MIN Y, MAX Y, RESPECTIVELY).
C                           USPLOD WILL CALCULATE EACH AXIS RANGE IF THE
C                           MIN AND MAX OF THAT RANGE ARE SET TO 0.0.
C                ICHAR  - CHARACTER STRING OF LENGTH M TO DEFINE
C                           THE SYMBOLS TO BE USED TO PLOT THE M
C                           FUNCTIONS. (INPUT)
C                IOPT   - OPTION INDICATING NUMBER OF PRINTER
C                           COLUMNS AVAILABLE. (INPUT)
C                             IOPT=0, 80 COLUMNS.
C                             IOPT=1, 129 COLUMNS.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           WARNING ERROR (WITH FIX)
C                           IER=65 IMPLIES M IS LESS THAN 1.
C                             ONE FUNCTION IS PLOTTED.
C                           IER=66 IMPLIES M IS GREATER THAN 10.
C                             ONLY 10 FUNCTIONS ARE PLOTTED.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,USPKD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE PLOT TITLE AND AXIS LABELS ARE AUTOMATICALLY
C                CENTERED.
C            2.  FOR MULTIPLE PLOTS, THE CHARACTER M IS USED IN THE
C                EVENT THE SAME PRINT POSITION IS SHARED BY TWO OR MORE
C                FUNCTIONS.
C            3.  OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO. SEE THE UGETIO DOCUMENT FOR DETAILS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USPLOD (X,Y,IY,N,M,INC,ITITLE,NTITLE,IXLABL,NXLABL,
     *                   IYLABL,NYLABL,RANGE,ICHAR,IOPT,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITITLE(1),IXLABL(1),IYLABL(1),ICHAR(1)
      INTEGER            IY,N,M,INC,NTITLE,NXLABL,NYLABL,IOPT,IER
      DOUBLE PRECISION   X(1),Y(IY,1),RANGE(4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            BLNK4,I,I1,II,IPSW,IT,ITEMP,IWORK(101),IX,IXA,
     *                   IXD,IXDP1,IXS,IYA,IYD,IYDP1,IYS,IYY,J,JI,JM1,
     *                   JTEMP,JTITLE(144),K,KTEMP,L,LL,MID,MP2,M2,NIN,
     *                   NCHMTB,NLEN,NLOG,NMIN,NOUT,PLUS,STAR,SYM4(10)
      DOUBLE PRECISION   A,H,HTRY(6),SFX,SFY,TEMP,X1,X2,XLEN,XMAX,XMIN,
     *                   XNOT(6),Y1,Y2,YMAX,YMIN,YNOT
      DATA               MP2/1HM/,IYS/1HI/,IXS/1H-/,STAR/1H*/,PLUS/1H+/
      DATA               HTRY /1.0D0,2.0D0,4.0D0,5.0D0,8.0D0,10.0D0/
      DATA               BLNK4 /1H /
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      DO 5 I=1,144
    5 JTITLE(I) = BLNK4
C                                  INITIALIZE LOCAL VARIABLE M2 = M
      M2 = M
      NLEN = MIN0(10,M2)
      NLEN = MAX0(1,NLEN)
      CALL USPKD(ICHAR,NLEN,SYM4,NCHMTB)
      NLEN = MIN0(72,NTITLE)
      IF (NTITLE.LT.0) NLEN = 0
      MID = 37-NLEN/2
      CALL USPKD(ITITLE,NLEN,JTITLE(MID),NCHMTB)
      NLEN = MIN0(36,NXLABL)
      IF (NXLABL.LT.0) NLEN = 0
      MID = 91-NLEN/2
      CALL USPKD(IXLABL,NLEN,JTITLE(MID),NCHMTB)
      NLEN = MIN0(36,NYLABL)
      IF (NYLABL.LT.0) NLEN = 0
      MID = 127-NLEN/2
      CALL USPKD(IYLABL,NLEN,JTITLE(MID),NCHMTB)
      IYD = 50
      IXD = 100
      IF (IOPT.EQ.0) IXD = 60
      IER = 0
      IF (M2.GT.0) GO TO 10
C                                  SET WARNING ERROR
      IER = 65
      M2 = 1
      GO TO 15
   10 IF (M2.LE.10) GO TO 15
C                                  SET WARNING ERROR
      IER = 66
      M2 = 10
C                                  PRINT TITLE
   15 IF (IOPT.EQ.0) WRITE (NOUT,20) (JTITLE(I),I=1,72)
      IF (IOPT.NE.0) WRITE (NOUT,25) (JTITLE(I),I=1,72)
   20 FORMAT (1H1, 4X, 72A1//)
   25 FORMAT (1H1, 27X, 72A1//)
      IXDP1 = IXD+1
      IYDP1 = IYD+1
C                                  FIND MIN,MAX OF X
      IF (RANGE(1).EQ.RANGE(2)) GO TO 30
      XMIN = RANGE(1)
      XMAX = RANGE(2)
      GO TO 60
   30 XMIN = X(1)
      XMAX = X(1)
      DO 40 I=1,N,INC
         IF (X(I).GE.XMIN) GO TO 35
         XMIN = X(I)
         GO TO 40
   35    IF (X(I).GT.XMAX) XMAX = X(I)
   40 CONTINUE
      IF (XMIN.NE.XMAX) GO TO 45
      XMIN = XMIN-1.0D0
      XMAX = XMAX+1.0D0
C
   45 XLEN = IXD/10.0D0
      H = (XMAX-XMIN)/(XLEN-1.0D0)
      TEMP = DLOG10(H)
      NLOG = TEMP
      IF (TEMP.LT.0.0D0) NLOG = NLOG-1
      DO 50 I=1,6
         H = HTRY(I)*10.0D0**NLOG
         TEMP = XMIN/H+0.01D0
         NMIN = TEMP
         IF (TEMP.LT.0.0D0) NMIN = NMIN-1
         A = NMIN*H
         IF ((XMAX-A).LE.(XLEN+0.01D0)*H) GO TO 55
   50 CONTINUE
   55 XMIN = A
      XMAX = A+XLEN*H
C                                  FIND MIN,MAX OF Y
   60 CONTINUE
      IF (RANGE(3).EQ.RANGE(4)) GO TO 65
      YMIN = RANGE(3)
      YMAX = RANGE(4)
      GO TO 100
   65 YMIN = Y(1,1)
      YMAX = Y(1,1)
      DO 80 I=1,M2
         Y1 = Y(1,I)
         Y2 = Y(1,I)
         DO 75 J=1,N,INC
            IF (Y(J,I).GE.Y1) GO TO 70
            Y1 = Y(J,I)
            GO TO 75
   70       IF (Y(J,I).GT.Y2) Y2 = Y(J,I)
   75    CONTINUE
         IF (Y1.LT.YMIN) YMIN = Y1
         IF (Y2.GT.YMAX) YMAX = Y2
   80 CONTINUE
      IF (YMIN.NE.YMAX) GO TO 85
      YMIN = YMIN-1.0D0
      YMAX = YMAX+1.0D0
   85 H = (YMAX-YMIN)/9.0D0
      TEMP = DLOG10(H)
      NLOG = TEMP
      IF (TEMP.LT.0.0D0) NLOG = NLOG-1
      DO 90 I=1,6
         H = HTRY(I)*10.0D0**NLOG
         TEMP = YMIN/H+0.01D0
         NMIN = TEMP
         IF (TEMP.LT.0.0D0) NMIN = NMIN-1
         A = NMIN*H
         IF ((YMAX-A).LE.10.01D0*H) GO TO 95
   90 CONTINUE
   95 YMIN = A
      YMAX = A+10.0D0*H
C                                  SCALE FACTOR FOR X-AXIS
  100 SFX = IXD/(XMAX-XMIN)
C                                  SCALE FACTOR FOR Y-AXIS
      SFY = IYD/(YMIN-YMAX)
C                                  SHIFT FACTOR FOR X
      IXA = 1.5D0-SFX*XMIN
C                                  SHIFT FACTOR FOR Y
      IYA = 1.5D0-SFY*YMAX
      DO 240 J=1,51
C                                  CLEAR PRINT IMAGE
         DO 105 I=1,IXDP1
  105    IWORK(I) = BLNK4
C                                  HORIZONTAL BORDERS
         IF (J.NE.1 .AND. J.NE.51) GO TO 120
         DO 115 I=1,IXDP1
            IF (MOD(I-1,10).NE.0) GO TO 110
            IWORK(I) = PLUS
            GO TO 115
  110       IWORK(I) = STAR
  115    CONTINUE
  120    CONTINUE
C                                  VERTICAL   BORDERS
         IF (J.EQ.1 .OR. J.EQ.51) GO TO 130
         IF (MOD(J-1,5).NE.0) GO TO 125
         IWORK(1) = PLUS
         IWORK(IXDP1) = PLUS
         GO TO 130
  125    IWORK(1) = STAR
         IWORK(IXDP1) = STAR
  130    CONTINUE
C                                  Y-AXIS
         IF (IXA.LE.1 .OR. IXA.GT.IXD) GO TO 135
         IF (J.EQ.1 .OR. J.EQ.51) GO TO 135
         IWORK(IXA) = IYS
C                                  X-AXIS
  135    IF (IYA.LE.1 .OR. IYA.GT.IYD) GO TO 145
         IF (J.NE.IYA) GO TO 145
         DO 140 I=2,IXD
            IWORK(I) = IXS
  140    CONTINUE
C                                  NOW COMPUTE PRINT IMAGE COORDINATES
  145    CONTINUE
         DO 180 I=1,N,INC
            ITEMP = I
            IX = 1.5D0+SFX*(X(ITEMP)-XMIN)
            IF (IX.GE.1 .AND. IX.LE.IXDP1) GO TO 150
            GO TO 180
  150       CONTINUE
            DO 175 JI=1,M2
               JTEMP = JI
               IYY = 1.5D0+SFY*(Y(ITEMP,JTEMP)-YMAX)
               IF (IYY.EQ.J) GO TO 155
               GO TO 175
  155          IF (JTEMP.EQ.1) GO TO 170
               IT = IWORK(IX)
               IF (IT.EQ.MP2) GO TO 175
               JM1 = JTEMP-1
               DO 160 K=1,JM1
                  KTEMP = K
                  IF (IT.NE.SYM4(KTEMP)) GO TO 160
                  IF (IT.NE.SYM4(JTEMP)) GO TO 165
  160          CONTINUE
               GO TO 170
C                                  IDENTICAL POINTS BECOME M
  165          IWORK(IX) = MP2
               GO TO 175
  170          IWORK(IX) = SYM4(JTEMP)
  175       CONTINUE
  180    CONTINUE
         Y1 = 1.0D0/SFY
         Y2 = YMAX-Y1
         IPSW = 0
         IF (MOD(J-1,5).NE.0) GO TO 185
         IPSW = IPSW+1
         YNOT = J*Y1+Y2
         IF (DABS(YNOT).LT.DABS(Y1)) YNOT = 0.0D0
  185    IF (J.LT.8 .OR. J.GT.43) GO TO 190
         IPSW = IPSW+2
         I1 = J+101
  190    IF (IPSW.NE.0) GO TO 200
         WRITE (NOUT,195) (IWORK(II),II=1,IXDP1)
  195    FORMAT (14X, 101A1)
         GO TO 235
  200    IF (IPSW-2) 205, 215, 225
  205    WRITE (NOUT,210) YNOT, (IWORK(II),II=1,IXDP1)
  210    FORMAT (3X, D10.2, 1X, 101A1)
         GO TO 235
  215    WRITE (NOUT,220) JTITLE(I1), (IWORK(II),II=1,IXDP1)
  220    FORMAT (1X, A1, 12X, 101A1)
         GO TO 235
  225    WRITE (NOUT,230) JTITLE(I1), YNOT, (IWORK(II),II=1,IXDP1)
  230    FORMAT (1X, A1, 1X, D10.2, 1X, 101A1)
  235    CONTINUE
  240 CONTINUE
      X1 = 1.0D0/SFX
      X2 = XMIN-X1
      L = 1
      DO 245 I=1,IXDP1,20
         XNOT(L) = I*X1+X2
         IF (DABS(XNOT(L)).LT.X1) XNOT(L) = 0.0D0
         L = L+1
  245 CONTINUE
C                                  PRINT X-AXIS VALUES
      LL = L-1
      WRITE (NOUT,250) (XNOT(I),I=1,LL)
  250 FORMAT (9X, D10.2, 5(10X, D10.2))
      IF (IOPT.EQ.0) WRITE (NOUT,255) (JTITLE(I),I=73,108)
      IF (IOPT.NE.0) WRITE (NOUT,260) (JTITLE(I),I=73,108)
  255 FORMAT (//5X, 36(A1, 1X))
  260 FORMAT (//28X, 36(A1, 1X))
      IF (IER.EQ.0) GO TO 9005
      CALL UERTST(IER,6HUSPLOD)
 9005 RETURN
      END

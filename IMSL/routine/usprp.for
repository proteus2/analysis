C   IMSL ROUTINE NAME   - USPRP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - PROBABILITY PLOT
C
C   USAGE               - CALL USPRP (X,N,N1,N2,IDIST,IOPT,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH N2-N1+1 CONTAINING THE DATA.
C                           THE DATA SET IS POSSIBLY A CENSORED ONE
C                           FROM A COMPLETE SAMPLE OF SIZE N. (INPUT)
C                N      - TOTAL NUMBER OF OBSERVATIONS IN UNCENSORED
C                           SAMPLE.  N MUST BE GREATER THAN OR EQUAL
C                           TO N2-N1+1.  IF THERE HAS BEEN NO CENSORING
C                           THEN  N1=1 AND N2=N. (INPUT)
C                N1     - ON INPUT, N1 IS THE RANK NUMBER OF THE
C                           SMALLEST OBSERVATION IN X IF RANKED IN
C                           THE COMPLETE SAMPLE.  (THE NUMBER OF
C                           OBSERVATIONS THAT HAVE BEEN CENSORED
C                           FROM BELOW IS N1-1.)
C                           IF, BECAUSE PROPERTIES OF SAMPLE IN X
C                           DO NOT MATCH THOSE OF THE DISTRIBUTION
C                           SPECIFIED, IT IS NECESSARY FOR USPRP TO
C                           DELETE THE SMALLEST K ITEMS IN X, N1 ON
C                           OUTPUT IS K PLUS THE INPUT VALUE OF N1.
C                           (INPUT/OUTPUT)
C                N2     - ON INPUT, N2 IS THE RANK NUMBER OF THE
C                           LARGEST OBSERVATION IN X IF RANKED IN
C                           THE COMPLETE SAMPLE.  (THE NUMBER OF
C                           OBSERVATIONS THAT HAVE BEEN CENSORED
C                           FROM ABOVE IS N-N2.)
C                           IF, BECAUSE PROPERTIES OF SAMPLE IN X
C                           DO NOT MATCH THOSE OF THE DISTRIBUTION
C                           SPECIFIED, IT IS NECESSARY FOR USPRP TO
C                           DELETE THE LARGEST L ITEMS IN X, N2 ON
C                           OUTPUT IS L PLUS THE INPUT VALUE OF N2.
C                           (INPUT/OUTPUT)
C                IDIST  - PARAMETER TO INDICATE DISTRIBUTION. (INPUT)
C                           IDIST=1, NORMAL DISTRIBUTION.
C                           IDIST=2, LOGNORMAL DISTRIBUTION.
C                           IDIST=3, HALF-NORMAL DISTRIBUTION.
C                           IDIST=4, EXPONENTIAL DISTRIBUTION.
C                           IDIST=5, WEIBULL DISTRIBUTION.
C                           IDIST=6, EXTREME VALUE DISTRIBUTION.
C                IOPT   - OPTION INDICATING NUMBER OF PRINTER COLUMNS
C                           AVAILABLE. (INPUT)
C                           IOPT=0, 80 COLUMNS.
C                           IOPT=1, 129 COLUMNS.
C                WK     - WORK VECTOR OF LENGTH 2N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING WITH FIX ERROR
C                           IER = 67 INDICATES THAT IT WAS NECESSARY
C                             TO DELETE SOME ITEMS FROM THE PLOTTING
C                             BECAUSE THOSE ITEMS DID NOT SATISFY
C                             PROPERTIES OF THE DISTRIBUTION.
C                         TERMINAL ERROR
C                           IER = 131 INDICATES THAT N1 OR N2 IS
C                             SPECIFIED INCORRECTLY.
C                           IER = 132 INDICATES THAT THE SAMPLE SIZE
C                             IS LESS THAN 2.
C                           IER = 133 INDICATES IDIST IS SPECIFIED
C                             INCORRECTLY.
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO,VSRTA
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO. SEE THE UGETIO DOCUMENT FOR DETAILS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USPRP(X,N,N1,N2,IDIST,IOPT,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,N1,N2,IDIST,IOPT,IER
      REAL               X(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            BLNK,I,I1,IDCML,II,IPCT(18),IPSW,ITITLE(12),
     *                   IWORK(101),IX,IXD,IXDP1,IXLAB(9),IXS,IYA,IYD,
     *                   IYDP1,IYY,J,KER,LOGY,MARK(9),MINUS,N11,NDTA,
     *                   NDTAM1,NIN,NOUT,NP1,NS,PLUS,STAR,SYMB
      REAL               ALG,DENOM,P,SFX,SFY,TRI,WKIPN,XDEN,XLAB(9),
     *                   XMAX,XMIN,XN,XNOT(9),Y1,Y2,YDIF,YLAB(9),YLOG,
     *                   YMAX,YMIN,YNOT
      REAL               HTRY(6)
      DATA               HTRY/1.,2.,4.,5.,8.,10./
      DATA               STAR/1H*/,PLUS/1H+/,BLNK/1H /,SYMB/1HX/
      DATA               XNOT/.01,.05,.10,.25,.50,.75,.90,.95,.99/
      DATA               IPCT/1H0,1H1,1H0,1H5,1H1,1H0,1H2,1H5,1H5,1H0,
     *                   1H7,1H5,1H9,1H0,1H9,1H5,1H9,1H9/
      DATA               IDCML/1H./
      DATA               ITITLE/1HO,1HB,1HS,1HE,1HR,1HV,1HA,1HT,1HI,1HO,
     *                   1HN,1HS/,MINUS/1H-/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      CALL UGETIO(1,NIN,NOUT)
      IF (N.GE.N2.AND.N2.GE.N1.AND.N1.GE.1) GO TO 5
      IER = 131
      GO TO 9000
    5 NS = N2 - N1 + 1
      IF (NS.GT.1.AND.N.GT.1) GO TO 10
      IER = 132
      GO TO 9000
   10 XN = N
      NP1 = N + 1
      NDTA = N + N1
      NDTAM1 = NDTA - 1
C                                  LOGY INDICATES WHETHER THE VERTICAL
C                                  AXIS HAS BEEN TRANSFORMED.
      LOGY = 0
      IF (IDIST.EQ.2 .OR. IDIST.EQ.5) LOGY = 1
C                                  TRANSFER DATA AND SORT
      DO 15 I = 1, NS
        WK(I+NDTAM1) = X(I)
   15 CONTINUE
      CALL VSRTA(WK(NDTA),NS)
      IF (IDIST.NE.1) GO TO 35
C                                  NORMAL PLOT OPTION
      XDEN = 3.0*XN + 1.0
      DO 20 I = N1, N2
        P = (3.0*I-1.0)/XDEN
        CALL MDNRIS(P,WK(I),KER)
   20 CONTINUE
      DO 25 I = 1, 9
        CALL MDNRIS(XNOT(I),XLAB(I),KER)
   25 CONTINUE
      IF (N1.NE.1 .OR. N2.NE.N) GO TO 30
      WRITE (NOUT,360)
      GO TO 160
   30 WRITE (NOUT,390)
      GO TO 160
   35 IF (IDIST.NE.4) GO TO 55
C                                  EXPONENTIAL PLOT OPTION
      DO 40 I = N1, N2
        P = (I-.5)/XN
        WK(I) = -ALOG(1.-P)
   40 CONTINUE
      DO 45 I = 1, 9
        XLAB(I) = -ALOG(1.-XNOT(I))
   45 CONTINUE
      IF (N1.NE.1 .OR. N2.NE.N) GO TO 50
      WRITE (NOUT,365)
      GO TO 160
   50 WRITE (NOUT,395)
      GO TO 160
   55 IF (IDIST.NE.6) GO TO 75
C                                  EXTREME VALUE PLOT OPTION
      DO 60 I = N1, N2
        P = (I-.5)/XN
        ALG = -ALOG(1.0-P)
        WK(I) = ALOG(ALG)
   60 CONTINUE
      DO 65 I = 1, 9
        ALG = -ALOG(1.0-XNOT(I))
        XLAB(I) = ALOG(ALG)
   65 CONTINUE
      IF (N1.NE.1 .OR. N2.NE.N) GO TO 70
      WRITE (NOUT,370)
      GO TO 160
   70 WRITE (NOUT,400)
      GO TO 160
   75 IF (IDIST.NE.2) GO TO 105
C                                  LOGNORMAL PLOT OPTION
      XDEN = 3.0*XN + 1.0
      N11 = N1
      DO 85 I = N1, N2
        P = (3*I-1.0)/XDEN
        CALL MDNRIS(P,WK(I),KER)
        WKIPN = WK(I+N)
        IF (WKIPN.LE.0.0) GO TO 80
        WK(I+N) = ALOG(WKIPN)
        GO TO 85
   80   N11 = N11 + 1
   85 CONTINUE
      DO 90 I = 1, 9
        CALL MDNRIS(XNOT(I),XLAB(I),KER)
   90 CONTINUE
      IF (N1.EQ.N11) GO TO 95
      IER = 67
      N1 = N11
      NS = N2 - N1 + 1
      IF (NS.GT.1) GO TO 95
      IER = 132
      GO TO 9000
   95 IF (N1.NE.1 .OR. N2.NE.N) GO TO 100
      WRITE (NOUT,375)
      GO TO 160
  100 WRITE (NOUT,405)
      GO TO 160
  105 IF (IDIST.NE.5) GO TO 135
C                                  WEIBULL PLOT OPTION
      N11 = N1
      DO 115 I = N1, N2
        P = (I-.5)/XN
        ALG = -ALOG(1.0-P)
        WK(I) = ALOG(ALG)
        WKIPN = WK(I+N)
        IF (WKIPN.LE.0.0) GO TO 110
        WK(I+N) = ALOG(WKIPN)
        GO TO 115
  110   N11 = N11 + 1
  115 CONTINUE
      DO 120 I = 1, 9
        ALG = -ALOG(1.0-XNOT(I))
        XLAB(I) = ALOG(ALG)
  120 CONTINUE
      IF (N1.EQ.N11) GO TO 125
      IER = 67
      N1 = N11
      NS = N2 - N1 + 1
      IF (NS.GT.1) GO TO 125
      IER = 132
      GO TO 9000
  125 IF (N1.NE.1 .OR. N2.NE.N) GO TO 130
      WRITE (NOUT,380)
      GO TO 160
  130 WRITE (NOUT,410)
      GO TO 160
  135 IF (IDIST.EQ.3) GO TO 140
      IER = 133
      GO TO 9000
C                                  HALF-NORMAL PLOT OPTION
  140 XDEN = 6.0*XN + 1.0
      TRI = 3.0*XN - 1.0
      DO 145 I = N1, N2
        P = (3*I+TRI)/XDEN
        CALL MDNRIS(P,WK(I),KER)
  145 CONTINUE
      DO 150 I = 1, 9
        P = 0.5 + XNOT(I)/2.0
        CALL MDNRIS(P,XLAB(I),KER)
  150 CONTINUE
      IF (N1.NE.1 .OR. N2.NE.N) GO TO 155
      WRITE (NOUT,385)
      GO TO 160
  155 WRITE (NOUT,415)
C
C                                  DETERMINE SCALE AND SHIFT FACTORS
C
  160 CONTINUE
      IYD = 50
      IXD = 100
      IF (IOPT.EQ.0) IXD = 60
      IXDP1 = IXD + 1
      IYDP1 = IYD + 1
C                                  FIND MIN,MAX FOR X-AXIS
C                                  (N.B. THE VALUES OF THE INPUT
C                                   VECTOR X GO ON THE Y-AXIS)
      XMIN = WK(N1)
      XMAX = WK(N2)
      IF (XMIN.GT.XLAB(1)) XMIN = XLAB(1)
      IF (XMAX.LT.XLAB(9)) XMAX = XLAB(9)
C                                  FIND MIN,MAX FOR Y-AXIS
      YMIN = WK(N+N1)
      YMAX = WK(N+N2)
C                                  SCALE FACTOR FOR X-AXIS
      SFX = IXD/(XMAX-XMIN)
C                                  SCALE FACTOR FOR Y-AXIS
      IF (YMIN.NE.YMAX) GO TO 165
      YMIN = YMIN - 1.0
      YMAX = YMAX + 1.0
  165 IF (LOGY.EQ.1) GO TO 180
C                                  TRY TO GET *PRETTY* NUMBERS FOR
C                                  LABELS.
      H = (YMAX-YMIN)/9.0
      TEMP = ALOG10(H)
      NLOG = TEMP
      IF (TEMP.LT.0.0) NLOG = NLOG - 1
      DO 170 I = 1, 6
        H = HTRY(I)*10.0**NLOG
        TEMP = YMIN/H + 0.01
        NMIN = TEMP
        IF (TEMP.LT.0.0) NMIN = NMIN - 1
        A = NMIN*H
        IF ((YMAX-A).LE.10.01*H) GO TO 175
  170 CONTINUE
  175 YMIN = A
      YMAX = A + 10.0*H
      SFY = IYD/(YMIN-YMAX)
      GO TO 185
C                                  A LOGRITHMIC SCALE IS USED FOR THE
C                                  VERTICAL AXIS.  DO NOT TRY TO GET
C                                  *PRETTY* NUMBERS FOR LABELS.
  180 SFY = IYD/(YMIN-YMAX)
C                                  ORIGIN POINT FOR Y-AXIS
  185 IYA = 1.5 - SFY*YMAX
C                                  FIX LABEL OF VERTICAL AXIS FOR
C                                  LOG PLOTS.
      IF (LOGY.EQ.0) GO TO 195
      YLAB(1) = EXP(WK(N+N1))
      YLAB(9) = EXP(WK(N+N2))
      YLAB(8) = (YLAB(9)-YLAB(1))/2. + YLAB(1)
      YLAB(6) = (YLAB(8)-YLAB(1))/2. + YLAB(1)
      YLAB(7) = (YLAB(8)-YLAB(6))/2. + YLAB(6)
      YDIF = (YLAB(6)-YLAB(1))/4.
      YLAB(5) = YLAB(6) - YDIF
      YLAB(4) = YLAB(5) - YDIF
      YLAB(3) = YLAB(4) - YDIF
      YLAB(2) = (YLAB(3)-YLAB(1))/2. + YLAB(1)
      MARK(1) = 51
      MARK(9) = 1
      DENOM = WK(N+N2) - WK(N1+N)
      DO 190 I = 2, 8
        MARK(I) = (WK(N+N2)-ALOG(YLAB(I)))/DENOM*51
  190 CONTINUE
  195 CONTINUE
C
C
      DO 305 J = 1, 51
C                                  VERTICAL   BORDERS
        IWORK(1) = STAR
        IWORK(IXDP1) = STAR
        IF (J.EQ.IYA.AND.LOGY.NE.1) GO TO 215
        IF (J.NE.1.AND.J.NE.51) GO TO 205
C                                  FIX HORIZONTAL BORDERS
        DO 200 I = 2, IXD
          IWORK(I) = STAR
  200   CONTINUE
        GO TO 225
C                                  CLEAR PRINT IMAGE
  205   DO 210 I = 2, IXD
          IWORK(I) = BLNK
  210   CONTINUE
        GO TO 235
C                                  HORIZONTAL AXIS FOR Y=0
  215   DO 220 I = 2, IXD
          IWORK(I) = MINUS
  220   CONTINUE
C                                  MARK NOTATIONS ON BORDERS
  225   IF (J.NE.51) GO TO 235
        DO 230 I = 1, 9
          IX = 1.5 + SFX*(XLAB(I)-XMIN)
          IWORK(IX) = PLUS
          IXLAB(I) = IX
  230   CONTINUE
  235   CONTINUE
C                                  NOW COMPUTE PRINT IMAGE COORDINATES
        DO 250 I = N1, N2
          IX = 1.5 + SFX*(WK(I)-XMIN)
          IF (IX.GE.1.AND.IX.LE.IXDP1) GO TO 240
          GO TO 250
  240     CONTINUE
          IYY = 1.5 + SFY*(WK(I+N)-YMAX)
          IF (IYY.EQ.J) GO TO 245
          GO TO 250
  245     IWORK(IX) = SYMB
  250   CONTINUE
        Y1 = 1./SFY
        Y2 = YMAX - Y1
        IPSW = 0
        IF (LOGY.EQ.0) GO TO 260
        Y1 = 1./(IYD/(WK(N+N1)-WK(N+N2)))
        Y2 = WK(N+N2) - Y1
        DO 255 I = 1, 9
          IF (J.EQ.MARK(I)) GO TO 265
  255   CONTINUE
        GO TO 275
  260   IF (MOD(J-1,5).NE.0) GO TO 275
  265   IPSW = IPSW + 1
        IWORK(1) = PLUS
        IF (LOGY.EQ.1) GO TO 270
        YNOT = J*Y1 + Y2
        IF (ABS(YNOT).LT.ABS(Y1)) YNOT = 0.0
        GO TO 275
  270   YNOT = YLAB(I)
  275   IF (J.LT.20 .OR. J.GT.31) GO TO 280
        IPSW = IPSW + 2
        I1 = J - 19
  280   IF (IPSW.NE.0) GO TO 285
        WRITE (NOUT,335)(IWORK(II),II=1,IXDP1)
        GO TO 305
  285   IF (IPSW-2) 290, 295, 300
  290   WRITE (NOUT,340) YNOT, (IWORK(II),II=1,IXDP1)
        GO TO 305
  295   WRITE (NOUT,345) ITITLE(I1), (IWORK(II),II=1,IXDP1)
        GO TO 305
  300   WRITE (NOUT,350) ITITLE(I1), YNOT, (IWORK(II),II=1,IXDP1)
C
  305 CONTINUE
C                                  INITIALIZE AND PRINT X-AXIS VALUES
      DO 310 I = 1, IXDP1
        IWORK(I) = BLNK
  310 CONTINUE
C                                  DETERMINE WHERE VALUES ARE TO BE
C                                  PRINTED
      II = IXLAB(1)
      DO 320 I = 2, 8
        IF (II+3.LE.IXLAB(I)) GO TO 315
        IXLAB(I) = 0
        GO TO 320
  315   II = IXLAB(I)
  320 CONTINUE
      IF (IXLAB(8)+3.LE.IXLAB(9)) GO TO 325
      IXLAB(8) = 0
      IF (IXLAB(7)+3.LE.IXLAB(9)) GO TO 325
      IXLAB(7) = 0
  325 CONTINUE
      DO 330 I = 1, 9
        II = IXLAB(I)
        IF (II.LE.0 .OR. II.GE.IXD-1) GO TO 330
        IWORK(II) = IDCML
        IWORK(II+1) = IPCT(2*I-1)
        IWORK(II+2) = IPCT(2*I)
  330 CONTINUE
      WRITE (NOUT,335)(IWORK(II),II=1,IXDP1)
      WRITE (NOUT,355)
      IF (IER.EQ.0) GO TO 9005
C                                  FORMATS
  335 FORMAT (14X,101A1)
  340 FORMAT (3X,E10.2,1X,101A1)
  345 FORMAT (1X,A1,12X,101A1)
  350 FORMAT (1X,A1,1X,E10.2,1X,101A1)
  355 FORMAT (/14X,22HCUMULATIVE PROBABILITY)
C                                  FORMATS FOR HEADING
  360 FORMAT (1H1,13X,40HPROBABILITY PLOT FOR NORMAL DISTRIBUTION)
  365 FORMAT (1H1,13X,44HPROBABILITY PLOT FOR EXPONENTIAL DISTRIBUTIO,
     *  1HN)
  370 FORMAT (1H1,13X,44HPROBABILITY PLOT FOR EXTREME VALUE DISTRIBUT,
     *  3HION)
  375 FORMAT (1H1,13X,43HPROBABILITY PLOT FOR LOGNORMAL DISTRIBUTION)
  380 FORMAT (1H1,13X,41HPROBABILITY PLOT FOR WEIBULL DISTRIBUTION)
  385 FORMAT (1H1,13X,44HPROBABILITY PLOT FOR HALF-NORMAL DISTRIBUTIO,
     *  1HN)
  390 FORMAT (1H1,13X,40HPROBABILITY PLOT FOR NORMAL DISTRIBUTION,
     *  11H (CENSORED))
  395 FORMAT (1H1,13X,44HPROBABILITY PLOT FOR EXPONENTIAL DISTRIBUTIO,
     *  1HN,11H (CENSORED))
  400 FORMAT (1H1,13X,44HPROBABILITY PLOT FOR EXTREME VALUE DISTRIBUT,
     *  3HION,11H (CENSORED))
  405 FORMAT (1H1,13X,43HPROBABILITY PLOT FOR LOGNORMAL DISTRIBUTION,
     *  11H (CENSORED))
  410 FORMAT (1H1,13X,41HPROBABILITY PLOT FOR WEIBULL DISTRIBUTION,
     *  11H (CENSORED))
  415 FORMAT (1H1,13X,44HPROBABILITY PLOT FOR HALF-NORMAL DISTRIBUTIO,
     *  1HN,11H (CENSORED))
 9000 CONTINUE
      CALL UERTST(IER,6HUSPRP )
 9005 RETURN
      END

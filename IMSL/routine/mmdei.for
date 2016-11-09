C   IMSL ROUTINE NAME   - MMDEI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - EXPONENTIAL INTEGRALS
C
C   USAGE               - FUNCTION MMDEI (IOPT,ARG,IER)
C
C   ARGUMENTS    MMDEI  - OUTPUT VALUE OF THE INTEGRAL. MMDEI MUST BE
C                           TYPED APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE PRECISION/HARDWARE SECTION.)
C                IOPT   - INPUT OPTION.
C                         FOR IOPT = 1, THE INTEGRAL (FROM -INFINITY TO
C                           ARG) OF EXP(T)/T DT WILL BE EVALUATED IF ARG
C                           IS GREATER THAN 0. IF ARG IS LESS THAN 0.0,
C                           (-1)*THE INTEGRAL (FROM -ARG TO INFINITY)
C                           OF EXP(-T)/T DT WILL BE EVALUATED.
C                         FOR IOPT = 2, THE INTEGRAL (FROM ARG TO
C                           INFINITY) OF EXP(-T)/T DT WILL BE EVALUATED.
C                           ARG MUST BE GREATER THAN 0.
C                         FOR IOPT = 3, EXP(-ARG)*THE INTEGRAL (FROM
C                           -INFINITY TO ARG) OF EXP(T)/T DT WILL BE
C                           EVALUATED IF ARG IS GREATER THAN 0. IF ARG
C                           IS LESS THAN 0.0, EXP(-ARG)*(-1)*THE
C                           INTEGRAL (FROM -ARG TO INFINITY) OF
C                           EXP(-T)/T DT WILL BE EVALUATED.
C                ARG    - INPUT PARAMETER. SEE IOPT DESCRIPTION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IOPT WAS LESS THAN
C                             1 OR GREATER THAN 3. MMDEI IS SET TO
C                             MACHINE INFINITY.
C                           IER = 130 INDICATES THAT ARG WAS EQUAL TO
C                             0.0. MMDEI IS SET TO MACHINE INFINITY IF
C                             IOPT = 2 AND NEGATIVE MACHINE INFINITY IF
C                             IOPT = 1 OR 3.
C                           IER = 131 INDICATES THAT AN OVERFLOW WOULD
C                             HAVE OCCURRED. IF IOPT = 1, MMDEI IS SET
C                             TO MACHINE INFINITY.
C                           IER = 132 INDICATES THAT AN UNDERFLOW WOULD
C                             HAVE OCCURRED. MMDEI IS SET TO 0.0 IF
C                             IOPT = 1 OR 2.
C                         WARNING WITH FIX
C                           IER = 69 INDICATES THAT ARG WAS NEGATIVE
C                             FOR IOPT = 2. CALCULATION CONTINUES USING
C                             ABS(ARG).
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
      DOUBLE PRECISION FUNCTION MMDEI (IOPT,ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      DOUBLE PRECISION   A(6),B(6),C(8),D(8),E(8),F(8),P1(9),Q1(9),
     1                   P2(9),Q2(8),P3(10),Q3(9),P4(10),Q4(9),P0(6),
     2                   Q0(6),PX(9),QX(9),FRAC,SUMP,SUMQ,T,W,X,X0,
     3                   XX0,X01,X02,XMX0,Y,DEXP40,XINF,XMAX,XMIN,
     4                   SIX,TWELVE,TWO,THREE,ZERO,ONE,HALF,TWENT4,
     5                   FOUR,FORTY
      DOUBLE PRECISION   SUMB,SUMC,SUMD,SUME,SUMF
      DATA               DEXP40/Z4F34441A72F2E5D5/
      DATA               XINF/Z7FFFFFFFFFFFFFFF/,
     1                   XMAX/Z42B3DC12A2B5787F/,
     2                   XMIN/ZC2AF0C2BE49BC591/
C
C                                  LET X = ARG
C                                  ZERO OF EI(X)
C
      DATA               X0/Z405F5CA54AD2D7F1/,
     1                   X01/Z405F5CA54AD20000/,
     2                   X02/Z36D7F0F264C30000/
C
C                                  COEFFICIENTS FOR R(5,5) APPROXIMATION
C                                  USED FOR X GREATER THAN OR EQUAL TO
C                                  -1 AND LESS THAN 0
C
      DATA               A(1)/ZC093C467E37DB0C8/,A(2)/Z40C110E996FE3178/
      DATA               A(3)/Z40213DCCA0D570ED/,A(4)/Z3F629544F2457AA6/
      DATA               A(5)/Z3E569011E23A2DC8/,A(6)/Z3D44F80EBCEB91CF/
      DATA               B(1)/Z4110000000000000/,B(2)/Z406D07BAC4D04469/
      DATA               B(3)/Z40146C6D6E730856/,B(4)/Z3F22015DE897CE8E/
      DATA               B(5)/Z3E1FE0E4A6BA9AF7/,B(6)/Z3CDB33FAB422E1DA/
C
C                                  COEFFICIENTS FOR R(7,7) APPROXIMATION
C                                  USED FOR X GREATER THAN OR EQUAL TO
C                                  -4 AND LESS THAN -1
C
      DATA               C(1)/Z3B174B1AD3E2E266/,C(2)/Z40FFFFB4D38BC730/
      DATA               C(3)/Z41BD92AE16CB369D/,C(4)/Z422D97D311F74081/
      DATA               C(5)/Z4245ED8DCFDDF7E2/,C(6)/Z422A852C0E187B8C/
      DATA               C(7)/Z418D633281818C45/,C(8)/Z4066C0AFC843D717/
      DATA               D(1)/Z4110000000000000/,D(2)/Z41CD923363206BCC/
      DATA               D(3)/Z4238717FD76B3FA3/,D(4)/Z426AA52AC377B0AB/
      DATA               D(5)/Z4259BB2A01911995/,D(6)/Z421F7F4782BEC4CE/
      DATA               D(7)/Z413CBABC9E89B478/,D(8)/Z401743F10E4FE979/
C
C                                  COEFFICIENTS FOR R(7,7) APPROXIMATION
C                                  USED FOR X LESS THAN -4
C
      DATA               E(1)/ZC0FFFFFFFFFFF884/,E(3)/ZC31AB885D23D00FD/
      DATA               E(4)/ZC395C04F986D432F/,E(5)/ZC41818DA23430F9D/
      DATA               E(2)/ZC22267FCB0C305CF/,E(6)/ZC419B018D42BE7A2/
      DATA               E(7)/ZC383A13CE9D21BEA/,E(8)/ZC1EE62A6F33B49AD/
      DATA               F(1)/Z4110000000000000/,F(2)/Z422467FCB0C2ED0F/
      DATA               F(3)/Z431EE585685725C3/,F(4)/Z43C7645BA5C7422E/
      DATA               F(5)/Z44286113476C61CF/,F(6)/Z443FC42536094A4D/
      DATA               F(7)/Z442B8DC679374D02/,F(8)/Z4394A2394EA64382/
C
C                                  COEFFICIENTS FOR R(5,5) APPROXIMATION
C                                  FOR LN(X/X0), ABS(1-X/X0) LESS THAN
C                                  .1
C
      DATA               P0(1)/Z429889F12723EF63/,
     1                   P0(2)/Z4315444AE68C49FC/,
     2                   P0(3)/Z43103BD17D0451EF/,
     3                   P0(4)/Z424DDEF79E7D4E0A/,
     4                   P0(5)/Z41775E404E1E69C8/,
     5                   P0(6)/Z3FE4572849F33E48/
      DATA               Q0(1)/Z429889F12723EF63/,
     1                   Q0(2)/Z431A089A6FC56977/,
     2                   Q0(3)/Z431A1294596B9C59/,
     3                   Q0(4)/Z42B9BD89B4DBA4B9/,
     4                   Q0(5)/Z4222E753B9A8CCCD/,
     5                   Q0(6)/Z4120000000000000/
C
C                                  COEFFICIENTS FOR R(8,8) APPROXIMATION
C                                  IN CHEBYSHEV POLYNOMIAL FORM USED FOR
C                                  X GREATER THAN 0.0 AND LESS THAN 6.
C
      DATA               P1(1)/Z415882FAB1C5ABD4/,
     1                   P1(2)/Z42CE5036A1E5C68A/,
     2                   P1(3)/Z4437C0581653B38E/,
     3                   P1(4)/Z448FF963B6E088AD/,
     4                   P1(5)/Z46449068754DD3A4/,
     5                   P1(6)/ZC619B411D22D93E0/,
     6                   P1(7)/Z481509C124CB86FB/,
     7                   P1(8)/ZC777652F57173F62/,
     8                   P1(9)/Z48B2AF8E86763BDA/
      DATA               Q1(1)/Z4219A10000000000/,
     1                   Q1(2)/ZC35E89E5032216A0/,
     2                   Q1(3)/Z44A6C08CCD646C5F/,
     3                   Q1(4)/ZC5B6965493822757/,
     4                   Q1(5)/Z4687293B6689CA00/,
     5                   Q1(6)/ZC74521D75324E9B5/,
     6                   Q1(7)/Z4817ED17037935DF/,
     7                   Q1(8)/ZC8535962239D6B97/,
     8                   Q1(9)/Z484C45A46809BFBA/
C
C                                  COEFFICIENTS FOR R(8,8) APPROXIMATION
C                                  IN J-FRACTION FORM, USED FOR
C                                  X GREATER THAN OR EQUAL TO 6.0 AND
C                                  LESS THAN 12.
C
      DATA               P2(1)/ZC12782B3E2F8723E/,
     1                   P2(2)/ZC224C73F69990A5C/,
     2                   P2(3)/Z421745E4DC48E91A/,
     3                   P2(4)/Z417E50C81D034A2D/,
     4                   P2(5)/ZC21369CDD0DEFDB4/,
     5                   P2(6)/Z415E2F7108D8FF75/,
     6                   P2(7)/Z4142E579A88DBC12/,
     7                   P2(8)/Z415BB2DC3A7A3A0C/,
     8                   P2(9)/Z40FFBBB08BC6EDB0/
      DATA               Q2(1)/Z412A3CBE7574E731/,
     1                   Q2(2)/Z433C567BC5456873/,
     2                   Q2(3)/ZC18633E65738ECF9/,
     3                   Q2(4)/Z4313D478C9B97198/,
     4                   Q2(5)/Z42345109DF0AD122/,
     5                   Q2(6)/Z431555D7E9168485/,
     6                   Q2(7)/ZC2C7264C3364A115/,
     7                   Q2(8)/Z4112570CE54009ED/
C
C                                  COEFFICIENTS FOR R(9,9) APPROXIMATION
C                                  IN J-FRACTION FORM USED FOR X GREATER
C                                  THAN OR EQUAL TO 12. AND LESS THAN 24
C
      DATA               P3(1)/ZC11A5D10E04A8481/,
     1                   P3(2)/ZC21299D5F90F9F7E/,
     2                   P3(3)/ZC1A01A4AF4BE8826/,
     3                   P3(4)/ZC2150EB24A571E03/,
     4                   P3(5)/ZC0E9DA0F2A953185/,
     5                   P3(6)/ZC2213C72BD7090B2/,
     6                   P3(7)/Z4218F472D6C8760A/,
     7                   P3(8)/Z421A869816A026F0/,
     8                   P3(9)/ZC11D85792426B9F4/,
     9                   P3(10)/Z40FFFF8FC5509274/
      DATA               Q3(1)/Z4261EC8D9F6EB5A3/,
     1                   Q3(2)/Z424009BAA237006F/,
     2                   Q3(3)/Z423BF306D954BF52/,
     3                   Q3(4)/Z42FDE1C64454B302/,
     4                   Q3(5)/Z422C4B4C38765403/,
     5                   Q3(6)/Z434A8D519BCB98C4/,
     6                   Q3(7)/Z42C719B6E66332E6/,
     7                   Q3(8)/ZC1AEF80FCFB6DE13/,
     8                   Q3(9)/Z411006485C45F248/
C
C                                  COEFFICIENTS FOR R(9,9) APPROXIMATION
C                                  IN J-FRACTION FORM USED FOR X GREATER
C                                  THAN OR EQUAL TO 24.
C
      DATA               P4(1)/Z42AF56BBAE030517/,
     1                   P4(2)/ZC2DF20AF08369607/,
     2                   P4(3)/ZC21231E952F731CB/,
     3                   P4(4)/ZC21BFAD7A31F07FA/,
     4                   P4(5)/ZC17A1A87A4CC3FB6/,
     5                   P4(6)/ZC1F49212B2DC262E/,
     6                   P4(7)/ZC17116FA4792119F/,
     7                   P4(8)/ZC1500045A1392820/,
     8                   P4(9)/ZC130000000DC93A1/,
     9                   P4(10)/Z4110000000000032/
      DATA               Q4(1)/Z449B689903F6E0E3/,
     1                   Q4(2)/Z413F90786B65B672/,
     2                   Q4(3)/Z4289CA5703B59D8A/,
     3                   Q4(4)/Z42752DE165133BF6/,
     4                   Q4(5)/Z42467BB1FE61EB5D/,
     5                   Q4(6)/ZC1C04CE86F570E7F/,
     6                   Q4(7)/ZC17FE10486DF7DDC/,
     7                   Q4(8)/ZC12FFFFEE39110D8/,
     8                   Q4(9)/Z411FFFFFFFFF588A/
      DATA               SIX/6.D0/,TWELVE/12.D0/,THREE/3.D0/,TWO/2.D0/,
     *                   ONE/1.D0/,HALF/.5D0/,TWENT4/24.D0/,FOUR/4.D0/,
     *                   FORTY/40.D0/,ZERO/0.D0/
C                                  FIRST EXECUTABLE STATEMENT
      IEND = 8
      JEND = 9
      KEND = 6
      IENDM1 = IEND-1
      IENDP1 = IEND+1
      JENDP1 = JEND+1
      KENDM1 = KEND-1
      X = ARG
      IER = 0
      IF (IOPT.LT.1.OR.IOPT.GT.3) GO TO 125
      GO TO (5,100,5), IOPT
C                                  IOPT = 1 OR 3
    5 IF (X) 60,115,10
   10 IF (X.GE.TWELVE) GO TO 35
      IF (X.GE.SIX) GO TO 25
C                                  X GREATER THAN OR EQUAL TO 0 AND
C                                  LESS THAN 6 - RATIONAL APPROXIMATION
C                                  USED IS EXPRESSED IN TERMS OF
C                                  CHEBYSHEV POLYNOMIALS TO IMPROVE
C                                  CONDITIONING
      T = X+X
      T = T/THREE-TWO
      PX(1) = ZERO
      QX(1) = ZERO
      PX(2) = P1(1)
      QX(2) = Q1(1)
      DO 15 I=2,IEND
         PX(I+1) = T*PX(I)-PX(I-1)+P1(I)
         QX(I+1) = T*QX(I)-QX(I-1)+Q1(I)
   15 CONTINUE
      SUMP = HALF*T*PX(IENDP1)-PX(IEND)+P1(IENDP1)
      SUMQ = HALF*T*QX(IENDP1)-QX(IEND)+Q1(IENDP1)
      FRAC = SUMP/SUMQ
      XMX0 = (X-X01)-X02
      IF (DABS(XMX0).LT.0.037D0) GO TO 20
      XX0 = X/X0
      MMDEI = DLOG(XX0)+XMX0*FRAC
      IF (IOPT.EQ.3) MMDEI = DEXP(-X)*MMDEI
      GO TO 9005
C                                  EVALUATE APPROXIMATION FOR LN(X/X0)
C                                  FOR X CLOSE TO X0
   20 Y = XMX0/X0
      SUMP = ((((P0(6)*Y+P0(5))*Y+P0(4))*Y+P0(3))*Y+P0(2))*Y+P0(1)
      SUMQ = ((((Q0(6)*Y+Q0(5))*Y+Q0(4))*Y+Q0(3))*Y+Q0(2))*Y+Q0(1)
      MMDEI = (SUMP/(SUMQ*X0)+FRAC)*XMX0
      IF (IOPT.EQ.3) MMDEI = DEXP(-X)*MMDEI
      GO TO 9005
C                                  X GREATER THAN OR EQUAL TO 6 AND
C                                  LESS THAN 12
   25 FRAC = ZERO
      DO 30 I=1,IEND
         FRAC = Q2(I)/(P2(I)+X+FRAC)
   30 CONTINUE
      MMDEI = (P2(IENDP1)+FRAC)/X
      IF (IOPT.NE.3) MMDEI = MMDEI*DEXP(X)
      GO TO 9005
C                                  X GREATER THAN OR EQUAL TO 12 AND
C                                  LESS THAN 24
   35 IF (X.GE.TWENT4) GO TO 45
      FRAC = ZERO
      DO 40 I=1,JEND
         FRAC = Q3(I)/(P3(I)+X+FRAC)
   40 CONTINUE
      MMDEI = (P3(JENDP1)+FRAC)/X
      IF (IOPT.NE.3) MMDEI = MMDEI*DEXP(X)
      GO TO 9005
C                                  X GREATER THAN OR EQUAL TO 24
   45 IF ((X.GE.XMAX).AND.(IOPT.LT.3)) GO TO 110
      Y = ONE/X
      FRAC = ZERO
      DO 50 I=1,JEND
         FRAC = Q4(I)/(P4(I)+X+FRAC)
   50 CONTINUE
      FRAC = P4(JENDP1)+FRAC
      MMDEI = Y+Y*Y*FRAC
      IF (IOPT.EQ.3) GO TO 9005
      IF (X.GT.170.0D0) GO TO 55
      MMDEI = MMDEI*DEXP(X)
      GO TO 9005
C                                  CALCULATION REFORMULATED TO AVOID
C                                  PREMATURE OVERFLOW
   55 MMDEI = (MMDEI*DEXP(X-FORTY))*DEXP40
      GO TO 9005
C                                  ORIGINAL X WAS NEGATIVE.
   60 Y = -X
   65 W = ONE/Y
      IF (Y.GT.FOUR) GO TO 85
      IF (Y.GT.ONE) GO TO 75
C                                  -X GREATER THAN 0 AND LESS THAN OR
C                                  EQUAL TO 1
      SUMB = B(KEND)
      DO 70 I=1,KENDM1
         J = KEND-I
         SUMB = (SUMB*Y)+B(J)
   70 CONTINUE
      MMDEI = DLOG(Y)-(((((A(6)*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y+A(1))
     1/SUMB
      IF (IOPT.EQ.3) MMDEI = MMDEI*DEXP(Y)
      GO TO 95
C                                  -X GREATER THAN -1 AND LESS THAN OR
C                                  EQUAL TO 4
   75 SUMC = C(IEND)
      SUMD = D(IEND)
      DO 80 I=1,IENDM1
         J = IEND-I
         SUMC = (SUMC*W)+C(J)
         SUMD = (SUMD*W)+D(J)
   80 CONTINUE
      MMDEI = -SUMC/SUMD
      IF (IOPT.EQ.3) GO TO 9005
      MMDEI = MMDEI*DEXP(-Y)
      GO TO 95
C                                  -X GREATER THAN 4
   85 IF ((-DABS(X).LT.XMIN).AND.(IOPT.LT.3)) GO TO 105
      SUME = E(IEND)
      SUMF = F(IEND)
      DO 90 I=1,IENDM1
         J = IEND-I
         SUME = (SUME*W)+E(J)
         SUMF = (SUMF*W)+F(J)
   90 CONTINUE
      MMDEI = -W*(1.0D0+W*SUME/SUMF)
      IF (IOPT.EQ.3) GO TO 9005
      MMDEI = MMDEI*DEXP(-Y)
   95 IF (IOPT.EQ.2) MMDEI = -MMDEI
      IF (IER.EQ.69) GO TO 9000
      GO TO 9005
  100 Y = X
      IF (Y) 120,115,65
C                                  TERMINAL ERROR - ARG IS LESS THAN
C                                  XMIN CAUSING UNDERFLOW
  105 MMDEI = ZERO
      IER = 132
      GO TO 9000
C                                  TERMINAL ERROR - X IS GREATER THAN
C                                  XMAX CAUSING OVERFLOW
  110 MMDEI = XINF
      IER = 131
      GO TO 9000
C                                  TERMINAL ERROR - ARG = 0
  115 MMDEI = -XINF
      IF (IOPT.EQ.2) MMDEI = -MMDEI
      IER = 130
      GO TO 9000
C                                  WARNING WITH FIX - ARG IS LESS THAN
C                                  0.0 FOR IOPT = 2
  120 IER = 69
      GO TO 60
C                                  TERMINAL ERROR - IOPT IS OUT OF RANGE
  125 MMDEI = XINF
      IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMMDEI )
 9005 RETURN
      END

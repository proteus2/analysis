C   IMSL ROUTINE NAME   - CTPR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - COMPUTE EXACT PROBABILITIES FOR CONTINGENCY
C                           TABLES
C
C   USAGE               - CALL CTPR (A,IA,IR,IC,PRT,PRE,PCHEK,IWK)
C
C   ARGUMENTS    A      - INPUT IR BY IC MATRIX CONTAINING THE OBSERVED
C                           COUNTS OF THE CONTINGENCY TABLE.
C                IA     - INPUT ROW DIMENSION OF A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                IR     - INPUT NUMBER OF ROWS IN THE CONTINGENCY
C                           TABLE.
C                IC     - INPUT NUMBER OF COLUMNS IN THE
C                           CONTINGENCY TABLE.
C                PRT    - OUTPUT PROBABILITY OF OBTAINING THE
C                           OBSERVED TABLE, GIVEN THE MARGINAL
C                           TOTALS.
C                PRE    - OUTPUT PROBABILITY OF OBTAINING THE
C                           OBSERVED OR ANY MORE EXTREME TABLE,
C                           GIVEN THE MARGINAL TOTALS.
C                PCHEK  - OUTPUT PROBABILITY OF ANY OF THE TABLES
C                           CONSIDERED FOR THE GIVEN MARGINAL
C                           TOTALS.  PCHEK SHOULD BE 1.0.  THE
C                           EXTENT TO WHICH IT DEVIATES FROM
C                           1.0 INDICATES THE ACCURACY OF THE
C                           COMPUTATIONS.
C                IWK    - WORK VECTOR OF LENGTH (IR+2)*(IC+2).
C
C   REQD. IMSL ROUTINES - SINGLE/CTPR1
C                       - DOUBLE/CTPR1,VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CTPR   (A,IA,IR,IC,PRT,PRE,PCHEK,IWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IR,IC,IWK(1)
      REAL               A(IA,1),PRT,PRE,PCHEK
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ICM1,IDIF,IEND,IFLAG,II,IND,IP,IRCRCR,IRCRCS,
     *                   IRIC1,IRICR1,IRICRC,IRICR,IRIC,IRJ,IRM2,IRP1,
     *                   IRP2,IRPIC,ISS,IST,ISUM,ITOT,IWKIJ,IWKI,I,JJ,
     *                   JND,JP,JSUM,J,K,MIN,M
      REAL               EPS,ONEPS,ONE,PX,ZERO,ECON,ECON1
      DOUBLE PRECISION   DSUM1,DSUM2,DSUM3,FACLOG(201),DEM,CLG,TEMP1,
     *                   TEMP2,X,CTPR1
      EXTERNAL           CTPR1
      DATA               EPS/Z3C100000/
      DATA               ZERO /0.0/,ONE /1.0/
      DATA               ECON1/Z00100000/
C                                  FIRST EXECUTABLE STATEMENT
C                                    INITIALIZATIONS
      ONEPS = ONE-100.0*EPS
      ECON = ALOG(ECON1) + 10
      ICM1 = IC-1
      IRM2 = IR-2
      IRIC = IR*IC
      IRIC1 = IRIC+1
      IRICR = IRIC+IR
      IRICRC = IRICR+IC
      IRPIC = IR+IC
      IRICR1 = IRICR+1
      IRP1 = IR+1
      IRP2 = IR+2
      IRCRCR = IRICRC+IR
      IRCRCS = IRCRCR+IC+1
      DO 5 I=IRIC1,IRICRC
         IWK(I) = 0
    5 CONTINUE
      ITOT = 0
      IRJ = -IR
      DO 15 J=1,IC
         IRJ = IRJ+IR
         DO 10 I=1,IR
            IWKIJ = A(I,J)+.5
            IWK(IRJ+I) = IWKIJ
            IWK(IRICR+J) = IWK(IRICR+J)+IWKIJ
            IWK(IRIC+I) = IWK(IRIC+I)+IWKIJ
   10    CONTINUE
         ITOT = ITOT+IWK(IRICR+J)
   15 CONTINUE
      IEND = MIN0(200,ITOT)+1
      X = 0.0D0
      FACLOG(1) = 0.0D0
      DSUM1 = 0.0D0
      DO 20 I=2,IEND
         X = X+1.0D0
         DSUM1 = DSUM1+DLOG(X)
         FACLOG(I) = DSUM1
   20 CONTINUE
      DSUM1 = 0.0D0
      DSUM2 = 0.0D0
      IF (ITOT.GT.200) GO TO 25
      DSUM1 = -FACLOG(ITOT+1)
      GO TO 30
   25 DSUM1 = -CTPR1(ITOT)
   30 DO 40 I=IRIC1,IRICRC
         IWKI = IWK(I)
         IF (IWKI.GT.200) GO TO 35
         DSUM1 = DSUM1+FACLOG(IWKI+1)
         GO TO 40
   35    DSUM1 = DSUM1+CTPR1(IWKI)
   40 CONTINUE
      CLG = DSUM1
C                                  COMPUTE PROBABILITY OF GIVEN TABLE
      DO 50 I=1,IRIC
         IWKI = IWK(I)
         IF (IWKI.GT.200) GO TO 45
         DSUM2 = DSUM2+FACLOG(IWKI+1)
         GO TO 50
   45    DSUM2 = DSUM2+CTPR1(IWKI)
   50 CONTINUE
      DEM = DSUM2
      PRT = DEXP(CLG-DEM)
C                                  SET UP TO GENERATE OTHER
C                                    COMBINATIONS
      DO 55 I=1,IRIC
         IWK(I) = 0
   55 CONTINUE
      DSUM1 = 0.0D0
      DSUM2 = 0.0D0
   60 IFLAG = 0
      IWK(IRP2) = -1
      IST = 2
      DO 160 J=2,IC
         IST = IST+IR
         IEND = IST+IRM2
         K = 1
         DO 155 I=IST,IEND
            K = K+1
            IWK(I) = IWK(I)+1
            ISUM = 0
            JSUM = 0
            II = I-IR
            DO 65 M=J,IC
               II = II+IR
               ISUM = ISUM+IWK(II)
   65       CONTINUE
            IF (ISUM.GT.IWK(IRIC+K)) GO TO 145
            DO 70 II=I,IEND
               JSUM = JSUM+IWK(II)
   70       CONTINUE
            IF (JSUM.GT.IWK(IRICR+J)) GO TO 145
            IF (IFLAG.EQ.1) GO TO 165
            IP = K
            JP = J
C                                  FORM NEXT TABLE SATISFYING MARGINALS
            DO 75 II=1,IR
               IWK(IRICRC+II) = IWK(IRIC+II)
   75       CONTINUE
            ISS = 1
            DO 80 II=IRICR1,IRICRC
               IWK(ISS) = 0
               ISS = ISS+IR
               IWK(II+IRPIC) = IWK(II)
   80       CONTINUE
            IRJ = 0
            DO 90 JJ=2,IC
               IRJ = IRJ+IR
               DO 85 II=2,IR
                  IWKIJ = IWK(IRJ+II)
                  IWK(IRICRC+II) = IWK(IRICRC+II)-IWKIJ
                  IWK(IRCRCR+JJ) = IWK(IRCRCR+JJ)-IWKIJ
   85          CONTINUE
   90       CONTINUE
            DO 105 JJ=1,IC
               JND = IRCRCS-JJ
               ISS = IRIC-JJ*IR
               DO 100 II=1,IR
                  ISS = ISS+1
                  IND = IRICRC+II
                  MIN = MIN0(IWK(IND),IWK(JND))
                  IF (MIN.EQ.0) GO TO 95
                  IWK(ISS) = IWK(ISS)+MIN
                  IWK(IND) = IWK(IND)-MIN
                  IWK(JND) = IWK(JND)-MIN
   95             IF (IWK(JND).EQ.0) GO TO 105
  100          CONTINUE
  105       CONTINUE
C                                  COMPUTE PROBABILITY
            DSUM3 = 0.D0
            DO 115 II=1,IRIC
               IWKI = IWK(II)
               IF (IWKI.GT.200) GO TO 110
               DSUM3 = DSUM3+FACLOG(IWKI+1)
               GO TO 115
  110          DSUM3 = DSUM3+CTPR1(IWKI)
  115       CONTINUE
            DEM = DSUM3
            IF (CLG-DEM.GT.DBLE(ECON)) PX = DEXP(CLG-DEM)
            IF (CLG-DEM.LE.DBLE(ECON)) PX = 0.0D0
C                                  IDENTIFY THE EXACT TABLE
            ICK1 = 0
            IF (IR.EQ.2 .AND. IC.EQ.2) GO TO 125
            IBEG = 3
            IF (IR.EQ.2) IBEG = IRP2+1
            DO 120 IK=IBEG,IRIC
               IF (IK.EQ.IRP1 .OR. IK.EQ.IRP2) GO TO 120
               JJ = (IK-1)/IR
               II = IK-JJ*IR
               JJ = JJ+1
               IAI = A(II,JJ)
               IF (IWK(IK).NE.IAI) GO TO 130
  120       CONTINUE
  125       ICK1 = 1
            IA21 = A(2,1)+.01
            IF (IWK(2).NE.IA21) GO TO 130
C                                  UPDATE PCHEK AND PRE
            DSUM1 = DSUM1+PX
            DSUM2 = DSUM2+PX
            GO TO 135
  130       CONTINUE
            DSUM1 = DSUM1+PX
            IF (PRT.GE.ONEPS*PX) DSUM2 = DSUM2+PX
  135       IF (IWK(IRP1).LT.1 .OR. IWK(2).LT.1) GO TO 150
            IWK(1) = IWK(1)+1
            IWK(IRP2) = IWK(IRP2)+1
            TEMP1 = IWK(2)*IWK(IRP1)
            TEMP2 = IWK(1)*IWK(IRP2)
            IF (PX.EQ.0.0D0) GO TO 1
            PX = PX*TEMP1/TEMP2
            GO TO 2
    1       DEM = DEM - DLOG(TEMP1/TEMP2)
            IF(CLG-DEM.GT.DBLE(ECON)) PX = DEXP(CLG-DEM)
    2       CONTINUE
            IWK(IRP1) = IWK(IRP1)-1
            IWK(2) = IWK(2)-1
            IF (ICK1.EQ.0) GO TO 140
            IF (IWK(2).NE.IA21) GO TO 140
            DSUM1 = DSUM1+PX
            DSUM2 = DSUM2+PX
            GO TO 135
  140       CONTINUE
            DSUM1 = DSUM1+PX
            IF (PRT.GE.ONEPS*PX) DSUM2 = DSUM2+PX
            GO TO 135
  145       IP = K
            JP = J
  150       IFLAG = 1
  155    CONTINUE
  160 CONTINUE
      GO TO 185
C                                  RESET CELLS TO ZERO
  165 JP = JP-1
      IST = JP*IR+2
      IEND = IST+IP-2
      DO 170 I=IST,IEND
         IWK(I) = 0
  170 CONTINUE
      IST = 2-IR
      DO 180 J=1,JP
         IST = IST+IR
         IEND = IST+IRM2
         DO 175 I=IST,IEND
            IWK(I) = 0
  175    CONTINUE
  180 CONTINUE
      GO TO 60
  185 CONTINUE
      PCHEK = DSUM1
      PRE = DSUM2
      RETURN
      END

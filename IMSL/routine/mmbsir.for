C   IMSL ROUTINE NAME   - MMBSIR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - MODIFIED BESSEL FUNCTION OF THE FIRST KIND
C                           OF NONNEGATIVE REAL ORDER FOR REAL POSITIVE
C                           ARGUMENTS WITH EXPONENTIAL SCALING OPTION
C
C   USAGE               - CALL MMBSIR (ARG,ORDER,N,IOPT,B,IER)
C
C   ARGUMENTS    ARG    - INPUT ARGUMENT. ARG MUST BE TYPED APPRO-
C                           PRIATELY IN THE CALLING PROGRAM. (SEE THE
C                           PRECISION/HARDWARE SECTION.) ARG MUST BE
C                           GREATER THAN OR EQUAL TO ZERO. IF SCALING
C                           IS NOT DESIRED (IOPT = 1), THEN ARG MUST BE
C                           LESS THAN EXPARG. EXPARG IS APPROXIMATELY
C                           THE LARGEST ACCEPTABLE ARGUMENT FOR THE
C                           FORTRAN EXP FUNCTION. SEE THE PROGRAMMING
C                           NOTES IN THE MANUAL FOR THE EXACT VALUE OF
C                           EXPARG.
C                ORDER  - INPUT VALUE SPECIFYING THE DESIRED ORDER OF
C                           THE BESSEL FUNCTION. ORDER MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) ORDER MUST
C                           BE GREATER THAN OR EQUAL TO ZERO AND LESS
C                           THAN ONE.
C                N      - INPUT PARAMETER SPECIFYING THE NUMBER OF
C                           FUNCTION VALUES TO BE COMPUTED.
C                IOPT   - INPUT SCALING OPTION SWITCH.
C                           IF IOPT = 1 THEN THE BESSEL FUNCTIONS
C                           CALCULATED ARE UNSCALED.
C                           IF IOPT = 2 THEN THE BESSEL FUNCTIONS
C                           CALCULATED ARE SCALED BY EXP(-ARG).
C                B      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           COMPUTED FUNCTION VALUES. B MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) B(1) WILL
C                           CONTAIN THE COMPUTED VALUE FOR THE INPUT
C                           ORDER, B(2) WILL CONTAIN THE COMPUTED
C                           FUNCTION VALUE FOR ORDER + 1, B(3) FOR
C                           ORDER + 2, ETC.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ONE OR MORE OF
C                             THE INPUT ARGUMENTS, ARG, ORDER, N,
C                             OR IOPT, IS OUT OF RANGE. B(I), (I=1,N)
C                             IS SET TO MACHINE INFINITY.
C                           IER = 129+J INDICATES THAT N WAS SO MUCH
C                             GREATER THAN THE LARGEST INTEGER IN ARG
C                             THAT ONLY B(I), (I=1,J) WAS CALCULATED
C                             CORRECTLY. SEE THE PROGRAMMING NOTES
C                             FOR THE ACCURACY OF THE REMAINING B(I).
C
C   REQD. IMSL ROUTINES - H32,H36/MGAMAD=DGAMMA,UERTST,UGETIO
C                         H48,H60/MGAMA=GAMMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMBSIR (ARG,ORDER,N,IOPT,B,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IOPT,IER
      DOUBLE PRECISION   ARG,ORDER,B(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDUM,NSIG,ISUM,MAGX,NCALC,NBMX,ISTEP
      INTEGER            NEND,NSTART,NN,L,I,INTEN
      DOUBLE PRECISION   P,TEST,TEMPA,TEMPB,TEMPC,EXPARG,SUM,TOVER
      DOUBLE PRECISION   PLAST,POLD,PSAVE,PSAVEL,ENTEN,ENSIG,RTNSIG
      DOUBLE PRECISION   EN,EM,EMPAL,EMP2AL,ENMTEN,STEP,STEP2,XLARGE
      DOUBLE PRECISION   HALFX,ZSIG,SMALL,XINF,XNTEN
      DOUBLE PRECISION   DGAMMA
      DATA               SMALL/Z0010000000000000/
      DATA               XLARGE /1.D+04/
      DATA               XINF/.723700557733226D+76/
      DATA               ZSIG/ 16.0D0/
      DATA               EXPARG/ 174.0D0/
      DATA               XNTEN/ 75.D0/
C                                  FIRST EXECUTABLE STATEMENT
      NSIG = IFIX(SNGL(ZSIG))
      ISUM = -NSIG/4
      RTNSIG = 10.D0**ISUM
      IER = 0
      INTEN = XNTEN
      ENTEN = 10.D0**INTEN
      ENSIG = 10.0D0**NSIG
      ENMTEN = 4.D0*SMALL
      TEMPA = ARG
C                                  CHECK ARGUMENT RANGE
      IF((N.GT.0) .AND. (ARG.GE.0.D0).AND.(((IOPT.EQ.1) .AND.
     *  (ARG.LE.EXPARG)) .OR. ((IOPT.EQ.2) .AND. (ARG.LE.XLARGE))) .AND.
     *  (ORDER.GE.0.D0) .AND. (ORDER.LT.1.D0)) GO TO 10
      IER = 129
      NN = N
      IF(N.LE.0) NN = 1
      DO 5  I=1,NN
         B(I) = XINF
    5 CONTINUE
      GO TO 9000
C                                  USE 2-TERM ASCENDING SERIES FOR
C                                    SMALL ARG
   10 MAGX = IFIX(SNGL(TEMPA))
      NCALC = N
      IF(TEMPA.LT.RTNSIG) GO TO 115
C                                  INITIALIZE THE CALCULATION OF P*S
      NBMX = N-MAGX
      NN = MAGX+1
      ISTEP = NN + NN
      STEP = ORDER + ORDER
      EN = DBLE(FLOAT(ISTEP)) + STEP
      PLAST = 1.D0
      P = EN/TEMPA
C                                  CALCULATE GENERAL SIGNIFICANCE TEST
      TEST = ENSIG + ENSIG
      STEP = TEST * P
      IF(2*MAGX .GT. 5*NSIG) TEST = DSQRT(STEP)
      IF(2*MAGX .LE. 5*NSIG) TEST = TEST/1.585D0**MAGX
      IF(NBMX .LT. 3) GO TO 20
C                                  CALCULATE P*S UNTIL NN = N-1.
C                                    CHECK FOR POSSIBLE OVERFLOW.
      TOVER = ENTEN/ENSIG
      NSTART = MAGX + 2
      NEND = N - 1
      ISTEP = NSTART + NSTART
      STEP = ORDER + ORDER
      EN = DBLE(FLOAT(ISTEP)) - 2.D0+STEP
      DO 15  NN=NSTART,NEND
         EN = EN + 2.D0
         POLD = PLAST
         PLAST = P
         P = EN*PLAST/TEMPA+POLD
         IF(P .GT. TOVER) GO TO 25
   15 CONTINUE
      NN = NEND
      ISTEP = NN + NN
      STEP = ORDER + ORDER
      EN = DBLE(FLOAT(ISTEP)) + STEP
C                                  CALCULATE SPECIAL SIGNIFICANCE
C                                    TEST FOR NBMX .GT. 2.
      STEP = PLAST * ENSIG
      STEP2 = P + P
      STEP = DSQRT(STEP) * DSQRT(STEP2)
      TEST = DMAX1(TEST,STEP)
C                                  CALCULATE P*S UNTIL SIGNIFICANCE
C                                    TEST PASSES
   20 NN = NN + 1
      EN = EN + 2.D0
      POLD = PLAST
      PLAST = P
      P = EN*PLAST/TEMPA+POLD
      IF(P.LT.TEST) GO TO 20
      GO TO 45
C                                  TO AVOID OVERFLOW, DIVIDE P*S BY
C                                    TOVER. CALCULATE P*S UNTIL
C                                    ABS(P).GT.1.
   25 TOVER = ENTEN
      P = P/TOVER
      PLAST = PLAST/TOVER
      PSAVE = P
      PSAVEL = PLAST
      NSTART = NN + 1
   30 NN = NN + 1
      EN = EN + 2.D0
      POLD = PLAST
      PLAST = P
      P = EN*PLAST/TEMPA+POLD
      IF(P .LE. 1.D0) GO TO 30
      TEMPB = EN/TEMPA
C                                  CALCULATE BACKWARD TEST, AND FIND
C                                    NCALC, THE HIGHEST N SUCH THAT
C                                    THE TEST IS PASSED.
      TEST = POLD*PLAST*(.5D0-.5D0/(TEMPB*TEMPB))/ENSIG
      P = PLAST * TOVER
      NN = NN - 1
      EN = EN - 2.D0
      NEND = MIN0(N,NN)
      DO 35  NCALC=NSTART,NEND
         POLD = PSAVEL
         PSAVEL = PSAVE
         PSAVE = EN*PSAVEL/TEMPA+POLD
         IF(PSAVE*PSAVEL .GT. TEST) GO TO 40
   35 CONTINUE
      NCALC = NEND + 1
   40 NCALC = NCALC - 1
C                                  INITIALIZE THE BACKWARD RECURSION
C                                    AND THE NORMALIZATION SUM
   45 NN = NN + 1
      EN = EN + 2.D0
      TEMPB = 0.D0
      TEMPA = 1.D0/P
      EM = DBLE(FLOAT(NN)) - 1.D0
      EMPAL = EM + ORDER
      EMP2AL = (EM-1.D0) + (ORDER+ORDER)
      SUM = TEMPA*EMPAL*EMP2AL/EM
      NEND = NN - N
      IF(NEND) 70, 60, 50
C                                  RECUR BACKWARD VIA DIFFERENCE
C                                    EQUATION, CALCULATING (BUT
C                                    NOT STORING) B(N), UNTIL NN = N.
   50 DO 55  L=1,NEND
         NN = NN - 1
         EN = EN - 2.D0
         TEMPC = TEMPB
         TEMPB = TEMPA
         TEMPA = (EN*TEMPB)/ARG+TEMPC
         EM = EM - 1.D0
         EMP2AL = EMP2AL - 1.D0
         IF(NN .EQ. 1) GO TO 60
         IF(NN .EQ. 2) EMP2AL = 1.D0
         EMPAL = EMPAL - 1.D0
         SUM = (SUM+TEMPA*EMPAL)*EMP2AL/EM
   55 CONTINUE
C                                  STORE B(N)
   60 B(NN) = TEMPA
      IF(N .GT. 1) GO TO 65
      SUM = (SUM+SUM) + TEMPA
      GO TO 100
C                                  CALCULATE AND STORE B(N-1)
   65 NN = NN - 1
      EN = EN - 2.D0
      B(NN) = (EN*TEMPA)/ARG+TEMPB
      IF(NN .EQ. 1) GO TO 95
      EM = EM - 1.D0
      EMP2AL = EMP2AL - 1.D0
      IF(NN .EQ. 2) EMP2AL = 1.D0
      EMPAL = EMPAL - 1.D0
      SUM = (SUM+B(NN)*EMPAL)*EMP2AL/EM
      GO TO 80
C                                  NN.LT.N, SO STORE B(NN) AND SET
C                                    HIGHER ORDERS TO ZERO
   70 B(NN) = TEMPA
      NEND = -NEND
      DO 75  L=1,NEND
         B(NN+L) = 0.D0
   75 CONTINUE
   80 NEND = NN - 2
      IF(NEND .EQ. 0) GO TO 90
C                                  CALCULATE VIA DIFFERENCE EQUATION
C                                    AND STORE B(NN), UNTIL NN = 2
      DO 85  L=1,NEND
         NN = NN - 1
         EN = EN - 2.D0
         B(NN) = (EN*B(NN+1))/ARG+B(NN+2)
         EM = EM - 1.D0
         EMP2AL = EMP2AL - 1.D0
         IF(NN .EQ. 2) EMP2AL = 1.D0
         EMPAL = EMPAL - 1.D0
         SUM = (SUM+B(NN)*EMPAL)*EMP2AL/EM
   85 CONTINUE
C                                  CALCULATE B(1)
   90 B(1) = 2.D0*EMPAL*B(2)/ARG+B(3)
   95 SUM = (SUM+SUM) + B(1)
C                                  NORMALIZE. DIVIDE ALL B(NN) BY SUM
  100 IF(ORDER .EQ. 0.D0) GO TO 105
      STEP = 1.D0 + ORDER
      SUM = SUM*DGAMMA(STEP)*(ARG*.5D0)**(-ORDER)
  105 STEP = -ARG
      IF(IOPT .EQ. 1) SUM = SUM*DEXP(STEP)
      TEMPA = ENMTEN
      IF(SUM .GT. 1.D0) TEMPA = TEMPA*SUM
      DO 110  NN=1,N
         IF(B(NN) .LT. TEMPA) B(NN) = 0.D0
         B(NN) = B(NN)/SUM
  110 CONTINUE
      GO TO 135
C                                  TWO-TERM ASCENDING SERIES FOR SMALL A
  115 TEMPA = 1.D0
      EMPAL = 1.D0 + ORDER
      HALFX = 0.D0
      IF(ARG .GT. ENMTEN) HALFX = .5D0*ARG
      IF(ORDER .NE. 0.D0) TEMPA = HALFX**ORDER/DGAMMA(EMPAL)
      STEP = -ARG
      IF(IOPT .EQ. 2) TEMPA = TEMPA*DEXP(STEP)
      TEMPB = 0.D0
      IF((ARG+1.D0) .GT. 1.D0) TEMPB = HALFX*HALFX
      B(1) = TEMPA+TEMPA*TEMPB/EMPAL
      IF((ARG .NE. 0.D0) .AND. (B(1) .EQ. 0.D0)) NCALC = 0
      IF(N .EQ. 1) GO TO 135
      IF(ARG .GT. 0.D0) GO TO 125
      DO 120  NN=2,N
         B(NN) = 0.D0
  120 CONTINUE
      GO TO 9005
  125 TEMPC = HALFX
      TOVER = (ENMTEN+ENMTEN)/ARG
      IF(TEMPB .NE. 0.D0) TOVER = ENMTEN/TEMPB
      DO 130  NN=2,N
         TEMPA = TEMPA/EMPAL
         EMPAL = EMPAL + 1.D0
         TEMPA = TEMPA * TEMPC
         IF(TEMPA .LE. TOVER*EMPAL) TEMPA = 0.D0
         B(NN) = TEMPA+TEMPA*TEMPB/EMPAL
         IF((B(NN) .EQ. 0.D0) .AND. (NCALC .GT. NN)) NCALC = NN-1
  130 CONTINUE
  135 IF(NCALC .EQ. N) GO TO 9005
      IER = 129 + NCALC
 9000 CONTINUE
      CALL UERTST(IER,6HMMBSIR)
 9005 RETURN
      END

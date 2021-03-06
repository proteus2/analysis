C   IMSL ROUTINE NAME   - MMBSJN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - BESSEL FUNCTION OF THE FIRST KIND OF
C                           NONNEGATIVE INTEGER ORDER FOR
C                           REAL ARGUMENTS
C
C   USAGE               - CALL MMBSJN (ARG,N,B,IER)
C
C   ARGUMENTS    ARG    - INPUT ARGUMENT. THE ABSOLUTE VALUE OF ARG MUST
C                           BE LESS THAN OR EQUAL TO 100000. ARG MUST BE
C                           TYPED APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE THE PRECISION/HARDWARE SECTION.)
C                N      - INPUT PARAMETER SPECIFYING THE NUMBER OF
C                           FUNCTION VALUES TO BE COMPUTED.
C                B      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           COMPUTED FUNCTION VALUES. B MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE THE PRECISION/HARDWARE SECTION.)
C                           B(1) WILL CONTAIN THE COMPUTED VALUE FOR
C                           ORDER ZERO, B(2) WILL CONTAIN THE COMPUTED
C                           VALUE FOR ORDER 1, B(3) FOR ORDER 2, ETC.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER ARG OR N IS
C                             OUT OF RANGE. B(I), (I=1,N) IS SET TO
C                             MACHINE INFINITY.
C                           IER = 129 + J INDICATES THAT B(I), (I=1,J)
C                             ARE COMPUTED TO MACHINE PRECISION, BUT
C                             PRECISION IS LOST FOR B(I), (I=J+1,N.)
C                             SEE THE PROGRAMMING NOTES.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMBSJN (ARG,N,B,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      DOUBLE PRECISION   ARG,B(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,LARGEX,MAGX,NCALC,NN,NBMX,M,NSTART,NEND
      DOUBLE PRECISION   TEST,TEMPA,TEMPB,TEMPC,EXPARG,P
      DOUBLE PRECISION   RSIGN,SUM,TOVER,PLAST,POLD,PSAVE,PSAVEL
      DOUBLE PRECISION   DSIG,RTEN,TMPA4,SMALLX,XINF
      DATA               DSIG/16.0D0/
      DATA               RTEN/75.D0/
      DATA               EXPARG/174.6730895011062D0/
      DATA               XINF/.723700557733226D+76/
      DATA               LARGEX /100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      TEMPA = DABS(ARG)
      MAGX = IFIX(SNGL(TEMPA))
      IF(N.GT.0 .AND. MAGX.LE.LARGEX) GO TO 10
C                                  ERROR RETURN -- ARG,N IS OUT OF RANGE
      IER = 129
      B(1) = XINF
      IF(N.LT.2) GO TO 9000
      DO 5 L=2,N
         B(L) = XINF
    5 CONTINUE
      GO TO 9000
   10 RSIGN = 1.D0
      NCALC = N
C                                  USE 2-TERM ASCENDING SERIES FOR
C                                    SMALL ARG
      TMPA4 = TEMPA**4.D0
      SMALLX = .1D0**DSIG
      IF(TMPA4.GE.SMALLX) GO TO 20
C                                  TWO-TERM ASCENDING SERIES FOR
C                                    SMALL ARG
      TEMPA = 1.D0
      TEMPB = -.25D0*ARG*ARG*RSIGN
      B(1) = 1.D0+TEMPB
      IF(N.EQ.1) GO TO 9005
      DO 15 NN=2,N
         TEMPA = TEMPA*ARG/DBLE(FLOAT(2*NN-2))
         B(NN) = TEMPA*(1.D0+TEMPB/DBLE(FLOAT(NN)))
   15 CONTINUE
      GO TO 9005
C                                  INITIALIZE THE CALCULATION OF P*S
   20 NBMX = N-MAGX
      NN = MAGX+1
      PLAST = 1.D0
      P = DBLE(FLOAT(2*NN))/TEMPA
C                                  CALCULATE GENERAL SIGNIFICANCE TEST
      TEST = 2.D0*1.D1**DSIG
      M = 0
      IF(NBMX.LT.3) GO TO 30
C                                  CALCULATE P*S UNTIL NN=N-1.
C                                    CHECK FOR POSSIBLE OVERFLOW.
      TOVER = 1.D1**(RTEN-DSIG)
      NSTART = MAGX+2
      NEND = N-1
      DO 25 NN=NSTART,NEND
         POLD = PLAST
         PLAST = P
         P = DBLE(FLOAT(2*NN))*PLAST/TEMPA-RSIGN*POLD
         IF(P-TOVER) 25, 25, 35
   25 CONTINUE
      NN = NEND
C                                  CALCULATE SPECIAL SIGNIFICANCE TEST
C                                    FOR NBMX.GT.2.
C
      TEST = DMAX1(TEST,DSQRT(PLAST*1.D1**DSIG)*DSQRT(2.D0*P))
C
C                                  CALCULATE P*S UNTIL SIGNIFICANCE
C                                    TEST PASSES
   30 NN = NN+1
      POLD = PLAST
      PLAST = P
      P = DBLE(FLOAT(2*NN))*PLAST/TEMPA-RSIGN*POLD
      IF(P.LT.TEST) GO TO 30
      IF(M.EQ.1) GO TO 55
C                                  FOR J*S, A STRONG VARIANT OF THE TEST
C                                    IS NECESSARY. CALCULATE IT, AND
C                                    CALCULATE P*S UNTIL THIS TEST IS
C                                    PASSED.
      M = 1
      TEMPB = P/PLAST
      TEMPC = DBLE(FLOAT(NN+1))/TEMPA
      IF(TEMPB+1.D0/TEMPB.GT.2.D0*TEMPC)TEMPB=TEMPC+DSQRT(TEMPC**2-1.D0)
      TEST = TEST/DSQRT(TEMPB-1.D0/TEMPB)
      IF(P-TEST) 30, 55, 55
C                                  TO AVOID OVERFLOW, DIVIDE P*S BY
C                                    TOVER.  CALCULATE P*S UNTIL
C                                    ABS(P).GT.1.
   35 TOVER = 1.D1**RTEN
      P = P/TOVER
      PLAST = PLAST/TOVER
      PSAVE = P
      PSAVEL = PLAST
      NSTART = NN+1
   40 NN = NN+1
      POLD = PLAST
      PLAST = P
      P = DBLE(FLOAT(2*NN))*PLAST/TEMPA-RSIGN*POLD
      IF(P.LE.1.D0) GO TO 40
      TEMPB = DBLE(FLOAT(2*NN))/TEMPA
      TEMPC = .5D0*TEMPB
      TEMPB = PLAST/POLD
      IF(TEMPB+1.D0/TEMPB.GT.2.D0*TEMPC)TEMPB=TEMPC+DSQRT(TEMPC**2-1.D0)
C
C                                  CALCULATE BACKWARD TEST, AND FIND
C                                    NCALC, THE HIGHEST NN SUCH THAT THE
C                                    TEST IS PASSED.
      TEST = .5D0*POLD*PLAST*(1.D0-1.D0/TEMPB**2)/1.D1**DSIG
      P = PLAST*TOVER
      NN = NN-1
      NEND = MIN0(N,NN)
      DO 45 NCALC=NSTART,NEND
         POLD = PSAVEL
         PSAVEL = PSAVE
         PSAVE = DBLE(FLOAT(2*NN))*PSAVEL/TEMPA-RSIGN*POLD
         IF(PSAVE*PSAVEL-TEST) 45, 45, 50
   45 CONTINUE
      NCALC = NEND+1
   50 NCALC = NCALC-1
C                                  THE SUM B(1)+2B(3)+2B(5)... IS USED
C                                    TO NORMALIZE. M, THE COEFFICIENT OF
C                                    B(NN), IS INITIALIZED TO 2 OR 0.
   55 NN = NN+1
      M = 2*NN-4*(NN/2)
C                                  INITIALIZE THE BACKWARD RECURSION AND
C                                    THE NORMALIZATION SUM
      TEMPB = 0.D0
      TEMPA = 1.D0/P
      SUM = DBLE(FLOAT(M))*TEMPA
      NEND = NN-N
      IF(NEND) 80, 70, 60
C                                  RECUR BACKWARD VIA DIFFERENCE
C                                    EQUATION, CALCULATING (BUT NOT
C                                    STORING) B(NN), UNTIL NN=N.
   60 DO 65 L=1,NEND
         NN = NN-1
         TEMPC = TEMPB
         TEMPB = TEMPA
         TEMPA = (DBLE(FLOAT(2*NN))*TEMPB)/ARG-RSIGN*TEMPC
         M = 2-M
         SUM = SUM+DBLE(FLOAT(M))*TEMPA
   65 CONTINUE
C                                  STORE B(NN)
   70 B(NN) = TEMPA
      IF(N.GT.1) GO TO 75
C                                  N=1.  SINCE 2*TEMPA IS ADDED TO THE
C                                    SUM, TEMPA MUST BE SUBTRACTED
      SUM = SUM-TEMPA
      GO TO 110
C                                  CALCULATE AND STORE B(NN-1)
   75 NN = NN-1
      B(NN) = (DBLE(FLOAT(2*NN))*TEMPA)/ARG-RSIGN*TEMPB
      IF(NN.EQ.1) GO TO 105
      M = 2-M
      SUM = SUM+DBLE(FLOAT(M))*B(NN)
      GO TO 90
C                                  NN.LT.N, SO STORE B(NN) AND SET
C                                  HIGHER ORDERS TO ZERO
   80 B(NN) = TEMPA
      NEND = -NEND
      DO 85 L=1,NEND
         ITEMP = NN+L
         B(ITEMP) = 0.0D0
   85 CONTINUE
   90 NEND = NN-2
      IF(NEND.EQ.0) GO TO 100
C                                  CALCULATE VIA DIFFERENCE EQUATION AND
C                                    STORE B(NN), UNTIL NN=2
      DO 95 L=1,NEND
         NN = NN-1
         B(NN) = (DBLE(FLOAT(2*NN))*B(NN+1))/ARG-RSIGN*B(NN+2)
         M = 2-M
         SUM = SUM+DBLE(FLOAT(M))*B(NN)
   95 CONTINUE
C                                  CALCULATE B(1)
  100 B(1) = 2.D0*B(2)/ARG-RSIGN*B(3)
  105 SUM = SUM+B(1)
C                                  NORMALIZE--IF IZE=1, DIVIDE SUM BY
C                                    COSH(ARG). DIVIDE ALL B(NN) BY SUM.
  110 CONTINUE
      DO 115 NN=1,N
  115 B(NN) = B(NN)/SUM
      IF(NCALC.EQ.N) GO TO 9005
      IER = 129+NCALC
 9000 CONTINUE
      CALL UERTST(IER,6HMMBSJN)
 9005 RETURN
      END

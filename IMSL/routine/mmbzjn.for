C   IMSL ROUTINE NAME   - MMBZJN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - BESSEL FUNCTION OF THE FIRST KIND OF
C                           NONNEGATIVE INTEGER ORDER FOR
C                           COMPLEX ARGUMENTS
C
C   USAGE               - CALL MMBZJN (X,Y,N,BR,BI,IER)
C
C   ARGUMENTS    X      - INPUT REAL PART OF THE COMPLEX ARGUMENT. X
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.) SEE THE DESCRIPTION OF IER BELOW
C                           FOR RESTRICTIONS ON THE SIZE OF X.
C                Y      - INPUT IMAGINARY PART OF THE COMPLEX ARGUMENT.
C                           Y MUST BE LESS THAN EXPARG IN ABSOLUTE
C                           VALUE. EXPARG IS APPROXIMATELY THE LARGEST
C                           ACCEPTABLE ARGUMENT FOR THE FORTRAN EXP
C                           FUNCTION. SEE THE PROGRAMMING NOTES IN THE
C                           MANUAL FOR THE EXACT VALUE OF EXPARG.
C                           Y MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.) SEE THE DESCRIPTION OF IER BELOW
C                           FOR FURTHER RESTRICTIONS ON THE SIZE OF Y.
C                N      - INPUT PARAMETER SPECIFYING THE NUMBER OF
C                           FUNCTION VALUES TO BE COMPUTED.
C                BR     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           REAL PART OF THE COMPUTED FUNCTION VALUES.
C                           BR MUST BE TYPED APPROPRIATELY IN THE
C                           CALLING PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.) BR(1) WILL CONTAIN THE REAL PART
C                           OF THE COMPUTED VALUE FOR ORDER ZERO, BR(2)
C                           WILL CONTAIN THE REAL PART OF THE COMPUTED
C                           VALUE FOR ORDER 1, BR(3) FOR ORDER 2, ETC.
C                BI     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           IMAGINARY PART OF THE COMPUTED FUNCTIONS.
C                           BI MUST BE TYPED APPROPRIATELY IN THE
C                           CALLING PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER X, Y,
C                             OR N IS OUT OF RANGE. BR(I), (I=1,N), IS
C                             SET TO MACHINE INFINITY. THE VALUE,
C                             SQRT(X*X + Y*Y), MUST BE LESS THAN
C                             OR EQUAL TO 10000.
C                           IER = 129 + J INDICATES THAT COMPUTED
C                             FUNCTION VALUES FOR I = 1, J ARE COMPUTED
C                             TO MACHINE PRECISION, BUT PRECISION IS
C                             LOST FOR COMPUTED FUNCTION VALUES
C                             IN THE RANGE I = J+1, N.
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
      SUBROUTINE MMBZJN (X,Y,N,BR,BI,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      DOUBLE PRECISION   X,Y,BR(N),BI(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LARGEZ,MAGZ,NCALC,NBMZ,NN,M,NSTART,NEND,
     *                   IPOS,MRECUR,K,L,MLAST,IMRECR,ITEMP,JMAG,I
      DOUBLE PRECISION   PR,PI,PLASTR,PLASTI,POLDR,POLDI,PSAVER,PSAVEI,
     *                   EXPARG,TEST,TOVER,TEMPAR,TEMPAI,TEMPBR,TEMPBI,
     *                   TEMPCR,TEMPCI,RSIGN,SUMR,SUMI,ZINVR,ZINVI,DSIG,
     *                   RTEN,XINF,TMPAR4,SMALLZ,SRPRI
C
C                                  MACHINE DEPENDENT CONSTANTS
C
      DATA               XINF/.723700557733226D+76/
      DATA               DSIG/16.0D0/
      DATA               RTEN/75.D0/
      DATA               LARGEZ /10000/
      DATA               EXPARG/174.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      TEMPAR = DSQRT(X*X + Y*Y)
      MAGZ = IFIX(SNGL(TEMPAR))
      IF(N.GT.0 .AND. MAGZ.LE.LARGEZ .AND.((DABS(Y).LE.EXPARG)))GO TO 10
C                                  ERROR RETURN - Z OR N IS OUT OF RANGE
      NN = N
      IF(N.LE.0) NN = 1
      DO 5  I=1,NN
         BR(I) = XINF
    5 CONTINUE
      IER = 129
      GO TO 9000
   10 RSIGN = 1.D0
      NCALC = N
C                                  USE 2-TERM ASCENDING SERIES FOR
C                                    SMALL Z
      TMPAR4 = TEMPAR**4
      SMALLZ = .1D0**DSIG
      IF(TMPAR4 .LT. SMALLZ) GO TO 175
C                                  INITIALIZE THE CALCULTION OF THE P*S
      NBMZ = N-MAGZ
      NN = MAGZ+1
      IF(DABS(X) .LT. DABS(Y)) GO TO 15
      ZINVR = 1.D0/(X+Y*Y/X)
      ZINVI = -Y*ZINVR/X
      GO TO 20
   15 ZINVI = -1.D0/(Y+X*X/Y)
      ZINVR = -X*ZINVI/Y
   20 PLASTR = 1.D0
      PLASTI = 0.D0
      PR = RSIGN*DBLE(FLOAT(2*NN))*ZINVR
      PI = RSIGN*DBLE(FLOAT(2*NN))*ZINVI
      TEST = 2.D0*1.D1**DSIG
      M = 0
      IF(NBMZ .LT. 3) GO TO 30
C                                  CALCULATE P*S UNTIL NN=N-1. CHECK FOR
C                                    POSSIBLE OVERFLOW
      TOVER = 1.D1**(RTEN-DSIG)
      NSTART = MAGZ+2
      NEND = N-1
      DO 25  NN=NSTART,NEND
         POLDR = PLASTR
         POLDI = PLASTI
         PLASTR = PR
         PLASTI = PI
         PR=RSIGN*(DBLE(FLOAT(2*NN))*(PLASTR*ZINVR-PLASTI*ZINVI)-POLDR)
         PI=RSIGN*(DBLE(FLOAT(2*NN))*(PLASTI*ZINVR+PLASTR*ZINVI)-POLDI)
         SRPRI = DSQRT(PR**2 + PI**2)
         IF(SRPRI - TOVER) 25, 25, 35
   25 CONTINUE
      NN = NEND
C                                  CALCULATE SPECIAL SIGNIFICANCE TEST
C                                    FOR NBMZ .GT. 2
      TEMPBI = DMAX1(DABS(PR),DABS(PI))
      TEMPBI = TEMPBI*DSQRT(2.D0*1.D1**DSIG*DSQRT(((PR/TEMPBI)**2+(PI/
     *         TEMPBI)**2)*((PLASTR/TEMPBI)**2+(PLASTI/TEMPBI)**2)))
      TEST = DMAX1(TEST,TEMPBI)
C                                  CALCULATE P*S UNTIL SIGNIFICANCE
C                                    TEST IS PASSED.
   30 NN = NN+1
      POLDR = PLASTR
      POLDI = PLASTI
      PLASTR = PR
      PLASTI = PI
      PR = RSIGN*(DBLE(FLOAT(2*NN))*(PLASTR*ZINVR-PLASTI*ZINVI)-POLDR)
      PI = RSIGN*(DBLE(FLOAT(2*NN))*(PLASTI*ZINVR+PLASTR*ZINVI)-POLDI)
      IF((PR/TEST)**2+(PI/TEST)**2.LT.1.D0) GO TO 30
      IF(M.EQ.1) GO TO 55
C                                  CALCULATE STRICT VARIANT OF
C                                    SIGNIFICANCE TEST. AND CALCULATE
C                                    P*S UNTIL THIS TEST IS PASSED.
      M = 1
      TEMPBI = DMAX1(DABS(PR),DABS(PI))
      TEMPBR = DSQRT(((PR/TEMPBI)**2+(PI/TEMPBI)**2)/((PLASTR/TEMPBI)**
     *         2+(PLASTI/TEMPBI)**2))
      TEMPBI = DBLE(FLOAT(NN+1))/TEMPAR
      IF(TEMPBR+1.D0/TEMPBR.GT.2.D0*TEMPBI) TEMPBR =
     *   TEMPBI+DSQRT(TEMPBI**2-1.D0)
      TEST = TEST/DSQRT(TEMPBR-1.D0/TEMPBR)
      IF((PR/TEST)**2+(PI/TEST)**2-1.D0) 30, 55, 55
   35 NSTART = NN+1
C                                  TO AVOID OVERFLOW, NORMALIZE P*S BY
C                                    DIVIDING BY TOVER. CALCULATE P*S
C                                    UNTIL UNNORMALIZED P WOULD OVERFLOW
      PR = PR/TOVER
      PI = PI/TOVER
      PLASTR = PLASTR/TOVER
      PLASTI = PLASTI/TOVER
      PSAVER = PR
      PSAVEI = PI
      TEMPCR = PLASTR
      TEMPCI = PLASTI
      TEST = 1.D1**(DSIG+DSIG)
   40 NN = NN+1
      POLDR = PLASTR
      POLDI = PLASTI
      PLASTR = PR
      PLASTI = PI
      PR = RSIGN*(DBLE(FLOAT(2*NN))*(PLASTR*ZINVR-PLASTI*ZINVI)-POLDR)
      PI = RSIGN*(DBLE(FLOAT(2*NN))*(PLASTI*ZINVR+PLASTR*ZINVI)-POLDI)
      TMPAR4 = PR*PR+PI*PI
      IF(TMPAR4 .LE. TEST) GO TO 40
C                                  CALCULATE BACKWARD TEST, AND FIND
C                                    NCALC, THE HIGHEST NN SUCH THAT THE
C                                    TEST IS PASSED.
      TEMPBR = DSQRT((PLASTR**2+PLASTI**2)/(POLDR**2+POLDI**2))
      TEMPBI = DBLE(FLOAT(NN))/TEMPAR
      IF(TEMPBR+1.D0/TEMPBR.GT.2.D0*TEMPBI) TEMPBR =
     *   TEMPBI+DSQRT(TEMPBI**2-1.D0)
      TEST = .5D0*(1.D0-1.D0/TEMPBR**2)/1.D1**DSIG
      TEST = ((PLASTR**2+PLASTI**2)*TEST)*((POLDR**2+POLDI**2)*TEST)
      PR = PLASTR*TOVER
      PI = PLASTI*TOVER
      NN = NN-1
      NEND = MIN0(N,NN)
      DO 45  NCALC=NSTART,NEND
         POLDR = TEMPCR
         POLDI = TEMPCI
         TEMPCR = PSAVER
         TEMPCI = PSAVEI
         PSAVER = RSIGN*(DBLE(FLOAT(2*NN))*(TEMPCR*ZINVR-TEMPCI*ZINVI)
     *            - POLDR)
         PSAVEI = RSIGN*(DBLE(FLOAT(2*NN))*(TEMPCI*ZINVR+TEMPCR*ZINVI)
     *            - POLDI)
         IF((PSAVER**2+PSAVEI**2)*(TEMPCR**2+TEMPCI**2)-TEST)45,45,50
   45 CONTINUE
      NCALC = NEND+1
   50 NCALC = NCALC-1
C                                  THE COEFFICIENT OF B(NN) IN THE
C                                    NORMALIZATION SUM IS
C                                    M*SQRT(-1)**JMAG, WHERE M=-2,0,
C                                    OR 2, AND IMAG IS 0 OR 1. CALCULATE
C                                    RECURSION RULES FOR M AND JMAG, AND
C                                    INITIALIZE THEM.
   55 NN = NN+1
      TEMPBR = Y
      IPOS = 0
      IF(TEMPBR) 60, 65, 60
   60 IPOS = IFIX(SNGL(1.1D0*TEMPBR/DABS(TEMPBR)))
   65 MRECUR = 4*((2+IPOS)/2)-3-2*(IPOS)
      K = 2+IPOS
      L = NN-4*(NN/4)
      MLAST = 2+8*((K*L)/4)-4*((K*L)/2)
      IF(IPOS.EQ.0 .AND. (L.EQ.1 .OR. L.EQ.3)) MLAST = 0
      L = L+3-4*((L+3)/4)
      M = 2+8*((K*L)/4)-4*((K*L)/2)
      IF(IPOS.EQ.0 .AND. (L.EQ.1 .OR. L.EQ.3)) M = 0
      IMRECR = IPOS*IPOS
      JMAG = IMRECR*(L-2*(L/2))
C                                  INITIALIZE THE BACKWARD RECURSION AND
C                                    THE NORMALIZATION SUM.
      TEMPBR = 0.D0
      TEMPBI = 0.D0
      IF(DABS(PI).GT.DABS(PR)) GO TO 70
      TEMPAR = 1.D0/(PR+PI*(PI/PR))
      TEMPAI = -(PI*TEMPAR)/PR
      GO TO 75
   70 TEMPAI = -1.D0/(PI+PR*(PR/PI))
      TEMPAR = -(PR*TEMPAI)/PI
   75 IF(JMAG.NE.0) GO TO 80
      SUMR = DBLE(FLOAT(M))*TEMPAR
      SUMI = DBLE(FLOAT(M))*TEMPAI
      GO TO 85
   80 SUMR = -DBLE(FLOAT(M))*TEMPAI
      SUMI = DBLE(FLOAT(M))*TEMPAR
   85 NEND = NN-N
      IF(NEND) 120, 105, 90
C                                  RECUR BACKWARD VIA DIFFERENCE
C                                    EQUATION CALCULATION (BUT NOT
C                                    STORING) BR(NN) AND BI(NN) UNTIL
C                                    NN = N.
   90 DO 100  L=1,NEND
         NN = NN-1
         TEMPCR = TEMPBR
         TEMPCI = TEMPBI
         TEMPBR = TEMPAR
         TEMPBI = TEMPAI
         PR = DBLE(FLOAT(2*NN))*ZINVR
         PI = DBLE(FLOAT(2*NN))*ZINVI
         TEMPAR = PR*TEMPBR-PI*TEMPBI-RSIGN*TEMPCR
         TEMPAI = PR*TEMPBI+PI*TEMPBR-RSIGN*TEMPCI
         JMAG = (1-JMAG)*IMRECR
         K = MLAST
         MLAST = M
         M = K*MRECUR
         IF(JMAG.NE.0) GO TO 95
         SUMR = SUMR+DBLE(FLOAT(M))*TEMPAR
         SUMI = SUMI+DBLE(FLOAT(M))*TEMPAI
         GO TO 100
   95    SUMR = SUMR-DBLE(FLOAT(M))*TEMPAI
         SUMI = SUMI+DBLE(FLOAT(M))*TEMPAR
  100 CONTINUE
C                                  STORE BR(NN),BI(NN)
  105 BR(NN) = TEMPAR
      BI(NN) = TEMPAI
      IF(NN.GT.1) GO TO 110
C                                  N=1.  SINCE 2*TEMPAR AND 2*TEMPAI
C                                    WERE ADDED TO SUMR AND SUMI
C                                    RESPECTIVELY, WE MUST SUBTRACT
C                                    TEMPAR AND TEMPAI
      SUMR = SUMR-TEMPAR
      SUMI = SUMI-TEMPAI
      GO TO 155
C                                  CALCULATE AND STORE BR(NN-1),(NN-1)
  110 NN = NN-1
      PR = DBLE(FLOAT(2*NN))*ZINVR
      PI = DBLE(FLOAT(2*NN))*ZINVI
      BR(NN) = PR*TEMPAR-PI*TEMPAI-RSIGN*TEMPBR
      BI(NN) = PR*TEMPAI+PI*TEMPAR-RSIGN*TEMPBI
      IF(NN.EQ.1) GO TO 150
      JMAG = (1-JMAG)*IMRECR
      K = MLAST
      MLAST = M
      M = K*MRECUR
      IF(JMAG.NE.0) GO TO 115
      SUMR = SUMR+DBLE(FLOAT(M))*BR(NN)
      SUMI = SUMI+DBLE(FLOAT(M))*BI(NN)
      GO TO 130
  115 SUMR = SUMR-DBLE(FLOAT(M))*BI(NN)
      SUMI = SUMI+DBLE(FLOAT(M))*BR(NN)
      GO TO 130
C                                  NN .LT. N, SO STORE BR(NN), AND SET
C                                    HIGHER ORDERS ZERO
  120 BR(NN) = TEMPAR
      BI(NN) = TEMPAI
      NEND = -NEND
      DO 125  L=1,NEND
         ITEMP = NN+L
         BR(ITEMP) = 0.0D0
         BI(ITEMP) = 0.0D0
  125 CONTINUE
  130 NEND = NN-2
      IF(NEND.EQ.0) GO TO 145
C                                  CALCULATE VIA DIFFERENCE EQUATION AND
C                                    STORE BR(NN),BI(NN), UNTIL NN=2
      DO 140  L=1,NEND
         NN = NN-1
         PR = DBLE(FLOAT(2*NN))*ZINVR
         PI = DBLE(FLOAT(2*NN))*ZINVI
         BR(NN) = PR*BR(NN+1)-PI*BI(NN+1)-RSIGN*BR(NN+2)
         BI(NN) = PR*BI(NN+1)+PI*BR(NN+1)-RSIGN*BI(NN+2)
         JMAG = (1-JMAG)*IMRECR
         K = MLAST
         MLAST = M
         M = K*MRECUR
         IF(JMAG.NE.0) GO TO 135
         SUMR = SUMR+DBLE(FLOAT(M))*BR(NN)
         SUMI = SUMI+DBLE(FLOAT(M))*BI(NN)
         GO TO 140
  135    SUMR = SUMR-DBLE(FLOAT(M))*BI(NN)
         SUMI = SUMI+DBLE(FLOAT(M))*BR(NN)
  140 CONTINUE
C                                  CALCULUTE AND STORE BR(1), BI(1)
  145 BR(1) = 2.D0*(BR(2)*ZINVR-BI(2)*ZINVI)-RSIGN*BR(3)
      BI(1) = 2.D0*(BR(2)*ZINVI+BI(2)*ZINVR)-RSIGN*BI(3)
  150 SUMR = SUMR+BR(1)
      SUMI = SUMI+BI(1)
C                                  CALCULATE NORMALIZATION FACTOR,
C                                    TEMPAR + I*TEMPAI
  155 TEMPCR = DBLE(FLOAT(IPOS))*Y
      TEMPCI = DBLE(FLOAT(-IPOS))*X
      TEMPCR = DEXP(TEMPCR)
      TEMPBR = DCOS(TEMPCI)
      TEMPBI = DSIN(TEMPCI)
      IF(DABS(SUMR).LT.DABS(SUMI)) GO TO 160
      TEMPCI = SUMI/SUMR
      TEMPCR = (TEMPCR/SUMR)/(1.D0+TEMPCI*TEMPCI)
      TEMPAR = TEMPCR*(TEMPBR+TEMPBI*TEMPCI)
      TEMPAI = TEMPCR*(TEMPBI-TEMPBR*TEMPCI)
      GO TO 165
  160 TEMPCI = SUMR/SUMI
      TEMPCR = (TEMPCR/SUMI)/(1.D0+TEMPCI*TEMPCI)
      TEMPAR = TEMPCR*(TEMPBR*TEMPCI+TEMPBI)
      TEMPAI = TEMPCR*(TEMPBI*TEMPCI-TEMPBR)
C                                  NORMALIZE
  165 DO 170  NN=1,N
         TEMPBR = BR(NN)*TEMPAR-BI(NN)*TEMPAI
         BI(NN) = BR(NN)*TEMPAI+BI(NN)*TEMPAR
         BR(NN) = TEMPBR
  170 CONTINUE
      GO TO 185
C                                  TWO-TERM ASCENDING SERIES FOR SMALL Z
  175 TEMPAR = 1.D0
      TEMPAI = 0.D0
      TEMPCR = .25D0*(X*X-Y*Y)
      TEMPCI = .5D0*X*Y
      BR(1) = 1.D0-RSIGN*TEMPCR
      BI(1) = -RSIGN*TEMPCI
      IF(N.EQ.1) GO TO 185
      DO 180  NN=2,N
         TEMPBR = (TEMPAR*X-TEMPAI*Y)/DBLE(FLOAT(2*NN-2))
         TEMPAI = (TEMPAR*Y+TEMPAI*X)/DBLE(FLOAT(2*NN-2))
         TEMPAR = TEMPBR
         TEMPBR = DBLE(FLOAT(NN))
         BR(NN) = TEMPAR*(1.D0-RSIGN*TEMPCR/TEMPBR)+TEMPAI*TEMPCI/TEMPBR
         BI(NN) = TEMPAI*(1.D0-RSIGN*TEMPCR/TEMPBR)-TEMPAR*TEMPCI/TEMPBR
  180 CONTINUE
  185 CONTINUE
      IF(N.EQ.NCALC) GO TO 9005
      IER = 129+NCALC
 9000 CONTINUE
      CALL UERTST(IER,6HMMBZJN)
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - MMBSJR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - BESSEL FUNCTION OF THE FIRST KIND OF
C                           NONNEGATIVE REAL ORDER FOR REAL POSITIVE
C                           ARGUMENTS
C
C   USAGE               - CALL MMBSJR (ARG,ORDER,N,RJ,WK,IER)
C
C   ARGUMENTS    ARG    - INPUT ARGUMENT. ARG MUST BE GREATER THAN OR
C                           EQUAL TO ZERO. ARG MUST BE TYPED APPROPRI-
C                           ATELY IN THE CALLING PROGRAM. (SEE THE
C                           PRECISION/HARDWARE SECTION.)
C                ORDER  - INPUT VALUE SPECIFYING THE DESIRED ORDER OF
C                           THE BESSEL FUNCTION. ORDER MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) ORDER
C                           MUST BE GREATER THAN OR EQUAL TO ZERO AND
C                           LESS THAN ONE.
C                N      - INPUT PARAMETER SPECIFYING THE NUMBER OF
C                           FUNCTION VALUES TO BE COMPUTED.
C                RJ     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           COMPUTED FUNCTION VALUES. RJ MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) RJ(1) WILL
C                           CONTAIN THE COMPUTED VALUE FOR THE INPUT
C                           ORDER, RJ(2) WILL CONTAIN THE COMPUTED
C                           VALUE FOR ORDER + 1, RJ(3) FOR ORDER + 2,
C                           ETC.
C                WK     - WORK VECTOR OF LENGTH 2*N. WK MUST BE
C                           TYPED APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE THE PRECISION/HARDWARE SECTION.)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER ARG, ORDER,
C                             OR N IS OUT OF RANGE. RJ(I), (I= 1,N) IS
C                             SET TO MACHINE INFINITY.
C                           IER = 130 INDICATES THAT OVERFLOW WOULD
C                             HAVE OCCURED IN SOME COMPUTATIONS.  RJ(I),
C                             (I= 1,N) IS SET TO MACHINE INFINITY.
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
      SUBROUTINE MMBSJR (ARG,ORDER,N,RJ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      DOUBLE PRECISION   ORDER,ARG,WK(1),RJ(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IOVER,NM1,NN,NP1,NP2,MHIT,I,M,LIMIT,NU,NPNN,II
      DOUBLE PRECISION   EPSLON,GAMA1,SUM,D1,R,Y,T,P,Z,RLGSML
      DOUBLE PRECISION   S,DIV,TOL,RLAMDA,RN,EPSLJN,EPSABS,DIFF
      DOUBLE PRECISION   XINF,RM,RL,D,RNMAX,OVER,SUMCK,PP
      DOUBLE PRECISION   DGAMMA
      DATA               IOVER/Z7FFFFFFF/
      DATA               RLGSML/ -180.2182D0/
      DATA               XINF/.723700557733226D+76/
      DATA               D/16.0D0/
      DATA               EPSLON/ Z3410000000000000/
      DATA               TOL/ Z0D10000000000000/
C                                  FIRST EXECUTABLE STATEMENT
      OVER = DBLE(FLOAT(IOVER))
      IER = 0
      NM1 = N-1
      IF(ORDER.LT.0.D0 .OR. ORDER.GE.1.D0 .OR. ARG.LT.0.D0 .OR. N.LT.1)
     *   GO TO 5
      GO TO 15
    5 RJ(1) = XINF
      IER = 129
      IF(N.LT.2) GO TO 9000
      DO 10 I=2,N
         RJ(I) = XINF
   10 CONTINUE
      GO TO 9000
   15 CONTINUE
      DO 20 NP1=1,N
         WK(NP1) = 0.D0
   20 CONTINUE
C                                  MAKE USE OF THE MGAMMA FUNCTION
C                                    IN THE LIBRARY
      SUM = ORDER+1.D0
      GAMA1 = DGAMMA(SUM)
      IF(ARG.EQ.0.D0) GO TO 30
C                                  CHECKING FOR UNDERFLOW IN SUM
C
      SUMCK = ORDER * (DLOG(ARG) + DLOG(.5D0))
      IF(SUMCK.LT.RLGSML) GO TO 25
      SUM = DEXP(SUMCK)
      GO TO 30
   25 SUM = 0.D0
   30 SUM = SUM/GAMA1
      D1 = 2.3026D0*D+1.3863D0
      R = 0.D0
      IF(N.LE.1) GO TO 45
      RNMAX = NM1
      Y = .5D0*D1/RNMAX
      IF(Y.GT.10.D0) GO TO 35
      P = 5.7941D-5*Y-1.76148D-3
      P = Y*P+2.08645D-2
      P = Y*P-.129013D0
      P = Y*P+.85777D0
      T = Y*P+1.0125D0
      GO TO 40
   35 Z = DLOG(Y)-.775D0
      P = (.775D0-DLOG(Z))/(1.D0+Z)
      P = 1.D0/(1+P)
      T = Y*P/Z
   40 R = RNMAX*T
C                                  CHECKING FOR OVERFLOW IN D1/ARG
   45 IF(ARG.GE.1) GO TO 50
      IF(D1.GT.ARG*XINF) GO TO 95
   50 Y = .73576D0*(D1/ARG)
      IF(Y.GT.10.D0) GO TO 80
      PP = -9.7561D0+DLOG(Y)
      P = 0.D0
      II = RLGSML/PP
      IF(II.LE.1) GO TO 75
      IF(II.GT.5) GO TO 55
      GO TO (75, 70, 65, 60, 55), II
   55 P = 5.7941D-5*Y-1.76148D-3
   60 P = Y*P+2.08645D-2
   65 P = Y*P-.129013D0
   70 P = Y*P+.85777D0
   75 T = Y*P+1.0125D0
      GO TO 85
   80 Z = DLOG(Y)-.775D0
      P = (.775D0-DLOG(Z))/(1.D0+Z)
      P = 1.D0/(1+P)
      T = Y*P/Z
C                                  CHECKING FOR OVERFLOW IN S
   85 IF(ARG.LE.1) GO TO 90
      IF(1.3591D0*T.GT.XINF/ARG) GO TO 95
   90 S = 1.3591D0*ARG*T
      IF(S.LE.OVER) GO TO 105
   95 IER = 130
      DO 100 I=1,N
         RJ(I) = XINF
  100 CONTINUE
      GO TO 9000
  105 NU = 1+S
      IF(R.GT.S) NU = 1+R
  110 M = 0
      RL = 1.D0
      LIMIT = (NU/2)
  115 M = M+1
      RM = M
      RL = RL*(RM+ORDER)/(RM+1.D0)
      IF(M.LT.LIMIT) GO TO 115
      NN = M+M
      R = 0.D0
      S = 0.D0
  120 DIV = (2.D0*(ORDER+NN)/ARG-R)
      TOL = TOL*DSIGN(1.D0,DIV)
      IF(DABS(DIV).LT.DABS(TOL)) DIV = TOL
      R = 1.D0/DIV
      RLAMDA = 0.D0
      IF(MOD(NN,2).EQ.1) GO TO 125
      RN = NN
      RL = RL*(RN+2.D0)/(RN+(ORDER+ORDER))
      RLAMDA = RL*(RN+ORDER)
  125 CONTINUE
      S = R*(RLAMDA+S)
      NPNN = N+NN
      IF(NN.LT.N) WK(NPNN) = R
      NN = NN-1
      IF(NN.GE.1) GO TO 120
      RJ(1) = SUM/(1.D0+S)
      IF(NM1.LE.0) GO TO 145
      N1 = MIN0(45,NM1)
      DO 130 NP1=1,N1
         NP2 = NP1+1
         NPNN = N+NP1
         RJ(NP2) = WK(NPNN)*RJ(NP1)
  130 CONTINUE
      IF(N1.EQ.NM1) GO TO 145
      N1 = N1+1
      DO 140 NP1=N1,NM1
         NP2 = NP1+1
         NPNN = N+NP1
C                                   CHECK FOR UNDERFLOW IN PRODUCT
         RJ(NP2) = 0.0D0
         IF(DABS(WK(NPNN)) .GE. 1.0D0) GO TO 135
         IF(DABS(RJ(NP1)) .LE. RLGSML/DABS(WK(NPNN))) GO TO 140
  135    RJ(NP2) = WK(NPNN)*RJ(NP1)
  140 CONTINUE
  145 CONTINUE
      DO 155 NP1=1,N
         DIFF = RJ(NP1)-WK(NP1)
         EPSABS = DABS(DIFF)
         EPSLJN = EPSLON*DABS(RJ(NP1))
         IF(EPSABS.LE.EPSLJN) GO TO 155
         DO 150 MHIT=1,N
            WK(MHIT) = RJ(MHIT)
  150    CONTINUE
         NU = NU+5
         GO TO 110
  155 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMBSJR)
 9005 RETURN
      END
 
R; T=0.06/0.55 00:25:35

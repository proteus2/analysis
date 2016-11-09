C   IMSL ROUTINE NAME   - FTML
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - MAXIMUM LIKELIHOOD ESTIMATION OF
C                           AUTOREGRESSIVE AND MOVING AVERAGE PARAMETERS
C                           IN AN ARIMA (BOX-JENKINS) STOCHASTIC MODEL
C
C   USAGE               - CALL FTML (X,IND,ARPS,PMAS,PMAC,WNV,GR,A,IER)
C
C   ARGUMENTS    X      - INPUT TIME SERIES OF LENGTH IND(1).
C                           X IS DESTROYED ON OUTPUT.
C                IND    - INPUT/OUTPUT VECTOR OF LENGTH 8.  IND(I) CON-
C                           TAINS WHEN
C                             I=1, LENGTH OF TIME SERIES X.
C                             I=2, NUMBER (NON-NEGATIVE) OF AUTOREGRES-
C                               SIVE PARAMETERS IN THE DIFFERENCED FORM
C                               OF THE ARIMA MODEL.
C                             I=3, NUMBER (NON-NEGATIVE) OF MOVING AVER-
C                               AGE PARAMETERS.
C                             IND(2)+IND(3) MUST BE POSITIVE.
C                             I=4, NUMBER (NON-NEGATIVE) OF DIFFERENCING
C                               OPERATIONS REQUIRED TO MAKE X STATION-
C                               ARY.  IF IND(4)=0 THE MEAN IS REMOVED
C                               FROM X.
C                             I=5, INPUT MAXIMUM NUMBER OF ITERATIONS
C                               DESIRED.  ON OUTPUT, IND(5) CONTAINS THE
C                               NUMBER OF ITERATIONS PERFORMED.
C                             I=6, NON-NEGATIVE CONVERGENCE PARAMETER.
C                               CONVERGENCE IS ASSUMED IF IND(6) SIGNIF-
C                               ICANT DIGITS OF THE OBJECTIVE FUNCTION
C                               DO NOT CHANGE AFTER IND(8) CONSECUTIVE
C                               ITERATIONS.
C                             I=7, IND(7) NONZERO IMPLIES INITIAL ESTI-
C                               MATES OF ARPS AND PMAS IN THE DIFFER-
C                               ENCED FORM OF THE MODEL ARE INPUT.
C                               IND(7)=0 IMPLIES FTML CALLS FTARPS AND
C                               FTMA TO CALCULATE INITIAL ESTIMATES.
C                             I=8, POSITIVE CONVERGENCE PARAMETER WHOSE
C                               FUNCTION IS DESCRIBED UNDER IND(6).
C                ARPS   - INPUT/OUTPUT VECTOR OF LENGTH IND(2)+IND(4).
C                           ON INPUT, IF IND(7) IS NONZERO, THE FIRST
C                           IND(2) LOCATIONS SHOULD CONTAIN INITIAL ES-
C                           TIMATES OF THE AUTOREGRESSIVE PARAMETERS IN
C                           THE DIFFERENCED FORM OF THE MODEL.  ON OUT-
C                           PUT, THE IND(2)+IND(4) PARAMETER ESTIMATES
C                           ARE FOR THE UNDIFFERENCED FORM OF THE MODEL.
C                PMAS   - INPUT/OUTPUT VECTOR OF LENGTH IND(3).  ON IN-
C                           PUT, IF IND(7) IS NONZERO, PMAS SHOULD CON-
C                           TAIN INITIAL ESTIMATES OF THE MOVING AVERAGE
C                           PARAMETERS.  ON OUTPUT, THE PARAMETER ESTI-
C                           MATES ARE RETURNED.
C                PMAC   - OUTPUT ESTIMATE OF OVERALL MOVING AVERAGE CON-
C                           STANT IN THE UNDIFFERENCED MODEL.
C                WNV    - OUTPUT ESTIMATE OF THE WHITE NOISE VARIANCE.
C                GR     - WORK AREA OF LENGTH 2*(IND(2)+IND(3)).
C                A      - WORK AREA OF LENGTH THE MAXIMUM OF
C                           1. (IND(2)+6)*IND(2)+IND(3)+1
C                           2. IND(2)+(IND(3)+1)*(3*IND(3)+24)/2
C                           3. 2*IND(1)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES AT LEAST ONE ELEMENT OF
C                             IND WAS OUT OF RANGE
C                           IER=130 INDICATES AN ERROR OCCURRED IN IMSL
C                             ROUTINE FTARPS.
C                         WARNING WITH FIX
C                           IER=67 INDICATES AN ERROR OCCURRED IN IMSL
C                             ROUTINE FTMA.  INITIAL PMAS ESTIMATES
C                             ARE SET TO ZERO.
C                           IER=68 INDICATES THAT ALL IND(5) ITERATIONS
C                             WERE PERFORMED.  CONVERGENCE IS ASSUMED
C                             AND PROGRAM CONTINUES CALCULATIONS.
C
C   REQD. IMSL ROUTINES - SINGLE/FTARPS,FTAUTO,FTMA,FTMA1,LEQT1F,
C                           LUDATN,LUELMN,UERTST,UGETIO,VABMXF,
C                           VBLA=SNRM2,ZSPOW,ZSPWA,ZSPWB,ZSPWC,ZSPWD,
C                           ZSPWE,ZSPWF,ZSPWG
C                       - DOUBLE/FTARPS,FTAUTO,FTMA,FTMA1,LEQT1F,
C                           LUDATN,LUELMN,UERTST,UGETIO,VABMXF,
C                           VBLA=DNRM2,VXADD,VXMUL,VXSTO,ZSPOW,ZSPWA,
C                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ESTIMATES OF THE RESIDUALS OR ONE-STEP FORECASTING
C                ERRORS ARE CONTAINED IN THE FIRST IND(1)-IND(4)
C                LOCATIONS OF WORK VECTOR A.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTML (X,IND,ARPS,PMAS,PMAC,WNV,GR,A,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,IND(8)
      REAL               X(1),ARPS(1),PMAS(1),PMAC,WNV,GR(1),A(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I0,I3,ICOUNT,ID,IP,IPN,IPQ,IQ,IRS,J,K,KIP,L,
     *                   LP2,N
      REAL               ALP,BK,DEL,EPS,FOURTH,G,HUNTH,ONE,SINIT,TEN,
     *                   TENTH,TOL,XBAR,Z0,ZERO
      DOUBLE PRECISION   C,HALF,S1,S2,S3,SSAVE,T,T1,TEMP,W
      DOUBLE PRECISION   DZERO
      DATA               I0 /0/,I3 /3/,Z0 /0.0/,TEN /10.0/,ZERO
     *                   /0.0/,FOURTH /0.25/,TENTH /0.1/,ONE
     *                   /1.0/,HUNTH /0.01/
      DATA               HALF /0.5D0/,TOL /1.E-6/,DZERO /0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  PARAMETER DEFINITIONS
      N = IND(1)
      IP = IND(2)
      IQ = IND(3)
      ID = IND(4)
      IPQ = IP+IQ
      IRS = MAX0(IP,IQ)+1
      EPS = TEN**(-IND(6))
C                                  ERROR CHECKS
      IF (IND(2).GE.0 .AND. IND(3).GE.0 .AND. IND(4).GE.0 .AND.
     *IND(6).GE.0 .AND. IND(8).GT.0 .AND. IND(2)+IND(3).GT.0) GO TO 5
      IER = 129
      GO TO 9000
    5 IF (ID.EQ.0) GO TO 20
      DO 15 K=1,ID
         N = N-1
         DO 10 I=1,N
            X(I) = X(I+1)-X(I)
   10    CONTINUE
   15 CONTINUE
   20 CALL FTAUTO(X,N,IPQ,I0,I3,XBAR,A(1),A(2),GR,GR,GR)
      DO 25 I=1,N
         X(I) = X(I)-XBAR
   25 CONTINUE
C             CALCULATE INITIAL PARAMETER ESTIMATES
      IF (IND(7).NE.0) GO TO 40
      IF (IP.EQ.0) GO TO 30
      CALL FTARPS(A,Z0,IP,IQ,ARPS,PMAC,A(2+IPQ),IER)
      IF (IER.EQ.0) GO TO 30
      IER = 130
      GO TO 9000
   30 IF (IQ.EQ.0) GO TO 40
      CALL FTMA(A,ARPS,IP,IQ,PMAS,WNV,A(2+IPQ),IER)
      IF (IER.EQ.0) GO TO 40
      IER = 67
      DO 35 I=1,IQ
         PMAS(I) = ZERO
   35 CONTINUE                                                          FTML1610
C                      INITIALIZE ITERATION
   40 DEL = TENTH
      ICOUNT = 0
      BK = ONE
      DO 45 I=1,IRS
         A(I) = ZERO
   45 CONTINUE
      DO 50 I=1,IPQ
         GR(I+IPQ) = ZERO
   50 CONTINUE
      ALP = ONE
C              MODIFIED STEEPEST DESCENT ALGORTIHM
   55 T = ALP
      L = 3
      GO TO 115
   60 S3 = W
      IF (ICOUNT.GT.0) GO TO 65
      SSAVE = S3
      SINIT = S3
      S1 = S3+S3
      S2 = S1
   65 ICOUNT = ICOUNT+1
C                      FIND BEST DESCENT POINT
      K = 3
      IF (SNGL(S2-S3).LT.ZERO) K = 2
      IF (K.EQ.2 .AND. SNGL(S1-S2).LT.ZERO) K = 1
      IF (K.EQ.3 .AND. SNGL(S1-S3).LT.ZERO) K = 1
      IF (K.GT.1) GO TO 70
      ALP = 0.0
      GO TO 90
   70 IF (K.EQ.3) GO TO 80
      T = BK
      L = 2
      GO TO 115
   75 S1 = S2
      ALP = BK
      GO TO 90
   80 T = ALP
      L = 1
      GO TO 115
   85 S1 = S3
C            CONVERGENCE CHECKS
   90 IF (ALP.LT.FOURTH*BK) DEL = DEL*FOURTH
      IF (DEL.LT.TOL) GO TO 210
      IF (DEL.GT.HUNTH .OR. ICOUNT.LT.5) GO TO 95
      IF (MOD(ICOUNT,IND(8)).NE.0) GO TO 95
      IF (DABS(S1-SSAVE).LT.EPS*S1) GO TO 210
      SSAVE = S1
   95 IF (ICOUNT.LE.IND(5)) GO TO 100
      IER = 68
      GO TO 210
C              SCALE GRADIENT TO STEPSIZE
  100 C = DZERO
      G = ZERO
      DO 105 I=1,IPQ
         C = C+DBLE(GR(I+IPQ))**2
         G = AMAX1(G,ABS(GR(I+IPQ)))
  105 CONTINUE
      C = DSQRT(C)
      BK = DEL/G
      T = BK
      L = 4
      GO TO 115
  110 S2 = W
      ALP = HALF*C**2*DBLE(BK)**2/(S2-S1+BK*C**2)
      GO TO 55
C             PERTURB PARAMETERS VIA L
  115 IF (IP.EQ.0) GO TO 125
      DO 120 I=1,IP
         GR(I) = ARPS(I)-T*GR(I+IPQ)
         IF (L.LT.3) ARPS(I) = GR(I)
  120 CONTINUE
      IF (IQ.EQ.0) GO TO 135
  125 DO 130 I=1,IQ
         GR(I+IP) = PMAS(I)-T*GR(I+IP+IPQ)
         IF (L.LT.3) PMAS(I) = GR(I+IP)
  130 CONTINUE
C                       CALCULATE RESIDUALS OF ARIMA MODEL
  135 IF (L.EQ.1) GO TO 165
      W = DZERO
      DO 160 I=IRS,N
         TEMP = X(I)
         IF (IP.EQ.0) GO TO 145
         DO 140 J=1,IP
            TEMP = TEMP-DBLE(GR(J))*X(I-J)
  140    CONTINUE
         IF (IQ.EQ.0) GO TO 155
  145    DO 150 J=1,IQ
            TEMP = TEMP+DBLE(GR(J+IP))*A(I-J)
  150    CONTINUE
  155    W = W+TEMP*TEMP
         A(I) = TEMP
  160 CONTINUE
      IF (L.EQ.3) GO TO 60
      IF (L.EQ.4) GO TO 110
C                   CALCULATE GRADIENT
  165 DO 170 I=1,IRS
         A(I+N) = ZERO
  170 CONTINUE
      IF (IQ.NE.0) GO TO 185
      DO 180 K=1,IP
         TEMP = DZERO
         DO 175 I=IRS,N
            TEMP = TEMP+DBLE(A(I))*X(I-K)
  175    CONTINUE
         GR(K+IPQ) = -TEMP-TEMP
  180 CONTINUE
      GO TO 205
  185 DO 200 K=1,IPQ
         KIP = IP-K
         TEMP = DZERO
         DO 195 I=IRS,N
            IF (K.LE.IP) T1 = -X(I-K)
            IF (K.GT.IP) T1 = A(I+KIP)
            IPN = I+N
            DO 190 J=1,IQ
               T1 = T1+DBLE(GR(J+IP))*A(IPN-J)
  190       CONTINUE
            A(IPN) = T1
            TEMP = TEMP+T1*A(I)
  195    CONTINUE
         GR(K+IPQ) = TEMP+TEMP
  200 CONTINUE
  205 IF (L.EQ.1) GO TO 85
      IF (L.EQ.2) GO TO 75
C                         UNDIFFERENCE ARPS
  210 PMAC = ZERO
      IF (IP.EQ.0) GO TO 220
      DO 215 I=1,IP
         PMAC = PMAC+ARPS(I)
  215 CONTINUE
  220 PMAC = XBAR*(ONE-PMAC)
      IND(5) = ICOUNT-1
      GR(1) = SINIT/N
      WNV = S1/N
      IF (ID.EQ.0) GO TO 245
      IF (IP.GT.0) GO TO 225
      ARPS(1) = ONE
      IF (ID.EQ.1) GO TO 245
      ID = ID-1
      IP = 1
  225 KIP = IP-1
      DO 240 K=1,ID
         L = K+KIP
         ARPS(L+1) = -ARPS(L)
         IF (L.EQ.1) GO TO 235
         LP2 = L+2
         DO 230 J=2,L
            I = LP2-J
            ARPS(I) = ARPS(I)-ARPS(I-1)
  230    CONTINUE
  235    ARPS(1) = ONE+ARPS(1)
  240 CONTINUE
  245 IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'FTML  ')
 9005 CONTINUE
      RETURN
      END

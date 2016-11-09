C   IMSL ROUTINE NAME   - FTCP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NON-SEASONAL ARIMA (BOX-JENKINS) STOCHASTIC
C                           MODEL ANALYSIS FOR A SINGLE TIME SERIES
C                           WITH FULL PARAMETER ITERATION AND MAXIMUM
C                           LIKELIHOOD ESTIMATION
C
C   USAGE               - CALL FTCP (X,IND,DSEED,ALPHA,ARPS,PMAS,PMAC,
C                           WNV,FCST,SIM,WK,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH IND(1) CONTAINING THE
C                           TIME SERIES.  X IS DESTROYED ON OUTPUT.
C                IND    - INPUT AND OUTPUT VECTOR OF LENGTH 10.
C                           IND IS DESTROYED ON OUTPUT EXCEPT AS
C                           INDICATED.
C                         IND(1) CONTAINS THE INPUT LENGTH OF TIME
C                           SERIES. IND(1) MUST BE GREATER THAN
C                           IND(6)+IND(7)+IND(8)+20.
C                         IND(2) CONTAINS, ON INPUT, THE MINIMUM
C                           NUMBER OF AUTOREGRESSIVE PARAMETERS
C                           IN THE DIFFERENCED FORM OF THE MODEL.
C                           IND(2) MUST BE GREATER THAN OR EQUAL TO
C                           ZERO. COMMON VALUES FOR IND(2) ARE 0, 1,
C                           OR 2.
C                           ON OUTPUT, IND(2) CONTAINS THE NUMBER OF
C                           AUTOREGRESSIVE PARAMETERS IN THE
C                           DIFFERENCED MODEL SELECTED FOR FITTING.
C                         IND(3) CONTAINS, ON INPUT, THE MINIMUM
C                           NUMBER OF MOVING AVERAGE PARAMETERS
C                           IN THE MODEL. IND(3) MUST BE GREATER THAN
C                           OR EQUAL TO ZERO. COMMON VALUES FOR IND(3)
C                           ARE 0, 1, OR 2.
C                           ON OUTPUT, IND(3) CONTAINS THE NUMBER OF
C                           MOVING AVERAGE PARAMETERS IN THE COMPUTED
C                           MODEL.
C                         NOTE THAT IND(2)+IND(3) MUST BE GREATER
C                           THAN ZERO.
C                         IND(4) CONTAINS, ON INPUT, THE MINIMUM
C                           NUMBER OF DIFFERENCING OPERATIONS ON THE
C                           TIME SERIES. IND(4) MUST BE GREATER THAN
C                           OR EQUAL TO ZERO. COMMON VALUES FOR IND(4)
C                           ARE 0, 1, OR 2.
C                           ON OUTPUT, IND(4) CONTAINS THE NUMBER OF
C                           DIFFERENCING OPERATIONS PERFORMED ON THE
C                           TIME SERIES.
C                         IND(5) CONTAINS THE INPUT MAXIMUM NUMBER OF
C                           ITERATIONS DESIRED TO CALCULATE MAXIMUM
C                           LIKELIHOOD ESTIMATES IN IMSL ROUTINE FTML.
C                           A COMMON VALUE FOR IND(5) IS 25.
C                         IND(6) CONTAINS THE INPUT MAXIMUM NUMBER OF
C                           AUTOREGRESSIVE PARAMETERS DESIRED IN THE
C                           DIFFERENCED FORM OF THE MODEL. IND(6) MUST
C                           BE GREATER THAN OR EQUAL TO IND(2). A COMMON
C                           VALUE FOR IND(6) IS 0, 1, OR 2.
C                         IND(7) CONTAINS THE INPUT MAXIMUM NUMBER OF
C                           MOVING AVERAGE PARAMETERS DESIRED IN THE
C                           MODEL. IND(7) MUST BE GREATER THAN OR
C                           EQUAL TO IND(3). A COMMON VALUE FOR IND(3)
C                           IS 0, 1, OR 2.
C                         IND(8) CONTAINS THE INPUT MAXIMUM NUMBER OF
C                           DIFFERENCING OPERATIONS TO BE PERFORMED ON
C                           THE TIME SERIES. IND(8) MUST BE GREATER THAN
C                           OR EQUAL TO IND(4). A COMMON VALUE FOR
C                           IND(8) IS 0, 1, OR 2.
C                         IND(9) CONTAINS THE INPUT POSITIVE FORECASTING
C                           PARAMETER. FORECASTS UP TO IND(9) STEPS IN
C                           ADVANCE ARE CALCULATED. IND(9) MUST BE
C                           GREATER THAN ZERO. A COMMON CHOICE FOR THE
C                           VALUE OF IND(9) IS THE LENGTH OF INTEREST
C                           IN THE FUTURE.
C                         IND(10) CONTAINS THE INPUT SIMULATION OPTION.
C                           IF IND(10) IS LESS THAN OR EQUAL TO ZERO,
C                             THEN SIMULATIONS ARE NOT DESIRED.
C                           IF IND(10) IS GREATER THAN ZERO,
C                             THEN IND(10) SIMULATIONS OF THE FUTURE UP
C                             TO IND(9) STEPS IN ADVANCE ARE DESIRED.
C                DSEED  - INPUT/OUTPUT. DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL. DSEED IS USED
C                           IN SIMULATING THE TIME SERIES, SO IT SHOULD
C                           BE CHOSEN RANDOMLY.
C                ALPHA  - INPUT AND OUTPUT VECTOR OF LENGTH 2.
C                         ALPHA(1) CONTAINS THE INPUT MINIMUM SIGNIFI-
C                           CANCE LEVEL, FOR MODEL SELECTION, IN THE
C                           EXCLUSIVE RANGE (0.0,1.0).
C                           A COMMON VALUE FOR ALPHA(1) IS .01 OR .05.
C                           LARGER VALUES OF ALPHA(1) CAUSE MODELS WITH
C                           MORE TERMS TO BE FITTED.
C                           ON OUTPUT, THE ESTIMATED SIGNIFICANCE LEVEL
C                           OY THE FITTED MODEL IS RETURNED.
C                           SEE THE ALGORITHM SECTION IN THE MANUAL
C                           DOCUMENT FOR FURTHER DETAILS.
C                         ALPHA(2) CONTAINS THE INPUT VALUE IN THE
C                           EXCLUSIVE RANGE (0.0,1.0) USED FOR COMPUTING
C                           100(1.-ALPHA(2)) PERCENT PROBABILITY
C                           LIMITS FOR THE FORECASTS. A COMMON VALUE
C                           FOR ALPHA(2) IS ANY CHOICE IN THE INTERVAL
C                           (.01,.5).
C                ARPS   - OUTPUT VECTOR OF LENGTH IND(6)+IND(8). THE
C                           FIRST IND(2)+IND(4) LOCATIONS CONTAIN THE
C                           AUTOREGRESSIVE PARAMETER ESTIMATES OF THE
C                           UNDIFFERENCED FORM OF THE MODEL.
C                PMAS   - OUTPUT VECTOR OF LENGTH 2*IND(7). THE FIRST
C                           IND(3) LOCATIONS CONTAIN THE MOVING AVERAGE
C                           PARAMETER ESTIMATES OF THE MODEL.
C                PMAC   - OUTPUT ESTIMATE OF THE OVERALL MOVING AVERAGE
C                           PARAMETER.
C                WNV    - OUTPUT ESTIMATE OF THE WHITE NOISE VARIANCE.
C                FCST   - OUTPUT MATRIX OF DIMENSION 3 BY IND(9).
C                         FOR LEAD TIMES J=1,2,...,IND(9),
C                         FCST(1,J) CONTAINS THE WEIGHTS FOR THE
C                           WEIGHTED SUM OF THE SHOCKS THAT GIVE THE
C                           FORECAST ERROR,
C                         FCST(2,J) CONTAINS THE FORECASTS,
C                         FCST(3,J) CONTAINS THE CORRESPONDING
C                           DEVIATIONS FROM EACH FORECAST FOR THE
C                           100*(1.-ALPHA(2)) PERCENT PROBABILITY
C                           LIMITS.
C                SIM    - OUTPUT VECTOR OF LENGTH IND(9)*IND(10) DEFINED
C                           ONLY FOR IND(10) GREATER THAN ZERO.
C                           SIM(J+(I-1)*IND(9)), FOR LEAD TIMES
C                           J=1,2,...,IND(9), CONTAINS THE RESULTS OF
C                           THE I-TH SIMULATION, FOR I=1,2,...,IND(10).
C                WK     - WORK AREA OF LENGTH THE MAXIMUM OF
C                           A. IND(1)+(IND(1)/10)+3*IND(6)+3*IND(7)
C                           B. 3*IND(6)+3*IND(7)+IND(8)+15+M, WHERE M
C                              IS THE MAXIMUM OF
C                              (1) IND(6)*(IND(6)+6)+IND(7)+1
C                              (2) 2*IND(1)
C                              (3) IND(6)+(IND(7)+1)*(3*IND(7)+24)/2
C                           C. 2*IND(6)+2*IND(7)+IND(8)+IND(9)+15
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES ONE OF THE INPUT
C                             PARAMETERS OF IND OR ALPHA WAS OUT OF
C                             RANGE.
C                           IER=130 INDICATES NO MODEL WAS TESTED THAT
C                             PASSED THE ALPHA(1) SIGNIFICANCE LEVEL.
C                         WARNING (WITH FIX)
C                           IER=67 INDICATES THAT ALL IND(5) ITERATIONS
C                             WERE PERFORMED IN FTML. CONVERGENCE IS
C                             ASSUMED AND THE PROGRAM CONTINUES
C                             CALCULATIONS.
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/FTARPS,FTAUTO,FTCAST,FTGEN,
C                           FTMA,FTMA1,FTML,GGNML,GGUBS,LEQT1F,
C                           LUDATN,LUELMN,MDCH,MDNOR,MDNRIS,MERFI,
C                           MERRC=ERFC,MGAMAD=DGAMMA,UERSET,UERTST,
C                           UGETIO,VABMXF,VBLA=SNRM2,ZSPOW,ZSPWA,
C                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG
C                       - SINGLE(H36,H48,H60)/FTARPS,FTAUTO,FTCAST,
C                           FTGEN,FTMA,FTMA1,FTML,GGNML,GGUBS,
C                           LEQT1F,LUDATN,LUELMN,MDCH,MDNOR,MDNRIS,
C                           MERFI,MERRC=ERFC,MGAMA=GAMMA,UERSET,UERTST,
C                           UGETIO,VABMXF,VBLA=SNRM2,ZSPOW,ZSPWA,
C                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG
C                       - DOUBLE/FTARPS,FTAUTO,FTCAST,FTGEN,FTMA,
C                           FTMA1,FTML,GGNML,GGUBS,LEQT1F,LUDATN,
C                           LUELMN,MDCH,MDNOR,MDNRIS,MERFI,
C                           MERRC=ERFC,MGAMAD=DGAMMA,UERSET,UERTST,
C                           UGETIO,VABMXF,VBLA=DNRM2,VXADD,VXMUL,VXSTO,
C                           ZSPOW,ZSPWA,ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,
C                           ZSPWG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE USER SHOULD EXPECT BETTER MODELS FOR LARGE
C                VALUES OF IND(1), ON THE ORDER OF 100.
C            2.  TO START FORECASTING AND SIMULATION PROCEDURES,
C                THE LAST IND(2)+IND(4) (AS ASSIGNED ON OUTPUT)
C                POINTS OF X ARE USED BY FTCP.
C            3.  IF IER=130 IS OBSERVED, ALPHA(1) CONTAINS THE
C                SIGNIFICANCE LEVEL OF THE LAST MODEL TESTED.
C                REDUCING OR INCREASING THE NUMBER OF PARAMETERS
C                OR DIFFERENCING AND A REDUCED ALPHA(1) ARE
C                POSSIBLE REMEDIES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTCP (X,IND,DSEED,ALPHA,ARPS,PMAS,PMAC,WNV,FCST,SIM,
     *                   WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,IND(10)
      REAL               X(1),ALPHA(2),ARPS(1),PMAS(1),PMAC,WNV,
     *                     FCST(3,1),SIM(1),WK(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I0,I3,ID,IDMIN,IDR,IEX,IP,IPD,IPQ,IQ,J,K,KP1,
     *                   L,LEVEL,LEVOLD,LS,LSIM,LSMIPD,LSP1,M,MPI,N,
     *                   NFOR,NMLSIM,NSIM
      REAL               CHI,T,XM
      REAL               ABAR,TNV,ZBAR
      DOUBLE PRECISION   TEMP,XBAR
      DATA               I0,I3/0,3/
C                                  FIRST EXECUTABLE STATEMENT
       IER=0
C                                  USEFUL CONSTANTS
      N = IND(1)
      IP = IND(2)
      IQ = IND(3)
      ID = IND(4)
      IDMIN = ID
      IPQ = IND(6)+IND(7)
      LSIM = IPQ+IND(8)+15
      LS = LSIM+1
      LSP1 = LS+1
      K = LSP1+IPQ
      KP1 = K+1
      NFOR = IND(9)
      NSIM = IND(10)
      CHI = 1.0
C                                  ERROR CHECKS
      IF (IP.GE.0 .AND. IQ.GE.0 .AND. ID.GE.0 .AND. IP+IQ.GT.0 .AND.
     *IND(6).GE.IP .AND. IND(7).GE.IQ .AND. IND(8).GE.ID .AND.
     *NFOR.GT.0 .AND. N.GE.LSIM+6 .AND. ALPHA(1).GT.0.0 .AND.
     *ALPHA(1).LT.1.0 .AND. ALPHA(2).GT.0.0 .AND. ALPHA(2).LT.1.0) GO
     *TO 5
      IER = 129
      GO TO 9000
C                                  STORE END OF SERIES
    5 LEVEL = 0
      CALL UERSET(LEVEL,LEVOLD)
      NMLSIM = N-LSIM
      DO 10 I=1,LSIM
         WK(I) = X(NMLSIM+I)
   10 CONTINUE
C                                  PREPROCESSING TIME SERIES
      IF (ID.EQ.0) GO TO 30
   15 DO 25 M=1,IDMIN
         N = N-1
         DO 20 I=1,N
            X(I) = X(I+1)-X(I)
   20    CONTINUE
   25 CONTINUE
   30 IDMIN = 1
C                                  PARAMETER ESTIMATES
      CALL FTAUTO(X,N,IPQ,I0,I3,ZBAR,WK(LS),WK(LSP1),WK,WK,WK)
      XBAR = ZBAR
      DO 35 I=1,N
         X(I) = X(I)-ZBAR
   35 CONTINUE
C                                  AUTOREGRESSIVE PARAMETERS
      IF (IP.EQ.0) GO TO 45
   40 CALL FTARPS(WK(LS),ZBAR,IP,IQ,ARPS,PMAC,WK(K),IER)
      IF (IER.GT.0) GO TO 90
C                                  MOVING AVERAGE PARAMETERS
      IF (IQ.EQ.0) GO TO 50
   45 CALL FTMA(WK(LS),ARPS,IP,IQ,PMAS,TNV,WK(K),IER)
      IF (IER.GT.0) GO TO 90
C                                  ONE-STEP RESIDUALS
   50 L = MAX0(IP,IQ)+1
      IDR = N/10+IP+IQ
      M = KP1+IDR
      DO 55 I=1,L
         WK(M+I) = 0.0
   55 CONTINUE
      DO 80 I=L,N
         MPI = M+I
         TEMP = X(I)
         IF (IP.EQ.0) GO TO 65
         DO 60 J=1,IP
            TEMP = TEMP-DBLE(ARPS(J))*X(I-J)
   60    CONTINUE
         IF (IQ.EQ.0) GO TO 75
   65    DO 70 J=1,IQ
            TEMP = TEMP+DBLE(PMAS(J))*WK(MPI-J)
   70    CONTINUE
   75    WK(MPI) = TEMP
   80 CONTINUE
C                                  RESIDUAL TEST OF BOX AND JENKINS
      CALL FTAUTO(WK(M+1),N,IDR,I0,I3,ABAR,WK(K),WK(KP1),WK,WK,WK)
      T = 0.0
      DO 85 I=1,IDR
         T = T+WK(K+I)**2
   85 CONTINUE
C                                  CHI-SQUARED STATISTIC
      T = (T*(N-ID))/WK(K)**2
      XM = IDR-IP-IQ
      CALL MDCH(T,XM,CHI,IEX)
      IF (1.0-CHI.GE.ALPHA(1)) GO TO 95
C                                  INCREMENT AND/OR RESET IQ,IP, OR ID
C                                  IF TEST NOT PASSED AT ALPHA(1) LEVEL
   90 IQ = IQ+1
      IF (IQ.LE.IND(7)) GO TO 45
      IQ = IND(3)
      IP = IP+1
      IF (IP.LE.IND(6)) GO TO 40
      IP = IND(2)
      ID = ID+1
      IF (ID.LE.IND(8)) GO TO 15
      IER = 130
      GO TO 115
C                                  MAXIMUM LIKELIHOOD PARAMETERS
   95 IND(1) = N
      IND(2) = IP
      IND(3) = IQ
      IND(4) = 0
      IND(6) = 5
      IND(7) = 1
      IND(8) = MAX0(4,IP+IQ**2)
      K = LS+IPQ+IPQ
      IF (K/2*2.EQ.K) K = K+1
      CALL FTML(X,IND,ARPS,PMAS,PMAC,WNV,WK(LS),WK(K),IER)
C                                  MOVING AVERAGE CONSTANT
      PMAC = 0.0
      IF (IP.EQ.0) GO TO 105
      DO 100 I=1,IP
         WK(I+LSIM) = ARPS(I)
         PMAC = PMAC+ARPS(I)
  100 CONTINUE
  105 PMAC = XBAR*(1.0-PMAC)
C                                  FORECASTING THE FUTURE
      IND(1) = LSIM
      IND(4) = ID
      IND(5) = NFOR
      CALL FTCAST(WK,WK(LSIM+1),PMAS,PMAC,ALPHA(2),IND,ARPS,FCST,TNV,
     *IEX)
      IF (NSIM.LE.0) GO TO 115
C                                  SIMULATING THE FUTURE
      IPD = IP+ID
      LSMIPD = LS-IPD
      DO 110 I=1,NSIM
         J = 1+(I-1)*NFOR
         CALL FTGEN(ARPS,PMAS,PMAC,WK(LSMIPD),WNV,DSEED,IPD,IQ,NFOR,
     *   SIM(J),WK(LS))
  110 CONTINUE
  115 ALPHA(1) = 1.0-CHI
      CALL UERSET(LEVOLD,LEVEL)
      IF (IER.EQ.0) GO TO 9005
      IF (IER.EQ.68) IER = 67
 9000 CONTINUE
      CALL UERTST(IER,'FTCP  ')
 9005 RETURN
      END

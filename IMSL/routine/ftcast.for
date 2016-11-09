C   IMSL ROUTINE NAME   - FTCAST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TIME SERIES FORECASTS AND PROBABILITY LIMITS
C                           USING AN ARIMA (BOX-JENKINS) MODEL
C
C   USAGE               - CALL FTCAST (Z,ARPS,PMAS,PMAC,ALPHA,LV,DARPS,
C                           FCST,WNV,IER)
C
C   ARGUMENTS    Z      - INPUT VECTOR OF LENGTH LV(1) CONTAINING THE
C                           TIME SERIES.
C                ARPS   - INPUT VECTOR OF LENGTH LV(2) CONTAINING
C                           ESTIMATES OF THE AUTOREGRESSIVE PARAMETERS.
C                PMAS   - INPUT VECTOR OF LENGTH 2*LV(3) CONTAINING
C                           ESTIMATES OF THE MOVING AVERAGE PARAMETERS
C                           IN THE FIRST LV(3) LOCATIONS.  THE REMAINING
C                           LOCATIONS ARE WORK STORAGE.
C                PMAC   - INPUT.  ESTIMATE OF OVERALL MOVING AVERAGE
C                           CONSTANT.
C                ALPHA  - INPUT.  A VALUE IN THE EXCLUSIVE INTERVAL
C                           (0,1) USED FOR COMPUTING 100(1-ALPHA) PER
C                           CENT PROBABILITY LIMITS FOR THE FORECASTS.
C                           THE VALUE 0.05 IS A COMMON CHOICE.
C                LV     - INPUT VECTOR OF LENGTH 5. LV(I) CONTAINS, WHEN
C                           I = 1, LENGTH OF TIME SERIES Z.
C                           I = 2, NUMBER OF AUTOREGRESSIVE PARAMETERS
C                             IN THE MODEL.
C                           I = 3, NUMBER OF MOVING AVERAGE PARAMETERS
C                             IN THE MODEL.
C                           I = 4, NUMBER OF DIFFERENCING OPERATIONS
C                             REQUIRED TO OBTAIN THE SERIES USED IN
C                             FITTING THE ARIMA MODEL.
C                           I = 5, MAXIMUM LEAD TIME DESIRED FOR A
C                             FORECAST.
C                DARPS  - OUTPUT VECTOR OF LENGTH LV(2)+LV(4) CONTAINING
C                           THE CONSTANTS, CORRESPONDING TO THOSE
C                           IN ARPS, FOR THE UNDIFFERENCED FORM OF
C                           THE MODEL.
C                FCST   - OUTPUT MATRIX OF DIMENSION 3 BY LV(5).
C                           FCST(I,J), FOR LEAD TIMES J=1,2,3,...,LV(5),
C                           CONTAINS WHEN
C                             I = 1, THE WEIGHTS FOR THE WEIGHTED SUM
C                               OF SHOCKS THAT GIVES THE FORECAST ERROR.
C                             I = 2, THE FORECASTS.
C                             I = 3, THE CORRESPONDING DEVIATIONS FROM
C                               EACH FORECAST FOR THE PROBABILITY
C                               LIMITS.
C                WNV    - OUTPUT.  ESTIMATE OF WHITE NOISE VARIANCE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES PARAMETER ALPHA WAS NOT IN
C                             THE EXCLUSIVE INTERVAL (0,1).
C                           IER=130 INDICATES ONE OR MORE OF LV(2),
C                             LV(3), OR LV(4) WERE LESS THAN ZERO OR
C                             LV(5) WAS LESS THAN ONE.
C                           IER=131 INDICATES LV(1) IS LESS THAN OR
C                             EQUAL TO LV(2)+LV(3)+LV(4).
C
C   REQD. IMSL ROUTINES - SINGLE/MDNRIS,MERFI,UERTST,UGETIO
C                       - DOUBLE/MDNRIS,MERFI,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE ESTIMATES CONTAINED IN THE INPUT PARAMETERS ARPS,
C                PMAC, AND PMAS, MAY BE COMPUTED BY USING IMSL ROUTINES
C                FTARPS AND FTMPS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTCAST (Z,ARPS,PMAS,PMAC,ALPHA,LV,DARPS,FCST,WNV,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LV(5),IER
      REAL               Z(1),ARPS(1),PMAS(1),PMAC,DARPS(1),FCST(3,1),
     1                   WNV,ALPHA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,LM,J,L,LV4,LV2,LL,K,KI,KK,IPD,IQP,IQ22,LV1,
     1                   LV24P,IQQ,IQ,IQ2,K1,MQL,LV3,KM1,LV2M,L1,
     2                   N1,LV3P,IJ,LV5,M
      REAL               TA,X
      REAL               ONEN,ZERO,ONE,ARP
      DOUBLE PRECISION   TEMP,TEMP1,TEMP2,S,DZERO
      DATA               DZERO/0.0D0/
      DATA               ONEN/-1.0/,ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CHECK FOR ERRORS IN LV.
      IF(LV(1) .LE. LV(2)+LV(3)+LV(4)) GO TO 205
      IF(LV(2).LT.0 .OR. LV(3).LT.0 .OR. LV(4).LT.0 .OR. LV(5).LT.1)
     1  GO TO 200
      LV2M = LV(2)-1
      IPD = LV(2)+LV(4)
      LV1 = LV(1)
      LV2 = LV(2)
      LV3 = LV(3)
      LV4 = LV(4)
      LV5 = LV(5)
      LV24P = IPD+1
C                                  COMPUTE THE CONSTANTS FOR THE
C                                  DIFFERENCE EQUATION FORM OF THE MODEL
C                                  CORRESPONDING TO THOSE IN ARPS
      IF(LV2 .NE. 0) GO TO 20
      IF(LV4 .EQ. 0) GO TO 30
      LM = LV4-1
      DARPS(1) = -ONEN
      IF(LM .EQ. 0) GO TO 30
      ONE = ONEN
      DO 15 I=1,LM
         L = I
         IF(I .EQ. 1) GO TO 10
         DO 5 J=2,I
            DARPS(L) = DARPS(L)-DARPS(L-1)
    5    L = L-1
   10    DARPS(I+1) = ONE
         ONE = -ONE
   15 DARPS(1) = DARPS(1)-ONEN
      DARPS(LV4) = ONEN
      IF(DARPS(LM) .LT. ZERO) DARPS(LV4) = -DARPS(LV4)
      GO TO 30
   20 DO 22 I=1,LV2
   22 DARPS(I) = ARPS(I)
      IF(LV4 .EQ. 0) GO TO 30
      LL = LV2+1
      ARP = ARPS(LV2)
      DO 27 I=1,LV4
         DARPS(LL) = -ARP
         L = LL-1
         M = LV2M+I
         IF(M .EQ. 1) GO TO 26
         DO 24 J=2,M
            DARPS(L) = DARPS(L)-DARPS(L-1)
   24    L = L-1
   26    DARPS(1) = DARPS(1)-ONEN
         ARP = -ARP
   27 LL = LL+1
      DARPS(LL-1) =  ARP
   30 K = MIN0(LV(5),LV(3))
      K1 = MAX0(LV(5),LV(3))
      MQL = K
      DO 35 I=1,LV5
   35 FCST(1,I) = ZERO
C                                  COMPUTE THE WEIGHTS FOR THE WEIGHTED
C                                  SUM OF SHOCKS THAT GIVES THE FORECAST
C                                  ERROR
      IF(K .EQ. 0) GO TO 45
      DO 40 I=1,K
   40 FCST(1,I) = -PMAS(I)
   45 K1 = LV24P
      KI = MIN0(IPD,LV5)
      IF(KI .EQ. 0) GO TO 52
      DO 50 I=1,KI
   50 FCST(1,I) = FCST(1,I)+DARPS(I)
   52 K = MIN0(K1,LV(5))
      K1 = LV5
      IF(K .LT. 2) GO TO 80
      DO 60 I=2,K
         L = I-1
         TEMP = DZERO
         DO 55 J=1,L
            L1 = I-J
   55    TEMP = TEMP+(DBLE(DARPS(J))*DBLE(FCST(1,L1)))
   60 FCST(1,I) = FCST(1,I)+TEMP
      IF(LV24P .GE. LV5) GO TO 80
      KK = LV24P+1
      DO 75 I=KK,LV5
         TEMP = DZERO
         DO 70 J=1,IPD
            L1 = I-J
            IF(L1 .NE. 0) GO TO 65
            TEMP = TEMP+DBLE(DARPS(J))
            GO TO 70
   65       TEMP = TEMP+(DBLE(DARPS(J))*DBLE(FCST(1,L1)))
   70    CONTINUE
   75 FCST(1,I) = FCST(1,I)+TEMP
   80 TEMP = DZERO
      K = K-1
      IQ = LV(3)-1
      IQQ = LV(3)
      IQP = IQQ+1
      IQ22 = IQQ+IQQ
      IF(IQQ .EQ. 0) GO TO 90
      DO 85 I=IQP,IQ22
   85 PMAS(I) = ZERO
C                                  COMPUTE THE ESTIMATE OF WHITE NOISE
C                                  VARIANCE
   90 IF(LV1 .LT. LV24P) GO TO 125
      DO 120 I=LV24P,LV1
         TEMP1 = DZERO
         TEMP2 = DZERO
         IF(IPD .EQ. 0) GO TO 100
         DO 95 J=1,IPD
   95    TEMP1 = TEMP1+(DBLE(DARPS(J))*DBLE(Z(I-J)))
  100    IF(IQQ .LT. 1) GO TO 115
         IQ2 = IQ22
         DO 105 J=1,IQQ
  105    TEMP2 = TEMP2+(DBLE(PMAS(J))*DBLE(PMAS(J+IQQ)))
         IF(IQ .LT. 1) GO TO 115
         DO 110 J=1,IQ
            IQ2 = IQ2-1
  110    PMAS(IQ2+1) = PMAS(IQ2)
  115    PMAS(IQP) = Z(I)-PMAC-TEMP1+TEMP2
         TEMP = TEMP+DBLE(PMAS(IQP))**2
  120 CONTINUE
  125 WNV = TEMP/LV1
      TA = 1.0-ALPHA/2.0
      CALL MDNRIS(TA,X,IER)
      IF(IER .GT. 127) GO TO 195
      S = X*SQRT(WNV)
C                                  COMPUTE THE CORRESPONDING DEVIATIONS
C                                  FROM EACH FORECAST FOR THE
C                                  PROBABILITY LIMITS
      FCST(3,1) = S
      IF(LV5 .EQ. 1) GO TO 140
      DO 135 I=2,LV5
         TEMP = DZERO
         L = I-1
         DO 130 J=1,L
  130    TEMP = TEMP+(DBLE(FCST(1,J))*DBLE(FCST(1,J)))
  135 FCST(3,I) = S*DSQRT(1.0D0+TEMP)
  140 K = KI
      IF(K .EQ. 0) GO TO 155
C                                  COMPUTE THE FORECASTS
      DO 150 I=1,K
         N1 = LV1+I
         TEMP = DZERO
         IF(I .GT. IPD) GO TO 150
         DO 145 J=I,IPD
  145    TEMP = TEMP+(DBLE(DARPS(J))*DBLE(Z(N1-J)))
  150 FCST(2,I) = PMAC+TEMP
  155 K = K+1
      IF(K .GT. K1) GO TO 165
      DO 160 I=K,K1
  160 FCST(2,I) = PMAC
  165 IF(MQL .EQ. 0) GO TO 180
      LV3P = LV3+1
      DO 175 I=1,MQL
         IF(I .GT. LV(3)) GO TO 180
         TEMP = DZERO
         DO 170 J=I,LV3
  170    TEMP = TEMP+DBLE(PMAS(J))*DBLE(PMAS(J-I+LV3P))
         FCST(2,I) = FCST(2,I)-TEMP
  175 CONTINUE
  180 IF(LV5 .EQ. 1) GO TO 9005
      DO 190 I=2,LV5
         KM1 = I-1
         TEMP = DZERO
         IF(IPD .EQ. 0) GO TO 9005
         IF(KM1 .GT. IPD) KM1 = IPD
         DO 185 J=1,KM1
            IJ = I-J
            TEMP = TEMP+(DBLE(DARPS(J))*DBLE(FCST(2,IJ)))
  185    CONTINUE
         FCST(2,I) = FCST(2,I)+TEMP
  190 CONTINUE
      GO TO 9005
  195 IER = 129
      GO TO 9000
  200 IER = 130
      GO TO 9000
  205 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,6HFTCAST)
 9005 CONTINUE
      RETURN
      END
 
R; T=0.08/0.84 00:00:15

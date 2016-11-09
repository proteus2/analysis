C   IMSL ROUTINE NAME   - FTFREQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - SINGLE OR MULTICHANNEL TIME SERIES ANALYSIS
C                           IN THE TIME AND FREQUENCY DOMAINS
C
C   USAGE               - CALL FTFREQ (X,IND,XIND,XYMV,ACV,FREQ,PS,
C                           XCOV,XSPECT,AMPHAS,XFER,COHER,IER)
C
C   ARGUMENTS    X      - INPUT/OUTPUT VECTOR OF LENGTH 2*IND(2). THE
C                           FIRST IND(Q) LOCATIONS CONTAIN A BASE TIME
C                           SERIES X, AND THE SECOND IND(2) LOCATIONS CO
C                           CONTAIN A TIME SERIES Y TO BE CROSSED WITH
C                           X.  X MAY BE MODIFIED ON OUTPUT.
C                IND    - INPUT VECTOR OF LENGTH 6. IND(I) CONTAINS
C                           WHEN
C                           I=1, A CONTROL PARAMETER.  IND(1) = 0
C                             IMPLIES CALCULATIONS FOR X SERIES ONLY.
C                             NONZERO IND(1) IMPLIES CALCULATIONS FOR
C                             Y SERIES AND CROSS-CALCULATIONS ONLY -
C                             ASSUMES X CALCULATIONS DONE AND
C                             LEFT IN X, XYMV, ACV, AND PS.
C                           I=2, THE LENGTH OF TIME SERIES X AND Y.
C                             IND(2) MUST BE GREATER THAN 3.
C                           I=3, CONTROL PARAMETER.  IF IND(3) IS
C                             NEGATIVE -- XYMV,ACV,XCOV
C                             ZERO     -- AND FREQ,PS,XSPECT
C                             POSITIVE -- AND AMPHAS,XFER,COHER
C                             ARE CALCULATED. (XSPECT, NOT IF IND(1)=0).
C                             FOR NONZERO IND(5) AND IND(1), IND(3)
C                             MUST BE NEGATIVE.
C                           I=4, THE NUMBER OF LAGS.
C                             IND(4) MUST BE GREATER THAN 2 AND
C                             LESS THAN IND(2).
C                             (CALL IND(4) M, AND LET MP1 = M+1)
C                           I=5, OPTION PARAMETER. NONZERO IND(5)
C                             IMPLIES PREWHITENING IS DESIRED.
C                             OTHERWISE, SET IND(5)=0.
C                           I=6, OPTION PARAMETER. NONZERO IND(6)
C                             IMPLIES DETRENDING IS DESIRED.
C                             OTHERWISE, SET IND(6)=0.
C                XIND   - INPUT VECTOR OF LENGTH 2. XIND(I) CONTAINS
C                           WHEN
C                           I=1, CONSTANT TIME INTERVAL (DECIMAL) IN
C                             SOME UNIT. XIND(1) MUST BE GREATER THAN
C                             0.0.
C                           I=2, CONSTANT IN EXCLUSIVE INTERVAL (-1,1)
C                             USED FOR PREWHITENING WHEN IND(5)
C                             IS NONZERO. XIND(2) MUST BE IN THE
C                             EXCLUSIVE INTERVAL (-1.0,1.0).
C                XYMV   - OUTPUT VECTOR OF LENGTH 6.
C                         XYMV(I) CONTAINS WHEN
C                           I=1, THE MEAN OF TIME SERIES X
C                           I=2, THE VARIANCE OF SERIES X
C                           I=3, THE MEAN OF TIME SERIES Y
C                           I=4, THE VARIANCE OF SERIES Y
C                           I=5 AND 6, WORK AREA
C                ACV    - OUTPUT VECTOR OF LENGTH 2*M. ACV(I) AND
C                           ACV(I+M) CONTAIN THE AUTOCOVARIANCES OF
C                           SERIES X AND Y, RESPECTIVELY, AT TIME
C                           LAG I, I=1,2,..,M.  XYMV(2) AND XYMV(4)
C                           CONTAIN THE AUTOCOVARIANCES AT LAG 0.
C                FREQ   - OUTPUT VECTOR OF LENGTH MP1 CONTAINING THE
C                           FREQUENCIES AT WHICH THE SPECTRAL QUANTITIES
C                           ARE CALCULATED IN CYCLES PER UNIT OF TIME.
C                PS     - OUTPUT VECTOR OF LENGTH 2*MP1. PS(I) AND
C                           PS(I+MP1) CONTAIN THE POWER SPECTRUM
C                           ESTIMATES OF X AND Y, RESPECTIVELY, AT
C                           FREQUENCIES FREQ(I), I=1,..,MP1.
C                XCOV   - OUTPUT VECTOR OF LENGTH M+MP1. XCOV(I)
C                           CONTAINS THE CROSS-COVARIANCE BETWEEN SERIES
C                           X AND Y AT TIME LAG I-MP1, I=1,2,..,M+MP1.
C                           THUS, (I-MP1) = -M,..,0,..,M .
C                XSPECT - OUTPUT VECTOR OF LENGTH 2*MP1, THE CROSS-
C                           SPECTRUM. XSPECT(I) AND XSPECT(I+MP1)
C                           CONTAIN THE COSPECTRUM (REAL PART) AND
C                           QUADRATURE SPECTRUM (IMAGINARY PART),
C                           RESPECTIVELY, AT FREQ(I), I=1,2,..,MP1 .
C                AMPHAS - OUTPUT VECTOR OF LENGTH 2*MP1. AMPHAS(I) AND
C                           AMPHAS(I+MP1) CONTAIN THE AMPLITUDE AND
C                           PHASE, RESPECTIVELY, OF THE CROSS-SPECTRUM
C                           AT FREQ(I),I=1,2,..,MP1. THE PHASE IS GIVEN
C                           IN FRACTIONS OF A CIRCLE, BETWEEN 0 AND 1.
C                XFER   - OUTPUT VECTOR OF LENGTH 2*MP1. XFER(I) AND
C                           XFER(I+MP1) CONTAIN THE AMPLITUDE OF THE
C                           TRANSFER FUNCTIONS FROM X TO Y AND Y TO X,
C                           RESPECTIVELY, AT FREQ(I), I=1,2,...,MP1.
C                           XFER(I) IS SET TO MACHINE INFINITY IF THE
C                           CORRESPONDING POWER SPECTRAL ESTIMATE IS NOT
C                           POSITIVE. (SEE IER=67)
C                COHER  - OUTPUT VECTOR OF LENGTH MP1. COHER(I) CONTAINS
C                           THE COHERENCE SQUARE AT FREQ(I), I=1,2,...,
C                           MP1.  COHER(I) IS SET TO ONE IF IER=67 OR IF
C                           THE CALCULATED ESTIMATE IS GREATER THAN 1.0
C                           (SEE IER=68).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES ONE OF IND OUT OF RANGE.
C                           IER=130 INDICATES ONE OF XIND OUT OF RANGE.
C                         WARNING (WITH FIX)
C                           IER=67 INDICATES A POWER SPECTRAL ESTIMATE
C                             WAS NOT POSITIVE. IT IS SET TO ZERO.
C                             THE CORRESPONDING XFER AND COHER ARE SET
C                             TO MACHINE INFINITY AND ONE, RESPECTIVELY.
C                           IER=68 INDICATES A COHERENCE ESTIMATE WAS
C                             GREATER THAN 1.0. IT IS SET TO 1.0.
C
C   REQD. IMSL ROUTINES - SINGLE/FTAUTO,UERTST,UGETIO
C                       - DOUBLE/FTAUTO,UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      SUPPOSE THE USER HAS A BASE SERIES X AND ONE OR MORE
C                CROSS SERIES, Y1 AND Y2 (TWO FOR EXAMPLE). THE USER
C                SHOULD FIRST SET ALL INPUT PARAMETERS AS FOR NONZERO
C                IND(1). THEN THE FIRST CALL TO FTFREQ SHOULD BE WITH
C                IND(1)=0 TO CALCULATE XYMV, ACV, AND PS (IF DESIRED)
C                FOR X ALONE. NEXT, IND(1) IS SET TO A NONZERO VALUE AND
C                SERIES Y1 IS PLACED IN X(I+IND(2)), I=1,2,...,IND(2)
C                AND FTFREQ IS CALLED AGAIN TO CALCULATE XYMV, ACV, AND
C                PS FOR Y1 AND ALSO THE CROSS-QUANTITIES. NOW REPLACE
C                Y1 BY Y2 IN THE CALLING SEQUENCE TO REPEAT THE ABOVE
C                CALCULATIONS. THIS PROCEDURE MAY BE REPEATED FOR AS
C                MANY CROSS-SERIES AS DESIRED. NOTE THAT OUTPUT
C                PARAMETERS FOR SERIES X SHOULD BE USED AS INPUT
C                PARAMETERS WHEN IND(1) IS NONZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTFREQ (X,IND,XIND,XYMV,ACV,FREQ,PS,XCOV,XSPECT,
     1                   AMPHAS,XFER,COHER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IND(6),IER
      REAL               X(1),XIND(2),XYMV(6),ACV(1),FREQ(1),PS(1),
     1                   XCOV(1),XSPECT(1),AMPHAS(1),XFER(1),COHER(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I0,I1,I2,IERX,IT1,ITEMP,ITM,ITMP,ITMPP1,ITN,
     1                   ITNP1,J,K,L,M,ML,MLM1,MLP1,N,NM1
      REAL               PI,XN,XML,WK,Z(1),XBAR,DEN,D,XI,XJ,E,PISXML,
     1                   PIPPI,SCAL1,SCAL2,SINF
      DOUBLE PRECISION   X54,X46,X23,SAVE,SAVEX,HEMP,HEMPX,AMBDA,
     1                   C,F,G,C1,C2,C3,TEMP,TEMP1,TEMP2,X1S8,XP5
      DATA               X1S8/.8333333333333333D-01/
      DATA               XP5/.5D0/,X54/.54D0/,X46/.46D0/,X23/.23D0/
      DATA               PI/3.141593/
      DATA               SINF/Z7FFFFFFF/
      DATA               I0/0/,I1/1/,I2/2/,IERX/0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  ERROR CHECKS
      IF (IND(2).GE.4.AND.IND(4).LT.IND(2).AND.IND(4).GE.3) GO TO 10
    5 IER = 129
      GO TO 9000
   10 IF (XIND(1).GT.0.0) GO TO 20
   15 IER = 130
      GO TO 9000
   20 IF (IND(5).EQ.0) GO TO 25
      IF (XIND(2).LE.-1.0.OR.XIND(2).GE.1.0) GO TO 15
      IF (IND(1).NE.0.AND.IND(3).GE.0) GO TO 5
C                                  USEFUL CONSTANTS
   25 N = IND(2)
      XN = N
      ML = IND(4)
      XML = ML
      MLP1 = ML+1
      IF (IND(1).EQ.0) IERX = 0
      MLM1 = ML-1
      IF (IND(1).NE.0) GO TO 30
      ITM = 0
      ITN = 0
      ITMP = 0
      IT1 = 0
      GO TO 35
   30 ITM = IND(4)
      ITN = IND(2)
      ITMP = ITM+1
      IT1 = 2
C                                  PREWHITENING
   35 IF (IND(5).EQ.0) GO TO 45
      NM1 = N-1
      ITNP1 = ITN+1
      SAVE = XIND(2)
      WK = SAVE*DBLE(X(N+ITN))-DBLE(X(1+ITN))
      DO 40 I=1,NM1
         X(I+ITN) = DBLE(X(I+ITNP1))-SAVE*DBLE(X(I+ITN))
   40 CONTINUE
      X(N+ITN) = WK
C                                  CALCULATE MEAN AND VARIANCE
   45 CALL FTAUTO (X(1+ITN),N,I0,I0,I1,SCAL1,SCAL2,Z,Z,Z,Z)
      XYMV(1+IT1) = SCAL1
      XYMV(2+IT1) = SCAL2
      IF (IND(6).NE.0) XYMV(2+IT1) = XYMV(2+IT1)+XYMV(1+IT1)*XYMV(1+IT1)
      IF (IND(6).NE.0) GO TO 55
      XBAR = XYMV(IT1+1)
      DO 50 I=1,N
         X(I+ITN) = X(I+ITN)-XBAR
   50 CONTINUE
C                                  COVARIANCES
   55 SCAL1 = 0.0
      SCAL2 = 0.0
      CALL FTAUTO (X(1+ITN),N,ML,I0,I2,SCAL1,SCAL2,ACV(1+ITM),Z,Z,Z)
C                                  SCALE COV
      SAVE = N
      DO 60 I=1,ML
         TEMP = SAVE/(SAVE-I)
         ACV(I+ITM) = DBLE(ACV(I+ITM))*TEMP
   60 CONTINUE
C                                  DETRENDING
      IF (IND(6).EQ.0) GO TO 75
      K = (N+2)/3
      ITNP1 = ITN+N-K
      TEMP = 0.0D0
      DO 65 I=1,K
         TEMP = TEMP+DBLE(X(I+ITNP1))-DBLE(X(I+ITN))
   65 CONTINUE
      ITEMP = 5 + IT1/2
      XYMV(ITEMP) = TEMP/(DBLE(XIND(1))*(K*(N-K)))
      XYMV(2+IT1) = DBLE(XYMV(2+IT1))-DBLE(XYMV(1+IT1))**2-(DBLE(XN)*
     1 DBLE(XIND(1)))**2*DBLE(XYMV(ITEMP))**2*X1S8
      HEMP = N
      TEMP = DBLE(XIND(1))**2
      TEMP1 = DBLE(XYMV(1+IT1))**2
      TEMP2 = DBLE(XYMV(ITEMP))**2*X1S8
      DO 70 I=1,ML
         SAVE = I
         AMBDA = TEMP*(HEMP*(HEMP-SAVE-SAVE)-2.D0*SAVE**2)
         ACV(I+ITM) = DBLE(ACV(I+ITM))-TEMP1-AMBDA*TEMP2
   70 CONTINUE
   75 IF (IND(3).LT.0) GO TO 120
      IF (IND(1).NE.0) GO TO 85
      DEN = 1.0/((XML+XML)*XIND(1))
      DO 80 I=1,MLP1
         FREQ(I) = (I-1)*DEN
   80 CONTINUE
C                                  SPECTRUM  - DOUBLE
   85 C = (XIND(1)+XIND(1))/PI
      D = PI/XML
      TEMP1 = DBLE(XYMV(2+IT1))*XP5
      TEMP2 = DBLE(ACV(ML+ITM))*XP5
      DO 95 I=1,MLP1
         XI = I-1
         XJ = D*XI
         TEMP = TEMP1
         DO 90 J=1,MLM1
            E = XJ*J
            TEMP = TEMP+DBLE(ACV(J+ITM))*DBLE(COS(E))
   90    CONTINUE
         TEMP = TEMP+TEMP2*DBLE(COS(XI*PI))
         PS(I+ITMP) = TEMP*C
   95 CONTINUE
C                                  SMOOTH SPECTRUM
      SAVE = PS(1+ITMP)
      PS(1+ITMP) = X54*SAVE+X46*DBLE(PS(2+ITMP))
      ITMPP1 = ITMP+1
      DO 100 I=2,ML
         HEMP = PS(I+ITMP)
         PS(I+ITMP) = X23*(SAVE+DBLE(PS(I+ITMPP1))) + X54*HEMP
         SAVE = HEMP
  100 CONTINUE
      K = MLP1+ITMP
      PS(K) = X54*DBLE(PS(K))+X46*SAVE
C                                  RECOLORING
      IF (IND(5).EQ.0) GO TO 110
      SAVE = XIND(2)
      PISXML = PI/XML
      DO 105 I=1,MLP1
         XI = I-1
         C = 1.D0+SAVE*(SAVE-2.D0*DBLE(COS(XI*PISXML)))
         PS(I+ITMP) = DBLE(PS(I+ITMP))/C
  105 CONTINUE
C                                  CHECK SIGN OF SPECTRAL ESTIMATE
  110 DO 115 I=1,MLP1
         IF (PS(I+ITMP).GT.0.0) GO TO 115
         PS(I+ITMP) = 0.0
         IF (IND(1).EQ.0) IERX = 1
         IER = 67
  115 CONTINUE
  120 IF (IND(1).EQ.0) GO TO 200
C                                  MULTICHANNEL CALCULATIONS
C                                  CROSS-COVARIANCE
C                                  NEGATIVE TAU
      DO 130 I=1,MLP1
         K = MLP1-I
         L = N-K
         TEMP = 0.0D0
         DO 125 J=1,L
            TEMP = TEMP+DBLE(X(J+K))*DBLE(X(J+N))
  125    CONTINUE
         XCOV(I) = TEMP/(N-K)
  130 CONTINUE
C                                  POSITIVE TAU
      DO 140 I=1,ML
         K = N+I
         L = N-I
         TEMP = 0.0D0
         DO 135 J=1,L
            TEMP = TEMP+DBLE(X(J))*DBLE(X(J+K))
  135    CONTINUE
         XCOV(I+MLP1) = TEMP/L
  140 CONTINUE
C                                  DETREND
      IF (IND(6).EQ.0) GO TO 150
      SAVE = XYMV(1)
      SAVEX = XYMV(2)
      HEMP = XYMV(3)
      HEMPX = XYMV(4)
      F = XYMV(5)
      G = XYMV(6)
      C1 = SAVE*HEMP
      C2 = (SAVE*G-SAVEX*F)*XP5
      C3 = F*G*X1S8
      HEMP = N
      TEMP2 = DBLE(XIND(1))**2
      DO 145 I=1,ML
         SAVE = I
         AMBDA = TEMP2*(HEMP*(HEMP-SAVE-SAVE)-2.D0*SAVE**2)
         K = MLP1+I
         M = MLP1-I
         TEMP = -C1 - AMBDA*C3
         TEMP1 = SAVE * C2
         XCOV(K) = DBLE(XCOV(K)) + TEMP - TEMP1
         XCOV(M) = DBLE(XCOV(M)) + TEMP + TEMP1
  145 CONTINUE
      XCOV(MLP1) = DBLE(XCOV(MLP1))-C1-C3*(DBLE(XIND(1))*HEMP)**2
C                                  CROSS-SPECTRUM
  150 IF (IND(3).LT.0) GO TO 200
      C = XIND(1)/PI
      D = PI/XML
      DO 160 I=1,MLP1
         XI = I-1
         XJ = D*XI
         TEMP1 = XCOV(MLP1)
         TEMP2 = 0.0D0
         DO 155 J=1,MLM1
            K = MLP1+J
            M = MLP1-J
            E = XJ*J
            F = DBLE(XCOV(K))
            G = DBLE(XCOV(M))
            TEMP1 = TEMP1+(F+G)*DBLE(COS(E))
            TEMP2 = TEMP2+(F-G)*DBLE(SIN(E))
  155    CONTINUE
         E = XI*PI
         F = DBLE(XCOV(ML+MLP1))
         G = DBLE(XCOV(1))
         XSPECT(I) = C*(TEMP1+(F+G)*DBLE(COS(E))*XP5)
         XSPECT(I+MLP1) = C*(TEMP2+(F-G)*DBLE(SIN(E))*XP5)
  160 CONTINUE
C                                  SMOOTH X-SPECTRUM
      SAVE = XSPECT(1)
      SAVEX = XSPECT(1+ITMP)
      XSPECT(1) = X54*SAVE+X46*DBLE(XSPECT(2))
      XSPECT(1+ITMP) = X54*SAVEX+X46*DBLE(XSPECT(2+ITMP))
      DO 165 I=2,ML
         K = I+ITMP
         HEMP = XSPECT(I)
         HEMPX = XSPECT(K)
         XSPECT(I) = X23*SAVE+X54*HEMP+X23*DBLE(XSPECT(I+1))
         XSPECT(K) = X23*SAVEX+X54*HEMPX+X23*DBLE(XSPECT(K+1))
         SAVE = HEMP
         SAVEX = HEMPX
  165 CONTINUE
      XSPECT(MLP1) = X54*DBLE(XSPECT(MLP1))+X46*SAVE
      XSPECT(MLP1+ITMP) = X54*DBLE(XSPECT(MLP1+ITMP))+X46*SAVEX
      IF (IND(3).EQ.0) GO TO 200
C                                  AMPLITUDE AND PHASE
      PIPPI = PI+PI
      DO 175 I=1,MLP1
         K = I+ITMP
         AMPHAS(I) = DSQRT(DBLE(XSPECT(I))**2+DBLE(XSPECT(K))**2)
         IF (XSPECT(I).NE.0.0) GO TO 170
         IF (XSPECT(K).GE.0.0) AMPHAS(K) = 0.25
         IF (XSPECT(K).LT.0.0) AMPHAS(K) = 0.75
         GO TO 175
  170    AMPHAS(K) = ATAN(XSPECT(K)/XSPECT(I))
         IF (XSPECT(I).LT.0.0) AMPHAS(K) = AMPHAS(K)+PI
         AMPHAS(K) = AMPHAS(K)/PIPPI
         IF (AMPHAS(K).LT.0.0) AMPHAS(K) = 1.0+AMPHAS(K)
  175 CONTINUE
C                                  TRANSFER FUNCTION
      IF (IERX.EQ.1) IER = 67
      DO 195 I=1,MLP1
         K = I+ITMP
         IF(PS(I).EQ.0.0) GO TO 180
         XFER(I) = AMPHAS(I)/PS(I)
         IF(PS(K).EQ.0.0) GO TO 185
         XFER(K) = AMPHAS(I)/PS(K)
         COHER(I) = AMPHAS(I)*AMPHAS(I)/(PS(I)*PS(K))
         GO TO 190
  180    XFER(I) = SINF
  185    XFER(K) = SINF
         IF(PS(K).NE.0.0) XFER(K) = AMPHAS(I)/PS(K)
         COHER(I) = SINF
  190    IF(COHER(I).LE.1.0) GO TO 195
         IF(IER.NE.67) IER = 68
         COHER(I) = 1.0
  195 CONTINUE
  200 IF(IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'FTFREQ')
 9005 RETURN
      END

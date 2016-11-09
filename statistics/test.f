      program test

      real ffff, gggg
      data ffff/3.402823E38/
      data gggg/2.3841857910156E-7/

      print*, 3.402823E38
      print*, 2.3841857910156E-7
      print*,ffff
      print*,gggg

      end


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
C                           FIRST IND(2) LOCATIONS CONTAIN A BASE TIME
C                           SERIES X, AND THE SECOND IND(2) LOCATIONS
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
      DATA               SINF/3.402823E38/
C      DATA               SINF/268435455./
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


C   IMSL ROUTINE NAME   - FTAUTO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - MEAN, VARIANCE, AUTOCOVARIANCES,
C                           AUTOCORRELATIONS, AND PARTIAL
C                           AUTOCORRELATIONS FOR A STATIONARY
C                           TIME SERIES.
C
C   USAGE               - CALL FTAUTO (W,LW,K,L,ISW,AMEAN,VAR,ACV,AC,
C                           PACV,WKAREA)
C
C   ARGUMENTS    W      - INPUT VECTOR OF LENGTH LW CONTAINING THE TIME
C                           SERIES.
C                LW     - INPUT LENGTH OF W.
C                K      - INPUT NUMBER OF AUTOCOVARIANCES AND
C                           AUTOCORRELATIONS TO BE COMPUTED.
C                L      - INPUT NUMBER OF PARTIAL AUTOCORRELATIONS TO
C                           BE COMPUTED.  L MUST BE LESS THAN OR EQUAL
C                           TO K.
C                ISW    - INPUT CONTROL PARAMETER FOR DETERMINING TASK
C                           TO BE PERFORMED. IF
C                             ISW = 1  FIND MEAN AND VARIANCE.
C                             ISW = 2  FIND AUTOCOVARIANCE.
C                             ISW = 3  FIND MEAN, VARIANCE, AND
C                                        AUTOCOVARIANCES.
C                             ISW = 4  FIND AUTOCOVARIANCES AND
C                                        AUTOCORRELATIONS.
C                             ISW = 5  FIND MEAN, VARIANCE, AUTO-
C                                        COVARIANCES, AND AUTOCORRELAT-
C                                        IONS.
C                             ISW = 6  FIND AUTOCOVARIANCES, AUTOCORREL-
C                                        ATIONS, AND PARTIAL AUTOCO-
C                                        RRELATIONS.
C                             ISW = 7  FIND MEAN, VARIANCE, AUTOCOVAR-
C                                        IANCES, AUTOCORRELATIONS, AND
C                                        PARTIAL AUTOCORRELATIONS.
C                AMEAN  - OUTPUT FOR ISW = 1,3,5, AND 7.
C                           INPUT FOR ISW = 2,4, AND 6. MEAN VALUE OF
C                           THE TIME SERIES W.
C                VAR    - OUTPUT FOR ISW = 1,3,5, AND 7.
C                           INPUT FOR ISW = 2,4, AND 6. VARIANCE OF
C                           TIME SERIES W.
C                ACV    - VECTOR OF LENGTH K.
C                           OUTPUT FOR ISW = 2,3,4,5,6 AND 7.
C                           AUTOCOVARIANCES FOR TIME SERIES W.  ACV(I)
C                           CORRESPONDS TO A TIME LAG OF I TIME UNITS.
C                AC     - VECTOR OF LENGTH K.
C                           OUTPUT FOR ISW = 4,5,6, AND 7.
C                           AUTOCORRELATIONS FOR TIME SERIES W.  AC(I)
C                           CORRESPONDS TO A TIME LAG OF I TIME UNITS.
C                PACV   - VECTOR OF LENGTH L.
C                           OUTPUT FOR ISW = 6 AND 7.
C                           PARTIAL AUTOCORRELATIONS OF TIME SERIES W.
C                           PACV(1) = AC(1).
C                WKAREA - WORK AREA VECTOR OF LENGTH L.
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IN SOME ROUTINES OF THIS CHAPTER ACV(I) CORRESPONDS
C                TO A TIME LAG OF (I-1) TIME UNITS RATHER THAN I TIME
C                UNITS. THUS, IN THESE SUBROUTINES, ACV(1) IS THE SAME
C                AS THE VARIANCE VAR. IN THE CALLING PROGRAM TO FTAUTO,
C                IF THE USER WISHES THE VARIANCE TO BE THE FIRST ENTRY
C                IN HIS AUTOCOVARIANCE ARRAY THE FOLLOWING CALL CAN BE
C                MADE
C                  CALL FTAUTO(W,LW,K,L,ISW,AMEAN,ACV(1),ACV(2),AC,PACV,
C                 *            WKAREA)
C                THE USER SHOULD NOTE THAT IN THIS CASE, ACV MUST BE
C                DIMENSIONED K+1 IN THE MAIN PROGRAM.
C            2.  IF THE TIME SERIES W IS CONSTANT, THEN ANY OF ACV, AC,
C                AND PACV THAT ARE OUTPUT, ACCORDING TO THE ISW SETTING,
C                ARE SET TO ZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTAUTO  (W,LW,K,L,ISW,AMEAN,VAR,ACV,AC,PACV,WKAREA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LW,K,L,ISW
      REAL               W(LW),ACV(1),AC(1),PACV(1),WKAREA(1),AMEAN,VAR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,KEND,J1,K0,J2,J1MK,IFLAG,IM
      REAL               TEMP2,ZERO
      DOUBLE PRECISION   TEMP,TEMP1,ONE,DZERO
      DATA               ZERO/0.0/,ONE/1.0D0/,DZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (ISW.EQ.1) GO TO 40
      IFLAG = 0
      IF (K.LE.0) GO TO 15
      DO 5 I=1,K
    5 ACV(I) = ZERO
      IF (ISW.LT.4) GO TO 25
      DO 10 I=1,K
   10 AC(I) = ZERO
   15 IF (L.LE.0) GO TO 25
      DO 20 I=1,L
   20 PACV(I) = ZERO
   25 IF ((ISW/2)*2.NE.ISW) VAR = ZERO
      DO 30 I=2,LW
         IF (W(I).NE.W(I-1)) GO TO 35
   30 CONTINUE
      IFLAG = 1
      AMEAN = W(1)
   35 IF (IFLAG.EQ.1) GO TO 9005
   40 IM = (ISW/2)*2-ISW
C                                  COMPUTE THE MEAN
      IF (IM.EQ.0) GO TO 55
      TEMP = DZERO
      DO 45 I=1,LW
         TEMP = TEMP+DBLE(W(I))
   45 CONTINUE
      AMEAN = TEMP/LW
C                                  COMPUTE THE VARIANCE
      TEMP = DZERO
      DO 50 I=1,LW
         TEMP = TEMP+(DBLE(W(I)-AMEAN)*DBLE(W(I)-AMEAN))
   50 CONTINUE
      VAR = TEMP/LW
      IF (ISW.EQ.1) GO TO 9005
C                                  COMPUTE AUTOCOVARIANCES
   55 DO 65 J=1,K
         KEND = LW-J
         IF (KEND .LE. 0) GO TO 65
         TEMP = DZERO
         DO 60 I=1,KEND
            TEMP = TEMP+(DBLE(W(I)-AMEAN)*DBLE(W(I+J)-AMEAN))
   60    CONTINUE
         ACV(J) = TEMP/LW
   65 CONTINUE
      IF (ISW.LT.4) GO TO 9005
C                                  COMPUTE AUTOCORRELATIONS
      DO 70 J=1,K
         AC(J) = ACV(J)/VAR
   70 CONTINUE
      IF (ISW.LT.6) GO TO 9005
C                                  COMPUTE PARTIAL AUTOCOVARIANCE
      PACV(1) = AC(1)
      DO 90 J=2,L
         J1 = J-1
         WKAREA(J1) = PACV(J1)
         J2 = (J1)/2
         IF (J.EQ.2) GO TO 80
         DO 75 K0=1,J2
            J1MK = J1-K0
            TEMP2 = WKAREA(K0)-PACV(J1)*WKAREA(J1MK)
            WKAREA(J1MK) = WKAREA(J1MK)-PACV(J1)*WKAREA(K0)
            WKAREA(K0) = TEMP2
   75    CONTINUE
   80    CONTINUE
         TEMP = DZERO
         TEMP1 = DZERO
         DO 85 I=1,J1
            TEMP = TEMP+(DBLE(AC(J-I))*DBLE(WKAREA(I)))
            TEMP1 = TEMP1+(DBLE(AC(I))*DBLE(WKAREA(I)))
   85    CONTINUE
         PACV(J) = (AC(J)-TEMP)/(ONE-TEMP1)
   90 CONTINUE
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - UERTST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
C
C   USAGE               - CALL UERTST (IER,NAME)
C
C   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
C                           IER = I+J WHERE
C                             I = 128 IMPLIES TERMINAL ERROR MESSAGE,
C                             I =  64 IMPLIES WARNING WITH FIX MESSAGE,
C                             I =  32 IMPLIES WARNING MESSAGE.
C                             J = ERROR CODE RELEVANT TO CALLING
C                                 ROUTINE.
C                NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING
C                           THE NAME OF THE CALLING ROUTINE. (INPUT)
C
C   REQD. IMSL ROUTINES - UGETIO,USPKD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
C                TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
C                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
C                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
C                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
C                UGETIO AS FOLLOWS..
C                                NIN = 0
C                                NOUT = NEW OUTPUT UNIT NUMBER
C                                CALL UGETIO(3,NIN,NOUT)
C                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UERTST (IER,NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      INTEGER            NAME(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEQ,IEQDF,IOUNIT,LEVEL,LEVOLD,NAMEQ(6),
     *                   NAMSET(6),NAMUPK(6),NIN,NMTB
      DATA               NAMSET/1HU,1HE,1HR,1HS,1HE,1HT/
      DATA               NAMEQ/6*1H /
      DATA               LEVEL/4/,IEQDF/0/,IEQ/1H=/
C                                  UNPACK NAME INTO NAMUPK
C                                  FIRST EXECUTABLE STATEMENT
      CALL USPKD (NAME,6,NAMUPK,NMTB)
C                                  GET OUTPUT UNIT NUMBER
      CALL UGETIO(1,NIN,IOUNIT)
C                                  CHECK IER
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
C                                  PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
C                                  PRINT WARNING WITH FIX MESSAGE
C      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
C      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
C                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
      GO TO 30
   15 CONTINUE
C                                  CHECK FOR UERSET CALL
      DO 20 I=1,6
         IF (NAMUPK(I).NE.NAMSET(I)) GO TO 25
   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
C                                  PRINT NON-DEFINED MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK
   30 IEQDF = 0
      RETURN
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   40 FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
C
C                                  SAVE P FOR P = R CASE
C                                    P IS THE PAGE NAMUPK
C                                    R IS THE ROUTINE NAMUPK
   55 IEQDF = 1
      DO 60 I=1,6
   60 NAMEQ(I) = NAMUPK(I)
   65 RETURN
      END


C   IMSL ROUTINE NAME   - UGETIO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
C                           VALUES FOR INPUT AND OUTPUT UNIT
C                           IDENTIFIERS.
C
C   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
C
C   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
C                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
C                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
C                           AND NOUT, RESPECTIVELY.
C                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
C                           RESET FOR SUBSEQUENT USE.
C                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
C                           RESET FOR SUBSEQUENT USE.
C                NIN    - INPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
C                NOUT   - OUTPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
C                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
C                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
C                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
C                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/5/,NOUTD/6/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - VXADD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EXTENDED PRECISION ADD
C
C   USAGE               - CALL VXADD (A,ACC)
C
C   ARGUMENTS    A      - DOUBLE PRECISION NUMBER TO BE ADDED TO THE
C                           ACCUMULATOR. (INPUT)
C                ACC    - ACCUMULATOR. (INPUT AND OUTPUT)
C                           ACC IS A DOUBLE PRECISION VECTOR OF LENGTH
C                           2. ON OUTPUT, ACC CONTAINS THE SUM OF
C                           INPUT ACC AND A.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      VXADD ADDS THE DOUBLE PRECISION NUMBER A TO THE
C                EXTENDED PRECISION ACCUMULATOR, ACC. THE SUBROUTINE
C                ASSUMES THAT AN EXTENDED PRECISION NUMBER IS ALREADY IN
C                THE ACCUMULATOR. THEREFORE, BEFORE THE FIRST CALL TO
C                VXADD, ACC(1) AND ACC(2) MUST BE SET TO ZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VXADD(A,ACC)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   A,ACC(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   X,Y,Z,ZZ
C                                  FIRST EXECUTABLE STATEMENT
      X = ACC(1)
      Y = A
      IF (DABS(ACC(1)).GE.DABS(A)) GO TO 1
      X = A
      Y = ACC(1)
C                                  COMPUTE Z+ZZ = ACC(1)+A EXACTLY
    1 Z = X+Y
      ZZ = (X-Z)+Y
C                                  COMPUTE ZZ+ACC(2) USING DOUBLE
C                                    PRECISION ARITHMETIC
      ZZ = ZZ+ACC(2)
C                                  COMPUTE ACC(1)+ACC(2) = Z+ZZ EXACTLY
      ACC(1) = Z+ZZ
      ACC(2) = (Z-ACC(1))+ZZ
      RETURN
      END


C   IMSL ROUTINE NAME   - VXMUL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EXTENDED PRECISION MULTIPLY
C
C   USAGE               - CALL VXMUL (A,B,ACC)
C
C   ARGUMENTS    A      - INPUT DOUBLE PRECISION NUMBER
C                B      - INPUT DOUBLE PRECISION NUMBER
C                ACC    - ACCUMULATOR. (INPUT AND OUTPUT)
C                           ACC IS A DOUBLE PRECISION VECTOR OF LENGTH
C                           2.  ON OUTPUT, ACC CONTAINS THE SUM OF
C                           INPUT ACC AND A*B.
C
C   REQD. IMSL ROUTINES - VXADD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      VXMUL ADDS THE PRODUCT A*B TO THE EXTENDED PRECISION
C                ACCUMULATOR, ACC. THE SUBROUTINE ASSUMES THAT AN
C                EXTENDED PRECISION NUMBER IS ALREADY IN THE
C                ACCUMULATOR.  THEREFORE, BEFORE THE FIRST CALL TO
C                VXMUL, ACC(1) AND ACC(2) MUST BE SET TO ZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VXMUL (A,B,ACC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   A,B,ACC(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   X,HA,TA,HB,TB
      INTEGER            IX(2),I
      LOGICAL*1          LX(8),LI(4)
      EQUIVALENCE        (X,LX(1),IX(1)),(I,LI(1))
      DATA               I/0/
C                                  SPLIT A = HA+TA
C                                        B = HB+TB
C                                  FIRST EXECUTABLE STATEMENT
      X = A
      LI(4) = LX(5)
      IX(2) = 0
      I = (I/16)*16
      LX(5) = LI(4)
      HA=X
      TA=A-HA
      X = B
      LI(4) = LX(5)
      IX(2) = 0
      I = (I/16)*16
      LX(5) = LI(4)
      HB = X
      TB = B-HB
C                                  COMPUTE HA*HB,HA*TB,TA*HB, AND TA*TB
C                                    AND CALL VXADD TO ACCUMULATE THE
C                                    SUM
      X = TA*TB
      CALL VXADD(X,ACC)
      X = HA*TB
      CALL VXADD(X,ACC)
      X = TA*HB
      CALL VXADD(X,ACC)
      X = HA*HB
      CALL VXADD(X,ACC)
      RETURN
      END


C   IMSL ROUTINE NAME   - VXSTO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DOUBLE PRECISION STORE.
C
C   USAGE               - CALL VXSTO(ACC,D)
C
C   ARGUMENTS    ACC    - ACCUMULATOR. (INPUT)
C                           ACC IS A DOUBLE PRECISION VECTOR OF LENGTH
C                           2. ACC IS ASSUMED TO BE THE RESULT OF
C                           CALLING VXADD OR VXMUL TO PERFORM EXTENDED
C                           PRECISION OPERATIONS.
C                D      - DOUBLE PRECISION SCALAR. (OUTPUT)
C                           ON OUTPUT, D CONTAINS A DOUBLE PRECISION
C                           APPROXIMATION TO THE VALUE OF THE EXTENDED
C                           PRECISION ACCUMULATOR.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VXSTO (ACC,D)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   ACC(2),D
C                                  FIRST EXECUTABLE STATEMENT
      D = ACC(1)+ACC(2)
      RETURN
      END


C   IMSL ROUTINE NAME   - USPKD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE
C                           CHARACTER STRING ARGUMENTS
C
C   USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C
C   ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)
C                NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.
C                UNPAKD - INTEGER ARRAY TO RECEIVE THE UNPACKED
C                         REPRESENTATION OF THE STRING. (OUTPUT)
C                NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE
C
C   REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO AN INTEGER ARRAY
C                IN (A1) FORMAT.
C            2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF
C                THAT ARE IGNORED.
C
C-----------------------------------------------------------------------
      SUBROUTINE USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,NCHARS,NCHMTB
C
      LOGICAL*1          UNPAKD(1),PACKED(1),LBYTE,LBLANK
      INTEGER*2          IBYTE,IBLANK
      EQUIVALENCE (LBYTE,IBYTE)
      DATA               LBLANK /1H /
      DATA               IBYTE /1H /
      DATA               IBLANK /1H /
C                                  INITIALIZE NCHMTB
      NCHMTB = 0
C                                  RETURN IF NCHARS IS LE ZERO
      IF(NCHARS.LE.0) RETURN
C                                  SET NC=NUMBER OF CHARS TO BE DECODED
      NC = MIN0 (129,NCHARS)
      NWORDS = NC*4
      J = 1
      DO 110 I = 1,NWORDS,4
      UNPAKD(I) = PACKED(J)
      UNPAKD(I+1) = LBLANK
      UNPAKD(I+2) = LBLANK
      UNPAKD(I+3) = LBLANK
  110 J = J+1
C                                  CHECK UNPAKD ARRAY AND SET NCHMTB
C                                  BASED ON TRAILING BLANKS FOUND
      DO 200 N = 1,NWORDS,4
         NN = NWORDS - N - 2
         LBYTE = UNPAKD(NN)
         IF(IBYTE .NE. IBLANK) GO TO 210
  200 CONTINUE
      NN = 0
  210 NCHMTB = (NN + 3) / 4
      RETURN
      END


C   IMSL ROUTINE NAME   - OFROTA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ORTHOGONAL ROTATION OF A FACTOR LOADING MATRIX
C                           USING A GENERALIZED ORTHOMAX CRITERION,
C                           INCLUDING QUARTIMAX, VARIMAX, AND EQUAMAX
C
C   USAGE               - CALL OFROTA (A,IA,NV,NF,NORM,II,MAXIT,W,
C                           EPS,DELTA,B,IB,T,IT,F,WK,IER)
C
C   ARGUMENTS    A      - INPUT NV BY NF UNROTATED FACTOR LOADING
C                           MATRIX.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NV     - INPUT NUMBER OF VARIABLES.
C                NF     - INPUT NUMBER OF FACTORS.
C                NORM   - INPUT OPTION PARAMETER.  NORM=0 INDICATES NO
C                           ROW NORMALIZATION OF A IS REQUIRED.
C                           OTHERWISE, ROW NORMALIZATION IS PERFORMED.
C                II     - INPUT OPTION PARAMETER.  II=0 INDICATES AN
C                           IMAGE ANALYSIS IS NOT BEING PERFORMED.
C                           OTHERWISE, AN IMAGE ANALYSIS IS ASSUMED, AND
C                           T BECOMES THE IMAGE TRANSFORMATION MATRIX.
C                MAXIT  - INPUT MAXIMUM NUMBER OF ITERATIONS ALLOWED
C                           FOR ROTATION.  MAXIT=30 IS TYPICAL.
C                W      - INPUT CONSTANT FOR ROTATION (SEE REMARKS).
C                EPS    - INPUT CONVERGENCE CONSTANT FOR ROTATION
C                           (ANGLE).  EPS=0.0001 IS TYPICAL.
C                DELTA  - INPUT CONVERGENCE CONSTANT FOR ROTATION
C                           (CRITERION FUNCTION). DELTA=.001 IS TYPICAL.
C                B      - OUTPUT NV BY NF ORTHOGONALLY ROTATED FACTOR
C                           LOADING MATRIX.
C                IB     - INPUT ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                T      - OUTPUT NF BY NF TRANSFORMATION MATRIX.
C                           IF II IS NON-ZERO, T CONTAINS THE
C                           IMAGE TRANSFORMATION MATRIX.
C                IT     - INPUT ROW DIMENSION OF MATRIX T EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                F      - OUTPUT VECTOR OF LENGTH NF.  F(I) CONTAINS
C                           THE VARIANCE ACCOUNTED FOR BY FACTOR I.
C                WK     - WORK VECTOR. THE LENGTH OF WK IS EQUAL TO
C                           THE MAXIMUM OF (NV,(NF*(NF+1))/2).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT AT LEAST ONE OF
C                             NV, NF, IA, IT, OR IB WAS SPECIFIED
C                             INCORRECTLY.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66 INDICATES CONVERGENCE DID NOT
C                             OCCUR IN MAXIT ITERATIONS. CONVERGENCE
C                             WAS ASSUMED AND CALCULATIONS CONTINUED.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1P,LUDECP,LUELMP,OFIMA3,UERTST,
C                           UGETIO,VIPRFF,VMULFM,VTPROF
C                       - DOUBLE/LEQT1P,LUDECP,LUELMP,OFIMA3,UERTST,
C                           UGETIO,VIPRFF,VMULFM,VTPROF,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      W IS A PARAMETER DETERMINING THE KIND OF SOLUTION TO
C                BE COMPUTED. W MAY BE SET AS FOLLOWS.
C
C                W = 0.0 IS THE QUARTIMAX METHOD, WHICH ATTEMPTS TO
C                GET EACH VARIABLE TO LOAD HIGHLY ON ONLY ONE (OR A
C                FEW) FACTOR(S).
C
C                W = 1.0 IS THE VARIMAX METHOD, WHICH ATTEMPTS TO LOAD
C                HIGHLY A RELATIVELY LOW NUMBER OF VARIABLES ON EACH
C                FACTOR. VARIMAX IS MOST WIDELY USED.
C
C                W = NF/2.0 IS THE EQUAMAX METHOD, WHICH IS A COM-
C                PROMISE OF THE ABOVE TWO.
C
C                W CAN BE ANY REAL NUMBER, BUT BEST VALUES LIE IN THE
C                CLOSED INTERVAL (1.0, 5.0*NF). GENERALLY, THE LARGER
C                W IS, THE MORE EQUAL IS THE DISPERSION OF THE VARIANCE
C                ACCOUNTED FOR ACROSS THE FACTORS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFROTA  (A,IA,NV,NF,NORM,II,MAXIT,W,EPS,
     1                   DELTA,B,IB,T,IT,F,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NV,NF,NORM,II,MAXIT,IB,IT,IER
      REAL               A(IA,NF),B(IB,NF),F(NF),T(IT,NF),WK(1),W,EPS,
     *                   DELTA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NFF,I,J,MV,NC,NCOUNT,NFM1,JP1,K
      REAL               AS,BB,BS,COSP,EPS4,FOURTH,ONE,HOLD,PHI,SINP,
     *                   TT,TVV,TWO,U,V,VV,VVV,WSNV,ZERO
      DOUBLE PRECISION   TEMP,SS,DD,SAVE,DNV
      DATA               ZERO,ONE,TWO,FOURTH/0.0E0,1.0E0,2.0E0,0.25E0/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (NF.LE.NV.AND.IA.GE.NV.AND.IB.GE.NV.AND.IT.GE.NF) GO TO 5
      IER = 129
      GO TO 9000
    5 IER = 0
      DNV = NV*NV
      NFF = ((NF-1)*NF)/2
      EPS4 = EPS*FOURTH
      WSNV = W/NV
      DO 15 I=1,NF
         DO 10 J=1,NF
            T(I,J) = ZERO
   10    CONTINUE
         T(I,I) = ONE
   15 CONTINUE
      IF (NORM.EQ.0) GO TO 35
C                                  ROW NORMALIZATION PERFORMED
      DO 30 I=1,NV
         TEMP = 0.0D0
         DO 20 J=1,NF
            TEMP = TEMP+DBLE(A(I,J))**2
   20    CONTINUE
         HOLD = DSQRT(TEMP)
         WK(I) = HOLD
         HOLD = ONE/HOLD
         DO 25 J=1,NF
            B(I,J) = A(I,J)*HOLD
   25    CONTINUE
   30 CONTINUE
      GO TO 50
   35 DO 45 I=1,NV
         DO 40 J=1,NF
            B(I,J) = A(I,J)
   40    CONTINUE
   45 CONTINUE
   50 MV = 1
      NC = 0
      NCOUNT = 0
      VVV = ZERO
C                                  BEGIN ORTHOMAX ITERATION
   55 MV = MV+1
      VV = VVV
C                                  CALCULATE ROTATION CRITERION
      TEMP = 0.0D0
      DO 65 J=1,NF
         SS = 0.0D0
         DD = 0.0D0
         DO 60 I=1,NV
            SAVE = DBLE(B(I,J))**2
            DD = DD+SAVE
            SS = SS+SAVE*SAVE
   60    CONTINUE
         TEMP = TEMP+((NV*SS)-(W*DD*DD))/DNV
   65 CONTINUE
      VVV = TEMP
      IF (NF.LE.1) GO TO 115
      IF (MV.LE.MAXIT) GO TO 70
      IER = 66
      GO TO 115
   70 TVV = VVV-VV
      IF (TVV.GT.DELTA*VV) GO TO 75
      NC = NC+1
      IF (NC.GE.2) GO TO 115
   75 NFM1 = NF-1
      DO 110 J=1,NFM1
         JP1 = J+1
         DO 105 K=JP1,NF
C                                  CALCULATE RATIO OF KAISER TAN(4*PHI)
            AS = ZERO
            BS = ZERO
            TT = ZERO
            BB = ZERO
            DO 80 I=1,NV
               U = (B(I,J)+B(I,K))*(B(I,J)-B(I,K))
               V = TWO*B(I,J)*B(I,K)
               AS = AS+U
               BS = BS+V
               BB = BB+(U+V)*(U-V)
               TT = TT+U*V
   80       CONTINUE
            TT = TT+TT
            TT = TT-TWO*AS*BS*WSNV
            BB = BB-(AS+BS)*(AS-BS)*WSNV
            IF (ABS(TT)+ABS(BB).GT.EPS) GO TO 90
   85       NCOUNT = NCOUNT+1
            IF (NCOUNT.LT.NFF) GO TO 105
C                                  COMPLETE CYCLE WITHOUT ROTATION
            GO TO 115
   90       PHI = FOURTH*ATAN2(TT,BB)
C                                  IS ANGLE ESSENTIALLY ZERO
            IF (ABS(PHI).LT.EPS4) GO TO 85
            COSP = COS(PHI)
            SINP = SIN(PHI)
            NCOUNT = 0
C                                  ROTATE MATRICES BY ANGLE PHI
            DO 95 I=1,NV
               SAVE = B(I,J)*COSP+B(I,K)*SINP
               B(I,K) = -B(I,J)*SINP+B(I,K)*COSP
               B(I,J) = SAVE
   95       CONTINUE
            DO 100 I=1,NF
               SAVE = T(I,J)*COSP+T(I,K)*SINP
               T(I,K) = -T(I,J)*SINP+T(I,K)*COSP
               T(I,J) = SAVE
  100       CONTINUE
  105    CONTINUE
  110 CONTINUE
      GO TO 55
  115 IF (NORM.EQ.0) GO TO 130
C                                  RESCALING
      DO 125 I=1,NV
         HOLD = WK(I)
         DO 120 J=1,NF
            B(I,J) = B(I,J)*HOLD
  120    CONTINUE
  125 CONTINUE
C                                  VARIANCE ACCOUNTED FOR BY FACTORS
  130 DO 140 I=1,NF
         TEMP = 0.0D0
         DO 135 K=1,NV
            TEMP = TEMP+DBLE(B(K,I))**2
  135    CONTINUE
         F(I) = TEMP
  140 CONTINUE
      IF (II.NE.0) CALL OFIMA3 (A,IA,B,IB,NV,NF,NF,T,IT,WK,IER)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFROTA)
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - OFIMA3
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LEAST SQUARES SOLUTION TO THE MATRIX
C                           EQUATION AT = B.
C
C   USAGE               - CALL OFIMA3 (A,IA,B,IB,NV,NS,NF,T,IT,WK,IER)
C
C   ARGUMENTS    A      - INPUT MATRIX OF DIMENSION NV BY NF OF
C                           COLUMN RANK NF.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                B      - INPUT NV BY NS MATRIX.
C                IB     - INPUT ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NV     - INPUT NUMBER OF ROWS IN MATRICES A AND B.
C                NS     - INPUT NUMBER OF COLUMNS IN MATRICES T AND B.
C                NF     - INPUT NUMBER OF ROWS IN MATRIX T AND
C                           NUMBER OF COLUMNS IN MATRIX A.
C                T      - OUTPUT MATRIX OF DIMENSION NF BY NS CONTAINING
C                           THE LEAST SQUARES EQUATION SOLUTION.
C                IT     - INPUT ROW DIMENSION OF MATRIX T EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                WK     - WORK VECTOR OF LENGTH (NF+1)*NF/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE RANK OF A IS
C                             LESS THAN NF NUMERICALLY.
C                           IER = 130 INDICATES AT LEAST ONE OF IA, IB,
C                             OR IT WAS SPECIFIED INCORRECTLY.
C
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1P,LUDECP,LUELMP,UERTST,UGETIO,
C                           VIPRFF,VMULFM,VTPROF
C                       - DOUBLE/LEQT1P,LUDECP,LUELMP,UERTST,UGETIO,
C                           VIPRFF,VMULFM,VTPROF,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFIMA3 (A,IA,B,IB,NV,NS,NF,T,IT,WK,IER)
C
      REAL               D1,D2
      REAL               A,B,T,WK
      DIMENSION          A(IA,NF),B(IB,NS),T(IT,NS),WK(1)
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (NV.LE.IA.AND.NV.LE.IB.AND.NF.LE.IT.AND.NF.LE.NV) GO TO 5
      IER = 130
      GO TO 9000
C                                  CALCULATE (A-TRANSPOSE) * A
    5 CALL VTPROF (A,NV,NF,IA,WK)
C                                  CALCULATE (A-TRANSPOSE) * B
      CALL VMULFM (A,B,NV,NF,NS,IA,IB,T,IT,IER)
C                                  SOLVE FOR T MATRIX
      CALL LEQT1P (WK,NS,NF,T,IT,0,D1,D2,IER)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFIMA3)
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - VTPROF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TRANSPOSE PRODUCT OF MATRIX (FULL STORAGE
C                           MODE)
C
C   USAGE               - CALL VTPROF (A,L,M,IA,ATA)
C
C   ARGUMENTS    A      - L BY M MATRIX IN FULL STORAGE MODE. (INPUT)
C                L      - MAXIMUM VALUE OF FIRST SUBSCRIPT OF A.
C                           (INPUT)
C                M      - MAXIMUM VALUE OF SECOND SUBSCRIPT OF A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF A EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                ATA    - M BY M MATRIX STORED IN SYMMETRIC MODE AS A
C                           VECTOR. (OUTPUT)
C
C   REQD. IMSL ROUTINES - SINGLE/VIPRFF
C                       - DOUBLE/VIPRFF,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VTPROF (A,L,M,IA,ATA)
C
      DIMENSION          A(IA,1),ATA(1)
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I=1,M
         K=I*(I-1)/2
            DO 5 J=1,I
            K=K+1
    5 CALL VIPRFF(A(1,I),A(1,J),L,1,1,ATA(K))
      RETURN
      END


C   IMSL ROUTINE NAME   - VIPRFF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - VECTOR INNER PRODUCT OF TWO VECTORS OR
C                           SUBSETS OF TWO VECTORS
C
C   USAGE               - CALL VIPRFF (X,Y,L,IX,IY,XYIP)
C
C   ARGUMENTS    X      - FIRST VECTOR (ADDRESS OF THE FIRST ELEMENT).
C                           (INPUT)
C                Y      - SECOND VECTOR (ADDRESS OF THE FIRST ELEMENT).
C                           (INPUT)
C                L      - NUMBER OF ELEMENTS OF VECTOR X OR Y INVOLVED
C                           IN THE INNER PRODUCT. (INPUT)
C                IX     - INCREMENT BETWEEN SUCCESSIVE ELEMENTS OF X.
C                           (INPUT)
C                IY     - INCREMENT BETWEEN SUCCESSIVE ELEMENTS OF Y.
C                           (INPUT)
C                XYIP   - INNER PRODUCT. (OUTPUT)  XYIP = X(1)*Y(1)+
C                           X(1+IX)*Y(1+IY)+...+X(1+L*IX)*Y(1+L*IY).
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SDOT
C                       - DOUBLE/VBLA=DDOT
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VIPRFF (X,Y,L,IX,IY,XYIP)
C                                 SPECIFICATIONS FOR ARGUMENTS
      REAL               X(1),XYIP,Y(1)
      INTEGER            IX,IY,L
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               SDOT
C                                  FIRST EXECUTABLE STATEMENT
      XYIP = SDOT(L,X,IX,Y,IY)
      RETURN
      END


C   IMSL ROUTINE NAME   - VBLA=SDOT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE SINGLE PRECISION DOT PRODUCT
C
C   USAGE               - FUNCTION SDOT (N,SX,INCX,SY,INCY)
C
C   ARGUMENTS    SDOT   - SUM FROM I=1 TO N OF X(I)*Y(I). (OUTPUT)
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                SX     - REAL VECTOR OF LENGTH MAX(N*IABS(INCX),1).
C                           (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF SX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           SX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           SX(1+(I-N)*INCX) IF INCX.LT.0.
C                SY     - REAL VECTOR OF LENGTH MAX(N*IABS(INCY),1).
C                           (INPUT)
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF SY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           SY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           SY(1+(I-N)*INCY) IF INCY.LT.0.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION SDOT (N,SX,INCX,SY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY
      REAL               SX(1),SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,MP1,NS,IX,IY
C                                  FIRST EXECUTABLE STATEMENT
      SDOT = 0.0E0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 5,15,35
    5 CONTINUE
C                                  CODE FOR UNEQUAL INCREMENTS OR
C                                    NONPOSITIVE INCREMENTS.
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I=1,N
         SDOT = SDOT+SX(IX)*SY(IY)
         IX = IX+INCX
         IY = IY+INCY
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 5.
   15 M = N-(N/5)*5
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         SDOT = SDOT+SX(I)*SY(I)
   20 CONTINUE
      IF (N.LT.5) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,5
         SDOT = SDOT+SX(I)*SY(I)+SX(I+1)*SY(I+1)+SX(I+2)*SY(I+2)+SX(I
     1   +3)*SY(I+3)+SX(I+4)*SY(I+4)
   30 CONTINUE
      RETURN
C                                  CODE FOR POSITIVE EQUAL INCREMENTS
C                                    .NE.1.
   35 CONTINUE
      NS = N*INCX
      DO 40 I=1,NS,INCX
         SDOT = SDOT+SX(I)*SY(I)
   40 CONTINUE
      RETURN
      END


C   IMSL ROUTINE NAME   - VMULFM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - MATRIX MULTIPLICATION OF THE TRANSPOSE OF
C                           MATRIX A BY MATRIX B (FULL STORAGE MODE)
C
C   USAGE               - CALL VMULFM (A,B,L,M,N,IA,IB,C,IC,IER)
C
C   ARGUMENTS    A      - L BY M MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                B      - L BY N MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                L      - NUMBER OF ROWS IN A AND B. (INPUT)
C                M      - NUMBER OF COLUMNS IN MATRIX A AND NUMBER
C                           OF ROWS IN MATRIX C. (INPUT)
C                N      - NUMBER OF COLUMNS IN B AND C. (INPUT)
C                IA     - ROW DIMENSION OF A EXACTLY AS SPECIFIED IN
C                           IN THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT OF THE CALLING
C                           PROGRAM. (INPUT)
C                C      - M BY N MATRIX CONTAINING THE PRODUCT
C                           C = (A-TRANSPOSE) * B. (OUTPUT)
C                IC     - ROW DIMENSION OF C EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         IER=0 IMPLIES NO ERROR.
C                         TERMINAL ERROR
C                           IER=129 INDICATES A,B, OR C WAS DIMENSIONED
C                             INCORRECTLY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VMULFM (A,B,L,M,N,IA,IB,C,IC,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L,M,N,IA,IB,IC,IER
      REAL               A(IA,1),B(IB,1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TEMP
C                                  FIRST EXECUTABLE STATEMENT
      IF (IA.GE.L .AND. IB.GE.L .AND. IC.GE.M) GO TO 5
C                                  TERMINAL ERROR
      IER = 129
      GO TO 9000
C                                  ROW INDICATOR
    5 IER = 0
      DO 20 I = 1,M
C                                  COLUMN INDICATOR
         DO 15 J = 1,N
            TEMP = 0.0D0
C                                  VECTOR DOT PRODUCT
            DO 10 K = 1,L
               TEMP = TEMP + DBLE(A(K,I))*DBLE(B(K,J))
   10       CONTINUE
            C(I,J) = TEMP
   15    CONTINUE
   20 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HVMULFM)
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - LEQT1P
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE DEFINITE
C                           MATRIX - SYMMETRIC STORAGE MODE - SPACE
C                           ECONOMIZER SOLUTION
C
C   USAGE               - CALL LEQT1P (A,M,N,B,IB,IDGT,D1,D2,IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING THE
C                           N BY N COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B. A IS A POSITIVE DEFINITE SYMMETRIC
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
C                         ON OUTPUT, A IS REPLACED BY THE LOWER
C                           TRIANGULAR MATRIX L WHERE A = L*L-TRANSPOSE.
C                           L IS STORED IN SYMMETRIC STORAGE MODE WITH
C                           THE DIAGONAL ELEMENTS OF L IN RECIPROCAL
C                           FORM.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
C                B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           THE RIGHT-HAND SIDES OF THE EQUATION
C                           AX = B.
C                         ON OUTPUT, THE N BY M SOLUTION MATRIX X
C                           REPLACES B.
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT
C                           TO IDGT DECIMAL DIGITS. (CURRENTLY NOT USED)
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE INPUT MATRIX
C                             A IS ALGORITHMICALLY NOT POSITIVE
C                             DEFINITE. (SEE THE CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - LUDECP,LUELMP,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQT1P (A,M,N,B,IB,IDGT,D1,D2,IER)
C
      DIMENSION          A(1),B(IB,1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
C                                  DECOMPOSE A
      CALL LUDECP (A,A,N,D1,D2,IER)
      IF (IER.NE.0) GO TO 9000
C                                  PERFORM ELIMINATION
      DO 5 I = 1,M
         CALL LUELMP (A,B(1,I),N,B(1,I))
    5 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HLEQT1P)
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - LUDECP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DECOMPOSITION OF A POSITIVE DEFINITE MATRIX -
C                           SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL LUDECP (A,UL,N,D1,D2,IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE N BY N POSITIVE DEFINITE SYMMETRIC
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
C                UL     - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE DECOMPOSED MATRIX L SUCH THAT A = L*
C                           L-TRANSPOSE. L IS STORED IN SYMMETRIC
C                           STORAGE MODE. THE DIAGONAL OF L CONTAINS THE
C                           RECIPROCALS OF THE ACTUAL DIAGONAL ELEMENTS.
C                N      - ORDER OF A. (INPUT)
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C                             (SEE THE CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUDECP (A,UL,N,D1,D2,IER)
C
      DIMENSION          A(1),UL(1)
      DATA               ZERO,ONE,FOUR,SIXTN,SIXTH/
     *                   0.0,1.,4.,16.,.0625/
C                                  FIRST EXECUTABLE STATEMENT
      D1=ONE
      D2=ZERO
      RN = ONE/(N*SIXTN)
      IP = 1
      IER=0
      DO 45 I = 1,N
         IQ = IP
         IR = 1
         DO 40 J = 1,I
            X = A(IP)
            IF (J .EQ. 1) GO TO 10
            DO 5  K=IQ,IP1
               X = X - UL(K) * UL(IR)
               IR = IR+1
    5       CONTINUE
   10       IF (I.NE.J) GO TO 30
            D1 = D1*X
            IF (A(IP) + X*RN .LE. A(IP)) GO TO 50
   15       IF (ABS(D1).LE.ONE) GO TO 20
            D1 = D1 * SIXTH
            D2 = D2 + FOUR
            GO TO 15
   20       IF (ABS(D1) .GE. SIXTH) GO TO 25
            D1 = D1 * SIXTN
            D2 = D2 - FOUR
            GO TO 20
   25       UL(IP) = ONE/SQRT(X)
            GO TO 35
   30       UL(IP) = X * UL(IR)
   35       IP1 = IP
            IP = IP+1
            IR = IR+1
   40    CONTINUE
   45 CONTINUE
      GO TO 9005
   50 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HLUDECP)
 9005 RETURN
      END


      SUBROUTINE LUELMP (A,B,N,X)
      DIMENSION      A(1),B(1),X(1)
      DATA           ZERO/0.0/
C                                   first executable statement
C                                   solution of lY = B
      IP = 1
      IW = 0
      DO 15 I =1,N
         T = B(I)
         IM1 = I - 1
         IF (IW .EQ. 0) GO TO 9
         IP = IP + IW - 1
         DO 5 K = IW,IM1
            T = T - A(IP) * X(K)
            IP = IP + 1
    5    CONTINUE
         GO TO 10
    9    IF (T .NE. ZERO) IW = I
         IP = IP + IM1
   10    X(I) = T * A(IP)
         IP = IP + 1
   15 CONTINUE
C                                               solution of UX = Y
      N1 = N + 1
      DO 30 I = 1,N
         II = N1 - I
         IP = IP - 1
         IS = IP
         IQ = II + 1
         T = X(II)
         IF (N .LT. IQ) GO TO 25
         KK = N
         DO 20 K = IQ,N
            T = T -A(IS) * X(KK)
            KK = KK - 1
            IS = IS - KK
   20    CONTINUE
   25    X(II) = T * A(IS)
   30 CONTINUE
      RETURN
      END


C   IMSL ROUTINE NAME   - LGINF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERALIZED INVERSE OF A REAL MATRIX
C
C   USAGE               - CALL LGINF (A,IA,M,N,TOL,AINV,IAINV,S,WK,IER)
C
C   ARGUMENTS    A      - M BY N MATRIX. (INPUT)  A IS DESTROYED.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                           A IS USED BY LGINF AS WORK STORAGE FOR AN
C                           N BY N MATRIX. THEREFORE, IA MUST BE
C                           GREATER THAN OR EQUAL TO MAX(M,N).
C                M      - NUMBER OF ROWS IN A. (INPUT)
C                N      - NUMBER OF COLUMNS IN A. (INPUT)
C                TOL    - TOLERANCE PARAMETER. (INPUT)
C                           IF TOL IS LESS THAN OR EQUAL TO ZERO ON
C                           INPUT, LGINF COMPUTES THE GENERALIZED
C                           INVERSE OF A.
C                           IF TOL IS GREATER THAN ZERO ON INPUT, LGINF
C                           COMPUTES THE GENERALIZED INVERSE OF A MATRIX
C                           CLOSE TO A, BUT HAVING CONDITION NUMBER LESS
C                           THAN 1.0/TOL.
C                AINV   - N BY M MATRIX. (OUTPUT)
C                           AINV CONTAINS THE GENERALIZED INVERSE OF A
C                           OR THE GENERALIZED INVERSE OF A MATRIX
C                           CLOSE TO A. SEE TOL ARGUMENT DESCRIPTION.
C                IAINV  - ROW DIMENSION OF MATRIX AINV EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                S      - VECTOR OF LENGTH N. S CONTAINS THE ORDERED
C                           SINGULAR VALUES OF A.  S(1) .GE. S(2),...,
C                           .GE. S(N) .GE. 0. (OUTPUT)
C                WK     - WORK VECTOR OF LENGTH 2N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT MATRIX A IS NOT
C                             FULL RANK OR VERY ILL-CONDITIONED. SMALL
C                             SINGULAR VALUES MAY NOT BE VERY ACCURATE.
C                           IER=34 INDICATES THAT EITHER N.LE.0 OR
C                             M.LE.0.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT CONVERGENCE WAS
C                             NOT OBTAINED BY LSVDB AND COMPUTATION
C                             WAS DISCONTINUED.
C
C   REQD. IMSL ROUTINES - LGING,LSVDB,LSVG1,LSVG2,VHS12,UERSET,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LGINF  (A,IA,M,N,TOL,AINV,IAINV,S,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IAINV,IER,M,N
      REAL               A(IA,N),S(N),TOL,WK(N,2),AINV(IAINV,M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L
      REAL               T,TOLL
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      DO 10 I=1,N
         DO 5 J=1,M
            AINV(I,J)=0.0
            IF (I.EQ.J) AINV(I,J)=1.0
    5    CONTINUE
   10 CONTINUE
      CALL LGING(A,IA,M,N,AINV,IAINV,M,S,WK,IER)
      IF (IER.GT.128) GO TO 9000
      L = MIN0(M,N)
C                                  U**T HAS BEEN PLACED IN AINV
      TOLL = S(1)*AMAX1(TOL,0.0)
C                                  COMPUTE Q**(+) U**T
      DO 30 I=1,L
         IF (S(I).LE.TOLL) GO TO 20
         DO 15 J=1,M
   15    AINV(I,J)=AINV(I,J)/S(I)
         GO TO 30
   20    DO 25 J=1,M
   25    AINV(I,J)=0.0
   30 CONTINUE
C                                  COMPUTE V Q**(+) U**T
      DO 50 J=1,M
         DO 40 I=1,N
            T = 0.0
            DO 35 K=1,L
   35       T = T+A(I,K)*AINV(K,J)
            WK(I,1)=T
   40    CONTINUE
         DO 45 I=1,N
   45    AINV(I,J)=WK(I,1)
   50 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLGINF )
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - LGING
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - NUCLEUS USED ONLY BY IMSL ROUTINE LGINF
C
C
C   REQD. IMSL ROUTINES - LSVDB,LSVG1,LSVG2,VHS12,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LGING (A,IA,M,N,B,IB,NB,S,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,M,N,IB,NB,IER
      REAL               A(IA,N),B(IB,1),S(N),WK(N,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JP1,K,L,LNM,MM,NN,NNP1,NS,NSP1,NM
      REAL               ZERO,ONE,T,F,FF
      DATA               ZERO /0.0/,ONE /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  BEGIN SPECIAL FOR ZERO ROWS AND
C                                    COLS. PACK THE NONZERO COLS TO THE
C                                    LEFT
      NN = N
      IER = 34
      IF (NN.LE.0 .OR. M.LE.0) GO TO 9005
      IER = 0
      J = NN
    5 CONTINUE
      DO 10 I=1,M
         IF (A(I,J).NE.ZERO) GO TO 25
   10 CONTINUE
C                                  COL J IS ZERO. EXCHANGE IT WITH COL
C                                    NN
      IF (J.EQ.NN) GO TO 20
      DO 15 I=1,M
   15 A(I,J) = A(I,NN)
   20 CONTINUE
      A(1,NN) = J
      NN = NN-1
   25 CONTINUE
      J = J-1
      IF (J.GE.1) GO TO 5
C                                  IF NN=0 THEN A IS ENTIRELY ZERO AND
C                                    SVD COMPUTATION CAN BE SKIPPED
      NS = 0
      IF (NN.EQ.0) GO TO 135
C                                  PACK NONZERO ROWS TO THE TOP QUIT
C                                    PACKING IF FIND N NONZERO ROWS
      DO 30 I=1,N
   30 WK(I,1) = I
      I = 1
      MM = M
   35 IF (I.GT.N .OR. I.GE.MM) GO TO 70
      IF (A(I,I).NE.ZERO) GO TO 45
      DO 40 J=1,NN
         IF (A(I,J).NE.ZERO) GO TO 45
   40 CONTINUE
      GO TO 50
   45 I = I+1
      GO TO 35
C                                  ROW I IS ZERO EXCHANGE ROWS I AND MM
C                                  AND RECORD THE PERMUTATION IN WK(I,1)
   50 WK(I,1) = MM
      DO 55 J=1,NN
   55 A(I,J) = A(MM,J)
      IF (MM.GT.NN) GO TO 65
      DO 60 J=1,NN
   60 A(MM,J) = ZERO
   65 CONTINUE
C                                  EXCHANGE IS FINISHED
      MM = MM-1
      GO TO 35
C
   70 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
C                                  BEGIN SVD ALGORITHM..
C                                  (1) REDUCE THE MATRIX TO UPPER
C                                    BIDIAGONAL FORM WITH HOUSEHOLDER
C                                    TRANSFORMATIONS.
C                                    H(N)...H(1)AQ(1)...Q(N-2) =
C                                    (D**T,0)**T WHERE D IS UPPER
C                                    BIDIAGONAL.
C                                  (2) APPLY H(N)...H(1) TO B. HERE
C                                    H(N)...H(1)*B REPLACES B IN
C                                    STORAGE.
C                                  (3) THE MATRIX PRODUCT W=
C                                    Q(1)...Q(N-2) OVERWRITES THE FIRST
C                                    N ROWS OF A IN STORAGE.
C                                  (4) AN SVD FOR D IS COMPUTED. HERE K
C                                    ROTATIONS RI AND PI ARE COMPUTED
C                                    SO THAT RK...R1*D*P1**(T)...PK**(T)
C                                    = DIAG(S1,...,SM) TO WORKING
C                                    ACCURACY. THE SI ARE NONNEGATIVE
C                                    AND NONINCREASING. HERE RK...R1*B
C                                    OVERWRITES B IN STORAGE WHILE
C                                    A*P1**(T)...PK**(T) OVERWRITES A
C                                    IN STORAGE.
C                                  (5) IT FOLLOWS THAT,WITH THE PROPER
C                                    DEFINITIONS, U**(T)*B OVERWRITES
C                                    B, WHILE V OVERWRITES THE FIRST N
C                                    ROW AND COLUMNS OF A.
      L = MIN0(MM,NN)
C                                  THE FOLLOWING LOOP REDUCES A TO
C                                    UPPER BIDIAGONAL AND ALSO APPLIES
C                                    THE PREMULTIPLYING TRANSFORMATIONS
C                                    TO B.
      DO 80 J=1,L
         IF (J.GE.MM) GO TO 75
         JP1 = MIN0(J+1,NN)
         CALL VHS12(1,J,J+1,MM,A(1,J),1,T,A(1,JP1),1,IA,NN-J)
         S(J) = T
   75    IF (J.GE.NN-1) GO TO 80
         CALL VHS12(1,J+1,J+2,NN,A(J,1),IA,WK(J,2),A(J+1,1),IA,1,MM-J)
   80 CONTINUE
C                                  CONSTRUCT THE FIRST N ROWS OF THE
C                                    MATRIX HL*...*H2*H1 = B
      DO 85 JJ=1,L
         J = L+1-JJ
         IF (J.GE.MM) GO TO 85
         T = S(J)
         CALL VHS12(2,J,J+1,MM,A(1,J),1,T,B,IB,1,N)
   85 CONTINUE
C                                  PERMUTE COLUMNS OF B ACCORDING TO THE
C                                  PERMUTATIONS MADE TO THE ROWS OF A
      LNM = MIN0(N,MM)
      IF (LNM.LT.1) GO TO 100
      DO 95 K=1,LNM
         J = WK(K,1)
         IF (J.EQ.K) GO TO 95
C                                  EXCHANGE COLUMNS J AND K
         DO 90 I=1,N
            T = B(I,K)
            B(I,K) = B(I,J)
            B(I,J) = T
   90    CONTINUE
   95 CONTINUE
  100 CONTINUE
C                                  COPY THE BIDIAGONAL MATRIX INTO THE
C                                    ARRAY S FOR LSVDB
      IF (L.EQ.1) GO TO 110
      DO 105 J=2,L
         S(J) = A(J,J)
         WK(J,1) = A(J-1,J)
  105 CONTINUE
  110 S(1) = A(1,1)
C
      NS = NN
      IF (MM.GE.NN) GO TO 115
      NS = MM+1
      S(NS) = ZERO
      WK(NS,1) = A(MM,MM+1)
  115 CONTINUE
C                                  CONSTRUCT THE EXPLICIT NN BY NN
C                                    PRODUCT MATRIX, W=Q1*Q2*...*QL*I
C                                    IN THE ARRAY A
      DO 130 K=1,NN
         I = NN+1-K
         IF (I.GT.MIN0(MM,NN-2)) GO TO 120
         CALL VHS12(2,I+1,I+2,NN,A(I,1),IA,WK(I,2),A(1,I+1),1,IA,NN-I)
  120    DO 125 J=1,NN
  125    A(I,J) = ZERO
         A(I,I) = ONE
  130 CONTINUE
C                                  COMPUTE THE SVD OF THE BIDIAGONAL
C                                    MATRIX
C
      LEVEL = 1
      CALL UERSET(LEVEL,LEVOLD)
      CALL LSVDB(S(1),WK(1,1),NS,A,IA,NN,B,IB,NB,IER)
C
      IF (IER.GT.128) GO TO 9005
      CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.NE.33) GO TO 135
      FF = 0.0
      NM = MIN0(N,M)
      IF (S(1).NE.ZERO) FF = S(NM)/S(1)
      F = 100.0+FF
      IF (F.EQ.100.0) GO TO 135
      IER = 0
  135 CONTINUE
      IF (NS.GE.NN) GO TO 145
      NSP1 = NS+1
      DO 140 J=NSP1,NN
  140 S(J) = ZERO
  145 CONTINUE
      IF (NN.EQ.N) GO TO 9005
      NNP1 = NN+1
C                                  MOVE RECORD OF PERMUTATIONS AND
C                                    STORE ZEROS
      DO 155 J=NNP1,N
         S(J) = A(1,J)
         IF (NN.LT.1) GO TO 155
         DO 150 I=1,NN
  150    A(I,J) = ZERO
  155 CONTINUE
C                                  PERMUTE ROWS AND SET ZERO SINGULAR
C                                    VALUES
      DO 165 K=NNP1,N
         I = S(K)
         S(K) = ZERO
         DO 160 J=1,N
            A(K,J) = A(I,J)
            A(I,J) = ZERO
  160    CONTINUE
         A(I,K) = ONE
  165 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
 9000 CONTINUE
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - VHS12
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REAL HOUSEHOLDER TRANSFORMATION -
C                           COMPUTATION AND APPLICATIONS.
C
C   USAGE               - CALL VHS12 (MODE,LP,L1,M,U,INCU,UP,C,INCC,
C                           ICV,NCV)
C
C   ARGUMENTS    MODE   - OPTION PARAMETER. (INPUT)
C                         IF MODE=1, THE SUBROUTINE COMPUTES A
C                           HOUSEHOLDER TRANSFORMATION AND IF NCV.GT.0,
C                           MULTIPLIES IT BY THE SET OF NCV VECTORS
C                           (EACH OF LENGTH M) STORED IN C. FOR A
C                           GIVEN VECTOR V OF LENGTH M AND TWO INTEGER
C                           INDICES LP AND L1 THAT SATISFY
C                           1 .LE. LP .LT. L1 .LE. M, THE SUBROUTINE
C                           DEFINES AN M BY M HOUSEHOLDER
C                           TRANSFORMATION Q WHICH SATIFIES QV=W WHERE
C                           W(I)=V(I) FOR I.LT.LP
C                           W(LP)=-SIG*SQRT(V(LP)**2+V(L1)**2+...
C                             +V(M)**2)
C                             SIG=1  IF V(LP).GE.0
C                             SIG=-1 IF V(LP).LT.0
C                           W(I)=V(I) FOR LP.LT.I.LT.L1
C                           W(I)=0    FOR I.GE.L1.
C                         IF MODE=2, THE SUBROUTINE ASSUMES THAT A
C                           HOUSEHOLDER TRANSFORMATION HAS ALREADY
C                           BEEN DEFINED BY A PREVIOUS CALL WITH
C                           MODE=1, AND IF NCV.GT.0, MULTIPLIES IT BY
C                           THE SET OF NCV VECTORS (EACH OF LENGTH
C                           M) STORED IN C.
C                LP     - PARAMETERS THAT DEFINE THE DESIRED
C                L1         HOUSEHOLDER TRANSFORMATION. (INPUT)
C                M          IF THE CONDITION 1.LE.LP.LT.L1.LE.M IS
C                           NOT SATISFIED, THE SUBROUTINE RETURNS TO
C                           THE CALLING PROGRAM WITHOUT PERFORMING
C                           ANY COMPUTATIONS.
C                U      - VECTOR OF M ELEMENTS. (INPUT, AND OUTPUT IF
C                           MODE=1)
C                           THE STORAGE INCREMENT BETWEEN ELEMENTS
C                           OF U IS INCU. (I.E., U(1+(J-1)*INCU),
C                           J=1,...,M). IF MODE=1, THE ARRAY V IS
C                           DEFINED AS V(J)=U(1+(J-1)*INCU),
C                           J=1,...,M.
C                         ON OUTPUT, U(1+(LP-1)*INCU) IS SET TO
C                           W(LP) (AS DEFINED ABOVE IN THE DESCRIPTION
C                           OF MODE=1).
C                INCU   - INCREMENT BETWEEN ELEMENTS OF U. (INPUT)
C                UP     - SCALAR SET TO V(LP)-W(LP) TO DEFINE THE
C                           HOUSEHOLDER TRANSFORMATION Q. (INPUT IF
C                           MODE=2, OUTPUT IF MODE=1)
C                C      - VECTOR OF NCV*M ELEMENTS. (INPUT/OUTPUT)
C                           IF NCV.LE.0, C IS NOT USED.
C                           IF NCV.GT.0, C CONTAINS NCV VECTORS
C                           OF LENGTH M WITH INCREMENT INCC BETWEEN
C                           ELEMENTS OF VECTORS AND INCREMENT ICV
C                           BETWEEN VECTORS. ELEMENT I OF VECTOR J IS
C                           DEFINED AS C(1+(I-1)*INCC+(J-1)*ICV),
C                           I=1,...,M AND J=1,...,NCV.
C                         ON OUTPUT, C CONTAINS THE SET OF NCV
C                           VECTORS RESULTING FROM MULTIPLYING
C                           THE GIVEN VECTORS BY Q.
C                INCC   - INCREMENT BETWEEN ELEMENTS OF VECTORS
C                           IN C. (INPUT)
C                ICV    - INCREMENT BETWEEN VECTORS IN C. (INPUT)
C                NCV    - NUMBER OF VECTORS STORED IN C. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF U IS A SINGLE SUBSCRIPTED ARRAY OR THE J-TH COLUMN
C                OF A MATRIX, THEN INCU=1. IF U IS THE I-TH ROW OF A
C                MATRIX THEN INCU IS THE ROW DIMENSION OF THE MATRIX
C                EXACTLY AS SPECIFIED IN THE CALLING PROGRAM.
C            2.  IF C IS A DOUBLE SUBSCRIPTED MATRIX AND THE VECTORS
C                ARE THE FIRST NCV COLUMNS OF C, THEN INCC=1 AND ICV
C                IS THE ROW DIMENSION OF C EXACTLY AS SPECIFIED IN
C                THE CALLING PROGRAM. IN THIS CASE C IS REPLACED
C                BY QC. IF THE VECTORS ARE SUCCESSIVE ROWS OF C
C                THEN INCC IS THE ROW DIMENSION OF C EXACTLY AS
C                SPECIFIED IN THE CALLING PROGRAM AND ICV=1. IN THIS
C                CASE C IS REPLACED BY CQ.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VHS12  (MODE,LP,L1,M,U,INCU,UP,C,INCC,ICV,NCV)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MODE,LP,L1,M,INCU,INCC,ICV,NCV
      REAL               U(1),UP,C(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IJ,ILP,IL1,IM,INCR,I2,I3,I4,J
      DOUBLE PRECISION   SM,B
      REAL               ONE,CL,CLINV,SM1
C                                  FIRST EXECUTABLE STATEMENT
      ONE = 1.
C
      IF (0.GE.LP.OR.LP.GE.L1.OR.L1.GT.M) GO TO 9005
      ILP = (LP-1)*INCU+1
      IL1 = (L1-1)*INCU+1
      IM = (M-1)*INCU+1
      CL = ABS(U(ILP))
      IF (MODE.EQ.2) GO TO 15
C                                  CONSTRUCT THE TRANSFORMATION.
      DO 5 IJ=IL1,IM,INCU
    5 CL = AMAX1(ABS(U(IJ)),CL)
      IF (CL.LE.0.0) GO TO 9005
      CLINV = ONE/CL
      SM = (DBLE(U(ILP))*CLINV)**2
      DO 10 IJ=IL1,IM,INCU
   10 SM = SM+(DBLE(U(IJ))*CLINV)**2
C                                  CONVERT DBLE. PREC. SM TO SNGL.
C                                    PREC. SM1
      SM1 = SM
      CL = CL*SQRT(SM1)
      IF (U(ILP).GT.0.0) CL = -CL
      UP = U(ILP)-CL
      U(ILP) = CL
      GO TO 20
C                                  APPLY THE TRANSFORMATION
C                                    I+U*(U**T)/B TO C.
   15 IF (CL.LE.0.0) GO TO 9005
   20 IF (NCV.LE.0) GO TO 9005
      B = DBLE(UP)*U(ILP)
C                                  B MUST BE NONPOSITIVE HERE. IF B =
C                                    0., RETURN.
      IF (B.GE.0.0) GO TO 9005
      B = ONE/B
      I2 = 1-ICV+INCC*(LP-1)
      INCR = INCC*(L1-LP)
      DO 35 J=1,NCV
         I2 = I2+ICV
         I3 = I2+INCR
         I4 = I3
         SM = C(I2)*DBLE(UP)
         DO 25 IJ=IL1,IM,INCU
            SM = SM+C(I3)*DBLE(U(IJ))
            I3 = I3+INCC
   25    CONTINUE
         IF (SM.EQ.0.0) GO TO 35
         SM = SM*B
         C(I2) = C(I2)+SM*DBLE(UP)
         DO 30 IJ=IL1,IM,INCU
            C(I4) = C(I4)+SM*DBLE(U(IJ))
            I4 = I4+INCC
   30    CONTINUE
   35 CONTINUE
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - UERSET
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SET MESSAGE LEVEL FOR IMSL ROUTINE UERTST
C
C   USAGE               - CALL UERSET (LEVEL,LEVOLD)
C
C   ARGUMENTS    LEVEL  - NEW VALUE FOR MESSAGE LEVEL. (INPUT)
C                           OUTPUT FROM IMSL ROUTINE UERTST IS
C                           CONTROLLED SELECTIVELY AS FOLLOWS,
C                             LEVEL = 4 CAUSES ALL MESSAGES TO BE
C                                       PRINTED,
C                             LEVEL = 3 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 32,
C                             LEVEL = 2 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 64,
C                             LEVEL = 1 MESSAGES ARE PRINTED IF IER IS
C                                       GREATER THAN 128,
C                             LEVEL = 0 ALL MESSAGE PRINTING IS
C                                       SUPPRESSED.
C                LEVOLD - PREVIOUS MESSAGE LEVEL. (OUTPUT)
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE UERSET (LEVEL,LEVOLD)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LEVEL,LEVOLD
C                                  FIRST EXECUTABLE STATEMENT
      LEVOLD = LEVEL
      CALL UERTST (LEVOLD,6HUERSET)
      RETURN
      END


C   IMSL ROUTINE NAME   - LSVDB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - SINGULAR VALUE DECOMPOSITION OF A BIDIAGONAL
C                           MATRIX.
C
C   USAGE               - CALL LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
C
C   ARGUMENTS    D      - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                         ON INPUT, D CONTAINS THE DIAGONAL ELEMENTS
C                           OF THE BIDIAGONAL MATRIX B. D(I)=B(I,I),
C                           I=1,...,N.
C                         ON OUTPUT, D CONTAINS THE N (NONNEGATIVE)
C                           SINGULAR VALUES OF B IN NONINCREASING
C                           ORDER.
C                E      - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                         ON INPUT, E CONTAINS THE SUPERDIAGONAL
C                           ELEMENTS OF B. E(1) IS ARBITRARY,
C                           E(I)=B(I-1,I), I=2,...,N.
C                         ON OUTPUT, THE CONTENTS OF E ARE MODIFIED
C                           BY THE SUBROUTINE.
C                N      - ORDER OF THE MATRIX B. (INPUT)
C                V      - NRV BY N MATRIX. (INPUT/OUTPUT)
C                           IF NRV.LE.0, V IS NOT USED. OTHERWISE,
C                           V IS REPLACED BY THE NRV BY N PRODUCT
C                           MATRIX V*VB. SEE REMARKS.
C                IV     - ROW DIMENSION OF MATRIX V EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NRV    - NUMBER OF ROWS OF V. (INPUT)
C                C      - N BY NCC MATRIX. (INPUT/OUTPUT)
C                           IF NCC.LE.0 C IS NOT USED. OTHERWISE, C
C                           IS REPLACED BY THE N BY NCC PRODUCT
C                           MATRIX UB**(T) * C. SEE REMARKS.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NCC    - NUMBER OF COLUMNS IN C. (INPUT)
C                IER    - ERROR PARAMETER. (INPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT MATRIX B IS NOT FULL
C                             RANK OR VERY ILL-CONDITIONED. SMALL
C                             SINGULAR VALUES MAY NOT BE VERY ACCURATE.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT CONVERGENCE WAS
C                             NOT ATTAINED AFTER 10*N QR SWEEPS.
C                             (CONVERGENCE USUALLY OCCURS IN ABOUT
C                             2*N SWEEPS).
C
C   REQD. IMSL ROUTINES - SINGLE/LSVG2,VHS12,UERTST,UGETIO,VBLA=SROTG
C                       - DOUBLE/LSVG2,VHS12,UERTST,UGETIO,VBLA=DROTG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      LSVDB COMPUTES THE SINGULAR VALUE DECOMPOSITION OF
C                AN N BY N BIDIAGONAL MATRIX
C                     B = UB * S * VB**(T)    WHERE
C                UB AND VB ARE N BY N ORTHOGONAL MATRICES AND
C                S IS DIAGONAL.
C                IF ARGUMENTS V AND C ARE N BY N IDENTITY MATRICES,
C                ON EXIT THEY ARE REPLACED BY VB AND UB**T,
C                RESPECTIVELY.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IV,NRV,IC,NCC,IER
      REAL               D(N),E(N),V(IV,1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,J,K,KK,L,LL,LP1,NQRS,N10
      LOGICAL            WNTV,HAVERS,FAIL
      REAL               DNORM,ZERO,ONE,TWO,CS,F,SQINF,FTEMP,G,H,HTEMP,
     *                   SN,T,X,Y,Z
      DATA               SQINF/0.8507057E38/
      DATA               ZERO/0.0/,ONE/1.0/,TWO/2.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N.LE.0) GO TO 9005
      N10 = 10*N
      WNTV = NRV.GT.0
      HAVERS = NCC.GT.0
      FAIL = .FALSE.
      NQRS = 0
      E(1) = ZERO
      DNORM = ZERO
      DO 5 J=1,N
    5 DNORM = AMAX1(ABS(D(J))+ABS(E(J)),DNORM)
      DO 100 KK=1,N
         K = N+1-KK
C                                  TEST FOR SPLITTING OR RANK
C                                    DEFICIENCIES FIRST MAKE TEST FOR
C                                    LAST DIAGONAL TERM, D(K), BEING
C                                    SMALL.
   10    IF (K.EQ.1) GO TO 25
         T = DNORM+D(K)
         IF (T.NE.DNORM) GO TO 25
C
C                                  SINCE D(K) IS SMALL WE WILL MAKE A
C                                    SPECIAL PASS TO TRANSFORM E(K) TO
C                                    ZERO.
         CS = ZERO
         SN = -ONE
         DO 20 II=2,K
            I = K+1-II
            F = -SN*E(I+1)
            E(I+1) = CS*E(I+1)
            T = D(I)
            FTEMP = F
            CALL SROTG (D(I),FTEMP,CS,SN)
C                                  TRANSFORMATION CONSTRUCTED TO ZERO
C                                    POSITION (I,K).
            IF (.NOT.WNTV) GO TO 20
            DO 15 J=1,NRV
   15       CALL LSVG2 (CS,SN,V(J,I),V(J,K))
C
C                                  ACCUMULATE RT. TRANSFORMATIONS IN V.
   20    CONTINUE
C                                  THE MATRIX IS NOW BIDIAGONAL, AND OF
C                                    LOWER ORDER SINCE E(K) .EQ. ZERO
   25    DO 30 LL=1,K
            L = K+1-LL
            T = DNORM+E(L)
            IF (T.EQ.DNORM) GO TO 50
            T = DNORM+D(L-1)
            IF (T.EQ.DNORM) GO TO 35
   30    CONTINUE
C                                  THIS LOOP CANT COMPLETE SINCE E(1) =
C                                    ZERO.
         GO TO 50
C                                  CANCELLATION OF E(L), L.GT.1.
   35    CS = ZERO
         SN = -ONE
         DO 45 I=L,K
            F = -SN*E(I)
            E(I) = CS*E(I)
            T = DNORM+F
            IF (T.EQ.DNORM) GO TO 50
            T = D(I)
            FTEMP = F
            CALL SROTG (D(I),FTEMP,CS,SN)
            IF (.NOT.HAVERS) GO TO 45
            DO 40 J=1,NCC
   40       CALL LSVG2 (CS,SN,C(I,J),C(L-1,J))
   45    CONTINUE
C                                  TEST FOR CONVERGENCE
   50    Z = D(K)
         IF (L.EQ.K) GO TO 85
C                                  SHIFT FROM BOTTOM 2 BY 2 MINOR OF
C                                    B**(T)*B.
         X = D(L)
         Y = D(K-1)
         G = E(K-1)
         H = E(K)
         F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
         G = ABS(F)
         IF (ABS(F) .LT. SQINF) G = SQRT(ONE+F**2)
         IF (F.LT.ZERO) GO TO 55
         T = F+G
         GO TO 60
   55    T = F-G
   60    F = ((X-Z)*(X+Z)+H*(Y/T-H))/X
C                                  NEXT QR SWEEP
         CS = ONE
         SN = ONE
         LP1 = L+1
         DO 80 I=LP1,K
            G = E(I)
            Y = D(I)
            H = SN*G
            G = CS*G
            HTEMP = H
            CALL SROTG (F,HTEMP,CS,SN)
            E(I-1) = F
            F = X*CS+G*SN
            G = -X*SN+G*CS
            H = Y*SN
            Y = Y*CS
            IF (.NOT.WNTV) GO TO 70
C                                  ACCUMULATE ROTATIONS (FROM THE
C                                    RIGHT) IN V
            DO 65 J=1,NRV
   65       CALL LSVG2 (CS,SN,V(J,I-1),V(J,I))
            HTEMP = H
   70       CALL SROTG (F,HTEMP,CS,SN)
            D(I-1) = F
            F = CS*G+SN*Y
            X = -SN*G+CS*Y
            IF (.NOT.HAVERS) GO TO 80
            DO 75 J=1,NCC
   75       CALL LSVG2 (CS,SN,C(I-1,J),C(I,J))
C
C                                  APPLY ROTATIONS FROM THE LEFT TO
C                                    RIGHT HAND SIDES IN C
   80    CONTINUE
         E(L) = ZERO
         E(K) = F
         D(K) = X
         NQRS = NQRS+1
         IF (NQRS.LE.N10) GO TO 10
C                                  RETURN TO TEST FOR SPLITTING.
         FAIL = .TRUE.
C                                  CUTOFF FOR CONVERGENCE FAILURE. NQRS
C                                    WILL BE 2*N USUALLY.
   85    IF (Z.GE.ZERO) GO TO 95
         D(K) = -Z
         IF (.NOT.WNTV) GO TO 95
         DO 90 J=1,NRV
   90    V(J,K) = -V(J,K)
   95    CONTINUE
C                                  CONVERGENCE. D(K) IS MADE
C                                    NONNEGATIVE
  100 CONTINUE
      IF (N.EQ.1) GO TO 140
      DO 105 I=2,N
         IF (D(I).GT.D(I-1)) GO TO 110
  105 CONTINUE
      GO TO 140
C                                  EVERY SINGULAR VALUE IS IN ORDER
  110 DO 135 I=2,N
         T = D(I-1)
         K = I-1
         DO 115 J=I,N
            IF (T.GE.D(J)) GO TO 115
            T = D(J)
            K = J
  115    CONTINUE
         IF (K.EQ.I-1) GO TO 135
         D(K) = D(I-1)
         D(I-1) = T
         IF (.NOT.HAVERS) GO TO 125
         DO 120 J=1,NCC
            T = C(I-1,J)
            C(I-1,J) = C(K,J)
  120    C(K,J) = T
  125    IF (.NOT.WNTV) GO TO 135
         DO 130 J=1,NRV
            T = V(J,I-1)
            V(J,I-1) = V(J,K)
  130    V(J,K) = T
  135 CONTINUE
C                                  END OF ORDERING ALGORITHM.
  140 IER = 129
      IF (FAIL) GO TO 9000
C                                  CHECK FOR POSSIBLE RANK DEFICIENCY
      IER = 33
      T = 0.0
      IF (D(1).NE.ZERO) T=D(N)/D(1)
      F=100.0+T
      IF (F.EQ.100.0) GO TO 9000
      IER = 0
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLSVDB )
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - VBLA=SROTG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - CONSTRUCT GIVENS PLANE ROTATION
C                           (SINGLE PRECISION)
C
C   USAGE               - CALL SROTG (SA,SB,SC,SS)
C
C   ARGUMENTS    SA     - FIRST ELEMENT OF VECTOR. (INPUT/OUTPUT)
C                         ON OUTPUT, R=(+/-)SQRT(SA**2 + SB**2)
C                           OVERWRITES SA.
C                SB     - SECOND ELEMENT OF VECTOR. (INPUT/OUTPUT)
C                         ON OUTPUT, Z OVERWRITES SB.
C                           Z IS DEFINED TO BE..
C                           SS,     IF ABS(SA).GT.ABS(SB)
C                           1.0/SC, IF ABS(SB).GE.ABS(SA) AND SC.NE.0.0
C                           1.,     IF SC.EQ.0.0.
C                SC     - ELEMENT OF OUTPUT TRANSFORMATION MATRIX.
C                           SEE REMARKS.
C                SS     - ELEMENT OF OUTPUT TRANSFORMATION MATRIX.
C                           SEE REMARKS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      SROTG CONSTRUCTS THE GIVENS TRANSFORMATION
C                    ( SC  SS )
C                G = (        ) ,  SC**2 + SS**2 = 1 ,
C                    (-SS  SC )
C                WHICH ZEROS THE SECOND ELEMENT OF (SA,SB)**T.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SROTG  (SA,SB,SC,SS)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               SA,SB,SC,SS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               R,U,V
C                                  FIRST EXECUTABLE STATEMENT
      IF (ABS(SA).LE.ABS(SB)) GO TO 5
C                                  HERE ABS(SA) .GT. ABS(SB)
      U = SA+SA
      V = SB/U
C                                  NOTE THAT U AND R HAVE THE SIGN OF
C                                    SA
      R = SQRT(.25+V**2)*U
C                                  NOTE THAT SC IS POSITIVE
      SC = SA/R
      SS = V*(SC+SC)
      SB = SS
      SA = R
      RETURN
C                                  HERE ABS(SA) .LE. ABS(SB)
    5 IF (SB.EQ.0.) GO TO 15
      U = SB+SB
      V = SA/U
C                                  NOTE THAT U AND R HAVE THE SIGN OF
C                                    SB (R IS IMMEDIATELY STORED IN SA)
      SA = SQRT(.25+V**2)*U
C                                  NOTE THAT SS IS POSITIVE
      SS = SB/SA
      SC = V*(SS+SS)
      IF (SC.EQ.0.) GO TO 10
      SB = 1./SC
      RETURN
   10 SB = 1.
      RETURN
C                                  HERE SA = SB = 0.
   15 SC = 1.
      SS = 0.
      SA = 0.
      SB = 0.
      RETURN
C
      END


C   IMSL ROUTINE NAME   - LSVG1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LSVG1  (A,B,CS,SN,SIG)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               A,B,CS,SN,SIG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               AA,BB
C                                  FIRST EXECUTABLE STATEMENT
      IF (ABS(A).LE.ABS(B)) GO TO 5
      AA = ABS(A+A)
      SIG = AA*SQRT(0.25+(B/AA)**2)
      CS = A/SIG
      SN = B/SIG
      RETURN
    5 IF (B.EQ.0.0) GO TO 10
      BB = ABS(B+B)
      SIG = BB*SQRT(0.25+(A/BB)**2)
      CS = A/SIG
      SN = B/SIG
      RETURN
   10 SIG = 0.0
      CS = 0.0
      SN = 1.0
      RETURN
      END


C   IMSL ROUTINE NAME   - LSVG2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LSVG2  (CS,SN,X,Y)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               CS,SN,X,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               XR
C                                  FIRST EXECUTABLE STATEMENT
      XR=CS*X+SN*Y
      Y=-SN*X+CS*Y
      X=XR
      RETURN
      END


C   IMSL ROUTINE NAME   - EIGRF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A REAL GENERAL MATRIX IN FULL STORAGE MODE
C
C   USAGE               - CALL EIGRF (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - THE INPUT REAL GENERAL MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                N      - THE INPUT ORDER OF THE MATRIX A.
C                IA     - THE INPUT ROW DIMENSION OF MATRIX A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IJOB   - THE INPUT OPTION PARAMETER. WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN WK(1). THE ROUTINES HAVE
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
C                           GREATER THAN 100).
C                W      - THE OUTPUT COMPLEX VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES OF A.
C                         NOTE - THE ROUTINE TREATS W AS A REAL VECTOR
C                           OF LENGTH 2*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                Z      - THE OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE W(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. IZ MUST BE GREATER
C                           THAN OR EQUAL TO N IF IJOB IS NOT EQUAL TO
C                           ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             (2+N)N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQRH3F FAILED
C                           TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                           J+1,J+2,...,N HAVE BEEN COMPUTED CORRECTLY.
C                           EIGENVALUES 1,...,J ARE SET TO ZERO.
C                           IF IJOB = 1 OR 2 EIGENVECTORS ARE SET TO
C                           ZERO. THE PERFORMANCE INDEX IS SET TO 1000.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS SET TO ZERO.
C
C   REQD. IMSL ROUTINES - EBALAF,EBBCKF,EHBCKF,EHESSF,EQRH3F,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGRF  (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IJOB,IZ,IER
      REAL               A(IA,1),WK(N,1),W(1),Z(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,IZ2,K,L,I,N1,N2,II,JJ,NP1,IIZ,NPI,JW,J,
     *                   IS,IG,IGZ,LW,LLZ,KKZ,LZ,KZ
      REAL               ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP,
     *                   ZERO,ONE,THOUS,AN,Z11
C      DATA               RDELP/2.3841857910156E-7/
      DATA               RDELP/8.78906E-03/
      DATA               ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IZ2 = IZ+IZ
      IF (IJOB .GE. 0 .AND. IJOB .LE. 3) GO TO 5
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 10
    5 IF (IJOB .EQ. 0) GO TO 16
   10 IF (IZ .GE. N) GO TO 15
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   15 IF (IJOB .EQ. 3) GO TO 95
C                                  PACK A INTO AN N BY N ARRAY
   16 K = 1
      L = 1
      DO 20 J=1,N
         DO 20 I=1,N
            A(K,L) = A(I,J)
C                                  SAVE INPUT A IF IJOB = 2
            IF (IJOB .EQ. 2) WK(I,J)=A(I,J)
            K = K+1
            IF (K .GT. IA) K = 1
            IF (K .EQ. 1) L = L+1
   20 CONTINUE
      N1 = 1
      IF (IJOB .EQ. 2) N1 = N+1
      N2 = N1+1
      IF (IJOB .EQ. 0) N2 = 1
C                                  BALANCE THE INPUT A
      CALL EBALAF (A,N,N,WK(1,N1),K,L)
      IF (IJOB .EQ. 0 .AND. L .EQ. 0) GO TO 35
C                                  IF L = 0, A IS ALREADY IN HESSENBERG
C                                    FORM
      CALL EHESSF (A,K,L,N,N,WK(1,N2))
      IF (IJOB .EQ. 0) GO TO 35
C                                  SET Z IDENTITY MATRIX
      II = 1
      JJ = 1
      NP1 = N+1
      DO 30 I=1,N
         DO 25 J=1,N
            Z(II) = ZERO
            II = II+1
   25    CONTINUE
         Z(JJ) = ONE
         JJ = JJ+NP1
   30 CONTINUE
      CALL EHBCKF (Z,A,WK(1,N2),N,N,N,K,L)
      IIZ = N
   35 IF (IJOB .EQ. 0) IIZ = 1
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z11 = Z(1)
      CALL EQRH3F (A,N,N,K,L,W(1),W(N+1),Z,IIZ,JER)
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z(1) = Z11
      IF (JER .GT. 128 .OR. IJOB .EQ. 0) GO TO 40
      CALL EBBCKF (WK(1,N1),Z,K,L,N,N,N)
C                                  CONVERT W (EIGENVALUES) TO COMPLEX
C                                    FORMAT
   40 DO 45 I=1,N
         NPI = N+I
         WK(I,N1) = W(NPI)
   45 CONTINUE
      JW = N+N
      J = N
      DO 50 I=1,N
         W(JW-1) = W(J)
         W(JW) = WK(J,N1)
         JW = JW-2
         J = J-1
   50 CONTINUE
      IF (IJOB .EQ. 0) GO TO 9000
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
      J = N
   60 IF (J .LT. 1) GO TO 85
      IF (W(J+J) .EQ. ZERO) GO TO 75
C                                  MOVE PAIR OF COMPLEX CONJUGATE
C                                    EIGENVECTORS
      IS = IZ2*(J-1)+1
      IG = N*(J-2)+1
      IGZ = IG+N
C                                  MOVE COMPLEX CONJUGATE EIGENVECTOR
      DO 65 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IGZ)
         IS = IS+2
         IG = IG+1
         IGZ = IGZ+1
   65 CONTINUE
C                                  MOVE COMPLEX EIGENVECTOR
      IS = IZ2*(J-2)+1
      IG = IS+IZ2
      DO 70 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IG+1)
         IS = IS+2
         IG = IG+2
   70 CONTINUE
      J = J-2
      GO TO 60
C                                  MOVE REAL EIGENVECTOR
   75 IS = IZ2*(J-1)+N+N
      IG = N*J
      DO 80 I=1,N
         Z(IS-1) = Z(IG)
         Z(IS) = ZERO
         IS = IS-2
         IG = IG-1
   80 CONTINUE
      J = J-1
      GO TO 60
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
C                                    NEXT, MOVE ORIGINAL MATRIX BACK
C                                    TO A
   85 IF (IJOB .LE. 1) GO TO 9000
      DO 90 I=1,N
         DO 90 J=1,N
            A(I,J) = WK(I,J)
   90 CONTINUE
      WK(1,1) = THOUS
      IF (JER .NE. 0) GO TO 9000
C                                  COMPUTE 1-NORM OF A
   95 ANORM = ZERO
      DO 105 J=1,N
         ASUM = ZERO
         DO 100 I=1,N
            ASUM = ASUM+ABS(A(I,J))
  100    CONTINUE
         ANORM = AMAX1(ANORM,ASUM)
  105 CONTINUE
      IF (ANORM .EQ. ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      LW = 1
      LLZ = 0
      KKZ = 0
      DO 120 J=1,N
         S = ZERO
         SUMZ = ZERO
         LZ = LLZ+1
         KZ = KKZ+1
         LW = J+J-1
         DO 115 L=1,N
            SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
            KZ = KKZ+1
            SUMR = -W(LW)*Z(LZ)+W(LW+1)*Z(LZ+1)
            SUMI = -W(LW)*Z(LZ+1)-W(LW+1)*Z(LZ)
            DO 110 K=1,N
               SUMR =SUMR+A(L,K)*Z(KZ)
               SUMI = SUMI+A(L,K)*Z(KZ+1)
               KZ = KZ+2
  110       CONTINUE
            S = S+CABS(CMPLX(SUMR,SUMI))
            LZ = LZ+2
  115    CONTINUE
         PI = AMAX1(PI,S/SUMZ)
         KKZ = KKZ+IZ2
         LLZ = LLZ+IZ2
  120 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1,1) = PI
 9000 CONTINUE
C      IF (IER .NE. 0) CALL UERTST (IER,6HEIGRF )
      IF (JER .EQ. 0) GO TO 9005
      IER = JER
C      CALL UERTST (IER,6HEIGRF )
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - EBALAF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBALAF (A,N,IA,D,K,L)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,K,L
      REAL               A(IA,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L1,K1,K1P1,K11,JJ,J,I,LL,NOCONV
      REAL               R,C,F,G,B,S,B2,ONE,ZERO,P95
      DATA               B/16.0/,B2/256.0/
      DATA               ZERO/0.0/,ONE/1.0/,P95/.95/
C                                  REDUCE NORM A BY DIAGONAL SIMILARITY
C                                  TRANSFORMATION STORED IN D
C                                  FIRST EXECUTABLE STATEMENT
      L1 = 1
      K1 = N
C                                  SEARCH FOR ROWS ISOLATING AN EIGEN-
C                                    VALUE AND PUSH THEM DOWN
    5 K1P1 = K1+1
      IF (K1.LT.1) GO TO 35
      K11=K1
      DO 30 JJ=1,K11
         J = K1P1-JJ
         R = ZERO
         DO 10 I=1,K1
            IF (I.EQ.J) GO TO 10
            R=R+ABS(A(J,I))
   10    CONTINUE
         IF (R.NE.ZERO) GO TO 30
         D(K1) = J
         IF (J.EQ.K1) GO TO 25
         DO 15 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,K1)
            A(I,K1) = F
   15    CONTINUE
         DO 20 I=L1,N
            F = A(J,I)
            A(J,I) = A(K1,I)
            A(K1,I) = F
   20    CONTINUE
   25    K1 = K1-1
         GO TO 5
   30 CONTINUE
C                                  SEARCH FOR COLUMNS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM LEFT
   35 IF (K1.LT.L1) GO TO 65
      LL = L1
      DO 60 J=LL,K1
         C = ZERO
         DO 40 I=L1,K1
            IF (I.EQ.J) GO TO 40
            C = C+ABS(A(I,J))
   40    CONTINUE
         IF (C.NE.ZERO) GO TO 60
         D(L1) = J
         IF (J.EQ.L1) GO TO 55
         DO 45 I=1,K1
            F = A(I,J)
            A(I,J) = A(I,L1)
            A(I,L1) = F
   45    CONTINUE
         DO 50  I=L1,N
            F = A(J,I)
            A(J,I) = A(L1,I)
            A(L1,I) = F
   50    CONTINUE
   55    L1 = L1+1
         GO TO 35
   60 CONTINUE
C                                  NOW BALANCE THE SUBMATRIX IN ROWS
C                                    L1 THROUGH K1
   65 K = L1
      L = K1
      IF (K1.LT.L1) GO TO 75
      DO 70  I=L1,K1
         D(I) = ONE
   70 CONTINUE
   75 NOCONV = 0
      IF (K1.LT.L1) GO TO 120
      DO 115 I=L1,K1
         C = ZERO
         R = ZERO
         DO 80 J=L1,K1
            IF (J.EQ.I) GO TO 80
            C = C+ABS(A(J,I))
            R = R+ABS(A(I,J))
   80    CONTINUE
         G = R/B
         F = ONE
         S = C+R
   85    IF (C.GE.G) GO TO 90
         F = F * B
         C = C*B2
         GO TO 85
   90    G = R*B
   95    IF (C.LT.G) GO TO 100
         F = F/B
         C = C/B2
         GO TO 95
C                                  NOW BALANCE
  100    IF ((C+R)/F.GE.P95*S) GO TO 115
         G = ONE/F
         D(I) = D(I)*F
         NOCONV = 1
         DO 105 J=L1,N
            A(I,J) = A(I,J)*G
  105    CONTINUE
         DO 110 J=1,K1
            A(J,I) = A(J,I)*F
  110    CONTINUE
  115 CONTINUE
  120 IF (NOCONV.EQ.1) GO TO 75
      RETURN
      END


C   IMSL ROUTINE NAME   - EHESSF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHESSF (A,K,L,N,IA,D)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,N,IA
      REAL               A(IA,N),D(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LA,KP1,M,I,MP,II,J,JJ
      REAL               F,G,H,SCALE,ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      LA = L - 1
      KP1 = K + 1
      IF (LA .LT. KP1) GO TO 50
      DO 45 M = KP1, LA
         H = ZERO
         D(M) = ZERO
         SCALE = ZERO
C                                  SCALE COLUMN
         DO 5 I = M, L
            SCALE = SCALE + ABS(A(I,M-1))
    5    CONTINUE
         IF (SCALE .EQ. ZERO ) GO TO 45
         MP = M + L
C                                  DO 10 I=L,M,-1
         DO 10 II = M, L
            I = MP - II
            D(I) = A(I,M-1) / SCALE
            H = H + D(I) * D(I)
   10    CONTINUE
         G = -SIGN(SQRT(H),D(M))
         H = H - D(M) * G
         D(M) = D(M) - G
C                                  FORM (I-(U*UT)/H) * A
         DO 25 J = M,N
            F = ZERO
C                                  DO 15 I=L,M,-1
            DO 15 II = M, L
               I = MP - II
               F = F + D(I) * A(I,J)
   15       CONTINUE
            F = F / H
            DO 20 I = M, L
               A(I,J) = A(I,J) - F * D(I)
   20       CONTINUE
   25    CONTINUE
C                                  FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
         DO 40 I = 1,L
            F = ZERO
C                                  DO 30 J=L,M,-1
            DO 30 JJ = M, L
               J = MP - JJ
               F = F + D(J) * A(I,J)
   30       CONTINUE
            F = F / H
            DO 35 J = M, L
               A(I,J) = A(I,J) - F * D(J)
   35       CONTINUE
   40    CONTINUE
         D(M) = SCALE * D(M)
         A(M,M-1) = SCALE * G
   45 CONTINUE
   50 RETURN
      END


C   IMSL ROUTINE NAME   - EHBCKF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHBCKF (Z,H,D,N,MM,IZH,K,L)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MM,IZH,K,L
      REAL               Z(IZH,1),H(IZH,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LM2,KI,LTEMP,M,MA,MP2,I,J
      REAL               G,ZERO
      DATA               ZERO/0.0/
C                                  ADAPTED FROM EISPACK ROUTINE ORTBAK
C                                  FIRST EXECUTABLE STATEMENT
      LM2=L-2
      IF(LM2.LT.K) GO TO 9005
      LTEMP=LM2+K
      DO 30 KI=K,LM2
         M=LTEMP-KI
         MA=M+1
         IF(H(MA,M).EQ.ZERO) GO TO 30
         MP2=M+2
         IF(MP2.GT.L) GO TO 10
         DO 5 I=MP2,L
            D(I)=H(I,M)
    5    CONTINUE
   10    IF(MA.GT.L) GO TO 30
         DO 25 J=1,MM
            G=ZERO
            DO 15 I=MA,L
               G=G+D(I)*Z(I,J)
   15       CONTINUE
C                                  DOUBLE DIVISION AVOIDS POSSIBLE
C                                  UNDERFLOW
            G = (G/D(MA))/H(MA,M)
            DO 20 I=MA,L
               Z(I,J)=Z(I,J)+G*D(I)
   20       CONTINUE
   25    CONTINUE
   30 CONTINUE
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - EQRH3F
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQRH3F (H,N,IH,K,L,WR,WI,Z,IZ,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IH,K,L,IZ,IER
      REAL               H(IH,N),WR(N),WI(N),Z(IZ,N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEN,ITS,IENM2,NPL,LL,LB,NAML,MM,M,MP2,KA,NA,
     *                   J,JJ
      REAL               T3(2),RDELP,P4,P5,P7,ZERO,ONE,T,X,Y,W,S,ZZ,R,P,
     *                   Q,RNORM,RA,SA,SCALE,VR,VI
      COMPLEX            Z3
      LOGICAL            NOTLAS
      EQUIVALENCE        (Z3,T3(1))
c      DATA               RDELP/2.3841857910156E-7/
      DATA               RDELP/8.78906E-03/

      DATA               P4 /0.4375/,P5 /0.5/,P7 /0.75/,ZERO /0.0/,ONE
     *                   /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  STORE ROOTS ISOLATED BY EBALAF
      RNORM = 0.0
      KA = 1
      DO 10 I=1,N
         DO 5 J=KA,N
    5    RNORM = RNORM+ABS(H(I,J))
         KA = I
         IF (I.GE.K .AND. I.LE.L) GO TO 10
         WR(I) = H(I,I)
         WI(I) = ZERO
   10 CONTINUE
      IEN = L
      T = ZERO
C                                  SEARCH FOR NEXT EIGENVALUES
   15 IF (IEN.LT.K) GO TO 145
      ITS = 0
      NA = IEN-1
      IENM2 = NA-1
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                  ELEMENT
   20 NPL = IEN+K
      DO 25 LL=K,IEN
         LB = NPL-LL
         IF (LB.EQ.K) GO TO 30
         S = ABS(H(LB-1,LB-1))+ABS(H(LB,LB))
         IF (S.EQ.0.0) S = RNORM
         IF (ABS(H(LB,LB-1)).LE.RDELP*S) GO TO 30
   25 CONTINUE
C
   30 X = H(IEN,IEN)
      IF (LB.EQ.IEN) GO TO 110
      Y = H(NA,NA)
      W = H(IEN,NA)*H(NA,IEN)
      IF (LB.EQ.NA) GO TO 115
      IF (ITS.EQ.30) GO TO 250
C                                  FORM SHIFT
      IF (ITS.NE.10 .AND. ITS.NE.20) GO TO 40
      T = T+X
      DO 35 I=K,IEN
         H(I,I) = H(I,I)-X
   35 CONTINUE
      S = ABS(H(IEN,NA))+ABS(H(NA,IENM2))
      X = P7*S
      Y = X
      W = -P4*S*S
   40 ITS = ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                  SUB-DIAGONAL ELEMENTS
      NAML = IENM2+LB
      DO 45 MM=LB,IENM2
         M = NAML-MM
         ZZ = H(M,M)
         R = X-ZZ
         S = Y-ZZ
         P = (R*S-W)/H(M+1,M)+H(M,M+1)
         Q = H(M+1,M+1)-ZZ-R-S
         R = H(M+2,M+1)
         S = ABS(P)+ABS(Q)+ABS(R)
         P = P/S
         Q = Q/S
         R = R/S
         IF (M.EQ.LB) GO TO 50
         IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)).LE.RDELP*ABS(P)*(ABS(H(M-1,
     *   M-1))+ABS(ZZ)+ABS(H(M+1,M+1)))) GO TO 50
   45 CONTINUE
   50 MP2 = M+2
      DO 55 I=MP2,IEN
         H(I,I-2) = ZERO
         IF (I.EQ.MP2) GO TO 55
         H(I,I-3) = ZERO
   55 CONTINUE
C                                  DOUBLE QR STEP INVOLVING ROWS
C                                  L TO EN AND COLUMNS M TO EN
      DO 105 KA=M,NA
         NOTLAS = KA.NE.NA
         IF (KA.EQ.M) GO TO 60
         P = H(KA,KA-1)
         Q = H(KA+1,KA-1)
         R = ZERO
         IF (NOTLAS) R = H(KA+2,KA-1)
         X = ABS(P)+ABS(Q)+ABS(R)
         IF (X.EQ.ZERO) GO TO 105
         P = P/X
         Q = Q/X
         R = R/X
   60    CONTINUE
         S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (KA.EQ.M) GO TO 65
         H(KA,KA-1) = -S*X
         GO TO 70
   65    IF (LB.NE.M) H(KA,KA-1) = -H(KA,KA-1)
   70    P = P+S
         X = P/S
         Y = Q/S
         ZZ = R/S
         Q = Q/P
         R = R/P
C                                  ROW MODIFICATION
         DO 80 J=KA,N
            P = H(KA,J)+Q*H(KA+1,J)
            IF (.NOT.NOTLAS) GO TO 75
            P = P+R*H(KA+2,J)
            H(KA+2,J) = H(KA+2,J)-P*ZZ
   75       H(KA+1,J) = H(KA+1,J)-P*Y
            H(KA,J) = H(KA,J)-P*X
   80    CONTINUE
         J = MIN0(IEN,KA+3)
C                                  COLUMN MODIFICATION
         DO 90 I=1,J
            P = X*H(I,KA)+Y*H(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 85
            P = P+ZZ*H(I,KA+2)
            H(I,KA+2) = H(I,KA+2)-P*R
   85       H(I,KA+1) = H(I,KA+1)-P*Q
            H(I,KA) = H(I,KA)-P
   90    CONTINUE
         IF (IZ.LT.N) GO TO 105
C                                  ACCUMULATE TRANSFORMATIONS
         DO 100 I=K,L
            P = X*Z(I,KA)+Y*Z(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 95
            P = P+ZZ*Z(I,KA+2)
            Z(I,KA+2) = Z(I,KA+2)-P*R
   95       Z(I,KA+1) = Z(I,KA+1)-P*Q
            Z(I,KA) = Z(I,KA)-P
  100    CONTINUE
  105 CONTINUE
      GO TO 20
C                                  ONE ROOT FOUND
  110 H(IEN,IEN) = X+T
      WR(IEN) = H(IEN,IEN)
      WI(IEN) = ZERO
      IEN = NA
      GO TO 15
C                                  TWO ROOTS FOUND
  115 P = (Y-X)*P5
      Q = P*P+W
      ZZ = SQRT(ABS(Q))
      H(IEN,IEN) = X+T
      X = H(IEN,IEN)
      H(NA,NA) = Y+T
      IF (Q.LT.ZERO) GO TO 135
C                                  REAL PAIR
      ZZ = P+SIGN(ZZ,P)
      WR(NA) = X+ZZ
      WR(IEN) = WR(NA)
      IF (ZZ.NE.ZERO) WR(IEN) = X-W/ZZ
      WI(NA) = ZERO
      WI(IEN) = ZERO
      X = H(IEN,NA)
C                                  EMPLOY SCALE FACTOR IN CASE X AND
C                                  ZZ ARE VERY SMALL
      SCALE = ABS(X) + ABS(ZZ)
      R = SCALE * SQRT( (X/SCALE)**2 + (ZZ/SCALE)**2 )
      P = X/R
      Q = ZZ/R
C                                  ROW MODIFICATION
      DO 120 J=NA,N
         ZZ = H(NA,J)
         H(NA,J) = Q*ZZ+P*H(IEN,J)
         H(IEN,J) = Q*H(IEN,J)-P*ZZ
  120 CONTINUE
C                                  COLUMN MODIFICATION
      DO 125 I=1,IEN
         ZZ = H(I,NA)
         H(I,NA) = Q*ZZ+P*H(I,IEN)
         H(I,IEN) = Q*H(I,IEN)-P*ZZ
  125 CONTINUE
      IF (IZ.LT.N) GO TO 140
C                                  ACCUMULATE TRANSFORMATIONS
      DO 130 I=K,L
         ZZ = Z(I,NA)
         Z(I,NA) = Q*ZZ+P*Z(I,IEN)
         Z(I,IEN) = Q*Z(I,IEN)-P*ZZ
  130 CONTINUE
      GO TO 140
C                                  COMPLEX PAIR
  135 WR(NA) = X+P
      WR(IEN) = X+P
      WI(NA) = ZZ
      WI(IEN) = -ZZ
  140 IEN = IENM2
      GO TO 15
C                                  ALL ROOTS FOUND, NOW
C                                  BACKSUBSTITUTE
  145 IF (IZ.LT.N) GO TO 9005
      IF (RNORM.EQ.ZERO) GO TO 9005
      DO 220 NN=1,N
         IEN = N+1-NN
         P = WR(IEN)
         Q = WI(IEN)
         NA = IEN-1
         IF (Q.GT.ZERO) GO TO 220
         IF (Q.LT.ZERO) GO TO 180
C                                  REAL VECTOR
         M = IEN
         H(IEN,IEN) = ONE
         IF (NA.EQ.0) GO TO 220
         DO 175 II=1,NA
            I = IEN-II
            W = H(I,I)-P
            R = H(I,IEN)
            IF (M.GT.NA) GO TO 155
            DO 150 J=M,NA
               R = R+H(I,J)*H(J,IEN)
  150       CONTINUE
  155       IF (WI(I).GE.ZERO) GO TO 160
            ZZ = W
            S = R
            GO TO 175
  160       M = I
            IF (WI(I).NE.ZERO) GO TO 165
            T = W
            IF (W.EQ.ZERO) T = RDELP*RNORM
            H(I,IEN) = -R/T
            GO TO 175
C                                  SOLVE REAL EQUATIONS
  165       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
            T = (X*S-ZZ*R)/Q
            H(I,IEN) = T
            IF (ABS(X).LE.ABS(ZZ)) GO TO 170
            H(I+1,IEN) = (-R-W*T)/X
            GO TO 175
  170       H(I+1,IEN) = (-S-Y*T)/ZZ
  175    CONTINUE
C                                  END REAL VECTOR
         GO TO 220
C                                  LAST VECTOR COMPONENT CHOSEN
C                                    IMAGINARY SO THAT EIGENVECTOR
C                                    MATRIX IS TRIANGULAR
  180    M = NA
C                                  COMPLEX VECTOR
         IF (ABS(H(IEN,NA)).LE.ABS(H(NA,IEN))) GO TO 185
         H(NA,NA) = Q/H(IEN,NA)
         H(NA,IEN) = -(H(IEN,IEN)-P)/H(IEN,NA)
         GO TO 190
  185    CONTINUE
         Z3 = CMPLX(ZERO,-H(NA,IEN))/CMPLX(H(NA,NA)-P,Q)
         H(NA,NA) = T3(1)
         H(NA,IEN) = T3(2)
  190    H(IEN,NA) = ZERO
         H(IEN,IEN) = ONE
         IENM2 = NA-1
         IF (IENM2.EQ.0) GO TO 220
         DO 215 II=1,IENM2
            I = NA-II
            W = H(I,I)-P
            RA = ZERO
            SA = H(I,IEN)
            DO 195 J=M,NA
               RA = RA+H(I,J)*H(J,NA)
               SA = SA+H(I,J)*H(J,IEN)
  195       CONTINUE
            IF (WI(I).GE.ZERO) GO TO 200
            ZZ = W
            R = RA
            S = SA
            GO TO 215
  200       M = I
            IF (WI(I).NE.ZERO) GO TO 205
            Z3 = CMPLX(-RA,-SA)/CMPLX(W,Q)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            GO TO 215
C                                  SOLVE COMPLEX EQUATIONS
  205       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
            VI = (WR(I)-P)*Q
            VI = VI+VI
            IF (VR.EQ.ZERO .AND. VI.EQ.ZERO) VR = RDELP*RNORM*(ABS(W)
     *      +ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
            Z3 = CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/CMPLX(VR,VI)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            IF (ABS(X).LE.ABS(ZZ)+ABS(Q)) GO TO 210
            H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,IEN))/X
            H(I+1,IEN) = (-SA-W*H(I,IEN)-Q*H(I,NA))/X
            GO TO 215
  210       CONTINUE
            Z3 = CMPLX(-R-Y*H(I,NA),-S-Y*H(I,IEN))/CMPLX(ZZ,Q)
            H(I+1,NA) = T3(1)
            H(I+1,IEN) = T3(2)
  215    CONTINUE
C                                  END COMPLEX VECTOR
  220 CONTINUE
C                                  END BACKSUBSTITUTION
C                                  VECTORS OF ISOLATED ROOTS
      DO 230 I=1,N
         IF (I.GE.K .AND. I.LE.L) GO TO 230
         DO 225 J=I,N
            Z(I,J) = H(I,J)
  225    CONTINUE
  230 CONTINUE
      IF (L.EQ.0) GO TO 9005
C                                  MULTIPLY BY TRANSFORMATION MATRIX
      DO 245 JJ=K,N
         J = N+K-JJ
         M = MIN0(J,L)
         DO 240 I=K,L
            ZZ = ZERO
            DO 235 KA=K,M
               ZZ = ZZ+Z(I,KA)*H(KA,J)
  235       CONTINUE
            Z(I,J) = ZZ
  240    CONTINUE
  245 CONTINUE
      GO TO 9005
C                                  NO CONVERGENCE AFTER 30 ITERATIONS
C                                  SET ERROR INDICATOR  TO THE INDEX
C                                  OF THE CURRENT EIGENVALUE
  250 IER = 128+IEN
      DO 255 I=1,IEN
         WR(I) = ZERO
         WI(I) = ZERO
  255 CONTINUE
      IF (IZ.LT.N) GO TO 9000
      DO 265 I=1,N
         DO 260 J=1,N
            Z(I,J) = ZERO
  260    CONTINUE
  265 CONTINUE
 9000 CONTINUE
C      CALL UERTST (IER,6HEQRH3F)
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - EBBCKF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBBCKF (D,Z,K,L,MM,N,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,MM,N,IZ
      REAL               Z(IZ,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,KM1,II,JJ,LP1
      REAL               S
C                                  COLUMN SCALE Z BY APPROPRIATE D VALUE
C                                  FIRST EXECUTABLE STATEMENT
      IF (L.EQ.0) GO TO 15
      DO 10 I=K,L
         S = D(I)
         DO 5 J=1,MM
            Z(I,J) = Z(I,J)*S
    5    CONTINUE
   10 CONTINUE
C                                  INTERCHANGE ROWS IF PERMUTATIONS
C                                    OCCURRED IN EBALAF
   15 IF (K.EQ.1) GO TO 30
      KM1 = K-1
      DO 25 I=1,KM1
         II = K-I
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 25
         DO 20 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   20    CONTINUE
   25 CONTINUE
   30 IF (L.EQ.N) GO TO 45
      LP1 = L+1
      DO 40 II=LP1,N
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 40
         DO 35 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   35    CONTINUE
   40 CONTINUE
   45 RETURN
      END


C   IMSL ROUTINE NAME   - EIGCH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A COMPLEX HERMITIAN MATRIX
C
C   USAGE               - CALL EIGCH (A,N,JOBN,D,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - INPUT COMPLEX HERMITIAN MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                         NOTE - THE ROUTINE TREATS A AS A REAL VECTOR.
C                           AN EQUIVALENCE STATEMENT MAY BE REQUIRED-
C                           SEE DOCUMENT EXAMPLE.
C                N      - INPUT ORDER OF THE MATRIX A AND MATRIX Z.
C                JOBN   - INPUT OPTION PARAMETER. IF JOBN.GE.10
C                         A IS ASSUMED TO BE IN FULL COMPLEX STORAGE
C                         MODE (MUST BE DIMENSIONED EXACTLY N BY N).
C                         IF JOBN.LT.10 THEN A IS ASSUMED TO BE IN
C                         HERMITIAN STORAGE MODE.  DEFINE
C                         IJOB=MOD(JOBN,10).  THEN WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY.
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN WK(1). THE ROUTINES HAVE
C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
C                           GREATER THAN 100).
C                D      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           EIGENVALUES OF A.
C                Z      - OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE D(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IZ MUST BE GREATER THAN
C                           OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST 3N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 3N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             N*N+4N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQRT2S
C                             FAILED TO CONVERGE ON EIGENVALUE J.
C                             EIGENVALUES J+1,J+2,...,N HAVE BEEN
C                             COMPUTED CORRECTLY.
C                           THE PERFORMANCE INDEX IS SET TO 1000.0.
C                         WARNING ERROR (WITH FIX)
C                         IN THE FOLLOWING, IJOB = MOD(JOBN,10).
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB IS SET TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS SET TO ZERO.
C                           IER = 68, INDICATES THAT MATRIX A IS NOT
C                             HERMITIAN BECAUSE SOME DIAGONAL ELEMENT(S)
C                             ARE NOT REAL. EIGCH SETS THE IMAGINARY
C                             PART OF THESE ELEMENTS TO ZERO AND
C                             PROCEEDS WITH THE COMPUTATIONS.
C
C   REQD. IMSL ROUTINES - EHBCKH,EHOUSH,EQRT2S,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGCH  (A,N,JOBN,D,Z,IZ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,JOBN,IZ,IER
      REAL               A(1),D(N),Z(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,K,I,NE,NTAU,NA,NI,NI2,IM1,J,IIZ,NZ,IIZ1,
     1                   IJOB,JR,IR,IJ,JI,NP1,
     2                   JZ,JZI,L,M,II,IL,KK,LZ,MZ,LK,KZ
      REAL               ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP,
     1                   ZERO,ONE,THOUS,AN,SIGNA
      DATA               RDELP/2.3841857910156E-7/
      DATA               ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IF (JOBN.LT.10) GO TO 15
C                                  CONVERT TO HERMETIAN STORAGE MODE
      JR = N + N - 2
      IJ = 2
      K = 2
      DO 10 J=1,N
         DO 5 I=1,J
            A(K-1) = A(IJ-1)
            A(K) = -A(IJ)
            K = K+2
            IJ = IJ + 2
    5    CONTINUE
         IJ = IJ + JR
         JR = JR - 2
   10 CONTINUE
   15 IJOB = MOD(JOBN,10)
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 20
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 25
   20 IF (IJOB.EQ.0) GO TO 45
   25 IF (IZ.GE.N) GO TO 30
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   30 K = 2
      DO 40 I=1,N
         IF (A(K).EQ.ZERO) GO TO 35
         A(K) = ZERO
C                                  WARNING ERROR - SOME DIAGONAL
C                                    ELEMENT(S) NOT REAL
         IER = 68
   35    K = K+I+I+2
   40 CONTINUE
      IF (IJOB.EQ.3) GO TO 110
   45 NE = 1
      NTAU = NE+N
      NA = NTAU+N+N
      NI = (N*(N+1))/2
      NI2 = NI+NI
      IF (IJOB.NE.2) GO TO 55
C                                  SAVE INPUT A IF IJOB = 2
      K = NA
      DO 50 I=1,NI2
         WK(K) = A(I)
         K = K+1
   50 CONTINUE
C                                  SEPARATE A INTO REAL AND IMAGINARY
C                                    PARTS
   55 IF (NI.LT.2) GO TO 70
      IM1 = 1
      DO 65 I=2,NI
         K = IM1+I
         PI = A(K)
         DO 60 J=1,IM1
            A(K) = A(K-1)
            K = K-1
   60    CONTINUE
         A(I) = PI
         IM1 = I
   65 CONTINUE
C                                  REDUCE HERMITIAN MATRIX TO A REAL
C                                    SYMMETRIC TRIDIAGONAL MATRIX
   70 CALL EHOUSH (A(1),A(NI+1),N,D,WK(NE),WK(NTAU))
      IIZ = 1
      IF (IJOB.NE.0) IIZ = IZ+IZ
      IF (IIZ.EQ.1) GO TO 85
C                                  SET Z TO AN IDENTITY MATRIX
      NZ = (IZ+IZ)*N
      DO 75 I=1,NZ
         Z(I) = ZERO
   75 CONTINUE
      K = 1
      IIZ1 = IIZ+1
      DO 80 I=1,N
         Z(K) = ONE
         K = K+IIZ1
   80 CONTINUE
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS
   85 CALL EQRT2S (D,WK(NE),N,Z(1),IIZ,JER)
      IF (IJOB.EQ.0) GO TO 9000
C                                  BACKTRANSFORM THE EIGENVECTORS
      CALL EHBCKH (A(1),A(NI+1),WK(NTAU),N,Z(1),Z(IZ+1),IIZ)
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
      JZ = 0
      DO 100 J=1,N
         JZI = JZ+IZ
         DO 90 I=1,N
            K = JZI+I
            WK(I) = Z(K)
   90    CONTINUE
         K = JZ+N
         L = K+N-1
         M = N
         DO 95 I=1,N
            Z(L) = Z(K)
            Z(L+1) = WK(M)
            K = K-1
            L = L-2
            M = M-1
   95    CONTINUE
         JZ = JZ+IZ+IZ
  100 CONTINUE
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
      IF (IJOB.NE.2) GO TO 9000
C                                  MOVE ORIGINAL MATRIX BACK TO A
      K = NA
      DO 105 I=1,NI2
         A(I) = WK(K)
         K = K+1
  105 CONTINUE
      WK(1) = THOUS
      IF (JER.NE.0) GO TO 9000
C                                  COMPUTE 1-NORM OF A
  110 ANORM = ZERO
      II = 1
      DO 120 I=1,N
         ASUM = ZERO
         IL = II
         KK = 2
         DO 115 L=1,N
            ASUM = ASUM+CABS(CMPLX(A(IL),A(IL+1)))
            IF (L.GE.I) KK = L+L
            IL = IL+KK
  115    CONTINUE
         ANORM = AMAX1(ANORM,ASUM)
         II = II+I+I
  120 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 135 I=1,N
         II = 1
         S = ZERO
         SUMZ = ZERO
         LZ = (IZ+IZ)*(I-1)+1
         LZ = IZ*(I-1)*2+1
         MZ = LZ
         DO 130 L=1,N
            LK = II
            KK = 2
            SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
            SUMR = -D(I)*Z(LZ)
            SUMI = -D(I)*Z(LZ+1)
            KZ = MZ
            DO 125 K=1,N
               SIGNA = ONE
               IF (K.GT.L) SIGNA = -ONE
               SUMR = SUMR+A(LK)*Z(KZ)-SIGNA*A(LK+1)*Z(KZ+1)
               SUMI = SUMI+A(LK)*Z(KZ+1)+SIGNA*A(LK+1)*Z(KZ)
               IF (K.GE.L) KK = K+K
               LK = LK+KK
               KZ = KZ+2
  125       CONTINUE
            S = S+CABS(CMPLX(SUMR,SUMI))
            LZ = LZ+2
            II = II+L+L
  130    CONTINUE
         IF (SUMZ.EQ.ZERO) GO TO 135
         PI = AMAX1(PI,S/SUMZ)
  135 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
      IF (JOBN.LT.10) GO TO 9000
C                                  CONVERT BACK TO FULL COMPLEX MODE
      NP1 = N + 1
      IJ = (N-1) * NP1
      IJ = IJ + IJ + 2
      K = N * NP1
      DO 145 JR=1,N
         J = N+1-JR
         DO 140 IR=1,J
            A(IJ-1) = A(K-1)
            A(IJ) = -A(K)
            K = K-2
            IJ = IJ - 2
  140    CONTINUE
         IJ = IJ - JR - JR
  145 CONTINUE
      JR = N + N
      II = 2
      JI = 2
      DO 155 I=1,N
         IJ = II
         DO 150 J=1,I
            A(IJ-1) = A(JI-1)
            A(IJ) = -A(JI)
            JI = JI+2
            IJ = IJ+JR
  150    CONTINUE
         JI = JI + JR - I - I
         II = II + 2
  155 CONTINUE
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST (IER,6HEIGCH )
      IF (JER.EQ.0) GO TO 9005
      IER = JER
      CALL UERTST (IER,6HEIGCH )
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - EHOUSH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCH
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHOUSH (AR,AI,N,D,E,TAU)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               AR(1),AI(1),D(1),E(1),TAU(2,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NM1,NN,I,NR,NRM1,L,INDX,J,JJ,INX1,INX2,JP1,KK,
     *                   IX,IM1
      REAL               RHO,TOLER,ZERO,ONE,T1,T2,TESTBB,VR,ROOT,DELTA,
     *                   RATIO,RDELP,Q1,Q2,X1,X2,TT1,TT2,BB
      DATA               ZERO/0.0/,ONE/1.0/
      DATA               RDELP/2.3841857910156E-7/
C                                  FIRST EXECUTABLE STATEMENT
      NM1=N-1
      TOLER=ZERO
      NN=(N*(N+1))/2
      DO 5 I=1,NN
         T1=ABS(AR(I))
         T2=ABS(AI(I))
         IF(T2.GT.T1) T1=T2
         IF (T1.GT.TOLER) TOLER=T1
    5 CONTINUE
      TESTBB=RDELP*TOLER
      IF (N.LE.2) GO TO 65
C                                  PERFORM N - 2 SIMILARITY
C                                    TRANSFORMATIONS
      DO 60 NR=2,NM1
         NRM1=NR-1
         VR=ZERO
         TAU(1,NR)=ZERO
         TAU(2,NR)=ZERO
         TAU(2,1)=ZERO
         DO 10 L=NR,N
            INDX=(L*(L-1))/2+NRM1
            VR=AR(INDX)**2+AI(INDX)**2+VR
   10    CONTINUE
         INDX=(NR*NRM1)/2+NRM1
         IF ((TESTBB)**2 .GE. VR) GO TO 60
         ROOT = CABS(CMPLX(AR(INDX),AI(INDX)))*SQRT(VR)
         IF(ROOT.NE.ZERO) GO TO 15
         AR(INDX)=SQRT(VR)
         DELTA=VR
         TAU(1,1)=-AR(INDX)
         GO TO 20
   15    DELTA=VR+ROOT
         RATIO=VR/ROOT
         TAU(1,1)=-RATIO*AR(INDX)
         TAU(2,1)= RATIO*AI(INDX)
         AR(INDX)=(RATIO+ONE)*AR(INDX)
         AI(INDX)=(RATIO+ONE)*AI(INDX)
C                                  THE MATRIX TO BE USED IN THE
C                                    SIMILARITY TRANSFORMATION HAS
C                                    BEEN DETERMINED. THE TRANSFOR-
C                                    MATION FOLLOWS
   20    DO 35 J=NR,N
            JJ=(J*(J-1))/2
            INDX=JJ+NRM1
            TAU(1,J)=AR(INDX)/DELTA
            TAU(2,J)=AI(INDX)/DELTA
            D(J)=ZERO
            E(J)=ZERO
            DO 25 L=NR,J
               INX1=(L*(L-1))/2+NRM1
               INX2=JJ+L
               D(J)= D(J)+AR(INX2)*AR(INX1)-AI(INX2)*AI(INX1)
               E(J)= E(J)+AR(INX2)*AI(INX1)+AI(INX2)*AR(INX1)
   25       CONTINUE
            JP1=J+1
            IF (JP1 .GT. N) GO TO 40
            DO 30 L=JP1,N
               KK=(L*(L-1))/2
               INX1=KK+NRM1
               INX2=KK+J
               D(J)=D(J)+AR(INX2)*AR(INX1)+AI(INX2)*AI(INX1)
               E(J)=E(J)+AR(INX2)*AI(INX1)-AI(INX2)*AR(INX1)
   30       CONTINUE
   35    CONTINUE
   40    RHO=ZERO
         DO 45 L=NR,N
            RHO=RHO+D(L)*TAU(1,L)+E(L)*TAU(2,L)
   45    CONTINUE
         IX=(NRM1*(NR-2))/2
         DO 55 I=NR,N
            IX=IX+I-1
            INX2=IX+NRM1
            DO 50 J=NR,I
               INX1=IX+J
               X1=TAU(1,I)*D(J)+TAU(2,I)*E(J)
               X2=TAU(2,I)*D(J)-TAU(1,I)*E(J)
               Q1=D(I)-RHO*AR(INX2)
               Q2=E(I)-RHO*AI(INX2)
               T1=Q1*TAU(1,J)+Q2*TAU(2,J)
               T2=Q2*TAU(1,J)-Q1*TAU(2,J)
               AR(INX1)=AR(INX1)-X1-T1
               AI(INX1)=AI(INX1)-X2-T2
   50       CONTINUE
   55    CONTINUE
         TAU(1,NR)=TAU(1,1)
         TAU(2,NR)=TAU(2,1)
   60 CONTINUE
C                                  THE MATRIX HAS BEEN REDUCED TO TRI-
C                                    DIAGONAL HERMITIAN FORM. THE SUB-
C                                    DIAGONAL HAS BEEN TEMPORARILY
C                                    STORED IN VECTOR TAU. STORE THE
C                                    DIAGONAL OF THE REDUCED MATRIX IN D
   65 INDX=0
      DO 70 I=1,N
         INDX=INDX+I
         D(I)=AR(INDX)
   70 CONTINUE
C                                  PERFORM THE DIAGONAL UNITARY SIMILA-
C                                    RITY TRANSFORMATION
      TAU(1,1)=ONE
      TAU(2,1)=ZERO
      E(1)=ZERO
      IF (N .EQ. 1) GO TO 85
      INDX=(N*NM1)/2+NM1
      TAU(1,N)=AR(INDX)
      TAU(2,N)=-AI(INDX)
C                                  CALCULATE SUBDIAGONAL E OF THE REAL
C                                    SYMMETRIC TRIDIAGONAL MATRIX. CAL-
C                                    CULATE TAU, THE DIAGONAL OF THE
C                                    DIAGONAL UNITARY MATRIX
      INDX=1
      DO 80 I=2,N
         INDX=INDX+I
         IM1=I-1
         BB= SQRT(TAU(1,I)**2+TAU(2,I)**2)
         E(I)=BB
         AI(INDX)=BB
         IF (TESTBB .LT. BB) GO TO 75
         TAU(1,I)=ONE
         TAU(2,I)=ZERO
         BB=ONE
   75    TT1=TAU(1,I)*TAU(1,IM1)-TAU(2,I)*TAU(2,IM1)
         TT2=TAU(1,I)*TAU(2,IM1)+TAU(2,I)*TAU(1,IM1)
         TAU(1,I)=TT1/BB
         TAU(2,I)=TT2/BB
   80 CONTINUE
   85 RETURN
      END


C   IMSL ROUTINE NAME   - EQRT2S
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A SYMMETRIC TRIDIAGONAL MATRIX USING THE
C                           QL METHOD.
C
C   USAGE               - CALL EQRT2S (D,E,N,Z,IZ,IER)
C
C   ARGUMENTS    D      - ON INPUT, THE VECTOR D OF LENGTH N CONTAINS
C                           THE DIAGONAL ELEMENTS OF THE SYMMETRIC
C                           TRIDIAGONAL MATRIX T.
C                           ON OUTPUT, D CONTAINS THE EIGENVALUES OF
C                           T IN ASCENDING ORDER.
C                E      - ON INPUT, THE VECTOR E OF LENGTH N CONTAINS
C                           THE SUB-DIAGONAL ELEMENTS OF T IN POSITION
C                           2,...,N. ON OUTPUT, E IS DESTROYED.
C                N      - ORDER OF TRIDIAGONAL MATRIX T.(INPUT)
C                Z      - ON INPUT, Z CONTAINS THE IDENTITY MATRIX OF
C                           ORDER N.
C                           ON OUTPUT, Z CONTAINS THE EIGENVECTORS
C                           OF T. THE EIGENVECTOR IN COLUMN J OF Z
C                           CORRESPONDS TO THE EIGENVALUE D(J).
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IF IZ IS LESS THAN N, THE
C                           EIGENVECTORS ARE NOT COMPUTED. IN THIS CASE
C                           Z IS NOT USED.
C                IER    - ERROR PARAMETER
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQRT2S FAILED
C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                             AND EIGENVECTORS 1,...,J-1 HAVE BEEN
C                             COMPUTED CORRECTLY, BUT THE EIGENVALUES
C                             ARE UNORDERED.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQRT2S (D,E,N,Z,IZ,IER)
C
      DIMENSION          D(1),E(1),Z(IZ,1)
      DATA               RDELP/2.3841857910156E-7/
      DATA               ZERO,ONE/0.0,1.0/
C                                  MOVE THE LAST N-1 ELEMENTS
C                                  OF E INTO THE FIRST N-1 LOCATIONS
C                                  FIRST EXECUTABLE STATEMENT
      IER  = 0
      IF (N .EQ. 1) GO TO 9005
      DO 5  I=2,N
         E(I-1) = E(I)
    5 CONTINUE
      E(N) = ZERO
      B = ZERO
      F = ZERO
      DO  60  L=1,N
         J = 0
         H = RDELP*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C                                  LOOK FOR SMALL SUB-DIAGONAL ELEMENT
         DO 10  M=L,N
            K=M
            IF (ABS(E(K)) .LE. B) GO TO 15
   10    CONTINUE
   15    M = K
         IF (M.EQ.L) GO TO 55
   20    IF (J .EQ. 30) GO TO 85
         J = J+1
         L1 = L+1
         G = D(L)
         P = (D(L1)-G)/(E(L)+E(L))
         R = ABS(P)
         IF (RDELP*ABS(P) .LT. 1.0) R = SQRT(P*P+ONE)
         D(L) = E(L)/(P+SIGN(R,P))
         H = G-D(L)
         DO 25 I = L1,N
            D(I) = D(I)-H
   25    CONTINUE
         F = F+H
C                                  QL TRANSFORMATION
         P = D(M)
         C = ONE
         S = ZERO
         MM1 = M-1
         MM1PL = MM1+L
         IF (L.GT.MM1) GO TO 50
         DO 45 II=L,MM1
            I = MM1PL-II
            G = C*E(I)
            H = C*P
            IF (ABS(P).LT.ABS(E(I))) GO TO 30
            C = E(I)/P
            R = SQRT(C*C+ONE)
            E(I+1) = S*P*R
            S = C/R
            C = ONE/R
            GO TO 35
   30       C = P/E(I)
            R = SQRT(C*C+ONE)
            E(I+1) = S*E(I)*R
            S = ONE/R
            C = C*S
   35       P = C*D(I)-S*G
            D(I+1) = H+S*(C*G+S*D(I))
            IF (IZ .LT. N) GO TO 45
C                                  FORM VECTOR
            DO 40 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I)+C*H
               Z(K,I) = C*Z(K,I)-S*H
   40       CONTINUE
   45    CONTINUE
   50    E(L) = S*P
         D(L) = C*P
         IF ( ABS(E(L)) .GT.B) GO TO 20
   55    D(L) = D(L) + F
   60 CONTINUE
C                                  ORDER EIGENVALUES AND EIGENVECTORS
      DO  80  I=1,N
         K = I
         P = D(I)
         IP1 = I+1
         IF (IP1.GT.N) GO TO 70
         DO 65  J=IP1,N
            IF (D(J) .GE. P) GO TO 65
            K = J
            P = D(J)
   65    CONTINUE
   70    IF (K.EQ.I) GO TO 80
         D(K) = D(I)
         D(I) = P
         IF (IZ .LT. N) GO TO 80
         DO 75 J = 1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
   75    CONTINUE
   80 CONTINUE
      GO TO 9005
   85 IER = 128+L
 9000 CONTINUE
      CALL UERTST(IER,6HEQRT2S)
 9005 RETURN
      END


C   IMSL ROUTINE NAME   - EHBCKH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCH
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHBCKH (AR,AI,TAU,N,ZR,ZI,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IZ
      REAL               AR(1),AI(1),TAU(2,1),ZR(IZ,1),ZI(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,K,NR,L,NRM1,INX1,INX2,K1
      REAL               DELTA,ZERO,ALPHA1,ALPHA2
      DATA               ZERO/0.0/
C                                  TRANSFORM THE EIGENVECTORS OF THE
C                                    REAL SYMMETRIC TRIDIAGONAL MATRIX
C                                    TO THOSE OF THE HERMITIAN TRIDIA-
C                                    GONAL MATRIX
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 J=1,N
         DO 5 K=1,N
            ZI(J,K)=-ZR(J,K)*TAU(2,J)
            ZR(J,K)=ZR(J,K)*TAU(1,J)
    5 CONTINUE
      IF (N .LE. 2) GO TO 30
C                                  RECOVER THE HOUSEHOLDER MATRICES IN
C                                    REVERSE ORDER
      DO 25 L=3,N
         NR=N-L+2
         NRM1=NR-1
         INX1=(NR*(NRM1))/2+NR
         INX2=INX1-1
         IF (AI(INX1) .EQ. ZERO) GO TO 25
         DELTA=AI(INX1)* SQRT(AR(INX2)**2+AI(INX2)**2)
         DO 20 J=1,N
            ALPHA1=ZERO
            ALPHA2=ZERO
            DO 10 K=NR,N
               K1=(K*(K-1))/2+NRM1
               ALPHA1=ALPHA1+AR(K1)*ZR(K,J)+AI(K1)*ZI(K,J)
               ALPHA2=ALPHA2-AI(K1)*ZR(K,J)+AR(K1)*ZI(K,J)
   10       CONTINUE
            ALPHA1=ALPHA1/DELTA
            ALPHA2=ALPHA2/DELTA
            DO 15 K=NR,N
               K1=(K*(K-1))/2+NRM1
               ZR(K,J)=ZR(K,J)-AR(K1)*ALPHA1+AI(K1)*ALPHA2
               ZI(K,J)=ZI(K,J)-AR(K1)*ALPHA2-AI(K1)*ALPHA1
   15       CONTINUE
   20    CONTINUE
   25 CONTINUE
   30 RETURN
      END


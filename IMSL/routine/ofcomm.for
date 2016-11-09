C   IMSL ROUTINE NAME   - OFCOMM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - COMPUTE AN UNROTATED FACTOR LOADING MATRIX
C                           ACCORDING TO A COMMON FACTOR MODEL BY
C                           UNWEIGHTED OR GENERALIZED LEAST SQUARES,
C                           OR BY MAXIMUM LIKELIHOOD PROCEDURES
C
C   USAGE               - CALL OFCOMM (R,NV,NF,IND,NT,IV,MAXIT,MAXTRY,
C                           EPS,EPSE,ALPHA,V,A,IA,RI,Y,S,G,IS,WK,IER)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE NV BY NV CORRELATION MATRIX IN
C                           SYMMETRIC STORAGE MODE.
C                NV     - INPUT NUMBER OF VARIABLES (AFTER DATA
C                           TRANSFORMATION, IF ANY WAS PERFORMED).
C                NF     - INPUT NUMBER OF FACTORS IN THE MODEL.
C                           NF MUST BE LESS THAN NV.
C                IND    - INPUT OPTION PARAMETER
C                           IND=1 IMPLIES UNWEIGHTED LEAST SQUARES
C                           IND=2 IMPLIES GENERALIZED LEAST SQUARES
C                           IND=3 IMPLIES MAXIMUM LIKELIHOOD
C                NT     - INPUT NUMBER OF SUBJECTS USED IN CALCULATING
C                           R.  NOT REQUIRED IF IND=1.
C                IV     - INPUT OPTION PARAMETER FOR INPUTTING INITIAL
C                           ESTIMATES OF THE UNIQUE VARIANCES OF THE
C                           MODEL IN VECTOR V.  IV=0 IMPLIES THE PROGRAM
C                           USES THE DIAGONAL OF THE INVERSE OF R AS AN
C                           INITIAL ESTIMATE.  OTHERWISE, USER ESTIMATES
C                           ARE INPUT IN VECTOR V.
C                MAXIT  - INPUT MAXIMUM NUMBER OF ITERATIONS ALLOWED
C                           FOR THE NEWTON-RAPHSON MINIMIZATION
C                           PROCEDURE.  MAXIT=30 IS TYPICAL.
C                MAXTRY - INPUT MAXIMUM NUMBER OF SEARCHES FOR A DESCENT
C                           POINT ALLOWED DURING ANY ONE ITERATION
C                           DESCRIBED UNDER MAXIT. MAXTRY=25 IS TYPICAL.
C                EPS    - INPUT CONVERGENCE CRITERION DESCRIBED IN THE
C                           REMARKS SECTION. EPS=0.005 IS TYPICAL.
C                EPSE   - INPUT CONVERGENCE CRITERION, DESCRIBED IN THE
C                           REMARKS SECTION, FOR USING EXACT OR APPROX-
C                           IMATE SECOND DERIVATIVES. EPSE=0.1 IS
C                           TYPICAL.
C                ALPHA  - INPUT/OUTPUT MODEL SIGNIFICANCE LEVEL
C                           PARAMETER, DEFINED ONLY IF IND IS 2 OR 3.
C                           ON INPUT, ALPHA CONTAINS THE DESIRED LEVEL
C                           TO BE ATTAINED, TYPICALLY 0.05.  ON OUTPUT,
C                           THE ATTAINED LEVEL IS RETURNED.
C                V      - INPUT/OUTPUT UNIQUE VARIANCE VECTOR OF LENGTH
C                           NV.  ON INPUT, FOR NONZERO IV, V CONTAINS
C                           INITIAL UNIQUE VARIANCE ESTIMATES.  ON
C                           OUTPUT, V CONTAINS THE UNIQUE VARIANCE
C                           ESTIMATES FOR THE MODEL.  IF A TERMINAL
C                           ERROR OCCURS AND IER EXCEEDS 130, THE
C                           LATEST ESTIMATES ARE RETURNED.
C                A      - OUTPUT NV BY NV MATRIX CONTAINING THE
C                           UNROTATED FACTOR LOADING MATRIX IN THE
C                           FIRST NF COLUMNS. THE REMAINING LOCATIONS
C                           ARE WORK STORAGE.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                RI     - OUTPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE INVERSE OF R.  DEFINED ONLY IF IV=0
C                           OR IF IND IS 2 OR 3.
C                Y      - OUTPUT VECTOR OF LENGTH NV CONTAINING THE
C                           COMMUNALITIES OF THE VARIABLES.
C                S      - OUTPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE NV BY NV NORMALIZED RESIDUAL CORRELATION
C                           MATRIX IN SYMMETRIC STORAGE MODE.
C                           (SEE OFRESI FOR NORMALIZATION DESCRIPTION).
C                G      - OUTPUT SCALING VECTOR OF LENGTH NV NEEDED
C                           FOR LATER USE IN THE ROTATION PROGRAMS.
C                IS     - INTEGER WORK VECTOR OF LENGTH 2*NV. ON OUTPUT,
C                           IS(1) CONTAINS THE NUMBER OF DEGREES OF
C                           FREEDOM FOR THE SIGNIFICANCE TEST STATISTIC
C                           IF IND IS 2 OR 3.  IS(2) CONTAINS THE
C                           NUMBER OF ITERATIONS USED OUT OF THE
C                           MAXIT AVAILABLE.
C                WK     - WORK VECTOR OF LENGTH 5*NV.  ON OUTPUT, WK(1)
C                           CONTAINS THE VALUE OF THE FUNCTION MINIMUM.
C                           WK(2) CONTAINS THE TUCKER RELIABILITY COEFF-
C                           ICIENT.  WK(3) CONTAINS THE CHI-SQUARED TEST
C                           STATISTIC.  WK(4) CONTAINS A VALUE OF EPS
C                           WHICH WOULD HAVE ALLOWED CONVERGENCE WHEN
C                           IER EXCEEDED 130 DURING A TERMINAL ERROR.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES AT LEAST ONE OF NV, NF,
C                             IA, OR IND WAS SPECIFIED INCORRECTLY.
C                           IER = 130 INDICATES R WAS NOT POSITIVE
C                             DEFINITE.
C                           IER = 131 INDICATES MAXIT WAS EXCEEDED.
C                           IER = 132 INDICATES MAXTRY WAS EXCEEDED.
C                           IER = 133 INDICATES EIGRS DID NOT CONVERGE.
C                             CHECK DATA FOR REDUNDANCY OR INCREASE
C                             THE PRECISION OF THE PROGRAM TO DOUBLE.
C                           IER = 134 INDICATES A SECOND DERIVATIVE
C                             MATRIX WAS NOT POSITIVE DEFINITE. CHECK
C                             THE DATA FOR REDUNDANCY AND CHECK THAT
C                             THE NUMBER OF DEGREES OF FREEDOM IS(1) IS
C                             POSITIVE.
C                         WARNING ERROR
C                           IER = 39 INDICATES THAT THE NUMBER OF
C                             DEGREES OF FREEDOM (IS(1)) OF THE CHI-
C                             SQUARED STATISTIC WAS NOT POSITIVE. THE
C                             TEST WAS SKIPPED. ANY RESULTS SHOULD BE
C                             SUSPECT.
C                           IER = 40 INDICATES MDCH COULD NOT CALCULATE
C                             THE PROBABILITY ALPHA ASSOCIATED WITH
C                             IS(1) DEGREES OF FREEDOM WITH VALUE
C                             WK(3). ONLY OUTPUT VARIABLE ALPHA IS
C                             AFFECTED.
C                           IER = 41 INDICATES SIGNIFICANCE LEVEL NOT
C                             ATTAINED IN MODEL WITH NF FACTORS.
C                           IER = 42 INDICATES THAT THE HEYWOOD CASE
C                             WAS ENCOUNTERED ON THE ITERATION THAT A
C                             TERMINAL ERROR WAS DETECTED.
C
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/EHOBKS,EHOUSS,EIGRS,EQRT2S,
C                           LEQ1S,LINV1P,LUDECP,LUELMP,MDCH,MDNOR,
C                           MERRC=ERFC,MGAMAD=DGAMMA,OFRESI,UERTST,
C                           UGETIO,VSRTU
C                       - SINGLE(H36,H48,H60)/EHOBKS,EHOUSS,EIGRS,
C                           EQRT2S,LEQ1S,LINV1P,LUDECP,LUELMP,MDCH,
C                           MDNOR,MERRC=ERFC,MGAMA=GAMMA,OFRESI,
C                           UERTST,UGETIO,VSRTU
C                       - DOUBLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,LEQ1S,
C                           LINV1P,LUDECP,LUELMP,MDCH,MDNOR,
C                           MERRC=ERFC,MGAMAD=DGAMMA,OFRESI,UERTST,
C                           UGETIO,VSRTUD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  EPS - IF THE LARGEST ABSOLUTE RELATIVE CORRECTION TO
C                THE UNIQUE STANDARD DEVIATION VECTOR BETWEEN TWO
C                ITERATIONS IS LESS THAN EPS, CONVERGENCE IS ASSUMED.
C                AS THE UNIQUE VARIANCES APPROACH ZERO, THE RELATIVE
C                CORRECTION IS WEIGHTED TOWARDS ABSOLUTE CORRECTION.
C            2.  EPSE - APPROXIMATE SECOND DERIVATIVE MATRICES ARE
C                USED IN THE NEWTON-RAPHSON PROCEDURE UNTIL THE LARGEST
C                CORRECTION DESCRIBED UNDER EPS IS LESS THAN EPSE. A
C                GREAT DEAL OF COMPUTATION IS THEREBY SAVED.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFCOMM (R,NV,NF,IND,NT,IV,MAXIT,MAXTRY,EPS,EPSE,
     1                   ALPHA,V,A,IA,RI,Y,S,G,IS,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NV,NF,IND,NT,IV,MAXIT,MAXTRY,IA,IS(1),IER
      REAL               V(1),R(1),RI(1),Y(1),S(1),A(IA,1),G(1),WK(1),
     *                   EPS,EPSE,ALPHA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NFP1,NVP1,NV3,NV4,NV5,NVV,MORE,IR,ITER,ITRV,
     *                   IJUMP,I,J,M,K,N,IHM,INDX,NDF
C                                  PARAMETERS CHSQ, XNDF, AND PROB MUST
C                                  BE SINGLE PRECISION FOR CALL TO MDCH
      REAL               CHSQ,PROB,XNDF
      REAL               C0,C1,C2,DB,DD,DDD,DET,DM,D1,D2,D3,D4,D5,D6,D7,
     *                   EP,EXCH,F,FM0,FT,F0,HALF,ONE,RERR,SIX,TEN,
     *                   TENTH,THREE,TWO,ZERO,OLDIF,T57
      DOUBLE PRECISION   SAVE,TEMP,SUM,SSS,VVV
      DATA               ZERO,TENTH,HALF,ONE,TWO,THREE,SIX,TEN,C1,C2,D1,
     *                   D2,D3,D4,D5,D6 /0.0E0,0.1E0,0.5E0,1.0E0,2.0E0,
     *                   3.0E0,6.0E0,10.0E0,9.99999E-1,1.0E-6,-2.0E-6,
     *                   -2.7631E1,-2.302585E0,1.0E20,1.0E30,5.0E-2/
      DATA               D7 /1.0E-7/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (NF.GT.0 .AND. NV.GT.NF .AND. IA.GE.NV .AND. IND.GE.1 .AND.
     *    IND.LE.3) GO TO 5
      IER = 129
      GO TO 9000
    5 NFP1 = NF+1
      NVP1 = NV+1
      NV3 = NV+NV
      NV4 = NV3+NV
      NV5 = NV4+NVP1
      NVV = ((NV+1)*NV)/2
C                                  BEGIN NEWTON-RAPHSON MINIMIZATION
C                                  MORE INDICATES CONVERGENCE STATUS
      MORE = 0
      F0 = D4
      OLDIF = D4
C                                  INITIALIZE THE SQUARE ROOT OF
C                                  THE UNIQUE VARIANCE VECTOR
      IF (IND.EQ.1 .AND. IV.NE.0) GO TO 25
      DO 10 I=1,NVV
         S(I) = R(I)
   10 CONTINUE
      CALL LINV1P (S,NV,RI,0,DM,DB,IER)
      IF (IER.EQ.0) GO TO 15
      IER = 130
      GO TO 9000
   15 DET = ALOG(DM)+DB*ALOG(TWO)
      IF (IV.NE.0) GO TO 25
      IR = 0
      DO 20 I=1,NV
         IR = IR+I
         V(I) = RI(IR)
   20 CONTINUE
   25 FT = ONE-NF/(TWO*NV)
      DO 30 I=1,NV
         DD = V(I)
         IF (IV.EQ.0) DD = FT/DD
         IF (IND.EQ.1) V(I) = SQRT(DD)
         IF (IND.NE.1) V(I) = ALOG(DD)
   30 CONTINUE
      ITER = 0
C                                  BEGIN FACTORING OF S MATRIX
   35 ITRY = 0
C                                  CALCULATE THE S MATRIX
      IJUMP = 0
   40 IF (IND.NE.1) GO TO 55
      IR = 0
      DO 50 I=1,NV
         WK(I) = V(I)*V(I)
         DO 45 J=1,I
            IR = IR+1
            S(IR) = R(IR)
   45    CONTINUE
         S(IR) = S(IR)-WK(I)
         IS(I) = NVP1-I
   50 CONTINUE
      GO TO 75
   55 DO 60 I=1,NV
         WK(I) = EXP(HALF*V(I))
         IS(I) = NVP1-I
   60 CONTINUE
      IR = 0
      DO 70 I=1,NV
         DD = WK(I)
         DO 65 J=1,I
            IR = IR+1
            S(IR) = R(IR)/(DD*WK(J))
   65    CONTINUE
   70 CONTINUE
C                                  CALCULATE THE EIGENPAIR OF S
   75 CALL EIGRS (S,NV,1,Y,A,IA,WK (NV5),IER)
      IF (IER.EQ.0) GO TO 80
      IER = 133
      GO TO 365
C                                  REVERSE THE EIGENPAIR ORDER
   80 M = NV/2
      DO 85 I=1,M
         K = NVP1-I
         EXCH = Y(K)
         Y(K) = Y(I)
         Y(I) = EXCH
         IF (IND.EQ.1) GO TO 85
C                                  EIGENVALUES OF S-INVERSE
         Y(I) = ONE/Y(I)
         Y(K) = ONE/Y(K)
   85 CONTINUE
      IF (M+M.NE.NV .AND. IND.NE.1) Y(M+1) = ONE/Y(M+1)
      CALL VSRTU (A,IA,NV,NV,0,IS,WK (NV5))
      IF (IJUMP.EQ.1) GO TO 310
C                                  FIRST DERIVATIVE VECTOR
      DO 95 I=1,NV
         TEMP = 0.0D0
         DO 90 M=NFP1,NV
            SSS = Y(M)
            IF (IND.EQ.2) SSS = SSS*(SSS-1.0D0)
            IF (IND.EQ.3) SSS = 1.0D0-1.0D0/SSS
            TEMP = TEMP+SSS*DBLE(A(I,M))**2
   90    CONTINUE
         IF (IND.EQ.1) G(I) = -V(I)*(TEMP+TEMP)
         IF (IND.NE.1) G(I) = TEMP
   95 CONTINUE
C                                  EVALUATE FUNCTION TO BE MINIMIZED
      TEMP = 0.0D0
      DO 100 M=NFP1,NV
         IF (IND.EQ.1) TEMP = TEMP+DBLE(Y(M))**2
         IF (IND.EQ.2) TEMP = TEMP+(Y(M)-1.0D0)**2
         IF (IND.EQ.3) TEMP = TEMP+1.0D0/Y(M)+ALOG(Y(M))
  100 CONTINUE
      IF (IND.NE.3) F = 0.5D0*TEMP
      IF (IND.EQ.3) F = TEMP+(NF-NV)
C                                  IS NEW FUNCTION VALUE BETTER
      T57 = F-F0
      IF (T57.GT.C2) GO TO 110
C                                  IF YES, SAVE UNIQUE STD DEV VECTOR
      F0 = F
      OLDIF = T57
      DO 105 I=1,NV
         WK(NV+I) = V(I)
  105 CONTINUE
      GO TO 135
C                                  IF NOT, APPROACH OLD VECTOR BY
C                                  HALVING UP TO MAXTRY TIMES
  110 ITRY = ITRY+1
      IF (ITRY.LE.MAXTRY) GO TO 115
      IER = 132
      GO TO 365
  115 IF (ABS(T57).GT.ABS(OLDIF)) GO TO 125
      OLDIF = T57
      F0 = F
      DO 120 I=1,NV
         WK(NV+I) = V(I)
  120 CONTINUE
      GO TO 135
  125 DO 130 I=1,NV
         V(I) = HALF*(V(I)+WK(NV+I))
  130 CONTINUE
      IJUMP = 0
      GO TO 40
C                                  CHECK THE NUMBER OF ITERATIONS
  135 ITER = ITER+1
      IF (ITER.LE.MAXIT) GO TO 140
      IER = 131
      GO TO 365
C                                  MORE=2 IMPLIES CONVERGENCE
  140 IF (MORE.EQ.2) GO TO 305
C                                  CALCULATE UNIQUE VARIANCES GIVEN
C                                  FACTORED MATRIX
C                                  MORE=0 IMPLIES APPROXIMATE
C                                  ALGORITHM USED FOR SECOND DERIVATIVES
      IF (MORE.EQ.0) GO TO 240
C                                  MORE=1 IMPLIES EXACT ALGORITHM
C                                  USED FOR SECOND DERIVATIVE MATRIX
      IF (IND.NE.1) GO TO 170
C                                  EXACT FOR IND=1
      IR = 0
      DO 165 I=1,NV
         VVV = V(I)
         SAVE = VVV*VVV
         VVV = 4.0D0*VVV
         DO 155 J=1,I
            IR = IR+1
            TEMP = 0.0D0
            DO 150 M=NFP1,NV
               SUM = 0.0D0
               SSS = Y(M)
               DO 145 N=1,NF
                  SUM = SUM+(DBLE(A(I,N))*A(J,N)*(SSS+Y(N)))/(SSS-Y(N))
  145          CONTINUE
               TEMP = TEMP+SUM*A(I,M)*A(J,M)
  150       CONTINUE
            S(IR) = VVV*V(J)*TEMP
  155    CONTINUE
         TEMP = S(IR)
         DO 160 M=NFP1,NV
            TEMP = TEMP+4.0D0*(SAVE-Y(M)*0.5D0)*DBLE(A(I,M))**2
  160    CONTINUE
         S(IR) = TEMP
  165 CONTINUE
      GO TO 210
C                                  EXACT FOR IND=2,3
  170 DD = -ONE
      IF (IND.EQ.3) GO TO 185
      DD = ONE
      DO 180 J=1,NV
         DDD = SQRT(Y(J))
         DO 175 I=1,NV
            A(I,J) = A(I,J)*DDD
  175    CONTINUE
  180 CONTINUE
  185 IR = 0
      DO 205 I=1,NV
         DO 200 J=1,I
            IR = IR+1
            IF (IND.EQ.2) SSS = DBLE(RI(IR))*WK(I)*WK(J)
            TEMP = 0.0D0
            DO 195 M=NFP1,NV
               VVV = Y(M)
               SUM = 0.0D0
               DO 190 N=1,NF
                  SUM = SUM+(DBLE(A(I,N))*A(J,N)*(-2.0D0+VVV+Y(N)))/
     *            (VVV-Y(N))
  190          CONTINUE
               IF (IND.EQ.2) SUM = SUM+SSS
               IF (IND.EQ.3 .AND. I.EQ.J) SUM = SUM+1.0D0
               TEMP = TEMP+SUM*A(I,M)*A(J,M)
  195       CONTINUE
            S(IR) = TEMP
  200    CONTINUE
         S(IR) = S(IR)+DD*G(I)
  205 CONTINUE
  210 DO 215 I=1,NV
         WK(NV4+I) = G(I)
  215 CONTINUE
C                                  ADJUST FOR HEYWOOD CASE
  220 IR = 0
      IHM = 0
      DO 230 I=1,NV
         IR = IR+I
         WK(NV3+I) = D5
         IF (S(IR).GT.D6) GO TO 230
         IHM = IHM+1
         WK(NV3+I) = ZERO
         IF (S(IR).GE.D7) WK(NV3+I) = G(I)/S(IR)
         DO 225 J=1,NV
            K = MIN0(I,J)
            M = MAX0(I,J)
            INDX = K+((M-1)*M)/2
            S(INDX) = ZERO
  225    CONTINUE
         S(IR) = ONE
  230 CONTINUE
C                                  CALCULATE CORRECTION TO UNIQUE
C                                  STANDARD DEVIATION VECTOR
      CALL LEQ1S (S,NV,G,1,NV,0,IS,WK (NV5),IER)
      IF (MORE.EQ.0) GO TO 270
      IF (IER.EQ.0) GO TO 275
      DO 235 I=1,NV
         G(I) = WK(NV4+I)
  235 CONTINUE
C                                  APPROXIMATE SECOND DERIVATIVES
  240 IR = 0
      DO 255 I=1,NV
         DO 250 J=1,I
            IR = IR+1
            SUM = 0.0D0
            DO 245 M=NFP1,NV
               SUM = SUM+DBLE(A(I,M))*A(J,M)
  245       CONTINUE
            S(IR) = SUM*SUM
  250    CONTINUE
  255 CONTINUE
      IF (IND.NE.1) GO TO 220
      IR = 0
      DO 265 I=1,NV
         VVV = 4.0D0*V(I)
         DO 260 J=1,I
            IR = IR+1
            S(IR) = S(IR)*VVV*V(J)
  260    CONTINUE
  265 CONTINUE
      GO TO 220
  270 IF (IER.EQ.0) GO TO 275
      IER = 134
      GO TO 365
C                                  ADJUST CORRECTIONS FOR HEYWOOD CASE
  275 IF (IHM.LE.0) GO TO 285
      DO 280 I=1,NV
         IF (WK(NV3+I).LE.D4) G(I) = WK(NV3+I)
  280 CONTINUE
C                                  CALCULATE RELATIVE CHANGE
  285 EP = ZERO
      DO 300 I=1,NV
         V(I) = WK(NV+I)-G(I)
         RERR = ABS(G(I))
         IF (IND.EQ.1) GO TO 290
         IF (V(I).GT.D1) V(I) = D1
         IF (V(I).LT.D2) V(I) = D2
C                                  RERR IS DECREASED IF V(I) .LT. 0.1
         IF (V(I).LE.D3) RERR = RERR*TEN*EXP(V(I))
         GO TO 295
  290    IF (V(I).GT.C1) V(I) = C1
         IF (V(I).LT.C2) V(I) = C2
         RERR = RERR/AMAX1(V(I),TENTH)
  295    EP = AMAX1(EP,RERR)
  300 CONTINUE
C                                  SET MORE, THE CONVERGENCE STATUS
      MORE = 0
      IF (EP.LT.EPSE) MORE = 1
      IF (EP.LT.EPS) MORE = 2
      GO TO 35
C                                  ALGORITHM HAS CONVERGED
C                                  FACTOR MATRIX ONE LAST TIME
  305 IJUMP = 1
      GO TO 40
C                                  SCALE LOADING MATRIX
C                                  SET SCALING VECTOR
  310 IS(2) = ITER
      IF (IND.NE.1) GO TO 330
      DO 320 J=1,NF
         DD = SQRT(Y(J))
         DO 315 I=1,NV
            A(I,J) = A(I,J)*DD
  315    CONTINUE
  320 CONTINUE
      DO 325 I=1,NV
         G(I) = ONE
         V(I) = WK(I)
  325 CONTINUE
      GO TO 350
  330 DO 340 J=1,NF
         DD = SQRT(ONE/Y(J)-ONE)
         DO 335 I=1,NV
            A(I,J) = A(I,J)*WK(I)*DD
  335    CONTINUE
  340 CONTINUE
      DO 345 I=1,NV
         G(I) = WK(I)
         V(I) = WK(I)**2
  345 CONTINUE
C                                  COMMUNALITIES, RESIDUAL CORRELATIONS
  350 CALL OFRESI (R,NV,NF,A,IA,Y,S,WK (NV5))
      WK(1) = F0
C                                  SIGNIFICANCE TESTING FOR MODEL
      IF (IND.EQ.1) GO TO 9005
      NDF = ((NV-NF)**2-NV-NF)/2
      IS(1) = NDF
      IF (NDF.GT.0) GO TO 355
      IER = 39
      GO TO 9000
  355 C0 = (NT-1)-(NV+NV+5)/SIX
      CHSQ = (C0-2*NF/THREE)*F0
      FM0 = -DET/(HALF*(NV-1)*NV)
      FM0 = (FM0 - F0/NDF)/(FM0 - ONE/(C0-2*NF/THREE))
      IF (FM0.GT.ONE) FM0 = ONE
      WK(2) = FM0
      WK(3) = CHSQ
C                                  NOTE MDCH IS SINGLE PRECISION
      XNDF = NDF
      CALL MDCH (CHSQ,XNDF,PROB,IER)
      IF (IER.EQ.0) GO TO 360
      IER = 40
      GO TO 9000
  360 VVV = ALPHA
      ALPHA = ONE-PROB
      IF (ALPHA.GT.VVV) GO TO 9005
      IER = 41
      GO TO 9000
C                                  TRAP ROUTINE FOR OUTPUTTING
C                                  UNIQUE VARIANCES FOR IER .GT. 130
  365 DO 375 I=1,NV
         DD = WK(NV+I)
         IF (IND.NE.1) GO TO 370
         IF (DD.LT.ZERO .OR. DD.GT.ONE) DD = ZERO
         V(I) = DD**2
         GO TO 375
  370    IF (DD.LT.D2 .OR. DD.GT.D1) DD = D2
         V(I) = EXP(DD)
  375 CONTINUE
      IS(1) = ((NV-NF)**2-NV-NF)/2
      IS(2) = ITER
      WK(1) = F0
      WK(4) = EP
 9000 CONTINUE
      CALL UERTST (IER,6HOFCOMM)
      IF (IHM.GT.0 .AND. IER.GT.129) CALL UERTST (42,6HOFCOMM)
 9005 RETURN
      END

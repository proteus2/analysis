C   IMSL ROUTINE NAME   - OFCOEF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - COMPUTE A MATRIX OF FACTOR SCORE COEFFICIENTS
C                           FOR INPUT TO IMSL ROUTINE OFSCOR
C
C   USAGE               - CALL OFCOEF (B,IB,NV,NF,IND,S,T,IT,C,IC,IS,
C                           WK,IER)
C
C   ARGUMENTS    B      - NV BY NF INPUT ROTATED FACTOR LOADING OR
C                           STRUCTURE MATRIX.  IF IND = 5, THEN B SHOULD
C                           CONTAIN THE UNROTATED IMAGE SCORE COEFF-
C                           ICIENT MATRIX CI CALCULATED BY IMSL
C                           ROUTINE OFIMAG.
C                IB     - INPUT ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NV     - INPUT NUMBER OF VARIABLES.
C                NF     - INPUT NUMBER OF FACTORS.
C                IND    - INPUT OPTION PARAMETER
C                           IND=1 IMPLIES REGRESSION METHOD
C                           IND=2 IMPLIES LEAST SQUARES METHOD
C                           IND=3 IMPLIES THE BARTLETT METHOD
C                           IND=4 IMPLIES ANDERSON AND RUBIN METHOD
C                           IND=5 IMPLIES IMAGE SCORE FOR IMAGE ANALYSIS
C                S      - IF IND IS 1 OR 4, INPUT VECTOR OF LENGTH
C                           (NV+1)*NV/2 CONTAINING A COPY OF THE NV BY
C                           NV CORRELATION MATRIX (CALL IT R) IN SYM-
C                           METRIC STORAGE MODE. (IF IND=4, S MUST BE
C                           OF LENGTH THE MAXIMUM OF (NV*NF,
C                           (NV+1)*NV/2). IN THIS CASE, S IS DESTROYED
C                           ON OUTPUT.) IF IND IS 2 OR 3, S IS A WORK
C                           VECTOR OF LENGTH (NF+1)*NF/2. IF IND=5, S
C                           IS NOT USED.
C                T      - IF IND=5, INPUT NF BY NF MATRIX CONTAINING
C                           THE IMAGE TRANSFORMATION MATRIX T FROM IMSL
C                           ROUTINE OFROTA, OFHARR, OFPROT, OR OFSCHN.
C                           IF IND=1, T IS NOT USED.  OTHERWISE, A WORK
C                           VECTOR OF LENGTH NF*NF IF IND=4, OR OF
C                           LENGTH (NF+1)*NF/2 IF IND IS 2 OR 3.
C                IT     - IF IND=5, INPUT ROW DIMENSION OF MATRIX T
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM. OTHERWISE,
C                           IT IS NOT USED.
C                C      - OUTPUT NV BY NF FACTOR SCORE COEFFICIENT
C                           MATRIX.
C                IC     - INPUT ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IS     - IF IND=4, INTEGER WORK VECTOR OF LENGTH NF.
C                           OTHERWISE, IS IS NOT USED.
C                WK     - IF IND IS 3 OR 4, A WORK VECTOR OF LENGTH
C                           NV + 2*NF.  OTHERWISE, WK IS NOT USED.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE CORRELATION
C                             MATRIX IS NOT POSITIVE DEFINITE NUMERI-
C                             CALLY.
C                           IER = 130 INDICATES THE ANDERSON AND RUBIN
C                             METHOD FAILED TO CALCULATE EIGENVALUES,
C                             S IS NOT POSITIVE DEFINITE OR B IS NOT
C                             OF FULL COLUMN RANK NF.
C                           IER = 131 INDICATES AT LEAST ONE OF NV, NF,
C                             IB, OR IC, WAS SPECIFIED INCORRECTLY.
C
C   REQD. IMSL ROUTINES - SINGLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,LEQT1P,
C                           LINV1P,LUDECP,LUELMP,UERTST,UGETIO,VIPRFF,
C                           VMULFF,VMULFS,VSRTU,VTPROF
C                       - DOUBLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,LEQT1P,
C                           LINV1P,LUDECP,LUELMP,UERTST,UGETIO,VIPRFF,
C                           VMULFF,VMULFS,VSRTUD,VTPROF,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFCOEF  (B,IB,NV,NF,IND,S,T,IT,C,IC,IS,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IB,NV,NF,IND,IT,IC,IS(1),IER
      REAL               B(IB,NF),T(1),C(IC,NF),WK(1),S(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NF2,NFP1,IR,K,N,I1,I2,INDX,M
      REAL               DD,D1,D2,EXCH
      DOUBLE PRECISION   TEMP,VVV
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (NV.GE.NF.AND.IB.GE.NV.AND.IC.GE.NV) GO TO 5
      IER = 131
      GO TO 9000
    5 IER = 0
      GO TO (15,30,35,35,10), IND
C                                  IMAGE SCORE METHOD
   10 CALL VMULFF (B,T,NV,NF,NF,IB,IT,C,IC,IER)
      GO TO 9005
C                                  REGRESSION METHOD
   15 DO 25 I=1,NV
         DO 20 J=1,NF
            C(I,J) = B(I,J)
   20    CONTINUE
   25 CONTINUE
      CALL LEQT1P (S,NF,NV,C,IC,0,D1,D2,IER)
      IF (IER.NE.0) GO TO 9000
      GO TO 9005
C                                  LEAST SQUARES METHOD
   30 CALL VTPROF (B,NV,NF,IB,S)
      GO TO 125
   35 NF2 = NF+NF
      NFP1 = NF+1
      DO 45 I=1,NV
         TEMP = 0.0D0
         DO 40 J=1,NF
            TEMP = TEMP+DBLE(B(I,J))**2
   40    CONTINUE
         WK(NF2+I) = 1.0D0/(1.0D0-TEMP)
   45 CONTINUE
      IF (IND.EQ.4) GO TO 65
C                                  BARTLETT METHOD
      IR = 0
      DO 60 I=1,NF
         DO 55 J=1,I
            IR = IR+1
            TEMP = 0.0D0
            DO 50 K=1,NV
               TEMP = TEMP+DBLE(B(K,I))*B(K,J)*WK(NF2+K)
   50       CONTINUE
            S(IR) = TEMP
   55    CONTINUE
   60 CONTINUE
      GO TO 125
C                                  ANDERSON - RUBIN METHOD
   65 IR = 0
      DO 75 I=1,NV
         DD = WK(NF2+I)
         DO 70 J=1,I
            IR = IR+1
            S(IR) = S(IR)*DD*WK(NF2+J)
   70    CONTINUE
   75 CONTINUE
      IR = 0
      DO 95 I=1,NF
         IS(I) = NFP1-I
         DO 90 J=1,I
            IR = IR+1
            TEMP = 0.0D0
            DO 85 K=1,NV
               VVV = B(K,J)
               DO 80 N=1,NV
                  I1 = MIN0(K,N)
                  I2 = MAX0(K,N)
                  INDX = I1+((I2-1)*I2)/2
                  TEMP = TEMP+DBLE(B(N,I))*VVV*S(INDX)
   80          CONTINUE
   85       CONTINUE
            T(IR) = TEMP
   90    CONTINUE
   95 CONTINUE
      CALL EIGRS (T,NF,1,WK,S,NF,WK(NFP1),IER)
      IF (IER.EQ.0) GO TO 100
      IER = 130
      GO TO 9000
C                                  REVERSE EIGENPAIR ORDER
  100 M = NF/2
      IF (NF.EQ.1) GO TO 110
      DO 105 I=1,M
         K = NFP1-I
         EXCH = WK(K)
         WK(K) = SQRT(SQRT(WK(I)))
         WK(I) = SQRT(SQRT(EXCH))
  105 CONTINUE
  110 IF (M+M.NE.NF) WK(M+1) = SQRT(SQRT(WK(M+1)))
      CALL VSRTU (S,NF,NF,NF,0,IS,WK(NF+1))
C                                  SCALE AND TRANSPOSE S
      IR = 0
      DO 120 I=1,NF
         DD = WK(I)
         K = 0
         DO 115 J=1,NF
            IR = IR+1
            T(I+K) = S(IR)*DD
            K = K+NF
  115    CONTINUE
  120 CONTINUE
      CALL VTPROF (T,NF,NF,NF,S)
  125 CALL LINV1P (S,NF,T,0,D1,D2,IER)
      IF (IER.NE.0) GO TO 9000
      CALL VMULFS (B,T,NV,NF,IB,C,IC)
      IF (IND.EQ.2) GO TO 9005
      DO 135 I=1,NV
         DD = WK(NF2+I)
         DO 130 J=1,NF
            C(I,J) = C(I,J)*DD
  130    CONTINUE
  135 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFCOEF)
 9005 RETURN
      END

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

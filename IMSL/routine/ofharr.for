C   IMSL ROUTINE NAME   - OFHARR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - TRANSFORMATION OF UNROTATED FACTOR LOADING
C                           MATRIX TO OBLIQUE AXES BY HARRIS-KAISER
C                           METHOD
C
C   USAGE               - CALL OFHARR (A,IA,NV,NF,NORM,II,MAXIT,W,EPS,
C                           DELTA,C,G,B,IB,T,IT,F,S,WK,IER)
C
C   ARGUMENTS    A      - INPUT NV BY NF UNROTATED FACTOR LOADING
C                           MATRIX.  ON OUTPUT, A CONTAINS THE NV BY NF
C                           PATTERN MATRIX.
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
C                MAXIT  - INPUT MAXIMUM NUMBER OF ITERATIONS ALLOWED FOR
C                           ROTATION BY OFROTA.  MAXIT=30 IS TYPICAL.
C                W      - INPUT CONSTANT FOR ROTATION (SEE REMARKS).
C                EPS    - INPUT CONVERGENCE CONSTANT FOR ROTATION
C                           (ANGLE).  EPS=0.0001 IS TYPICAL.
C                DELTA  - INPUT CONVERGENCE CONSTANT FOR ROTATION
C                           (CRITERION FUNCTION). DELTA=.001 IS TYPICAL.
C                C      - INPUT CONSTANT BETWEEN ZERO AND ONE USED
C                           IN THE ROTATION.  (SEE REMARKS).
C                G      - INPUT SCALING VECTOR OF LENGTH NV. (PREVIOUS
C                           OUTPUT FROM OFCOMM, OFIMAG, OR OFPRI.)
C                B      - OUTPUT NV BY NF STRUCTURE MATRIX.
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
C                S      - OUTPUT VECTOR OF LENGTH (NF+1)*NF/2 CONTAINING
C                           THE NF BY NF FACTOR CORRELATION MATRIX IN
C                           SYMMETRIC STORAGE MODE.
C                WK     - WORK VECTOR OF LENGTH 2*NF.  IF II PARAMETER
C                           IS NONZERO, WK IS OF LENGTH (NV*NF) + THE
C                           MAXIMUM OF (2*NF) AND ((NF+1)*NF/2).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT AT LEAST ONE OF IA,
C                             IB, IT, OR C WAS SPECIFIED INCORRECTLY.
C                           IER = 130 INDICATES THAT (ONLY) THE IMAGE
C                             TRANSFORMATION MATRIX T COULD NOT BE
C                             CALCULATED BECAUSE THE INPUT MATRIX
C                             A WAS NOT OF FULL COLUMN RANK.
C                         WARNING ERROR (WITH FIX)
C                           IER = 67 INDICATES THAT MAXIT WAS EXCEEDED.
C                             THE APPROXIMATE SOLUTION WAS USED.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1P,LUDECP,LUELMP,OFIMA3,OFROTA,
C                           UERTST,UGETIO,VIPRFF,VMULFM,VMULFS,VTPROF
C                       - DOUBLE/LEQT1P,LUDECP,LUELMP,OFIMA3,OFROTA,
C                           UERTST,UGETIO,VIPRFF,VMULFM,VMULFS,VTPROF,
C                           VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  INPUT ARGUMENTS W, EPS, DELTA, AND NORM ARE INPUT
C                ARGUMENTS TO IMSL ROUTINE OFROTA TO OBTAIN THE TRANS-
C                FORMATION MATRIX. W MAY BE SET AS FOLLOWS.
C
C                W = 0.0 IS THE QUARTIMAX METHOD, WHICH ATTEMPTS TO
C                GET EACH VARIABLE TO LOAD HIGHLY ON ONLY ONE (OR A
C                FEW) FACTOR(S).
C
C                W = 1.0 IS THE VARIMAX METHOD, WHICH ATTEMPTS TO LOAD
C                HIGHLY A RELATIVELY LOW NUMBER OF VARIABLES ON EACH
C                FACTOR. VARIMAX IS MOST WIDELY USED.
C
C                W = NF/2.0 IS THE EQUAMAX METHOD, WHICH IS A COMPRO-
C                MISE OF THE ABOVE TWO.
C
C                W CAN BE ANY REAL NUMBER, BUT BEST VALUES LIE IN THE
C                CLOSED INTERVAL (1.0, 5.0*NF). GENERALLY, THE LARGER
C                W IS, THE MORE EQUAL IS THE DISPERSION OF THE VARIANCE
C                ACCOUNTED FOR ACROSS THE FACTORS.
C
C            2.  INPUT ARGUMENT C MUST BE BETWEEN ZERO AND ONE. AS C
C                INCREASES TO ONE, THE SOLUTION BECOMES ORTHOGONAL.
C                FOR FACTORIALLY SIMPLE DATA, C=0.0 IS OFTEN BEST. FOR
C                MORE COMPLEX DATA, C=0.5 IS NEAR OPTIMAL. RARELY
C                SHOULD C BE GREATER THAN A HALF.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFHARR (A,IA,NV,NF,NORM,II,MAXIT,W,EPS,DELTA,
     *                   C,G,B,IB,T,IT,F,S,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NV,NF,NORM,II,MAXIT,IB,IT,IER
      REAL               A(IA,NF),G(NV),S(1),B(IB,NF),F(NF),T(IT,NF),
     *                   WK(1),C,DELTA,EPS,W
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NF1,IR,J,I,K,JER,IJ
      REAL               CEX,DD,ONE,TWO,ZERO
      DOUBLE PRECISION   TEMP,CS2
      DATA               ZERO,ONE,TWO/0.0E0,1.0E0,2.0E0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  ERROR CHECKS
      IF (IA.GE.NV.AND.IB.GE.NV.AND.IT.GE.NF.AND.C.GE.ZERO.AND.C.LE.ONE
     *   .AND.NV.GE.NF) GO TO 5
      IER = 129
      GO TO 9000
    5 IER = 0
      IF (II.EQ.0) GO TO 20
      NF1 = MAX0(NF+NF,((NF+1)*NF)/2)+1
      IR = NF1
      DO 15 J=1,NF
         DO 10 I=1,NV
            WK(IR) = A(I,J)
            IR = IR+1
   10    CONTINUE
   15 CONTINUE
   20 CS2 = 0.5D0*C
C                                  SCALE FACTOR LOADING MATRIX
      DO 30 J=1,NF
         TEMP = 0.0D0
         DO 25 I=1,NV
            A(I,J) = A(I,J)/G(I)
            TEMP = TEMP+DBLE(A(I,J))**2
   25    CONTINUE
         WK(J) = 1.0D0/DSQRT(TEMP)
         WK(J+NF) = TEMP**CS2
   30 CONTINUE
      DO 40 J=1,NF
         DD = WK(J)*WK(J+NF)
         DO 35 I=1,NV
            A(I,J) = A(I,J)*DD
   35    CONTINUE
   40 CONTINUE
C                                  ROTATE FACTOR LOADING MATRIX
      CALL OFROTA (A,IA,NV,NF,NORM,0,MAXIT,W,EPS,DELTA,B,IB,T,IT,F,S,IER
     1)
      IF (IER.EQ.66) IER = 67
      CEX = TWO*(C-ONE)
      DO 45 I=1,NF
         WK(I) = WK(I)**CEX
   45 CONTINUE
C                                  CORRELATION MATRIX OF FACTORS
      IR = 0
      DO 60 I=1,NF
         DO 55 J=1,I
            IR = IR+1
            TEMP = 0.0D0
            DO 50 K=1,NF
               TEMP = TEMP+DBLE(T(K,I))*T(K,J)*WK(K)
   50       CONTINUE
            S(IR) = TEMP
   55    CONTINUE
   60 CONTINUE
      IR = 0
      IJ = 0
      DO 70 I=1,NF
         IJ = IJ + I
         DD = ONE/SQRT(S(IJ))
         WK(I) = DD
         DO 65 J=1,I
            IR = IR+1
            S(IR) = S(IR)*DD*WK(J)
   65    CONTINUE
   70 CONTINUE
C                                  PATTERN MATRIX
      DO 80 I=1,NV
         DD = G(I)
         DO 75 J=1,NF
            A(I,J) = B(I,J)*DD/WK(J)
   75    CONTINUE
   80 CONTINUE
C                                  STRUCTURE MATRIX
      CALL VMULFS (A,S,NV,NF,IA,B,IB)
C                                  VARIANCE ACCOUNTED FOR BY FACTORS
      DO 90 J=1,NF
         TEMP = 0.0D0
         DO 85 I=1,NV
            TEMP = TEMP+DBLE(A(I,J))*B(I,J)
   85    CONTINUE
         F(J) = TEMP
   90 CONTINUE
      IF (II.EQ.0) GO TO 95
      JER = IER
      CALL OFIMA3 (WK(NF1),NV,B,IB,NV,NF,NF,T,IT,WK,IER)
      IF (JER.GT.IER) IER = JER
   95 IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFHARR)
 9005 RETURN
      END

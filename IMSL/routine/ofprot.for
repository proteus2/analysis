C   IMSL ROUTINE NAME   - OFPROT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - OBLIQUE TRANSFORMATION OF THE FACTOR LOADING
C                           MATRIX USING A TARGET MATRIX, INCLUDING
C                           PIVOT AND POWER VECTOR OPTIONS
C
C   USAGE               - CALL OFPROT (A,IA,NV,NF,IND,NORM,II,MAXIT,W,
C                           EPS,DELTA,F,X,IX,B,IB,T,IT,S,WK,IER)
C
C   ARGUMENTS    A      - INPUT NV BY NF UNROTATED FACTOR LOADING
C                           MATRIX.  ON OUTPUT, A CONTAINS THE NV BY NF
C                           FACTOR PATTERN MATRIX.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NV     - INPUT NUMBER OF VARIABLES.
C                NF     - INPUT NUMBER OF FACTORS.
C                IND    - INPUT OPTION PARAMETER.
C                           IND=5 IMPLIES THE PROMAX METHOD.
C                           IND=6 IMPLIES THE PIVOTAL PROMAX METHOD.
C                           IND=7 IMPLIES OBLIQUE PROCRUSTES METHOD.
C                             (TARGET MATRIX X IS INPUT FOR IND=7).
C                NORM   - IF IND IS 5 OR 6, INPUT OPTION PARAMETER.
C                           NORM=0 INDICATES NO ROW NORMALIZATION OF A
C                           IS REQUIRED.  OTHERWISE, ROW NORMALIZATION
C                           IS PERFORMED.  NORM IS NOT USED IF IND=7.
C                II     - INPUT OPTION PARAMETER.  II=0 INDICATES AN
C                           IMAGE ANALYSIS IS NOT BEING PERFORMED.
C                           OTHERWISE, AN IMAGE ANALYSIS IS ASSUMED, AND
C                           T BECOMES THE IMAGE TRANSFORMATION MATRIX.
C                MAXIT  - INPUT MAXIMUM NUMBER OF ITERATIONS ALLOWED
C                           FOR ROTATION BY OFROTA. MAXIT=30 IS TYPICAL.
C                W      - INPUT CONSTANT FOR ROTATION (SEE REMARKS).
C                           NOT REQUIRED IF IND = 7.
C                EPS    - INPUT CONVERGENCE CONSTANT FOR ROTATION
C                           (ANGLE).  EPS=0.0001 IS TYPICAL.
C                           (SEE REMARKS).  NOT REQUIRED IF IND=7.
C                DELTA  - INPUT CONVERGENCE CONSTANT FOR ROTATION
C                           (CRITERION FUNCTION). DELTA=.001 IS TYPICAL.
C                           (SEE REMARKS).  NOT REQUIRED IF IND=7.
C                F      - INPUT VECTOR OF LENGTH NF IF IND IS 5 OR 6.
C                           IF IND=5, F CONTAINS THE POWER VECTOR.
C                           IF IND=6, F CONTAINS THE PIVOT LOADINGS.
C                           (SEE REMARKS).  ON OUTPUT, F(I) CONTAINS
C                           THE VARIANCE ACCOUNTED FOR BY FACTOR I.
C                X      - NV BY NF TARGET MATRIX OF THE TRANSFORMATION.
C                           IF IND=7, X IS AN INPUT MATRIX.  IF IND IS 5
C                           OR 6, OUTPUT X CONTAINS THE TARGET USED.
C                IX     - INPUT ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                B      - OUTPUT NV BY NF PRIMARY STRUCTURE MATRIX.
C                IB     - INPUT ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                T      - OUTPUT NF BY NF TRANSFORMATION MATRIX.
C                           IF II IS NON-ZERO, T CONTAINS THE
C                           IMAGE TRANSFORMATION MATRIX.
C                IT     - INPUT ROW DIMENSION OF MATRIX T EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                S      - OUTPUT VECTOR OF LENGTH (NF+1)*NF/2 CONTAINING
C                           THE NF BY NF PRIMARY FACTOR CORRELATION
C                           MATRIX IN SYMMETRIC STORAGE MODE.
C                WK     - WORK VECTOR OF LENGTH (NF+1)*NF/2.  IF II IS
C                           NONZERO, WK SHOULD BE OF LENGTH (NV*NF) +
C                           (NF*(NF+1))/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT MATRIX A IS NOT
C                             OF FULL COLUMN RANK NF.
C                           IER = 130 INDICATES THAT MATRIX A IS NOT
C                             OF FULL COLUMN RANK NF. ONLY OUTPUT
C                             MATRIX T IS AFFECTED.
C                           IER = 131 INDICATES ONE OF IA, IX, IT, IB,
C                             NV, NF, OR IND WAS SPECIFIED INCORRECTLY.
C                         WARNING ERROR (WITH FIX)
C                           IER = 68 INDICATES THAT MAXIT WAS EXCEEDED.
C                             THE APPROXIMATE SOLUTION WAS USED.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1P,LINV1P,LUDECP,LUELMP,OFIMA3,
C                           OFROTA,UERTST,UGETIO,VIPRFF,VMULFF,VMULFM,
C                           VMULFS,VTPROF
C                       - DOUBLE/LEQT1P,LINV1P,LUDECP,LUELMP,OFIMA3,
C                           OFROTA,UERTST,UGETIO,VIPRFF,VMULFF,VMULFM,
C                           VMULFS,VTPROF,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  INPUT ARGUMENTS W, EPS, DELTA, AND NORM ARE REQUIRED
C                FOR IND=5 OR 6 ONLY. THEY ARE INPUT ARGUMENTS TO IMSL
C                ROUTINE OFROTA TO OBTAIN THE INITIAL ORTHOGONAL SOLU-
C                TION (WHICH MAY HAVE ALREADY BEEN EXAMINED TO DETER-
C                MINE INPUT VECTOR F FOR IND=6). W MAY BE SET AS
C                FOLLOWS.
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
C            2.  FOR IND=5, F(J) SHOULD BE GREATER THAN 1.0, TYPICALLY
C                4.0. GENERALLY, THE LARGER THE VALUE OF F(J), THE
C                MORE OBLIQUE THE SOLUTION WILL BE.
C                FOR IND=6, F(J) SHOULD BE IN THE INTERVAL (0.,1.).
C            3.  THE TARGET MATRIX IS AN HYPOTHESIZED ROTATED FACTOR
C                LOADING MATRIX BASED ON PRIOR KNOWLEDGE, WITH LOADINGS
C                CHOSEN TO ENHANCE INTERPRETABILITY. FOR IND=7, THE
C                ELEMENTS OF THE INPUT MATRIX X MUST BE LESS THAN ONE
C                IN ABSOLUTE VALUE. THE WEIGHTS X(I,J) SHOULD BE CHOSEN
C                BASED ON PRIOR KNOWLEDGE OF OR INTUITION REGARDING THE
C                STRUCTURE OF THE VARIABLES. A SIMPLE STRUCTURE SOLU-
C                TION WILL HAVE MOST OF THE WEIGHTS X(I,J) EITHER ZERO
C                OR LARGE IN MAGNITUDE. NOTE THAT THE TWO OPTIONS IND=5
C                AND IND=6 ATTEMPT TO ACHIEVE THIS SIMPLE STRUCTURE
C                BASED ON AN INITIAL ORTHOGONAL SOLUTION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFPROT  (A,IA,NV,NF,IND,NORM,II,MAXIT,W,EPS,DELTA,
     1                   F,X,IX,B,IB,T,IT,S,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NV,NF,IND,NORM,II,MAXIT,IX,IB,IT,IER
      REAL               A(IA,NF),X(IX,NF),F(NF),B(IB,NF),T(IT,NF),
     1                   S(1),WK(1),W,EPS,DELTA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NF1,IR,J,I,JER
      REAL               ABSX,CONS,DD,D1,D2,ZERO
      DOUBLE PRECISION   TEMP
      DATA               CONS/-89.0E0/,ZERO/0.0E0/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (IA.GE.NV.AND.IX.GE.NV.AND.IT.GE.NF.AND.IB.GE.NV.AND.NV.GE.NF.A
     1ND.IND.GE.5.AND.IND.LE.7) GO TO 5
      IER = 131
      GO TO 9000
    5 IF (II.EQ.0) GO TO 20
      NF1 = ((NF+1)*NF)/2+1
      IR = NF1
      DO 15 J=1,NF
         DO 10 I=1,NV
            WK(IR) = A(I,J)
            IR = IR+1
   10    CONTINUE
   15 CONTINUE
   20 IF (IND.EQ.7) GO TO 40
C                                  PRELIMINARY SOLUTION IN X
C                                  NOTE BOTH B AND WK ARE WORK ARRAYS
      CALL OFROTA (A,IA,NV,NF,NORM,0,MAXIT,W,EPS,DELTA,X,IX,T,IT,WK,B,JE
     1R)
      IF (JER.EQ.66) JER = 68
C                                  GENERATE TARGET MATRIX
      DO 35 J=1,NF
         DD = F(J)
         DO 30 I=1,NV
            IF (X(I,J).EQ.ZERO) GO TO 30
            ABSX = ABS(X(I,J))
            D1 = DD*ALOG(ABSX)
            IF (IND.EQ.6) D1 = D1/ABSX
            IF (D1.LT.CONS) GO TO 25
            IF (IND.EQ.5) X(I,J) = SIGN(ABSX**DD,X(I,J))
            IF (IND.EQ.6) X(I,J) = SIGN(ABSX**(DD/ABSX),X(I,J))
            GO TO 30
   25       X(I,J) = ZERO
   30    CONTINUE
   35 CONTINUE
   40 CALL OFIMA3 (A,IA,X,IX,NV,NF,NF,T,IT,WK,IER)
      IF (IER.NE.0) GO TO 9000
      DO 55 J=1,NF
         TEMP = 0.0D0
         DO 45 I=1,NF
            TEMP = TEMP+DBLE(T(I,J))**2
   45    CONTINUE
         DD = 1.0D0/DSQRT(TEMP)
         DO 50 I=1,NF
            T(I,J) = T(I,J)*DD
   50    CONTINUE
   55 CONTINUE
C                                  T IS THE LAMBDA MATRIX TRANSFORMATION
C                                  TO THE REFERENCE STRUCTURE SOLUTION
      CALL VMULFF (A,T,NV,NF,NF,IA,IT,B,IB,IER)
C                                  B IS THE REFERENCE STRUCTURE SOLUTION
      CALL VTPROF (T,NF,NF,IT,S)
C                                  S IS CORRELATION OF REFERENCE VECTORS
      CALL LINV1P (S,NF,WK,0,D1,D2,IER)
      IF (IER.EQ.0) GO TO 60
      IER = 130
      GO TO 9000
C                                  CORRELATION OF PRIMARY FACTORS
   60 IR = 0
      DO 70 I=1,NF
         F(I) = SQRT(WK(IR+I))
         DD = F(I)
         DO 65 J=1,I
            IR = IR+1
            S(IR) = WK(IR)/(DD*F(J))
   65    CONTINUE
   70 CONTINUE
C                                  PRIMARY PATTERN SOLUTION
      DO 80 J=1,NF
         DD = F(J)
         DO 75 I=1,NV
            A(I,J) = B(I,J)*DD
   75    CONTINUE
   80 CONTINUE
      CALL VMULFS (A,S,NV,NF,IA,B,IB)
C                                  PRIMARY STRUCTURE SOLUTION
C                                  VARIANCE ACCOUNTED FOR BY FACTORS
      DO 90 J=1,NF
         TEMP = 0.0D0
         DO 85 I=1,NV
            TEMP = TEMP+DBLE(A(I,J))*B(I,J)
   85    CONTINUE
         F(J) = TEMP
   90 CONTINUE
      IF (II.EQ.0) GO TO 95
      CALL OFIMA3 (WK(NF1),NV,B,IB,NV,NF,NF,T,IT,WK,IER)
   95 IF (IND.NE.7) IER = MAX0(IER,JER)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFPROT)
 9005 RETURN
      END

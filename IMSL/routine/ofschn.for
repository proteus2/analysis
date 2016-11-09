C   IMSL ROUTINE NAME   - OFSCHN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ORTHOGONAL TRANSFORMATION OF THE FACTOR
C                           LOADING MATRIX USING A TARGET MATRIX
C
C   USAGE               - CALL OFSCHN (A,IA,NV,NF,II,X,IX,B,IB,
C                           T,IT,S,F,IS,WK,IER)
C
C   ARGUMENTS    A      - INPUT NV BY NF UNROTATED FACTOR LOADING
C                           MATRIX.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NV     - INPUT NUMBER OF VARIABLES.
C                NF     - INPUT NUMBER OF FACTORS.
C                II     - INPUT OPTION PARAMETER.  II=0 INDICATES AN
C                           IMAGE ANALYSIS IS NOT BEING PERFORMED.
C                           OTHERWISE, AN IMAGE ANALYSIS IS ASSUMED, AND
C                           T BECOMES THE IMAGE TRANSFORMATION MATRIX.
C                X      - ON INPUT, THE NV BY NF TARGET MATRIX OF THE
C                           ROTATION.  ON OUTPUT, X CONTAINS THE ERROR
C                           MATRIX, I.E., THE TARGET MINUS THE B MATRIX.
C                           SEE REMARKS.
C                IX     - INPUT ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                B      - OUTPUT NV BY NF ORTHOGONALLY ROTATED FACTOR
C                           LOADING MATRIX.
C                IB     - INPUT ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                T      - OUTPUT NF BY NF TRANSFORMATION MATRIX.
C                           IF II IS NON-ZERO, T CONTAINS THE
C                           IMAGE TRANSFORMATION MATRIX.
C                IT     - INPUT ROW DIMENSION OF MATRIX T EXACTLY AS
C                         SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                         CALLING PROGRAM.
C                S      - WORK VECTOR OF LENGTH (NF+3)*NF/2.
C                F      - OUTPUT VECTOR OF LENGTH NF.  F(I) CONTAINS
C                           THE VARIANCE ACCOUNTED FOR BY FACTOR I.
C                IS     - INTEGER WORK VECTOR OF LENGTH NF.
C                WK     - WORK VECTOR OF LENGTH 2*NF*NF.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ONE OF THE INPUT
C                             MATRICES A OR X WAS NOT OF FULL COLUMN
C                             RANK NF.
C                           IER = 130 INDICATES THAT AT LEAST ONE OF
C                             NV, NF, IA, IX, IT, OR IB WAS SPECIFIED
C                             INCORRECTLY.
C
C   REQD. IMSL ROUTINES - SINGLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,LEQT1P,
C                           LUDECP,LUELMP,OFIMA3,UERTST,UGETIO,VIPRFF,
C                           VMULFF,VMULFM,VMULFP,VSRTU,VTPROF
C                       - DOUBLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,LEQT1P,
C                           LUDECP,LUELMP,OFIMA3,UERTST,UGETIO,VIPRFF,
C                           VMULFF,VMULFM,VMULFP,VSRTUD,VTPROF,VXADD,
C                           VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE TARGET MATRIX IS AN HYPOTHESIZED ROTATED FACTOR
C                LOADING MATRIX BASED ON PRIOR KNOWLEDGE, WITH LOADINGS
C                CHOSEN TO ENHANCE INTERPRETABILITY. THE ELEMENTS OF
C                THE INPUT MATRIX X MUST BE LESS THAN ONE IN ABSOLUTE
C                VALUE. THE WEIGHTS X(I,J) SHOULD BE CHOSEN BASED ON
C                PRIOR KNOWLEDGE OF OR INTUITION REGARDING THE STRUC-
C                TURE OF THE VARIABLES. A SIMPLE STRUCTURE SOLUTION
C                WILL HAVE MOST OF THE WEIGHTS X(I,J) EITHER ZERO OR
C                LARGE IN MAGNITUDE.
C            2.  THIS SUBROUTINE MAY ALSO BE USED TO REFINE A SOLUTION
C                OBTAINED BY BLIND TRANSFORMATION IN IMSL ROUTINE
C                OFROTA, WITH THE TARGET MATRIX CLOSELY RESEMBLING THE
C                BLIND SOLUTION, MODIFIED TO HAVE A SIMPLE STRUCTURE
C                AS DESCRIBED ABOVE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFSCHN  (A,IA,NV,NF,II,X,IX,B,IB,T,IT,S,F,IS,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NV,NF,II,IX,IB,IT,IS(NF),IER
      REAL               A(IA,NF),X(IX,NF),T(IT,NF),B(IB,NF),
     *                   F(NF),WK(1),S(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NFP1,NF1,NF2,IR,I,J,K,L
      REAL               ZERO
      DOUBLE PRECISION   TEMP
      DATA               ZERO/0.0E0/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (IA.GE.NV.AND.IX.GE.NV.AND.IT.GE.NF.AND.IB.GE.NV.AND.NV.GE.NF)
     1GO TO 5
      IER = 130
      GO TO 9000
    5 IER = 0
      NFP1 = NF+1
      NF1 = (NFP1*NF)/2+1
      NF2 = NF*NF+1
      CALL VMULFM (A,X,NV,NF,NF,IA,IX,T,IT,IER)
      IR = 0
C                                  TRANSPOSE OF T INTO WK
      DO 15 I=1,NF
         IS(I) = NFP1-I
         DO 10 J=1,NF
            IR = IR+1
            WK(IR) = T(I,J)
   10    CONTINUE
   15 CONTINUE
C                                  CALCULATE T*(T-TRANSPOSE)
      CALL VTPROF (WK,NF,NF,NF,S)
      CALL EIGRS (S,NF,1,F,WK,NF,S(NF1),IER)
      IF (IER.NE.0) GO TO 20
C                                  REVERSE EIGENVECTOR ORDER
      CALL VSRTU  (WK,NF,NF,NF,0,IS,S)
C                                  CALCULATE (T-TRANSPOSE)*T
      CALL VTPROF (T,NF,NF,IT,S)
      CALL EIGRS (S,NF,1,F,WK(NF2),NF,S(NF1),IER)
      IF (IER.EQ.0) GO TO 25
   20 IER = 129
      GO TO 9000
   25 DO 30 I=1,NF
         IS(I) = NFP1-I
   30 CONTINUE
C                                  REVERSE EIGENVECTOR ORDER
      CALL VSRTU  (WK(NF2),NF,NF,NF,0,IS,S)
      CALL VMULFF (T,WK(NF2),NF,NF,NF,IT,NF,B,IB,IER)
      CALL VMULFM (B,WK,NF,NF,NF,IB,NF,T,IT,IER)
      K = 0
      DO 40 J=1,NF
         IF (T(J,J).GT.ZERO) GO TO 40
         DO 35 I=1,NF
            L = K+I
            WK(L) = -WK(L)
   35    CONTINUE
   40 K = K+NF
C                                  TRANSFORMATION MATRIX
      CALL VMULFP (WK,WK(NF2),NF,NF,NF,NF,NF,T,IT,IER)
C                                  TRANSFORM LOADING MATRIX
      CALL VMULFF (A,T,NV,NF,NF,IA,IT,B,IB,IER)
C                                  VARIANCE ACCOUNTED FOR BY FACTORS
C                                  ERROR MATRIX = TARGET - OUTPUT B
      DO 50 J=1,NF
         TEMP = 0.0D0
         DO 45 I=1,NV
            TEMP = TEMP+DBLE(B(I,J))**2
            X(I,J) = X(I,J)-B(I,J)
   45    CONTINUE
         F(J) = TEMP
   50 CONTINUE
      IF (II.EQ.0) GO TO 9005
      CALL OFIMA3 (A,IA,B,IB,NV,NF,NF,T,IT,WK,IER)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFSCHN)
 9005 RETURN
      END

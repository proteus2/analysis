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

C   IMSL ROUTINE NAME   - EIGZF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           THE SYSTEM A*X=LAMBDA*B*X WHERE A AND B ARE
C                           REAL MATRICES
C
C   USAGE               - CALL EIGZF (A,IA,B,IB,N,IJOB,ALFA,BETA,Z,IZ,
C                           WK,IER)
C
C   ARGUMENTS    A      - THE INPUT REAL GENERAL MATRIX OF ORDER N.
C                           INPUT A IS DESTROYED IF IJOB IS EQUAL
C                           TO 0 OR 1.
C                IA     - THE INPUT ROW DIMENSION OF MATRIX A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                B      - THE INPUT REAL GENERAL MATRIX OF ORDER N.
C                           INPUT B IS DESTROYED IF IJOB IS EQUAL
C                           TO 0 OR 1.
C                IB     - THE INPUT ROW DIMENSION OF MATRIX B EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                N      - THE INPUT ORDER OF THE MATRICES A AND B.
C                IJOB   - INPUT OPTION PARAMETER. WHEN
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
C                ALFA   - OUTPUT VECTORS OF LENGTH N.
C                BETA       ALFA IS TYPE COMPLEX AND BETA IS TYPE REAL.
C                           IF A AND B WERE SIMULTANEOUSLY REDUCED
C                           TO TRIANGULAR FORM BY UNITARY EQUIVALENCES,
C                           ALFA AND BETA WOULD CONTAIN THE DIAGONAL
C                           ELEMENTS OF THE RESULTING MATRICES. (SEE
C                           MOLER-STEWART REFERENCE).
C                         THE J-TH EIGENVALUE IS THE COMPLEX NUMBER
C                           GIVEN BY ALFA(J)/BETA(J).
C                         NOTE - THE ROUTINE TREATS ALFA AS A REAL
C                           VECTOR OF LENGTH 2*N. AN APPROPRIATE
C                           EQUIVALENCE STATEMENT MAY BE REQUIRED.
C                           SEE DOCUMENT EXAMPLE.
C                Z      - THE OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE ALFA(J)/BETA(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IZ MUST BE GREATER THAN
C                           OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB AS FOLLOWS,
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             2*N*N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQZTF FAILED
C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                             J+1,J+2,...,N HAVE BEEN COMPUTED COR-
C                             RECTLY. EIGENVALUES 1,...,J MAY BE
C                             INACCURATE. IF IJOB = 1 OR 2 EIGENVECTORS
C                             MAY BE INACCURATE. THE PERFORMANCE INDEX
C                             IS SET TO 1000.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB IS RESET
C                             TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             0, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS RESET TO 0.
C
C   REQD. IMSL ROUTINES - EQZQF,EQZTF,EQZVF,UERTST,UGETIO,VHSH2C,
C                           VHSH2R,VHSH3R
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGZF  (A,IA,B,IB,N,IJOB,ALFA,BETA,Z,IZ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,N,IJOB,IZ,IER
      REAL               A(IA,N),B(IB,N),ALFA(1),WK(N,1),Z(1),BETA(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,KER,J,I,IIZ,NPI,JA,IZ2,N2,IS,IG,IGZ,LW,LZ,
     *                   KZ
      REAL               ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,Z11,
     1                   BSUM,REPS,ZERO,ONE,THOUS,BNORM,SUMS,SUMJ,
     2                   EPSA,EPSB
      DATA               REPS/Z3C100000/
      DATA               ZERO/0.0/,ONE/1.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      KER = 0
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 5
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      KER = 66
      IJOB = 1
      GO TO 10
    5 IF (IJOB.EQ.0) GO TO 20
   10 IF (IZ.GE.N) GO TO 15
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      KER = 67
      IJOB = 0
   15 IF (IJOB.EQ.3) GO TO 85
   20 IF (IJOB.NE.2) GO TO 30
C                                  SAVE INPUT A AND B IF IJOB = 2
      DO 25 J=1,N
      DO 25 I=1,N
         WK(I,J) = A(I,J)
         WK(I,J+N) = B(I,J)
   25 CONTINUE
   30 CONTINUE
      IIZ = N
      IF (IJOB.EQ.0) IIZ = 1
      IF (IJOB.EQ.0.AND.N.EQ.1) Z11 = Z(1)
      CALL EQZQF (A,IA,B,IB,N,Z,IIZ)
      CALL EQZTF (A,IA,B,IB,N,EPSA,EPSB,Z,IIZ,JER)
      CALL EQZVF (A,IA,B,IB,N,EPSA,EPSB,ALFA(1),ALFA(N+1),BETA,Z,IIZ)
      IF (IJOB.EQ.0.AND.N.EQ.1) Z(1) = Z11
      IF (IJOB.LE.1) GO TO 40
C                                  MOVE ORIGINAL MATRICES BACK
C                                    TO A AND B
      DO 35 I=1,N
      DO 35 J=1,N
         A(I,J) = WK(I,J)
         B(I,J) = WK(I,J+N)
   35 CONTINUE
   40 CONTINUE
C                                  CONVERT ALFA TO COMPLEX FORMAT
      DO 45 I=1,N
         NPI = N+I
         WK(I,1) = ALFA(NPI)
   45 CONTINUE
      JA = N+N
      J = N
      DO 50 I=1,N
         ALFA(JA-1) = ALFA(J)
         ALFA(JA) = WK(J,1)
         JA = JA-2
         J = J-1
   50 CONTINUE
      IF (IJOB.EQ.0) GO TO 115
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
      IZ2 = IZ+IZ
      N2 = N+N
      J = N
   55 IF (J.LT.1) GO TO 80
      IF (ALFA(J+J).EQ.ZERO) GO TO 70
C                                  MOVE PAIR OF COMPLEX CONJUGATE
C                                    EIGENVECTORS
      IS = IZ2*(J-1)+1
      IG = N*(J-2)+1
      IGZ = IG+N
C                                  MOVE COMPLEX CONJUGATE EIGENVECTOR
      DO 60 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IGZ)
         IS = IS+2
         IG = IG+1
         IGZ = IGZ+1
   60 CONTINUE
C                                  MOVE COMPLEX EIGENVECTOR
      IS = IZ2*(J-2)+1
      IG = IS+IZ2
      DO 65 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IG+1)
         IS = IS+2
         IG = IG+2
   65 CONTINUE
      J = J-2
      GO TO 55
C                                  MOVE REAL EIGENVECTOR
   70 IS = IZ2*(J-1)+N2
      IG = N*J
      DO 75 I=1,N
         Z(IS-1) = Z(IG)
         Z(IS) = ZERO
         IS = IS-2
         IG = IG-1
   75 CONTINUE
      J = J-1
      GO TO 55
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
   80 IF (IJOB.LE.1) GO TO 115
      WK(1,1) = THOUS
      IF (JER.NE.0) GO TO 115
C                                  COMPUTE MAX-NORM OF A AND B
   85 ANORM = ZERO
      BNORM = ZERO
      DO 95 I=1,N
         ASUM = ZERO
         BSUM = ZERO
         DO 90 J=1,N
            ASUM = ASUM+ABS(A(I,J))
            BSUM = BSUM+ABS(B(I,J))
   90    CONTINUE
         ANORM = AMAX1(ANORM,ASUM)
         BNORM = AMAX1(BNORM,BSUM)
   95 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
      IF (BNORM.EQ.ZERO) BNORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      LW = 1
      IZ2 = IZ+IZ
      DO 110 J=1,N
         S = ZERO
         SUMZ = ZERO
         LZ = IZ2*(J-1)+1
         DO 105 L=1,N
            SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
            KZ = IZ2*(J-1)+1
            SUMR = ZERO
            SUMS = ZERO
            SUMI = ZERO
            SUMJ = ZERO
            DO 100 K=1,N
               SUMR = SUMR+A(L,K)*Z(KZ)
               SUMI = SUMI+A(L,K)*Z(KZ+1)
               SUMS = SUMS+B(L,K)*Z(KZ)
               SUMJ = SUMJ+B(L,K)*Z(KZ+1)
               KZ = KZ+2
  100       CONTINUE
            SUMR = BETA(J)*SUMR-ALFA(LW)*SUMS+ALFA(LW+1)*SUMJ
            SUMI = BETA(J)*SUMI-ALFA(LW)*SUMJ-ALFA(LW+1)*SUMS
            S = AMAX1(S,CABS(CMPLX(SUMR,SUMI)))
            LZ = LZ+2
  105    CONTINUE
         SUMR = CABS(CMPLX(ALFA(LW),ALFA(LW+1)))*BNORM
         SUMR = SUMZ*(ABS(BETA(J))*ANORM+SUMR)
         IF(SUMR .NE. ZERO) PI = AMAX1(PI,S/SUMR)
         LW = LW+2
  110 CONTINUE
      PI = PI/REPS
      WK(1,1) = PI
  115 CONTINUE
      IER = MAX0(KER,JER)
      IF (KER.NE.0) CALL UERTST (KER,6HEIGZF )
      IF (JER.NE.0) CALL UERTST (JER,6HEIGZF )
      RETURN
      END

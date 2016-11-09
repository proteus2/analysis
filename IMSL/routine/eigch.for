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
      DATA               RDELP/Z3C100000/
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

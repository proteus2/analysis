C   IMSL ROUTINE NAME   - EIGRS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A REAL SYMMETRIC MATRIX
C
C   USAGE               - CALL EIGRS (A,N,JOBN,D,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - INPUT REAL SYMMETRIC MATRIX OF ORDER N,
C                           WHOSE EIGENVALUES AND EIGENVECTORS
C                           ARE TO BE COMPUTED. INPUT A IS
C                           DESTROYED IF IJOB IS EQUAL TO 0 OR 1.
C                N      - INPUT ORDER OF THE MATRIX A.
C                JOBN   - INPUT OPTION PARAMETER.  IF JOBN.GE.10
C                         A IS ASSUMED TO BE IN FULL STORAGE MODE
C                         (IN THIS CASE, A MUST BE DIMENSIONED EXACTLY
C                         N BY N IN THE CALLING PROGRAM).
C                         IF JOBN.LT.10 THEN A IS ASSUMED TO BE IN
C                         SYMMETRIC STORAGE MODE.  DEFINE
C                         IJOB=MOD(JOBN,10).  THEN WHEN
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
C                D      - OUTPUT VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES OF A.
C                Z      - OUTPUT N BY N MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE D(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             N(N+1)/2+N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQRT2S FAILED
C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                             AND EIGENVECTORS 1,...,J-1 HAVE BEEN
C                             COMPUTED CORRECTLY, BUT THE EIGENVALUES
C                             ARE UNORDERED. THE PERFORMANCE INDEX
C                             IS SET TO 1000.0
C                         WARNING ERROR (WITH FIX)
C                           IN THE FOLLOWING, IJOB = MOD(JOBN,10).
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS SET TO ZERO.
C
C   REQD. IMSL ROUTINES - EHOBKS,EHOUSS,EQRT2S,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGRS  (A,N,JOBN,D,Z,IZ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,JOBN,IZ,IER
      REAL               A(1),D(1),WK(1),Z(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IJOB,IR,JR,IJ,JI,NP1
      INTEGER            JER,NA,ND,IIZ,IBEG,IL,KK,LK,I,J,K,L
      REAL               ANORM,ASUM,PI,SUMZ,SUMR,AN,S,TEN,RDELP,ZERO,
     1                   ONE,THOUS
!     DATA               RDELP/Z3C100000/
      DATA               RDELP/Z'34000000'/
      DATA               ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IF (JOBN.LT.10) GO TO 15
C                                  CONVERT TO SYMMETRIC STORAGE MODE
      K = 1
      JI = N-1
      IJ = 1
      DO 10 J=1,N
         DO 5 I=1,J
            A(K) = A(IJ)
            IJ = IJ+1
            K = K+1
    5    CONTINUE
         IJ = IJ + JI
         JI = JI - 1
   10 CONTINUE
   15 IJOB = MOD(JOBN,10)
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 20
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 25
   20 IF (IJOB.EQ.0) GO TO 35
   25 IF (IZ.GE.N) GO TO 30
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   30 IF (IJOB.EQ.3) GO TO 75
   35 NA = (N*(N+1))/2
      IF (IJOB.NE.2) GO TO 45
      DO 40 I=1,NA
         WK(I) = A(I)
   40 CONTINUE
C                                  SAVE INPUT A IF IJOB = 2
   45 ND = 1
      IF (IJOB.EQ.2) ND = NA+1
C                                  REDUCE A TO SYMMETRIC TRIDIAGONAL
C                                    FORM
      CALL EHOUSS (A,N,D,WK(ND),WK(ND))
      IIZ = 1
      IF (IJOB.EQ.0) GO TO 60
      IIZ = IZ
C                                  SET Z TO THE IDENTITY MATRIX
      DO 55 I=1,N
         DO 50 J=1,N
            Z(I,J) = ZERO
   50    CONTINUE
         Z(I,I) = ONE
   55 CONTINUE
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS
   60 CALL EQRT2S (D,WK(ND),N,Z,IIZ,JER)
      IF (IJOB.EQ.0) GO TO 9000
      IF (JER.GT.128) GO TO 65
C                                  BACK TRANSFORM EIGENVECTORS
      CALL EHOBKS (A,N,1,N,Z,IZ)
   65 IF (IJOB.LE.1) GO TO 9000
C                                  MOVE INPUT MATRIX BACK TO A
      DO 70 I=1,NA
         A(I) = WK(I)
   70 CONTINUE
      WK(1) = THOUS
      IF (JER.NE.0) GO TO 9000
C                                  COMPUTE 1 - NORM OF A
   75 ANORM = ZERO
      IBEG = 1
      DO 85 I=1,N
         ASUM = ZERO
         IL = IBEG
         KK = 1
         DO 80 L=1,N
            ASUM = ASUM+ABS(A(IL))
            IF (L.GE.I) KK = L
            IL = IL+KK
   80    CONTINUE
         ANORM = AMAX1(ANORM,ASUM)
         IBEG = IBEG+I
   85 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 100 I=1,N
         IBEG = 1
         S = ZERO
         SUMZ = ZERO
         DO 95 L=1,N
            LK = IBEG
            KK = 1
            SUMZ = SUMZ+ABS(Z(L,I))
            SUMR = -D(I)*Z(L,I)
            DO 90 K=1,N
               SUMR = SUMR+A(LK)*Z(K,I)
               IF (K.GE.L) KK = K
               LK = LK+KK
   90       CONTINUE
            S = S+ABS(SUMR)
            IBEG = IBEG+L
   95    CONTINUE
         IF (SUMZ.EQ.ZERO) GO TO 100
         PI = AMAX1(PI,S/SUMZ)
  100 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
      IF (JOBN.LT.10) GO TO 9000
C                                  CONVERT BACK TO FULL STORAGE MODE
      NP1 = N+1
      IJ = (N-1)*NP1 + 2
      K = (N*(NP1))/2
      DO 110 JR=1,N
         J = NP1-JR
         DO 105 IR=1,J
            IJ = IJ-1
            A(IJ) = A(K)
            K = K-1
  105    CONTINUE
         IJ = IJ-JR
  110 CONTINUE
      JI = 0
      K = N-1
      DO 120 I=1,N
         IJ = I-N
         DO 115 J=1,I
            IJ = IJ+N
            JI = JI+1
            A(IJ) = A(JI)
  115    CONTINUE
         JI = JI + K
         K = K-1
  120 CONTINUE
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST (IER,6HEIGRS )
      IF (JER.EQ.0) GO TO 9005
      IER = JER
      CALL UERTST (IER,6HEIGRS )
 9005 RETURN
      END


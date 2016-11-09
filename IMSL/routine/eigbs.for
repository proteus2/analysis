C   IMSL ROUTINE NAME   - EIGBS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - FIND SOME EIGENVALUES AND (OPTIONALLY)
C                           EIGENVECTORS OF A REAL SYMMETRIC
C                           BAND MATRIX
C
C   USAGE               - CALL EIGBS (A,N,IA,IJOB,NC,M,D,Z,IZ,WORK,IER)
C
C   ARGUMENTS    A      - INPUT MATRIX WHOSE EIGENVALUES ARE TO
C                           BE DETERMINED.  A IS ASSUMED TO BE STORED
C                           IN BAND SYMMETRIC STORAGE MODE AND THEREFORE
C                           HAS DIMENSION N BY (NC+1). A IS DESTROYED
C                           BY EIGBS ON OUTPUT IF ABS(IJOB) = 1.
C                N      - INPUT ORDER OF THE MATRIX A
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IJOB   - INPUT OPTION PARAMETER. WHEN
C                           IJOB=-1, COMPUTE THE M SMALLEST EIGENVALUES.
C                           IJOB=-2, COMPUTE THE M SMALLEST EIGENVALUES
C                             AND CORRESPONDING EIGENVECTORS.
C                           IJOB=-3, COMPUTE THE M SMALLEST EIGENVALUES
C                             AND CORRESPONDING EIGENVECTORS AND THE
C                             PERFORMANCE INDEX.
C                           IJOB=1, COMPUTE THE M LARGEST EIGENVALUES
C                           IJOB=2, COMPUTE THE M LARGEST EIGENVALUES
C                             AND CORRESPONDING EIGENVECTORS.
C                           IJOB=3, COMPUTE THE M LARGEST EIGENVALUES
C                             AND CORRESPONDING EIGENVECTORS AND THE
C                             PERFORMANCE INDEX.
C                               IF THE PERFORMANCE INDEX IS COMPUTED, IT
C                             IS RETURNED IN WORK(1). THE ROUTINES HAVE
C                             PERFORMED (WELL, SATISFACTORILY, POORLY)
C                             IF WORK(1) IS (LESS THAN 1, BETWEEN 1 AND
C                             100, GREATER THAN 100).
C                NC     - INPUT NUMBER OF UPPER OR LOWER CODIAGONALS OF
C                           MATRIX A.
C                M      - INPUT NUMBER OF EIGENVALUES DESIRED
C                D      - OUTPUT VECTOR OF LENGTH AT LEAST M CONTAINING
C                           THE M LARGEST OR SMALLEST EIGENVALUES.
C                Z      - OUTPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           THE EIGENVECTORS. THE EIGENVECTOR CORRES-
C                           PONDING TO EIGENVALUE D(I) WILL BE PLACED
C                           IN COLUMN I OF Z. IF ABS(IJOB)=1, Z IS
C                           NOT USED.
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                WORK   - WORKSPACE VECTOR WITH DIMENSION=
C                           IF ABS(IJOB)=1, AT LEAST 3*N
C                           IF ABS(IJOB)=2 OR 3, AT LEAST N*(3*NC+6)
C                IER    - ERROR PARAMETER (OUTPUT)
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS OUT OF RANGE.
C                             IJOB EQUAL TO 1 (OR -1) IS USED.
C                           IER = 67, INDICATES ABS(IJOB) = 2 OR 3
C                             AND IZ IS LESS THAN THE ORDER OF MATRIX
C                             A. IJOB EQUAL TO 1 (OR -1) IS USED.
C                         TERMINAL ERROR
C                           IER = 129 IMPLIES THAT SOME EIGENVECTORS
C                             WERE NOT CALCULATED ACCEPTABLY. THE
C                             COLUMNS OF Z CORRESPONDING TO THOSE
C                             EIGENVECTORS ARE SET TO ZERO.
C
C   REQD. IMSL ROUTINES - EBNDR,EBNDV,EQRT1S,UERTST,UGETIO,VMULQF
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGBS  (A,N,IA,IJOB,NC,M,D,Z,IZ,WORK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IJOB,IZ,NC,M,IER
      REAL               A(IA,1),D(1),Z(IZ,1),WORK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,MB,I,LIM,IJ,L,K,IAJ,JER,N4,
     1                   KK,II,JBEG,JEND,JJ,KN,KSTRT
      REAL               E12,ASIGN,PI,ANORM,ZERO,ASUM,ONE,
     1                   S,SUMZ,SUMR,AN,TEN,RDELP,PIJ
      DATA               ZERO/0.0E0/,ONE/1.0E0/,TEN/10.0E0/
!     DATA               RDELP/Z3C100000/
      DATA               RDELP/Z'34000000'/
C                                  FIRST EXECUTABLE STATEMENT
      IAJ = IABS(IJOB)
      IER = 0
      JER = 0
      N4 = 4*N
      MB = NC+1
      IF (IAJ.GE.1.AND.IAJ.LE.3) GO TO 5
      IER = 66
      IAJ = 1
      GO TO 25
    5 IF (IAJ.EQ.1) GO TO 25
      IF (IZ.GE.N) GO TO 10
      IER = 67
      IAJ = 1
      GO TO 25
C                                  IF EIGENVECTORS TO BE COMPUTED
C                                  SAVE A, SINCE EBNDR DESTROYS A.
   10 DO 20 J=1,MB
         LIM = MB-J+1
         IJ = N4+LIM+(J-1)*N-1
         DO 15 I=LIM,N
            IJ = IJ+1
            WORK(IJ) = A(I,J)
   15    CONTINUE
   20 CONTINUE
C                                  REDUCE A TO SYMMETRIC TRIDIAGONAL FOR
   25 CALL EBNDR (A,IA,N,NC,0,WORK(1),WORK(N+1),WORK(2*N+1),Z,IZ)
      IF (IJOB.LT.0) GO TO 35
      DO 30 I=1,N
         WORK(I) = -WORK(I)
   30 CONTINUE
C                                  FIND M SMALLEST (OR LARGEST) EIGENVAL
   35 CALL EQRT1S (WORK(1),WORK(2*N+1),N,M,0,JER)
      E12 = 0.0
      ASIGN = 1.0
      IF (IJOB.LT.0) GO TO 40
      E12 = 2.0
      ASIGN = -1.0
   40 DO 45 I=1,M
         D(I) = ASIGN*WORK(I)
   45 CONTINUE
      IF (IAJ.EQ.1) GO TO 9000
C                                  USE INVERSE ITERATION TO FIND
C                                  EIGENVECTORS, IF DESIRED.
      CALL EBNDV (WORK(4*N+1),N,N,NC,E12,M,D,Z,IZ,WORK(IJ+1),WORK(3*N
     1+1),IER)
C                                  PERFORMANCE INDEX DESIRED
C                                    RETRIEVE INPUT MATRIX A
      L = NC
      IF (L.EQ.0) L = 1
      DO 55 K=1,L
         DO 50 KN=1,L
            A(K,KN) = ZERO
   50    CONTINUE
   55 CONTINUE
      DO 65 J=1,MB
         LIM = MB-J+1
         IJ = N4+LIM+(J-1)*N-1
         DO 60 I=LIM,N
            IJ = IJ+1
            A(I,J) = WORK(IJ)
   60    CONTINUE
   65 CONTINUE
      IF (IER.EQ.0) GO TO 70
      IER = 129
   70 IF (IAJ.NE.3) GO TO 9000
C                                  COMPUTE 1 NORM OF A
      ANORM = ZERO
      JBEG = MB
      DO 85 I=1,N
         JEND = JBEG+NC+NC+1
         JJ = JBEG
         ASUM = ZERO
         II = I
         DO 75 J=JBEG,JEND
            ASUM = ASUM+ABS(A(II,JJ))
            JJ = JJ+1
            IF (J.LT.MB) GO TO 75
            JJ = JJ-2
            IF (JJ.LT.1) GO TO 80
            II = II+1
            IF (II.GT.N) GO TO 80
   75    CONTINUE
   80    JBEG = JBEG-1
         ANORM = AMAX1(ANORM,ASUM)
         IF (JBEG.LE.1) JBEG = 1
   85 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 95 J=1,M
         CALL VMULQF (A,N,NC,IA,Z(1,J),1,N,WORK,N)
         SUMZ = ZERO
         SUMR = ZERO
         DO 90 I=1,N
            SUMR = SUMR+ABS(WORK(I)-D(J)*Z(I,J))
            SUMZ = SUMZ+ABS(Z(I,J))
   90    CONTINUE
         IF (SUMZ.EQ.ZERO) GO TO 95
         PIJ = SUMR/SUMZ
         IF (PIJ.GT.PI) PI = PIJ
   95 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WORK(1) = PI
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST (IER,6HEIGBS )
 9005 RETURN
      END


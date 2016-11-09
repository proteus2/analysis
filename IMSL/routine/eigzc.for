C   IMSL ROUTINE NAME   - EIGZC
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM/SINGLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           THE SYSTEM A*X=LAMBDA*B*X WHERE A AND B ARE
C                           COMPLEX MATRICES
C
C   USAGE               - CALL EIGZC (A,IA,B,IB,N,IJOB,EIGA,EIGB,Z,IZ,
C                           WK,INFER,IER)
C
C   ARGUMENTS    A      - THE INPUT COMPLEX GENERAL MATRIX OF ORDER N.
C                           INPUT A IS DESTROYED IF IJOB IS EQUAL
C                           TO 0 OR 1.
C                IA     - THE INPUT ROW DIMENSION OF MATRIX A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                B      - THE INPUT COMPLEX GENERAL MATRIX OF ORDER N.
C                           INPUT B IS DESTROYED IF IJOB IS EQUAL
C                           TO 0 OR 1.
C                IB     - THE INPUT ROW DIMENSION OF MATRIX B EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                N      - THE INPUT ORDER OF THE MATRICES A AND B.
C                IJOB   - THE INPUT OPTION PARAMETER. WHEN
C                           IJOB = 0, COMPUTE EIGENVALUES ONLY
C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
C                             VECTORS.
C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
C                             AND PERFORMANCE INDEX.
C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
C                           RETURNED IN THE REAL PART OF WK(1,1). THE
C                           ROUTINES HAVE PERFORMED (WELL,SATISFACTOR-
C                           ILY,POORLY) IF REAL(WK(1,1)) IS (LESS THAN
C                           1, BETWEEN 1 AND 100, GREATER THAN 100).
C                EIGA   - OUTPUT COMPLEX VECTORS OF LENGTH N.
C                EIGB       EIGA CONTAINS THE DIAGONAL OF THE TRIANGU-
C                           LARIZED A, AND EIGB CONTAINS THE DIAGONAL
C                           OF THE TRIANGULARIZED B.
C                         THE J-TH EIGENVALUE IS THE COMPLEX NUMBER
C                           GIVEN BY EIGA(J)/EIGB(J) WHEN EIGB(J) IS
C                           NOT EQUAL TO ZERO. IF EIGB(J) IS ZERO, THEN
C                           THE J-TH EIGENVALUE IS REGARDED AS BEING
C                           INFINITE.
C                Z      - THE OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE EIGA(J)/EIGB(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IZ MUST BE GREATER THAN
C                           OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
C                WK     - A COMPLEX MATRIX USED AS WORK SPACE.
C                           THE SIZE OF WK DEPENDS ON THE VALUE OF
C                           IJOB AS FOLLOWS,
C                           IJOB = 0, WK IS NOT USED.
C                           IJOB = 1, WK IS NOT USED.
C                           IJOB = 2, WK IS AT LEAST N BY 2*N.
C                           IJOB = 3, WK IS AT LEAST 1 BY 1.
C                INFER  - OUTPUT INTEGER VALUE. INFER = J INDICATES
C                           THAT ELZVC FAILED TO CONVERGE ON EIGENVALUE
C                           J WITHIN 30 ITERATIONS. EIGENVALUES J+1,
C                           J+2,...,N HAVE BEEN COMPUTED CORRECTLY.
C                           EIGENVALUES 1,2,...,J MAY BE INACCURATE. IF
C                           IJOB = 1 OR 2, EIGENVECTORS MAY BE
C                           INACCURATE. THE PERFORMANCE INDEX IS SET TO
C                           1000. J = 0 INDICATES THAT ELZVC CONVERGED
C                           FOR ALL N EIGENVALUES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, INDICATES THAT ELZVC FAILED
C                             TO CONVERGE ON EIGENVALUE J. SEE DESCRIP-
C                             TION OF INFER ABOVE.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB IS RESET
C                             TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             0, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS RESET TO 0.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - ELZHC,ELZVC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGZC  (A,IA,B,IB,N,IJOB,EIGA,EIGB,Z,IZ,WK,INFER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,N,IJOB,IZ,INFER,IER
      COMPLEX            A(IA,N),B(IB,N),EIGA(N),WK(N,1),Z(IZ,1),
     *                   EIGB(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JJOB,JER,KER,I,J,L,K
      REAL               ANORM,ASUM,PI,S,SUMR2,ZNORM,
     1                   BSUM,REPS,ZERO,ONE,THOUS,BNORM
      COMPLEX            CZERO,SUMR,SUMS
      DATA               CZERO/(0.0,0.0)/
      DATA               REPS/Z3C100000/
      DATA               ZERO/0.0/,ONE/1.0/,THOUS/1000.0/
C                                  FIRST EXECUTABLE STATEMENT
      JJOB = 0
C                                  INITIALIZE ERROR PARAMETERS
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
   15 IF (IJOB.EQ.3) GO TO 50
      JJOB = 1
   20 IF (IJOB.NE.2) GO TO 30
C                                  SAVE INPUT A AND B IF IJOB = 2
      DO 25 J=1,N
      DO 25 I=1,N
         WK(I,J) = A(I,J)
         WK(I,J+N) = B(I,J)
   25 CONTINUE
   30 CONTINUE
      CALL ELZHC (N,A,IA,B,IB,JJOB,Z,IZ)
      CALL ELZVC (N,A,IA,B,IB,Z,IZ,JJOB,EIGA,EIGB,INFER,JER)
   40 IF (IJOB.LE.1) GO TO 85
C                                  MOVE ORIGINAL MATRICES BACK
C                                    TO A AND B
      DO 45 J=1,N
      DO 45 I=1,N
         A(I,J) = WK(I,J)
         B(I,J) = WK(I,J+N)
   45 CONTINUE
      WK(1,1) = CMPLX(THOUS,ZERO)
      IF (JER.NE.0) GO TO 85
C                                  COMPUTE 1-NORM OF A AND B
   50 ANORM = ZERO
      BNORM = ZERO
      DO 60 J=1,N
         ASUM = ZERO
         BSUM = ZERO
         DO 55 I=1,N
            ASUM = ASUM+CABS(A(I,J))
            BSUM = BSUM+CABS(B(I,J))
   55    CONTINUE
         ANORM = AMAX1(ANORM,ASUM)
         BNORM = AMAX1(BNORM,BSUM)
   60 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
      IF (BNORM.EQ.ZERO) BNORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 80 J=1,N
         ZNORM = ZERO
         DO 65 I=1,N
            ZNORM = AMAX1(CABS(Z(I,J)),ZNORM)
   65    CONTINUE
         S = ZERO
         DO 75 L=1,N
            SUMR = CZERO
            SUMS = CZERO
            DO 70 K=1,N
               SUMR = SUMR+A(L,K)*Z(K,J)
               SUMS = SUMS+B(L,K)*Z(K,J)
   70       CONTINUE
            SUMR = EIGB(J)*SUMR-EIGA(J)*SUMS
            S = AMAX1(S,CABS(SUMR))
   75    CONTINUE
         SUMR2 = CABS(EIGA(J))*BNORM*ZNORM
         SUMR2 = CABS(EIGB(J))*ANORM*ZNORM+SUMR2
         PI = AMAX1(PI,S/SUMR2)
   80 CONTINUE
      PI = PI/REPS
      WK(1,1) = CMPLX(PI,ZERO)
   85 CONTINUE
      IER = MAX0(KER,JER)
 9000 CONTINUE
      IF (KER.NE.0) CALL UERTST (KER,6HEIGZC )
      IF (JER.NE.0) CALL UERTST (JER,6HEIGZC )
 9005 RETURN
      END
 
R; T=0.06/0.46 23:03:38
LL UERTST (KER,6HEIGZC )
      IF (JER.NE.0) CALL UERTST (JER,6HEIGZC )
 9005 RETURN
      END
 
R

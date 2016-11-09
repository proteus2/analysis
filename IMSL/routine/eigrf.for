C   IMSL ROUTINE NAME   - EIGRF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A REAL GENERAL MATRIX IN FULL STORAGE MODE
C
C   USAGE               - CALL EIGRF (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - THE INPUT REAL GENERAL MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                N      - THE INPUT ORDER OF THE MATRIX A.
C                IA     - THE INPUT ROW DIMENSION OF MATRIX A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IJOB   - THE INPUT OPTION PARAMETER. WHEN
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
C                W      - THE OUTPUT COMPLEX VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES OF A.
C                         NOTE - THE ROUTINE TREATS W AS A REAL VECTOR
C                           OF LENGTH 2*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                Z      - THE OUTPUT N BY N COMPLEX MATRIX CONTAINING
C                           THE EIGENVECTORS OF A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE W(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. IZ MUST BE GREATER
C                           THAN OR EQUAL TO N IF IJOB IS NOT EQUAL TO
C                           ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             (2+N)N.
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT EQRH3F FAILED
C                           TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                           J+1,J+2,...,N HAVE BEEN COMPUTED CORRECTLY.
C                           EIGENVALUES 1,...,J ARE SET TO ZERO.
C                           IF IJOB = 1 OR 2 EIGENVECTORS ARE SET TO
C                           ZERO. THE PERFORMANCE INDEX IS SET TO 1000.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS SET TO ZERO.
C
C   REQD. IMSL ROUTINES - EBALAF,EBBCKF,EHBCKF,EHESSF,EQRH3F,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGRF  (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IJOB,IZ,IER
      REAL               A(IA,1),WK(N,1),W(1),Z(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,IZ2,K,L,I,N1,N2,II,JJ,NP1,IIZ,NPI,JW,J,
     *                   IS,IG,IGZ,LW,LLZ,KKZ,LZ,KZ
      REAL               ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP,
     *                   ZERO,ONE,THOUS,AN,Z11
      DATA               RDELP/Z3C100000/
      DATA               ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IZ2 = IZ+IZ
      IF (IJOB .GE. 0 .AND. IJOB .LE. 3) GO TO 5
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 10
    5 IF (IJOB .EQ. 0) GO TO 16
   10 IF (IZ .GE. N) GO TO 15
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   15 IF (IJOB .EQ. 3) GO TO 95
C                                  PACK A INTO AN N BY N ARRAY
   16 K = 1
      L = 1
      DO 20 J=1,N
         DO 20 I=1,N
            A(K,L) = A(I,J)
C                                  SAVE INPUT A IF IJOB = 2
            IF (IJOB .EQ. 2) WK(I,J)=A(I,J)
            K = K+1
            IF (K .GT. IA) K = 1
            IF (K .EQ. 1) L = L+1
   20 CONTINUE
      N1 = 1
      IF (IJOB .EQ. 2) N1 = N+1
      N2 = N1+1
      IF (IJOB .EQ. 0) N2 = 1
C                                  BALANCE THE INPUT A
      CALL EBALAF (A,N,N,WK(1,N1),K,L)
      IF (IJOB .EQ. 0 .AND. L .EQ. 0) GO TO 35
C                                  IF L = 0, A IS ALREADY IN HESSENBERG
C                                    FORM
      CALL EHESSF (A,K,L,N,N,WK(1,N2))
      IF (IJOB .EQ. 0) GO TO 35
C                                  SET Z IDENTITY MATRIX
      II = 1
      JJ = 1
      NP1 = N+1
      DO 30 I=1,N
         DO 25 J=1,N
            Z(II) = ZERO
            II = II+1
   25    CONTINUE
         Z(JJ) = ONE
         JJ = JJ+NP1
   30 CONTINUE
      CALL EHBCKF (Z,A,WK(1,N2),N,N,N,K,L)
      IIZ = N
   35 IF (IJOB .EQ. 0) IIZ = 1
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z11 = Z(1)
      CALL EQRH3F (A,N,N,K,L,W(1),W(N+1),Z,IIZ,JER)
      IF(IJOB .EQ. 0 .AND. N .EQ. 1) Z(1) = Z11
      IF (JER .GT. 128 .OR. IJOB .EQ. 0) GO TO 40
      CALL EBBCKF (WK(1,N1),Z,K,L,N,N,N)
C                                  CONVERT W (EIGENVALUES) TO COMPLEX
C                                    FORMAT
   40 DO 45 I=1,N
         NPI = N+I
         WK(I,N1) = W(NPI)
   45 CONTINUE
      JW = N+N
      J = N
      DO 50 I=1,N
         W(JW-1) = W(J)
         W(JW) = WK(J,N1)
         JW = JW-2
         J = J-1
   50 CONTINUE
      IF (IJOB .EQ. 0) GO TO 9000
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
      J = N
   60 IF (J .LT. 1) GO TO 85
      IF (W(J+J) .EQ. ZERO) GO TO 75
C                                  MOVE PAIR OF COMPLEX CONJUGATE
C                                    EIGENVECTORS
      IS = IZ2*(J-1)+1
      IG = N*(J-2)+1
      IGZ = IG+N
C                                  MOVE COMPLEX CONJUGATE EIGENVECTOR
      DO 65 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IGZ)
         IS = IS+2
         IG = IG+1
         IGZ = IGZ+1
   65 CONTINUE
C                                  MOVE COMPLEX EIGENVECTOR
      IS = IZ2*(J-2)+1
      IG = IS+IZ2
      DO 70 I=1,N
         Z(IS) = Z(IG)
         Z(IS+1) = -Z(IG+1)
         IS = IS+2
         IG = IG+2
   70 CONTINUE
      J = J-2
      GO TO 60
C                                  MOVE REAL EIGENVECTOR
   75 IS = IZ2*(J-1)+N+N
      IG = N*J
      DO 80 I=1,N
         Z(IS-1) = Z(IG)
         Z(IS) = ZERO
         IS = IS-2
         IG = IG-1
   80 CONTINUE
      J = J-1
      GO TO 60
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
C                                    NEXT, MOVE ORIGINAL MATRIX BACK
C                                    TO A
   85 IF (IJOB .LE. 1) GO TO 9000
      DO 90 I=1,N
         DO 90 J=1,N
            A(I,J) = WK(I,J)
   90 CONTINUE
      WK(1,1) = THOUS
      IF (JER .NE. 0) GO TO 9000
C                                  COMPUTE 1-NORM OF A
   95 ANORM = ZERO
      DO 105 J=1,N
         ASUM = ZERO
         DO 100 I=1,N
            ASUM = ASUM+ABS(A(I,J))
  100    CONTINUE
         ANORM = AMAX1(ANORM,ASUM)
  105 CONTINUE
      IF (ANORM .EQ. ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      LW = 1
      LLZ = 0
      KKZ = 0
      DO 120 J=1,N
         S = ZERO
         SUMZ = ZERO
         LZ = LLZ+1
         KZ = KKZ+1
         LW = J+J-1
         DO 115 L=1,N
            SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
            KZ = KKZ+1
            SUMR = -W(LW)*Z(LZ)+W(LW+1)*Z(LZ+1)
            SUMI = -W(LW)*Z(LZ+1)-W(LW+1)*Z(LZ)
            DO 110 K=1,N
               SUMR =SUMR+A(L,K)*Z(KZ)
               SUMI = SUMI+A(L,K)*Z(KZ+1)
               KZ = KZ+2
  110       CONTINUE
            S = S+CABS(CMPLX(SUMR,SUMI))
            LZ = LZ+2
  115    CONTINUE
         PI = AMAX1(PI,S/SUMZ)
         KKZ = KKZ+IZ2
         LLZ = LLZ+IZ2
  120 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1,1) = PI
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST (IER,6HEIGRF )
      IF (JER .EQ. 0) GO TO 9005
      IER = JER
      CALL UERTST (IER,6HEIGRF )
 9005 RETURN
      END

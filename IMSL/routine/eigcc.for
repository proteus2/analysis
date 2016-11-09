C   IMSL ROUTINE NAME   - EIGCC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           A COMPLEX GENERAL MATRIX
C
C   USAGE               - CALL EIGCC (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - THE INPUT COMPLEX GENERAL MATRIX OF ORDER N
C                           WHOSE EIGENVALUES AND EIGENVECTORS ARE
C                           TO BE COMPUTED. INPUT A IS DESTROYED IF
C                           IJOB IS EQUAL TO 0 OR 1.
C                         NOTE - THE ROUTINE TREATS A AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                N      - THE INPUT ORDER OF THE MATRIX A AND MATRIX Z.
C                IA     - THE INPUT ROW DIMENSION OF THE MATRIX A
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
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
C                           THE EIGENVECTORS OF MATRIX A.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE W(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE THE DOCUMENT
C                           EXAMPLE.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IZ MUST BE GREATER THAN
C                           OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
C                           ON THE VALUE OF IJOB, WHEN
C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST 2N.
C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
C                             2(N*N+N).
C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 128+J, INDICATES THAT ELRH1C OR ELRH2C
C                             FAILED TO CONVERGE ON EIGENVALUE J.
C                           EIGENVALUES J+1,J+2,...,N HAVE BEEN
C                           COMPUTED CORRECTLY.
C                           THE PERFORMANCE INDEX IS SET TO 1000.0.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. IJOB IS SET TO ZERO.
C
C   REQD. IMSL ROUTINES - EBALAC,EBBCKC,EHESSC,ELRH1C,ELRH2C,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGCC  (A,N,IA,IJOB,W,Z,IZ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IJOB,IZ,IER
      REAL               A(1),W(1),Z(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,IAA,IZZ,NN,N1,N2,L,M,NM,KA1,K,I,INFER,NPI,
     *                   N2PI,JW,J,NP1,NM1,KA2,KA3,JB,JZ,JR,JI,KA,LW,
     *                   LZ,KZ
      REAL               ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,AN,RDELP,ZERO,
     *                   ONE,TEN,THOUS
      DATA               RDELP/Z3C100000/
      DATA               ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IAA = IA+IA
      IZZ = IZ+IZ
      IF (IJOB .GE. 0 .AND. IJOB .LE. 3) GO TO 5
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 10
    5 IF (IJOB .EQ. 0) GO TO 20
   10 IF (IZ .GE. N) GO TO 15
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   15 IF (IJOB .EQ. 3) GO TO 100
C                                  PACK A INTO AN N BY N ARRAY
C                                    SEPARATE A INTO REAL AND IMAGINARY
C                                    PARTS AND PACK INTO AN N BY N
C                                    MATRIX
   20 NN = N+N
      N1 = 1
      IF (IJOB .EQ. 2) N1 = NN*N+1
      N2 = N1+N
      IF (IJOB .EQ. 0) N2 = 1
      L = 1
      M = 1
      NM = N1-1
      KA1 = 2
      DO 45 J = 1,N
         K = KA1
         DO 30 I = 1,N
            IF (IJOB .NE. 2) GO TO 25
            WK(M) = A(K-1)
            WK(M+1) = A(K)
            M = M+2
   25       WK(NM+I) = A(K)
            A(L) = A(K-1)
            K = K+2
            L = L+1
   30    CONTINUE
         DO 40 I = 1,N
   35       A(L) = WK(NM+I)
            L = L+1
   40    CONTINUE
         KA1 = KA1+IAA
   45 CONTINUE
C                                  BALANCE THE INPUT A
      CALL EBALAC (A(1),A(N+1),N,NN,K,L,WK(N1))
C                                  REDUCES A TO HESSENBERG FORM
      CALL EHESSC (A(1),A(N+1),K,L,N,NN,WK(N2))
      IF (IJOB .NE. 0) GO TO 50
C                                  COMPUTE EIGENVALUES ONLY
      CALL ELRH1C (A(1),A(N+1),K,L,N,NN,W(1),W(N+1),INFER,JER)
      GO TO 55
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS
   50 CALL ELRH2C (A(1),A(N+1),K,L,N,NN,W(1),W(N+1),Z(1),Z(N+1),
     1   WK(N2),INFER,JER)
C                                  BACKTRANSFORM THE EIGENVECTORS
      CALL EBBCKC (Z(1),Z(N+1),N,NN,K,L,N,WK(N1))
   55 N2 = N2-1
      DO 60 I=1,N
         NPI = N+I
         N2PI = N2+I
         WK(N2PI) = W(NPI)
   60 CONTINUE
      JW = N+N
      J = N
      DO 65 I=1,N
         W(JW-1) = W(J)
         N2PJ = N2+J
         W(JW) = WK(N2PJ)
         JW = JW-2
         J = J-1
   65 CONTINUE
      IF (IJOB .EQ. 0) GO TO 9000
      IF (IJOB .NE. 2) GO TO 80
C                                  MOVE ORIGINAL MATRIX BACK TO A
      K = 1
      KA1 = 1
      DO 75 J = 1,N
         L = KA1
         DO 70 I = 1,N
            A(L) = WK(K)
            A(L+1) = WK(K+1)
            K = K+2
            L = L+2
   70    CONTINUE
         KA1 = KA1+IAA
   75 CONTINUE
C                                  CONVERT Z (EIGENVECTORS) TO COMPLEX
C                                    FORMAT Z(IZ,N)
   80 NP1 = N+1
      NM1 = N-1
      KA1 = NN*N-NM1
      KA2 = NN+IZZ*(N-1)
      KA3 = NN*N-N
      DO 95 JB = 1,N
         J = NP1-JB
         K = KA1
         DO 85 I = 1,N
            WK(I) = Z(K)
            K = K+1
   85    CONTINUE
         JZ = KA2
         JR = KA3
         JI = N
         DO 90 I = 1,N
            Z(JZ-1) = Z(JR)
            Z(JZ) = WK(JI)
            JZ = JZ-2
            JR = JR-1
            JI = JI-1
   90    CONTINUE
         KA1 = KA1-NN
         KA2 = KA2-IZZ
         KA3 = KA3-NN
   95 CONTINUE
C                                  Z IS NOW IN COMPLEX FORMAT Z(IZ,N).
      IF (IJOB .LE. 1) GO TO 9000
      WK(1) = THOUS
      IF (JER .NE. 0) GO TO 9000
C                                  COMPUTE 1-NORM OF A
  100 ANORM = ZERO
      KA1 = 1
      DO 110 J=1,N
         ASUM = ZERO
         KA = KA1
         DO 105 I=1,N
            ASUM = ASUM+CABS(CMPLX(A(KA),A(KA+1)))
            KA = KA+2
  105    CONTINUE
         ANORM = AMAX1(ANORM,ASUM)
         KA1 = KA1+IAA
  110 CONTINUE
      IF (ANORM .EQ. ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      LW = 1
      KA1 = 1
      DO 125 J=1,N
         S = ZERO
         SUMZ = ZERO
         LZ = KA1
         DO 120 L=1,N
            SUMR = ZERO
            SUMI = ZERO
            SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
            KA = L+L-1
            KZ = KA1
            DO 115 K=1,N
               SUMR =SUMR+A(KA)*Z(KZ)-A(KA+1)*Z(KZ+1)
               SUMI = SUMI+A(KA)*Z(KZ+1)+A(KA+1)*Z(KZ)
               KA = KA+IAA
               KZ = KZ+2
  115       CONTINUE
            SUMR = SUMR-W(LW)*Z(LZ)+W(LW+1)*Z(LZ+1)
            SUMI = SUMI-W(LW)*Z(LZ+1)-W(LW+1)*Z(LZ)
            S = S+CABS(CMPLX(SUMR,SUMI))
            LZ = LZ+2
  120    CONTINUE
         PI = AMAX1(PI,S/SUMZ)
         LW = LW+2
         KA1 = KA1+IZZ
  125 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST (IER,6HEIGCC )
      IF (JER .EQ. 0) GO TO 9005
      IER = JER+INFER
      CALL UERTST (IER,6HEIGCC )
 9005 RETURN
      END

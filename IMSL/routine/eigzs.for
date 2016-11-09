C   IMSL ROUTINE NAME   - EIGZS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
C                           THE SYSTEM A*X=LAMBDA*B*X WHERE A AND B ARE
C                           REAL SYMMETRIC MATRICES AND B IS POSITIVE
C                           DEFINITE
C
C   USAGE               - CALL EIGZS (A,B,N,IJOB,D,Z,IZ,WK,IER)
C
C   ARGUMENTS    A      - AN INPUT REAL SYMMETRIC MATRIX OF ORDER N,
C                           STORED IN SYMMETRIC STORAGE MODE.
C                B      - AN INPUT REAL SYMMETRIC MATRIX OF ORDER N,
C                           STORED IN SYMMETRIC STORAGE MODE.
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
C                D      - OUTPUT REAL VECTOR OF LENGTH N,
C                           CONTAINING THE EIGENVALUES.
C                Z      - THE OUTPUT N BY N REAL MATRIX CONTAINING
C                           THE EIGENVECTORS.
C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
C                           PONDS TO THE EIGENVALUE D(J).
C                           IF IJOB = 0, Z IS NOT USED.
C                IZ     - THE INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. IZ MUST BE GREATER THAN
C                           OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
C                WK     - WORK AREA OF LENGTH N*(N+2)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, INDICATES THAT B IS NOT
C                             POSITIVE DEFINITE.
C                           IER = 129+J, INDICATES THAT EQRT2S FAILED
C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES
C                             J+1,J+2,...,N HAVE BEEN COMPUTED COR-
C                             RECTLY. EIGENVALUES 1,...,J MAY BE
C                             INACCURATE. IF IJOB = 1 OR 2 EIGENVECTORS
C                             MAY BE INACCURATE. THE PERFORMANCE INDEX
C                             IS SET TO 1000.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
C                             IJOB IS GREATER THAN 3. CONTINUES AS IF
C                             IJOB IS 1.
C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
C                             0, AND IZ IS LESS THAN THE ORDER OF
C                             MATRIX A. CONTINUES AS IF IJOB IS 0.
C
C   REQD. IMSL ROUTINES - SINGLE/EIGRS,EHOBKS,EHOUSS,EQRT2S,EREDU,
C                           VMULSF,VNRMS1,VABSMS,UGETIO,UERTST
C                       - DOUBLE/EIGRS,EHOBKS,EHOUSS,EQRT2S,EREDU,
C                           VMULSF,VNRMS1,VABSMS,UGETIO,UERTST,
C                           VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIGZS (A,B,N,IJOB,D,Z,IZ,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IJOB,IZ,IER
      REAL               A(1),B(1),D(1),Z(IZ,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,III,II,IIJOB,IJOB1,I,JER,J,KI,K,NPI,NSYM
      REAL               ANORM,BNORM,ONE,PI,REPS,SN,THOUS,X,ZERO,ZJ
      DATA               REPS/Z3C100000/
      DATA               ZERO /0.0/,ONE /1.0/,THOUS /1000.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IIJOB = IJOB
      IF (IIJOB.GE.0 .AND. IIJOB.LE.3) GO TO 5
      IIJOB = 1
      IER = 66
    5 IF (IIJOB.EQ.0 .OR. IZ.GE.N) GO TO 10
      IIJOB = 0
      IER = 67
   10 IF (IIJOB.EQ.3) GO TO 35
      NSYM = N*(N+1)/2
      CALL EREDU (N,A,B,WK (NSYM+1),WK,JER)
      IER = MAX0(IER,JER)
      IF (IER.GE.129) GO TO 9000
      IJOB1 = MIN0(1,IIJOB)
      CALL EIGRS (WK (NSYM+1),N,IJOB1,D,Z,IZ,WK (2*NSYM+1),JER)
      IF (JER.GE.129) JER = JER+1
      IER = MAX0(IER,JER)
      IF (IER.GE.129) WK(1) = THOUS
      IF (IER.GE.129) GO TO 9000
      IF (IIJOB.EQ.0) GO TO 9000
      DO 30 J=1,N
         DO 25 II=1,N
            I = N+1-II
            I1 = I+1
            X = Z(I,J)
            IF (I.EQ.N) GO TO 20
            DO 15 K=I1,N
               KI = K*(K-1)/2+I
               X = X-WK(KI)*Z(K,J)
   15       CONTINUE
   20       III = I*(I+1)/2
            Z(I,J) = X/WK(III)
   25    CONTINUE
   30 CONTINUE
      IF (IIJOB.EQ.1) GO TO 9000
C                                  COMPUTE MAX-NORM OF A AND B
   35 CALL VNRMS1 (A,N,ANORM)
      CALL VNRMS1 (B,N,BNORM)
      IF (ANORM.EQ.ZERO) ANORM = ONE
      IF (BNORM.EQ.ZERO) BNORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 45 J=1,N
         CALL VMULSF (A,N,Z (1,J),1,N,WK (1),N)
         CALL VMULSF (B,N,Z (1,J),1,N,WK (N+1),N)
         ZJ = ZERO
         SN = ZERO
         DO 40 I=1,N
            ZJ = ZJ+ABS(Z(I,J))
            NPI = N+I
            SN = AMAX1(SN,ABS(WK(I)-D(J)*WK(NPI)))
   40    CONTINUE
         PI = AMAX1(PI,SN/(ANORM*ZJ+ABS(D(J))*BNORM*ZJ))
   45 CONTINUE
      PI = PI/REPS
      WK(1) = PI
 9000 IF (IER.GT.0) CALL UERTST (IER,6HEIGZS )
      RETURN
      END

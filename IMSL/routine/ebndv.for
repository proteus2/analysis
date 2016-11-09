C   IMSL ROUTINE NAME   - EBNDV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE EIGBS
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBNDV  (A,IA,N,NC,ARDR,M,W,Z,IZ,W1,W2,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,N,NC,M,IZ,IER
      REAL               A(IA,1),ARDR,W(M),Z(IZ,1),W1(1),W2(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,R,II,IJ,JJ,KJ,MB,M1,IJ1,ITS,KJ1,M21,
     1                   MAXJ,MAXK,GROUP,MBW
      REAL               U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,REPS
C
C                                  REPS IS A MACHINE DEPENDENT CONSTANT
C                                    SPECIFYING THE RELATIVE PRECISION
C                                    OF FLOATING POINT ARITHMETIC
!     DATA               REPS/Z3C100000/
      DATA               REPS/Z'34000000'/
C
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (M.EQ.0) GO TO 9005
      MBW = NC+1
      MB = MBW
      IF (ARDR.LT.0.0) MB = (MBW+1)/2
      M1 = MB-1
      M21 = M1+MB
      ORDER = 1.0-ABS(ARDR)
C                                  FIND VECTORS BY INVERSE ITERATION
      DO 180 R=1,M
         ITS = 1
         X1 = W(R)
         IF (R.NE.1) GO TO 20
C                                  COMPUTE NORM OF MATRIX
         NORM = 0.0
         DO 10 J=1,MB
            JJ = MB+1-J
            KJ = JJ+M1
            IJ = 1
            DO 5 I=JJ,N
               NORM = NORM+ABS(A(I,J))
               IF (ARDR.GE.0.0) GO TO 5
               NORM = NORM+ABS(A(IJ,KJ))
               IJ = IJ+1
    5       CONTINUE
   10    CONTINUE
         IF (ARDR.LT.0.0) NORM = 0.5*NORM
C                                  EPS2 IS THE CRITERION FOR GROUPING,
C                                  EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                                  ROOTS ARE MODIFIED BY EPS3,
C                                  EPS4 IS TAKEN VERY SMALL TO AVOID
C                                  OVERFLOW
         IF (NORM.EQ.0.0) NORM = 1.0
         EPS2 = 1.0E-3*NORM*ABS(ORDER)
         EPS3 = REPS*NORM
         UK = N
         UK = SQRT(UK)
         EPS4 = UK*EPS3
   15    GROUP = 0
         GO TO 25
C                                  LOOK FOR CLOSE OR COINCIDENT ROOTS
   20    IF (ABS(X1-X0).GE.EPS2) GO TO 15
         GROUP = GROUP+1
         IF (ORDER*(X1-X0).LE.0.0) X1 = X0+ORDER*EPS3
C                                  EXPAND MATRIX, SUBTRACT EIGENVALUE,
C                                  AND INITIALIZE VECTOR
   25    DO 55 I=1,N
            IJ = I+MIN0(0,I-M1)*N
            KJ = IJ+MB*N
            IJ1 = KJ+M1*N
            IF (M1.EQ.0) GO TO 50
            DO 45 J=1,M1
               IF (IJ.GT.M1) GO TO 30
               IF (IJ.GT.0) GO TO 35
               W1(IJ1) = 0.0
               IJ1 = IJ1+N
               GO TO 35
   30          W1(IJ) = A(I,J)
   35          IJ = IJ+N
               II = I+J
               IF (II.GT.N) GO TO 45
               JJ = MB-J
               IF (ARDR.GE.0.0) GO TO 40
               II = I
               JJ = MB+J
   40          W1(KJ) = A(II,JJ)
               KJ = KJ+N
   45       CONTINUE
   50       W1(IJ) = A(I,MB)-X1
            W2(I) = EPS4
            IF (ORDER.EQ.0.0) W2(I) = Z(I,R)
   55    CONTINUE
         IF (M1.EQ.0) GO TO 100
C                                  ELIMINATION WITH INTERCHANGES
         DO 95 I=1,N
            II = I+1
            MAXK = MIN0(I+M1-1,N)
            MAXJ = MIN0(N-I,M21-2)*N
            DO 65 K=I,MAXK
               KJ1 = K
               J = KJ1+N
               JJ = J+MAXJ
               DO 60 KJ=J,JJ,N
                  W1(KJ1) = W1(KJ)
                  KJ1 = KJ
   60          CONTINUE
               W1(KJ1) = 0.0
   65       CONTINUE
            IF (I.EQ.N) GO TO 95
            U = 0.0
            MAXK = MIN0(I+M1,N)
            MAXJ = MIN0(N-II,M21-2)*N
            DO 70 J=I,MAXK
               IF (ABS(W1(J)).LT.ABS(U)) GO TO 70
               U = W1(J)
               K = J
   70       CONTINUE
            J = I+N
            JJ = J+MAXJ
            IF (K.EQ.I) GO TO 80
            KJ = K
            DO 75 IJ=I,JJ,N
               V = W1(IJ)
               W1(IJ) = W1(KJ)
               W1(KJ) = V
               KJ = KJ+N
   75       CONTINUE
            IF (ORDER.NE.0.0) GO TO 80
            V = W2(I)
            W2(I) = W2(K)
            W2(K) = V
   80       IF (U.EQ.0.0) GO TO 95
            DO 90 K=II,MAXK
               V = W1(K)/U
               KJ = K
               DO 85 IJ=J,JJ,N
                  KJ = KJ+N
                  W1(KJ) = W1(KJ)-V*W1(IJ)
   85          CONTINUE
               IF (ORDER.EQ.0.0) W2(K) = W2(K)-V*W2(I)
   90       CONTINUE
   95    CONTINUE
C                                  BACK SUBSTITUTION
C                                  FOR I=N STEP -1 UNTIL 1 DO --
  100    DO 120 II=1,N
            I = N+1-II
            MAXJ = MIN0(II,M21)
            IF (MAXJ.EQ.1) GO TO 110
            IJ1 = I
            J = IJ1+N
            JJ = J+(MAXJ-2)*N
            DO 105 IJ=J,JJ,N
               IJ1 = IJ1+1
               W2(I) = W2(I)-W1(IJ)*W2(IJ1)
  105       CONTINUE
  110       V = W1(I)
            IF (ABS(V).GE.EPS3) GO TO 115
C                                  SET ERROR -- NEARLY SINGULAR
C                                    LINEAR SYSTEM
            IF (ORDER.EQ.0.0) IER = 130
            V = SIGN(EPS3,V)
  115       W2(I) = W2(I)/V
  120    CONTINUE
         XU = 1.0
         IF (ORDER.EQ.0.0) GO TO 170
C                                  ORTHOGONALIZE WITH RESPECT TO
C                                    PREVIOUS MEMBERS OF GROUP
         IF (GROUP.EQ.0) GO TO 140
         DO 135 JJ=1,GROUP
            J = R-GROUP-1+JJ
            XU = 0.0
            DO 125 I=1,N
  125       XU = XU+W2(I)*Z(I,J)
            DO 130 I=1,N
  130       W2(I) = W2(I)-XU*Z(I,J)
  135    CONTINUE
  140    NORM = 0.0
         DO 145 I=1,N
  145    NORM = NORM+ABS(W2(I))
         IF (NORM.GE.1.0E-1) GO TO 160
C                                  IN-LINE PROCEDURE FOR CHOOSING
C                                  A NEW STARTING VECTOR
         IF (ITS.GE.N) GO TO 155
         ITS = ITS+1
         XU = EPS4/(UK+1.0)
         W2(1) = EPS4
         DO 150 I=2,N
  150    W2(I) = XU
         W2(ITS) = W2(ITS)-EPS4*UK
         GO TO 100
C                                  SET ERROR -- NON-CONVERGED
C                                    EIGENVECTOR
  155    IER = 129
         XU = 0.0
         GO TO 170
C                                  NORMALIZE SO THAT SUM OF SQUARES IS
C                                  1 AND EXPAND TO FULL ORDER
  160    U = 0.0
         DO 165 I=1,N
  165    U = U+W2(I)**2
         XU = 1.0/SQRT(U)
  170    DO 175 I=1,N
  175    Z(I,R) = W2(I)*XU
         X0 = X1
  180 CONTINUE
 9000 CONTINUE
      IF(IER.NE.0) CALL UERTST (IER,6HEBNDV )
 9005 RETURN
      END


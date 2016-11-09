C   IMSL ROUTINE NAME   - LSVDF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - SINGULAR VALUE DECOMPOSITION OF A REAL
C                           MATRIX
C
C   USAGE               - CALL LSVDF (A,IA,M,N,B,IB,NB,S,WK,IER)
C
C   ARGUMENTS    A      - REAL M BY N MATRIX. (INPUT/OUTPUT)
C                         ON INPUT, A CONTAINS THE MATRIX TO BE
C                           DECOMPOSED.
C                         ON OUTPUT, A CONTAINS THE N BY N MATRIX V
C                           IN ITS FIRST N ROWS. SEE REMARKS. EITHER
C                           M.GE.N OR M.LT.N IS PERMITTED.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                           A IS USED BY LSVDF AS WORK STORAGE FOR AN
C                           N BY N MATRIX. THEREFORE, IA MUST BE
C                           GREATER THAN OR EQUAL TO MAX(M,N).
C                M      - NUMBER OF ROWS IN A. (INPUT)
C                N      - NUMBER OF COLUMNS IN A. (INPUT)
C                B      - M BY NB MATRIX. (INPUT/OUTPUT)
C                           B IS NOT USED IF NB.LE.0. OTHERWISE, B IS
C                           REPLACED BY THE MATRIX PRODUCT U**(T) * B.
C                           SEE REMARKS.
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NB     - NUMBER OF COLUMNS IN B. (INPUT)
C                           IF NB.LE.0, B IS NOT USED.
C                S      - VECTOR OF LENGTH N. (OUTPUT)
C                         ON OUTPUT, S CONTAINS THE ORDERED SINGULAR
C                           VALUES OF A.  S(1) .GE. S(2),...,
C                           .GE. S(N) .GE. 0.
C                WK     - WORK VECTOR OF LENGTH 2N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT MATRIX A IS NOT
C                             FULL RANK OR VERY ILL-CONDITIONED. SMALL
C                             SINGULAR VALUES MAY NOT BE VERY ACCURATE.
C                           IER=34 INDICATES THAT EITHER N.LE.0 OR
C                             M.LE.0.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT CONVERGENCE WAS
C                             NOT OBTAINED BY LSVDB AND COMPUTATION
C                             WAS DISCONTINUED.
C
C   REQD. IMSL ROUTINES - SINGLE/LSVDB,LSVG2,VHS12,UERSET,UERTST,UGETIO
C                           VBLA=SROTG
C                       - DOUBLE/LSVDB,LSVG2,VHS12,UERSET,UERTST,UGETIO
C                           VBLA=DROTG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  LSVDF COMPUTES THE SINGULAR VALUE DECOMPOSITION OF
C                A REAL M BY N MATRIX
C                     A = U * Q * V**(T)  WHERE
C                U IS AN M BY M ORTHOGONAL MATRIX,
C                V IS AN N BY N ORTHOGONAL MATRIX, AND
C                Q IS AN M BY N MATRIX WITH ALL ELEMENTS ZERO EXCEPT
C                     Q(I,I) = S(I) I=1,...,MIN(M,N).
C                V IS RETURNED IN THE FIRST N ROWS OF A.
C                U IS OBTAINED BY SETTING B TO THE M BY M IDENTITY
C                MATRIX, ON INPUT, AND SETTING NB=M. ON OUTPUT, B IS
C                REPLACED BY U**(T).
C            2.  THE NOTATION U**(T) AND V**(T) REPRESENTS U
C                TRANSPOSE AND V TRANSPOSE, RESPECTIVELY. Q**(+)
C                DENOTES THE GENERALIZED INVERSE OF Q.
C            3.  LSVDF IS USEFUL IN ANALYZING AND SOLVING THE LEAST
C                SQUARES PROBLEM A*X.APPR.B (WHERE .APPR. MEANS
C                APPROXIMATELY EQUALS). IN THIS CASE B IS A VECTOR OF
C                LENGTH M AND LSVDF IS CALLED WITH IB=M, NB=1. THE
C                SOLUTION IS X=V*Q**(+)*U**(T)*B. U**(T)*B REPLACES
C                B ON OUTPUT. THE SOLUTION X IS OBTAINED AS FOLLOWS...
C                (THE USER MAY WISH TO SET SMALL SINGULAR VALUES SUCH AS
C                S(I) TO ZERO JUST PRIOR TO COMPUTING X. SEE REFERENCE
C                FOR DETAILS.)
C
C                C                  COMPUTE Q**(+) * U**(T) * B
C                      L=MIN0(M,N)
C                      DO 10 I=1,L
C                         T=0.0
C                         IF (S(I).NE.0.0) T=B(I)/S(I)
C                         B(I)=T
C                   10 CONTINUE
C                C                  COMPUTE V * Q**(+) * U**(T) * B
C                   15 DO 25 I=1,N
C                         X(I)=0.0
C                         DO 20 J=1,L
C                   20    X(I)=X(I)+A(I,J)*B(J)
C                   25 CONTINUE
C                IF B IS SET TO THE J-TH COLUMN OF THE M BY M IDENTITY
C                MATRIX ON INPUT, X IS THE J=TH COLUMN OF A**(+) (THE
C                GENERALIZED INVERSE OF A).
C            4.  THE USER SHOULD BE AWARE OF SEVERAL PRACTICAL ASPECTS
C                OF THE SINGULAR VALUE ANALYSIS. ONE OF THESE IS THE
C                EFFECT OF THE UNCERTAINTY OF THE DATA AND THE NEED FOR
C                SCALING. SEE THE LAWSON-HANSON REFERENCE PAGES 180-198
C                FOR DETAILS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LSVDF  (A,IA,M,N,B,IB,NB,S,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,M,N,IB,NB,IER
      REAL               A(IA,N),B(IB,1),S(N),WK(N,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JP1,K,L,MM,NN,NNP1,NS,NSP1
      REAL               ZERO,ONE,T
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  BEGIN SPECIAL FOR ZERO ROWS AND
C                                    COLS. PACK THE NONZERO COLS TO THE
C                                    LEFT
      NN = N
      IER = 34
      IF (NN.LE.0.OR.M.LE.0) GO TO 9000
      IER = 0
      J = NN
    5 CONTINUE
      DO 10 I=1,M
         IF (A(I,J).NE.ZERO) GO TO 25
   10 CONTINUE
C                                  COL J IS ZERO. EXCHANGE IT WITH COL
C                                    N
      IF (J.EQ.NN) GO TO 20
      DO 15 I=1,M
   15 A(I,J) = A(I,NN)
   20 CONTINUE
      A(1,NN) = J
      NN = NN-1
   25 CONTINUE
      J = J-1
      IF (J.GE.1) GO TO 5
C                                  IF N=0 THEN A IS ENTIRELY ZERO AND
C                                    SVD COMPUTATION CAN BE SKIPPED
      NS = 0
      IF (NN.EQ.0) GO TO 120
C                                  PACK NONZERO ROWS TO THE TOP QUIT
C                                    PACKING IF FIND N NONZERO ROWS
      I = 1
      MM = M
   30 IF (I.GT.N.OR.I.GE.MM) GO TO 75
      IF (A(I,I).NE.ZERO) GO TO 40
      DO 35 J=1,NN
         IF (A(I,J).NE.ZERO) GO TO 40
   35 CONTINUE
      GO TO 45
   40 I = I+1
      GO TO 30
C                                  ROW I IS ZERO EXCHANGE ROWS I AND M
   45 IF (NB.LE.0) GO TO 55
      DO 50 J=1,NB
         T = B(I,J)
         B(I,J) = B(MM,J)
         B(MM,J) = T
   50 CONTINUE
   55 DO 60 J=1,NN
   60 A(I,J) = A(MM,J)
      IF (MM.GT.NN) GO TO 70
      DO 65 J=1,NN
   65 A(MM,J) = ZERO
   70 CONTINUE
C                                  EXCHANGE IS FINISHED
      MM = MM-1
      GO TO 30
C
   75 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
C                                  BEGIN SVD ALGORITHM..
C                                  (1) REDUCE THE MATRIX TO UPPER
C                                    BIDIAGONAL FORM WITH HOUSEHOLDER
C                                    TRANSFORMATIONS.
C                                    H(N)...H(1)AQ(1)...Q(N-2) =
C                                    (D**T,0)**T WHERE D IS UPPER
C                                    BIDIAGONAL.
C                                  (2) APPLY H(N)...H(1) TO B. HERE
C                                    H(N)...H(1)*B REPLACES B IN
C                                    STORAGE.
C                                  (3) THE MATRIX PRODUCT W=
C                                    Q(1)...Q(N-2) OVERWRITES THE FIRST
C                                    N ROWS OF A IN STORAGE.
C                                  (4) AN SVD FOR D IS COMPUTED. HERE K
C                                    ROTATIONS RI AND PI ARE COMPUTED
C                                    SO THAT RK...R1*D*P1**(T)...PK**(T)
C                                    = DIAG(S1,...,SM) TO WORKING
C                                    ACCURACY. THE SI ARE NONNEGATIVE
C                                    AND NONINCREASING. HERE RK...R1*B
C                                    OVERWRITES B IN STORAGE WHILE
C                                    A*P1**(T)...PK**(T) OVERWRITES A
C                                    IN STORAGE.
C                                  (5) IT FOLLOWS THAT,WITH THE PROPER
C                                    DEFINITIONS, U**(T)*B OVERWRITES
C                                    B, WHILE V OVERWRITES THE FIRST N
C                                    ROW AND COLUMNS OF A.
      L = MIN0(MM,NN)
C                                  THE FOLLOWING LOOP REDUCES A TO
C                                    UPPER BIDIAGONAL AND ALSO APPLIES
C                                    THE PREMULTIPLYING TRANSFORMATIONS
C                                    TO B.
      DO 85 J=1,L
         IF (J.GE.MM) GO TO 80
         JP1 = MIN0(J+1,NN)
         CALL VHS12 (1,J,J+1,MM,A(1,J),1,T,A(1,JP1),1,IA,NN-J)
         CALL VHS12 (2,J,J+1,MM,A(1,J),1,T,B,1,IB,NB)
   80    IF (J.GE.NN-1) GO TO 85
         CALL VHS12 (1,J+1,J+2,NN,A(J,1),IA,WK(J,2),A(J+1,1),IA,1,MM-J)
   85 CONTINUE
C                                  COPY THE BIDIAGONAL MATRIX INTO THE
C                                    ARRAY S FOR LSVDB
      IF (L.EQ.1) GO TO 95
      DO 90 J=2,L
         S(J) = A(J,J)
         WK(J,1) = A(J-1,J)
   90 CONTINUE
   95 S(1) = A(1,1)
C
      NS = NN
      IF (MM.GE.NN) GO TO 100
      NS = MM+1
      S(NS) = ZERO
      WK(NS,1) = A(MM,MM+1)
  100 CONTINUE
C                                  CONSTRUCT THE EXPLICIT N BY N
C                                    PRODUCT MATRIX, W=Q1*Q2*...*QL*I
C                                    IN THE ARRAY A
      DO 115 K=1,NN
         I = NN+1-K
         IF (I.GT.MIN0(MM,NN-2)) GO TO 105
         CALL VHS12 (2,I+1,I+2,NN,A(I,1),IA,WK(I,2),A(1,I+1),1,IA,NN-I)
  105    DO 110 J=1,NN
  110    A(I,J) = ZERO
         A(I,I) = ONE
  115 CONTINUE
C                                  COMPUTE THE SVD OF THE BIDIAGONAL
C                                    MATRIX
C
      LEVEL=1
      CALL UERSET(LEVEL,LEVOLD)
      CALL LSVDB (S(1),WK(1,1),NS,A,IA,NN,B,IB,NB,IER)
C                                  TEST FOR IER=33
C
      IF (IER.GT.128) GO TO 9000
      CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.NE.33) GO TO 120
      T=0.0
      NM=MIN0(M,N)
      IF (S(1).NE.ZERO) T=S(NM)/S(1)
      F=100.0+T
      IF (F.EQ.100.0) GO TO 120
      IER=0
  120 CONTINUE
      IF (NS .LT. MIN0(M,N)) IER = 33
      IF (NS.GE.NN) GO TO 130
      NSP1 = NS+1
      DO 125 J=NSP1,NN
  125 S(J) = ZERO
  130 CONTINUE
      IF (NN.EQ.N) GO TO 155
      NNP1 = NN+1
C                                  MOVE RECORD OF PERMUTATIONS AND
C                                    STORE ZEROS
      DO 140 J=NNP1,N
         S(J) = A(1,J)
         IF (NN.LT.1) GO TO 140
         DO 135 I=1,NN
  135    A(I,J) = ZERO
  140 CONTINUE
C                                  PERMUTE ROWS AND SET ZERO SINGULAR
C                                    VALUES
      DO 150 K=NNP1,N
         I = S(K)
         S(K) = ZERO
         DO 145 J=1,N
            A(K,J) = A(I,J)
  145    A(I,J) = ZERO
         A(I,K) = ONE
  150 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
  155 IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLSVDF )
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - LGING
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - NUCLEUS USED ONLY BY IMSL ROUTINE LGINF
C
C
C   REQD. IMSL ROUTINES - LSVDB,LSVG1,LSVG2,VHS12,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LGING (A,IA,M,N,B,IB,NB,S,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,M,N,IB,NB,IER
      REAL               A(IA,N),B(IB,1),S(N),WK(N,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JP1,K,L,LNM,MM,NN,NNP1,NS,NSP1,NM
      REAL               ZERO,ONE,T,F,FF
      DATA               ZERO /0.0/,ONE /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  BEGIN SPECIAL FOR ZERO ROWS AND
C                                    COLS. PACK THE NONZERO COLS TO THE
C                                    LEFT
      NN = N
      IER = 34
      IF (NN.LE.0 .OR. M.LE.0) GO TO 9005
      IER = 0
      J = NN
    5 CONTINUE
      DO 10 I=1,M
         IF (A(I,J).NE.ZERO) GO TO 25
   10 CONTINUE
C                                  COL J IS ZERO. EXCHANGE IT WITH COL
C                                    NN
      IF (J.EQ.NN) GO TO 20
      DO 15 I=1,M
   15 A(I,J) = A(I,NN)
   20 CONTINUE
      A(1,NN) = J
      NN = NN-1
   25 CONTINUE
      J = J-1
      IF (J.GE.1) GO TO 5
C                                  IF NN=0 THEN A IS ENTIRELY ZERO AND
C                                    SVD COMPUTATION CAN BE SKIPPED
      NS = 0
      IF (NN.EQ.0) GO TO 135
C                                  PACK NONZERO ROWS TO THE TOP QUIT
C                                    PACKING IF FIND N NONZERO ROWS
      DO 30 I=1,N
   30 WK(I,1) = I
      I = 1
      MM = M
   35 IF (I.GT.N .OR. I.GE.MM) GO TO 70
      IF (A(I,I).NE.ZERO) GO TO 45
      DO 40 J=1,NN
         IF (A(I,J).NE.ZERO) GO TO 45
   40 CONTINUE
      GO TO 50
   45 I = I+1
      GO TO 35
C                                  ROW I IS ZERO EXCHANGE ROWS I AND MM
C                                  AND RECORD THE PERMUTATION IN WK(I,1)
   50 WK(I,1) = MM
      DO 55 J=1,NN
   55 A(I,J) = A(MM,J)
      IF (MM.GT.NN) GO TO 65
      DO 60 J=1,NN
   60 A(MM,J) = ZERO
   65 CONTINUE
C                                  EXCHANGE IS FINISHED
      MM = MM-1
      GO TO 35
C
   70 CONTINUE
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
      DO 80 J=1,L
         IF (J.GE.MM) GO TO 75
         JP1 = MIN0(J+1,NN)
         CALL VHS12(1,J,J+1,MM,A(1,J),1,T,A(1,JP1),1,IA,NN-J)
         S(J) = T
   75    IF (J.GE.NN-1) GO TO 80
         CALL VHS12(1,J+1,J+2,NN,A(J,1),IA,WK(J,2),A(J+1,1),IA,1,MM-J)
   80 CONTINUE
C                                  CONSTRUCT THE FIRST N ROWS OF THE
C                                    MATRIX HL*...*H2*H1 = B
      DO 85 JJ=1,L
         J = L+1-JJ
         IF (J.GE.MM) GO TO 85
         T = S(J)
         CALL VHS12(2,J,J+1,MM,A(1,J),1,T,B,IB,1,N)
   85 CONTINUE
C                                  PERMUTE COLUMNS OF B ACCORDING TO THE
C                                  PERMUTATIONS MADE TO THE ROWS OF A
      LNM = MIN0(N,MM)
      IF (LNM.LT.1) GO TO 100
      DO 95 K=1,LNM
         J = WK(K,1)
         IF (J.EQ.K) GO TO 95
C                                  EXCHANGE COLUMNS J AND K
         DO 90 I=1,N
            T = B(I,K)
            B(I,K) = B(I,J)
            B(I,J) = T
   90    CONTINUE
   95 CONTINUE
  100 CONTINUE
C                                  COPY THE BIDIAGONAL MATRIX INTO THE
C                                    ARRAY S FOR LSVDB
      IF (L.EQ.1) GO TO 110
      DO 105 J=2,L
         S(J) = A(J,J)
         WK(J,1) = A(J-1,J)
  105 CONTINUE
  110 S(1) = A(1,1)
C
      NS = NN
      IF (MM.GE.NN) GO TO 115
      NS = MM+1
      S(NS) = ZERO
      WK(NS,1) = A(MM,MM+1)
  115 CONTINUE
C                                  CONSTRUCT THE EXPLICIT NN BY NN
C                                    PRODUCT MATRIX, W=Q1*Q2*...*QL*I
C                                    IN THE ARRAY A
      DO 130 K=1,NN
         I = NN+1-K
         IF (I.GT.MIN0(MM,NN-2)) GO TO 120
         CALL VHS12(2,I+1,I+2,NN,A(I,1),IA,WK(I,2),A(1,I+1),1,IA,NN-I)
  120    DO 125 J=1,NN
  125    A(I,J) = ZERO
         A(I,I) = ONE
  130 CONTINUE
C                                  COMPUTE THE SVD OF THE BIDIAGONAL
C                                    MATRIX
C
      LEVEL = 1
      CALL UERSET(LEVEL,LEVOLD)
      CALL LSVDB(S(1),WK(1,1),NS,A,IA,NN,B,IB,NB,IER)
C
      IF (IER.GT.128) GO TO 9005
      CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.NE.33) GO TO 135
      FF = 0.0
      NM = MIN0(N,M)
      IF (S(1).NE.ZERO) FF = S(NM)/S(1)
      F = 100.0+FF
      IF (F.EQ.100.0) GO TO 135
      IER = 0
  135 CONTINUE
      IF (NS.GE.NN) GO TO 145
      NSP1 = NS+1
      DO 140 J=NSP1,NN
  140 S(J) = ZERO
  145 CONTINUE
      IF (NN.EQ.N) GO TO 9005
      NNP1 = NN+1
C                                  MOVE RECORD OF PERMUTATIONS AND
C                                    STORE ZEROS
      DO 155 J=NNP1,N
         S(J) = A(1,J)
         IF (NN.LT.1) GO TO 155
         DO 150 I=1,NN
  150    A(I,J) = ZERO
  155 CONTINUE
C                                  PERMUTE ROWS AND SET ZERO SINGULAR
C                                    VALUES
      DO 165 K=NNP1,N
         I = S(K)
         S(K) = ZERO
         DO 160 J=1,N
            A(K,J) = A(I,J)
            A(I,J) = ZERO
  160    CONTINUE
         A(I,K) = ONE
  165 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
 9000 CONTINUE
 9005 RETURN
      END

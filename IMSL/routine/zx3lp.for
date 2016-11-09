C   IMSL ROUTINE NAME   - ZX3LP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - SOLVE THE LINEAR PROGRAMMING PROBLEM
C                           VIA THE REVISED SIMPLEX ALGORITHM -
C                           EASY TO USE VERSION
C
C   USAGE               - CALL ZX3LP (A,IA,B,C,N,M1,M2,S,PSOL,DSOL,RW,
C                           IW,IER)
C
C   ARGUMENTS    A      - MATRIX OF DIMENSION M1+M2+2 BY N
C                           CONTAINING THE COEFFICIENTS OF THE M1
C                           INEQUALITY CONSTRAINTS IN THE FIRST M1 ROWS
C                           FOLLOWED BY THE COEFFICIENTS OF THE M2
C                           EQUALITY CONSTRAINTS. (INPUT)
C                           THE LAST TWO ROWS OF A ARE USED ONLY
C                           AS WORKING STORAGE.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT) TWO ROWS OF A ARE REQUIRED
C                           FOR WORKING STORAGE, AND THEREFORE, IA MUST
C                           NOT BE LESS THAN M1+M2+2.
C                B      - VECTOR OF LENGTH M1+M2+2 CONTAINING THE
C                           RIGHT HAND SIDES OF THE INEQUALITY
C                           CONSTRAINTS IN ITS FIRST M1 LOCATIONS
C                           FOLLOWED BY THE M2 RIGHT HAND SIDES OF THE
C                           EQUALITY CONSTRAINTS. (INPUT)
C                           THE LAST TWO ELEMENTS OF B ARE USED
C                           AS WORKING STORAGE.
C                C      - VECTOR OF LENGTH N CONTAINING THE
C                           COEFFICIENTS OF THE OBJECTIVE FUNCTION.
C                           (INPUT)
C                N      - NUMBER OF UNKNOWNS IN THE MODEL. (INPUT)
C                M1     - NUMBER OF INEQUALITY CONSTRAINTS. (INPUT)
C                M2     - NUMBER OF EQUALITY CONSTRAINTS. (INPUT)
C                S      - VALUE OF THE OBJECTIVE FUNCTION. (OUTPUT)
C                PSOL   - VECTOR OF LENGTH N CONTAINING THE PRIMAL
C                           SOLUTION. (OUTPUT) PSOL IS ALSO USED AS
C                           WORK STORAGE AND THEREFORE MUST HAVE LENGTH
C                           AT LEAST MAX(N,M1+M2).
C                DSOL   - VECTOR OF LENGTH M1+M2+2 CONTAINING THE
C                           DUAL SOLUTION. (OUTPUT)
C                           THAT IS, DSOL(1),...,DSOL(M1+M2)
C                           CONTAIN THE SOLUTION TO THE
C                           PROBLEM MIN BT*Y SUBJECT TO AT*Y IS
C                           GREATER THAN OR EQUAL TO C AND Y GREATER
C                           THAN OR EQUAL TO 0 WHERE AT = A-TRANSPOSE
C                           AND BT = B-TRANSPOSE. WHEN THE PRIMAL
C                           PROBLEM HAS EQUALITY CONSTRAINTS, THE
C                           CORRESPONDING COMPONENTS OF THE DUAL
C                           SOLUTION ARE UNCONSTRAINED. DSOL(M1+M2+1)
C                           AND DSOL(M1+M2+2) ARE USED AS WORKING
C                           STORAGE.
C                RW     - WORK VECTOR OF LENGTH (M1+M2+2)*(M1+M2+2)+
C                           3*M1+2*M2+4.
C                IW     - WORK VECTOR OF LENGTH 2*M2+3*M1+4.
C                IER    - ERROR INDICATOR. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 130 INDICATES THAT IA IS LESS THAN
C                             M1+M2+2.
C                           IER = 131 INDICATES THAT THE COST CRITERION
C                             HAS UNBOUNDED VALUES.
C                           IER = 132 INDICATES THAT THE MAXIMUM NUMBER
C                             OF ITERATIONS WAS REACHED IN ZX0LP.
C                           IER = 133 INDICATES THAT NO FEASIBLE
C                             SOLUTION EXISTS.
C                         WARNING (WITH FIX)
C                           IER = 70 INDICATES THAT SOME ARTIFICIAL
C                             VARIABLES REMAINED IN THE SOLUTION BASIS
C                             AT A ZERO LEVEL AFTER PHASE 1.
C                             THIS CONDITION CAN BE CAUSED BY
C                             HAVING REDUNDANT CONSTRAINTS.
C                             NEVERTHELESS, A SOLUTION IS COMPUTED
C                             AND RETURNED IN PSOL AND DSOL.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZX0LP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ZX3LP IS INTENDED TO BE VERY EASY TO USE. THEREFORE,
C                THE CALLING PARAMETERS ARE VERY SIMPLE AND THERE ARE
C                NO OPTIONS. ZX0LP CAN BE USED DIRECTLY IN SITUATIONS
C                THAT REQUIRE MORE FLEXIBILITY THAN IS PROVIDED BY
C                ZX3LP.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZX3LP  (A,IA,B,C,N,M1,M2,S,PSOL,DSOL,RW,IW,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,N,M1,M2,IW(1),IER
      REAL               A(IA,1),B(1),C(1),S,PSOL(1),DSOL(1),RW(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IPHASE,ITMAX,I,IEND,IR,J,IBEG,INEXT,K,L,
     1                   LIC,M1P2,IDES,IDES2,IENDIW,ICOPI,II,JJ,
     2                   M12,M,JER,M1P1,M1P3,IEND1,IQ,M1MIQ,L1,IWK,
     3                   LICSV,JI
      REAL               ZERO,TEMP,ONE,EPS
      DATA               EPS/1.0E-5/
C                                  EPS IS USED IN TESTS FOR ZERO
C                                    IF ABS(T) .LE. EPS, THEN T IS
C                                    CONSIDERED TO BE ZERO
      DATA               ZERO/0.0/,ONE/1.0/
      DATA               ITMAX/10000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IEND = M1+M2
      M12 = IEND
C                                  TERMINAL ERROR - IA IS LESS THAN
C                                  M1+M2+2
      IF (IA.GE.M12+2) GO TO 5
      IER = 130
      GO TO 9000
    5 M1P1 = M1+1
      M1P2 = M1+2
      M1P3 = M1+3
      IEND1 = IEND+1
C                                  MOVE A AND B DOWN 2 ROWS
      DO 15 I=1,IEND
         K = IEND1-I
         IR = K+2
         B(IR) = B(K)
         PSOL(I) = I
         DO 10 J=1,N
            A(IR,J) = A(K,J)
   10    CONTINUE
   15 CONTINUE
      IR = IEND+2
C                                  CHECK EQUALITY CONSTRAINTS FOR
C                                  NEGATIVE RIGHT HAND SIDE.
      IF (M2.EQ.0) GO TO 30
      IBEG = M1+3
      DO 25 I=IBEG,IR
         IF (B(I).GE.ZERO) GO TO 25
         B(I) = -B(I)
         PSOL(I-2) = -PSOL(I-2)
         DO 20 J=1,N
            A(I,J) = -A(I,J)
   20    CONTINUE
   25 CONTINUE
C                                  RE-ORDER OTHER CONSTRAINTS SO
C                                  B(I) .GE. 0 I = 1,...,M1-IQ
   30 IQ = 0
      IF (M1.EQ.0) GO TO 60
      IEND = M1P2
      INEXT = IEND
   35 IF (B(IEND).GE.ZERO) GO TO 55
      IF (INEXT.EQ.IEND) GO TO 45
      TEMP = B(IEND)
      B(IEND) = B(INEXT)
      B(INEXT) = TEMP
      TEMP = PSOL(IEND-2)
      PSOL(IEND-2) = PSOL(INEXT-2)
      PSOL(INEXT-2) = TEMP
      DO 40 J=1,N
         TEMP = A(IEND,J)
         A(IEND,J) = A(INEXT,J)
         A(INEXT,J) = TEMP
   40 CONTINUE
   45 IQ = IQ+1
      PSOL(INEXT-2) = -PSOL(INEXT-2)
      B(INEXT) = -B(INEXT)
      DO 50 J=1,N
         A(INEXT,J) = -A(INEXT,J)
   50 CONTINUE
      INEXT = INEXT-1
   55 IEND = IEND-1
      IF (IEND.NE.2) GO TO 35
C                                  COMPUTE ROW 1 AND 2 OF A AND B
   60 DO 65 J=1,N
         A(2,J) = -C(J)
         A(1,J) = ZERO
      DO 65 I=2,IR
         A(1,J) = A(1,J)-A(I,J)
   65 CONTINUE
      B(1) = ZERO
      B(2) = ZERO
      DO 70 I=3,IR
         B(1) = B(1)-B(I)
   70 CONTINUE
      M = M12+1
      IF (IA.EQ.IR) GO TO 80
C                                  PACK A
      K = 0
      L = 0
      DO 75 J=1,N
      DO 75 I=1,IR
         K = MOD(K,IA)+1
         IF (K.EQ.1) L = L+1
         A(K,L) = A(I,J)
   75 CONTINUE
C                                  GET ICOLMS AND ROW
   80 LIC = IR+IQ
      M1MIQ = M1-IQ
      L1 = M1P1-IQ
      DO 95 I=1,LIC
         IF (I.GE.M1P3) GO TO 90
         IF (I.GT.L1) GO TO 85
         RW(I) = -ONE
         IW(I) = I+1
         GO TO 95
   85    RW(I) = ONE
         IW(I) = -I-1
         GO TO 95
   90    RW(I) = ZERO
         IW(I) = I-IQ
   95 CONTINUE
C                                  WORK STORAGE ASSIGNMENTS
C                                      ICOLMS(1) = IW(1)
C                                      IDES(1) = IW(IDES)
C                                      COPI(1,1) = RW(ICOPI)
C                                      ROW(1) = RW(1)
C                                      WA(1) = RW(IWK)
C                                      X(1) = DSOL(1)
      IW(M1P2) = 1
C                                  GET IDES
      IDES = LIC+1
      IW(IDES) = N+M1P2
      IDES2 = IDES+1
      IEND = IDES2+M1MIQ
      IENDIW = IDES2+M12
      K = N+1
      DO 100 I=IDES2,IENDIW
         IW(I) = K
         IF (I.EQ.IEND) K = N+M1P2
         K = K+1
  100 CONTINUE
C                                  GET COPI
      ICOPI = IDES
      IWK = IR*IR+ICOPI
      K = IWK-1
      DO 105 I=ICOPI,K
         RW(I) = ZERO
  105 CONTINUE
      L = ICOPI
      J = M1MIQ+1
      DO 110 I=1,J
         L = L+IR
         RW(L) = ONE
  110 CONTINUE
      J = 0
      DO 115 I=ICOPI,K,IR
         RW(I+J) = ONE
         J = J+1
  115 CONTINUE
      K = 1
C                                  SOLVE PHASE 1 PROBLEM
      IPHASE = 1
      CALL ZX0LP (IPHASE,A,B,IW,RW,K,M,N,ITMAX,LIC,IR,RW(ICOPI),
     1            IW(IDES),DSOL,RW(IWK),IER)
      IF (IER.NE.0) GO TO 185
      K = M1P2+N
      J = 2
C                                  CHECK PHASE 1 SOLUTION FOR ARTIFICIAL
C                                  VARIABLES
      DO 125 I=IDES2,IENDIW
         IF (IW(I).LE.K) GO TO 120
         IF (DSOL(J).GT.EPS) GO TO 180
C                                  ARTIFICIAL VARIABLES REMAIN IN THE
C                                  SOLUTION AT A ZERO LEVEL
         JER = 70
  120    J = J+1
  125 CONTINUE
C                                  INTERCHANGE FIRST TWO ROWS OF THE A
C                                  MATRIX AND NEGATE THE SECOND ROW
      K = (N-1)*IR+1
      DO 140 L=1,K,IR
         J = (L+IA-1)/IA
         I = L-(J-1)*IA
         TEMP = A(I,J)
         IF (I.LT.IA) GO TO 130
         II = 1
         JJ = J+1
         GO TO 135
  130    II = I+1
         JJ = J
  135    A(I,J) = A(II,JJ)
         A(II,JJ) = -TEMP
  140 CONTINUE
C                                  INTERCHANGE FIRST TWO ELEMENTS OF THE
C                                  RIGHT HAND SIDE AND NEGATE THE
C                                  SECOND ELEMENT
      B(2) = -B(1)
      B(1) = ZERO
C                                  INTERCHANGE FIRST TWO COLUMNS OF
C                                  COPI AND NEGATE THE SECOND COLUMN
      J = ICOPI+IR-1
      DO 145 I=ICOPI,J
         TEMP = RW(I)
         RW(I) = RW(I+IR)
         RW(I+IR) = -TEMP
  145 CONTINUE
      J = ICOPI+IR*IR-1
      DO 150 I=ICOPI,J,IR
         TEMP = RW(I)
         RW(I) = RW(I+1)
         RW(I+1) = TEMP
  150 CONTINUE
      K = 2
C                                  NEGATE ROW
      DO 155 I=1,LIC
         RW(I) = -RW(I)
  155 CONTINUE
C                                  ICOLMS(1) = 1
      IW(1) = 1
C                                  INTERCHANGE IDES(1) AND IDES(2)
      J = IW(IDES)
      IW(IDES) = IW(IDES2)
      IW(IDES2) = J
      IW(M1P2) = -2
      LICSV = LIC
C                                  REMOVE ARTIFICIAL VARIABLES FROM
C                                    TABLEAU, IF POSSIBLE
      IF (JER.EQ.0) LIC = M1P2
C                                  SOLVE PHASE 2 PROBLEM
      IPHASE = 2
      CALL ZX0LP (IPHASE,A,B,IW,RW,K,M,N,ITMAX,LIC,IR,RW(ICOPI),
     1            IW(IDES),DSOL,RW(IWK),IER)
      LIC = LICSV
      S = DSOL(1)
C                                  RE-ORDER THE PRIMAL SOLUTION
      DO 160 J=1,M12
         RW(J) = PSOL(J)
  160 CONTINUE
      DO 165 J=1,N
         PSOL(J) = ZERO
  165 CONTINUE
      DO 170 J=1,M
         K = IW(IDES+J)
         IF (K.GT.N) GO TO 170
         PSOL(K) = DSOL(J+1)
  170 CONTINUE
C                                  GET DUAL SOLUTION
      JI = LIC+1+IR
      TEMP = RW(JI)
      DO 175 I=1,M12
         J = ABS(RW(I))
         JI = JI+IR
         DSOL(J) = RW(JI)+TEMP
         IF (RW(I).LT.ZERO) DSOL(J) = -DSOL(J)
         IF (I.LE.M1) DSOL(J) = ABS(DSOL(J))
  175 CONTINUE
      GO TO 195
C                                  ARTIFICIAL VARIABLES ARE IN THE
C                                  SOLUTION HENCE NO FEASIBLE
C                                  SOLUTION EXISTS
  180 IER = 133
  185 IF (IER.EQ.130) IER = 133
      DO 190 J=1,M12
         RW(J) = PSOL(J)
  190 CONTINUE
C                                  RESTORE A AND B
  195 IF (IA.EQ.IR) GO TO 210
C                                  UNPACK A
      K = MOD(N*IR,IA)
      IF (K .EQ. 0) K = IA
      L = (N*IR+IA-1)/IA
      J = N
      DO 205 JJ=1,N
         I = IR
         DO 200 II=1,IR
            A(I,J) = A(K,L)
            I = I-1
            K = K-1
            IF (K.NE.0) GO TO 200
            K = IA
            L = L-1
  200    CONTINUE
         J = J-1
  205 CONTINUE
  210 DO 220 I=1,M12
         B(I) = B(I+2)
         DO 215 J=1,N
            A(I,J) = A(I+2,J)
  215    CONTINUE
  220 CONTINUE
C                                  PERMUTE ROWS OF A AND ELEMENTS OF B
C                                  ACCORDING TO PERMUTATIONS STORED IN
C                                  RW
      DO 230 I=1,M12
         J = ABS(RW(I))
         IF (J.EQ.I) GO TO 230
         TEMP = RW(I)
         RW(I) = RW(J)
         RW(J) = TEMP
         TEMP = B(I)
         B(I) = B(J)
         B(J) = TEMP
         DO 225 K=1,N
            TEMP = A(I,K)
            A(I,K) = A(J,K)
            A(J,K) = TEMP
  225    CONTINUE
  230 CONTINUE
      DO 240 I=1,M12
         IF (RW(I).GT.ZERO) GO TO 240
         B(I) = -B(I)
         DO 235 J=1,N
            A(I,J) = -A(I,J)
  235    CONTINUE
  240 CONTINUE
      IF (IER .EQ. 0) IER=JER
 9000 CONTINUE
      IF (JER.NE.0) CALL UERTST (JER,6HZX3LP )
      IF (IER.GT.JER) CALL UERTST (IER,6HZX3LP )
 9005 RETURN
      END

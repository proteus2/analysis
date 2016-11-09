C   IMSL ROUTINE NAME   - NAWRPU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - WILSONS ANOVA (1, 2, OR 3 WAY DESIGNS) WITH
C                           UNEQUAL REPLICATION
C
C   USAGE               - CALL NAWRPU (Y,N,NREPS,IR,NF,NDF,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR CONTAINING THE RESPONSE VALUES.
C                           Y IS OF LENGTH THE SUM OF THE ELEMENTS IN
C                           NREPS. SEE THE DESCRIPTION OF NREPS BELOW.
C                           SEE THE PROGRAMMING NOTES IN THE MANUAL
C                           DOCUMENT FOR THE ORDERING OF THE ELEMENTS
C                           IN Y.
C                         ON OUTPUT, Y WILL BE SORTED INTO ASCENDING
C                           ORDER.
C                N      - INPUT VECTOR OF LENGTH 4 INDICATING THE NUMBER
C                           OF FACTORS AND NUMBER OF LEVELS OF EACH
C                           FACTOR.
C                         N(1) CONTAINS THE NUMBER OF FACTORS.
C                           N(1) MUST BE ONE, TWO, OR THREE.
C                         N(2) CONTAINS THE NUMBER OF LEVELS OF FACTOR
C                           ONE. N(2) MUST BE GREATER THAN OR EQUAL TO
C                           TWO.
C                         N(3) CONTAINS THE NUMBER OF LEVELS OF FACTOR
C                           TWO. N(3) IS REQUIRED ONLY WHEN N(1) IS TWO
C                           OR THREE AND MUST BE GREATER THAN OR EQUAL
C                           TO TWO.
C                         N(4) CONTAINS THE NUMBER OF LEVELS OF FACTOR
C                           THREE. N(4) IS REQUIRED ONLY WHEN N(1) IS
C                           THREE AND MUST BE GREATER THAN OR EQUAL TO
C                           TO TWO.
C                NREPS  - INPUT VECTOR CONTAINING THE NUMBER OF
C                           REPLICATIONS PER CELL.
C                           IF N(1)=1, THEN NREPS IS OF LENGTH N(2).
C                           IF N(1)=2, THEN NREPS IS OF LENGTH N(2)*N(3)
C                           IF N(1)=3, THEN NREPS IS OF LENGTH
C                             N(2)*N(3)*N(4).
C                           EACH ELEMENT OF NREPS MUST BE GREATER
C                           THAN OR EQUAL TO ONE.
C                IR     - WORK VECTOR OF SAME LENGTH AS Y.
C                           THE ORIGINAL POSITION OF EACH ELEMENT IN THE
C                           VECTOR Y IS CONTAINED IN THE CORRESPONDING
C                           ELEMENT IN IR.
C                NF     - WORK VECTOR CONTAINING THE FREQUENCIES OF
C                           OBSERVATIONS BELOW THE MEDIAN AND NUMBER OF
C                           REPLICATIONS PER EFFECT. THE LENGTH OF NF
C                           IS DEPENDENT UPON N(1) AS FOLLOWS
C                           FOR N(1)=1, THE LENGTH OF NF IS N(2)
C                           FOR N(1)=2, THE LENGTH OF NF IS
C                             N(2)*N(3)+2*N(2)+2*N(3)
C                           FOR N(1)=3, THE LENGTH OF NF IS
C                             N(2)*N(3)*N(4)+2*N(2)+2*N(3)+2*N(4)+
C                             2*N(2)*N(3)+2*N(2)*N(4)+2*N(3)*N(4)
C                NDF    - OUTPUT VECTOR OF LENGTH 7 CONTAINING DEGREES
C                           OF FREEDOM FOR CHI-SQUARE TESTS.
C                         NDF(1) CONTAINS THE FACTOR ONE DEGREES OF
C                           FREEDOM
C                         NDF(2) CONTAINS THE FACTOR TWO DEGREES OF
C                           FREEDOM. NDF(2) IS DEFINED ONLY WHEN
C                           N(1) = 2 OR 3.
C                         NDF(3) CONTAINS THE FACTOR THREE DEGREES OF
C                           FREEDOM. NDF(3) IS DEFINED ONLY WHEN
C                           N(1) = 3.
C                         NDF(4) CONTAINS THE INTERACTION DEGREES OF
C                           FREEDOM FOR FACTORS ONE AND TWO. NDF(4) IS
C                           DEFINED ONLY WHEN N(1) = 2 OR 3.
C                         NDF(5) CONTAINS THE INTERACTION DEGREES OF
C                           FREEDOM FOR FACTORS ONE AND THREE. NDF(5) IS
C                           DEFINED ONLY WHEN N(1)=3.
C                         NDF(6) CONTAINS THE INTERACTION DEGREES OF
C                           FREEDOM FOR FACTORS TWO AND THREE. NDF(6) IS
C                           DEFINED ONLY WHEN N(1)=3.
C                         NDF(7) CONTAINS THE INTERACTION DEGREES OF
C                           FREEDOM FOR FACTORS ONE, TWO, AND THREE.
C                           NDF(7) IS DEFINED ONLY WHEN N(1)=3.
C                STAT   - OUTPUT VECTOR OF LENGTH 14 CONTAINING CHI-SQ
C                           TEST STATISTIC VALUES AND CORRESPONDING
C                           PROBABILITIES OF OBSERVING THESE VALUES OR
C                           GREATER IF THE NULL HYPOTHESES ARE TRUE.
C                         STAT(1) CONTAINS THE FACTOR ONE CHI-SQ
C                           STATISTIC
C                         STAT(2) CONTAINS THE FACTOR TWO CHI-SQ
C                           STATISTIC. STAT(2) IS DEFINED
C                           ONLY WHEN N(1) = 2 OR 3.
C                         STAT(3) CONTAINS THE FACTOR THREE CHI-SQ
C                           STATISTIC. STAT(3) IS DEFINED
C                           ONLY WHEN N(1) = 3.
C                         STAT(4) CONTAINS THE INTERACTION CHI-SQ
C                           STATISTIC FOR FACTORS ONE AND TWO. STAT(4)
C                           IS DEFINED ONLY WHEN N(1) = 2 OR 3.
C                         STAT(5) CONTAINS THE INTERACTION CHI-SQ
C                           STATISTIC FOR FACTORS ONE AND THREE. STAT(5)
C                           DEFINED ONLY WHEN N(1) = 3.
C                         STAT(6) CONTAINS THE INTERACTION CHI-SQ
C                           STATISTIC FOR FACTORS TWO AND THREE. STAT(6)
C                           IS  DEFINED ONLY WHEN N(1) = 3.
C                         STAT(7) CONTAINS THE INTERACTION CHI-SQ
C                           STATISTIC FOR FACTORS ONE, TWO, AND THREE.
C                           STAT(7) IS DEFINED ONLY WHEN N(1)=3.
C                         STAT(8) CONTAINS THE PROBABILITY
C                           CORRESPONDING TO STAT(1).
C                         STAT(9) CONTAINS THE PROBABILITY
C                           CORRESPONDING TO STAT(2). STAT(9) IS
C                           DEFINED ONLY WHEN N(1) = 2 OR 3.
C                         STAT(10) CONTAINS THE PROBABILITY
C                           CORRESPONDING TO STAT(3). STAT(10) IS
C                           DEFINED ONLY WHEN N(1) = 3.
C                         STAT(11) CONTAINS THE PROBABILITY
C                           CORRESPONDING TO STAT(4). STAT(11) IS
C                           DEFINED ONLY WHEN N(1) = 2 OR 3.
C                         STAT(12) CONTAINS THE PROBABILITY
C                           CORRESPONDING TO STAT(5). STAT(12) IS
C                           DEFINED ONLY WHEN N(1) = 3.
C                         STAT(13) CONTAINS THE PROBABILITY
C                           CORRESPONDING TO STAT(6). STAT(13) IS
C                           DEFINED ONLY WHEN N(1) = 3.
C                         STAT(14) CONTAINS THE PROBABILITY
C                           CORRESPONDING TO STAT(7). STAT(14) IS
C                           DEFINED ONLY WHEN N(1) = 3.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES N(1) IS NOT 1,2, OR 3
C                           IER=130 INDICATES N(2),N(3) OR N(4) IS LESS
C                             THAN 2
C                           IER=131 INDICATES SOME ELEMENT OF NREPS IS
C                             LESS THAN 1
C                         WARNING ERROR
C                           IER=36 INDICATES A NEGATIVE CHI-SQUARE VALUE
C                             FOR SOME INTERACTION TERM WAS CALCULATED
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,
C                           UERTST,UGETIO,VSRTR
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           UERTST,UGETIO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NAWRPU (Y,N,NREPS,IR,NF,NDF,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N(1),NREPS(1),IR(1),NF(1),NDF(1),IER
      REAL               Y(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IADD,J,K,L,LB,LC,LL,LR,MIDEL,NAA,NB,NBB,
     1                   NC,NFF,NFT,NR,NSUM,NT,NTOT,NUM
      REAL               CHI,CHITA,EXPC,FNDF,XNAA,XNBB,XNSUM,XREP,YMED
      DOUBLE PRECISION   CHIT,CHINT,DIF
C                                  FIRST EXECUTABLE STATEMENT
      NFF = N(1)
C                                  CHECK N(1)
      IF ((NFF .GE. 1) .AND. (NFF .LE. 3)) GO TO 5
      IER = 129
      GO TO 9000
C                                  CHECK N(2)
    5 NR = N(2)
      IF (NR .GE. 2) GO TO 10
      IER = 130
      GO TO 9000
C                                  CHECK N(3)
   10 IF (NFF .GE. 2) GO TO 20
      LR = 0
      NC = 1
      N(3) = 1
      LC = 0
   15 NB = 1
      NBB = 0
      N(4) = 1
      LB = 0
      GO TO 30
   20 LR = NR
      NC = N(3)
      LC = NC
      IF (NC .GE. 2) GO TO 25
      IER = 130
      GO TO 9000
   25 IF (NFF .NE. 3) GO TO 15
C                                  CHECK N(4)
      NBB = NR
      NB = N(4)
      LB = NB
      IF (NB .GE. 2) GO TO 30
      IER = 130
      GO TO 9000
   30 NTOT = NR * NC * NB
      NFT = NTOT + LR + LC + LB + NBB*LC + NR*LB + NC*LB
      NT = NFT + NFT - NTOT
C                                  INITIALIZE FREQUENCY ARRAY
      DO 35 I=1,NT
         NF(I) = 0
   35 CONTINUE
C                                  CHECK NREPS AND DISTRIBUTE IN NF
      NT = NFT
      LL = NTOT/NR
      K = 1
      I = 1
   40 IADD = NREPS(I)
      IF (IADD .GE. 1) GO TO 45
      IER = 131
      GO TO 9000
   45 IF (NFF .GT. 1) GO TO 80
   50 I = I + 1
      IF (I .LE. NTOT) GO TO 40
      IER = 0
C                                  ACCUMULATE ELEMENTS OF NREPS
      DO 55 I=2,NTOT
         NREPS(I) = NREPS(I) + NREPS(I-1)
   55 CONTINUE
      NSUM = NREPS(NTOT)
C                                  INITIALIZE IR
      DO 60 I=1,NSUM
         IR(I) = I
   60 CONTINUE
C                                  SORT DATA VECTOR
      CALL VSRTR (Y,NSUM,IR)
C                                  FIND MEDIAN
      MIDEL = NSUM/2 + 1
      YMED = Y(MIDEL)
C                                  COUNT Y BELOW MEDIAN
      IADD = 1
      NBB = 0
      NT = NTOT
      MIDEL = MIDEL - 1
      K = 2
      J = 1
   65 IF (Y(J) .EQ. YMED) GO TO 95
      NBB = NBB + 1
      LC = IR(J)
C                                  FIND CELL NUMBER
      DO 70 L=1,NTOT
         IF (LC .GT. NREPS(L)) GO TO 70
         I = L
         GO TO 75
   70 CONTINUE
   75 NF(I) = NF(I) + IADD
C                                  INCREMENT ROW NUMBER
      IF (NFF .EQ. 1) GO TO 90
   80 LR = I/LL
      IF (LR*LL .NE. I) LR = LR+1
      NUM = NT + LR
      NF(NUM) = NF(NUM) + IADD
C                                  INCREMENT COLUMN NUMBER
      NUM = I/NB
      IF (NUM*NB .NE. I) NUM=NUM+1
      LC = MOD(NUM,NC)
      IF (LC .EQ. 0) LC = NC
      NUM = NT + NR + LC
      NF(NUM) = NF(NUM) + IADD
      IF (NFF .EQ. 3) GO TO 85
      GO TO (50,90),K
C                                  INCREMENT BLOCK NUMBER
   85 LB = MOD(I,NB)
      IF (LB .EQ. 0) LB=NB
      NUM = NT + NR + NC + LB
      NF(NUM) = NF(NUM) + IADD
C                                  INCREMENT R X C NUMBER
      NUM = NT + NR + NC + NB + (LR-1)*NC + LC
      NF(NUM) = NF(NUM) + IADD
C                                  INCREMENT R X B NUMBER
      NUM = NT + NR + NC + NB + NR*NC + (LR-1)*NB + LB
      NF(NUM) = NF(NUM) + IADD
C                                  INCREMENT C X B NUMBER
      NUM = NT + NR + NC + NB + NR*NC + NR*NB + (LC-1)*NB + LB
      NF(NUM) = NF(NUM) + IADD
      IF (K .EQ. 1) GO TO 50
   90 J = J+1
      IF (J .LE. MIDEL) GO TO 65
   95 NAA = NSUM - NBB
C                                  CHANGE NREPS ARRAY TO INPUT MODE
      LC = 0
      DO 100 I=2,NTOT
         LC = NREPS(I-1) + LC
         NREPS(I) = NREPS(I) - LC
  100 CONTINUE
      STAT(3) = 0.0
C                                  DEFINE DEGREES OF FREEDOM
      NDF(1) = NR - 1
      IF (NFF .EQ. 1) GO TO 105
      NDF(2) = NC - 1
      NDF(4) = NDF(1) * NDF(2)
      IF (NFF .EQ. 2) GO TO 105
      NDF(3) = NB - 1
      NDF(5) = NDF(1) * NDF(3)
      NDF(6) = NDF(2) * NDF(3)
      NDF(7) = NDF(6) * NDF(1)
C                                  INITIALIZE VARIABLES FOR CHI-SQ CALC
  105 XNAA = NAA
      XNBB = NBB
      XNSUM = NSUM
      XNAA = XNAA/XNSUM
      XNBB = XNBB/XNSUM
      LL = 1
      LB = 1
      LC = NTOT
C                                  LOOP FOR CHI-SQ CALCULATIONS
      DO 190 I=1,8
         IF ((NFF .EQ. 2) .AND. (I .EQ. 4)) GO TO 190
         K = I - 1
         CHI = 0.0
         CHITA = 0.0
         GO TO (110,160,160,160,165,165,165,175), I
C                                  CHI-SQ STATISTIC CALCULATION
  110    DO 125 J=LB,LC
            IF (I .EQ. 1) GO TO 115
            NUM = NFT + J - NTOT
            XREP = NF(NUM)
            GO TO 120
  115       XREP = NREPS(J)
  120       XNSUM = NF(J)
            EXPC = XREP * XNAA
            DIF = XREP - XNSUM - EXPC
            CHITA = CHITA + DIF*DIF/EXPC
            EXPC = XREP * XNBB
            DIF = XNSUM - EXPC
            CHI = CHI + DIF*DIF/EXPC
  125    CONTINUE
         CHI = CHI + CHITA
         IF (I .GT. 1) GO TO 130
         CHIT = CHI
         IF (NFF .GT. 1) GO TO 190
         K = 1
         LL = 0
         GO TO 150
  130    GO TO (150,150,150,150,135,140,145,150), I
C                                  SUBTRACT MAIN EFFECT FROM INTERACTION
  135    CHI = CHI - STAT(1) - STAT(2)
         GO TO 150
  140    CHI = CHI - STAT(1) - STAT(3)
         GO TO 150
  145    CHI = CHI - STAT(2) - STAT(3)
  150    STAT(K) = CHI
C                                  CHI-SQ PROBABILITY
         FNDF = NDF(K)
         CALL MDCH(CHI,FNDF,CHITA,J)
  155    STAT(K+7) = 1. - CHITA
         GO TO 185
C                                  R,C, AND B CHI-SQ
  160    LB = LC + 1
         LC = LC + N(I)
         GO TO 110
C                                  2-FACTOR INTERACTIONS
  165    CHINT = CHIT - STAT(1) - STAT(2) - STAT(3)
         IF (NFF .EQ. 3) GO TO 170
         CHI = CHINT
         LL = 0
         IF(CHI .LT. 0.0) GO TO 180
         GO TO 150
  170    LB = LC + 1
         LC = LC + NTOT/N(9-I)
         GO TO 110
C                                  R X C X B INTERACTION
  175    CHI = CHINT - STAT(4) - STAT(5) - STAT(6)
         IF(CHI .LT. 0.0) GO TO 180
         GO TO 150
C                                  CHI-SQ LESS THAN ZERO
  180    STAT(K+7) = -1.0
         STAT(K) = 0.0
         IER = 36
         GO TO 9000
  185    IF (LL .EQ. 0) GO TO 9005
  190 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HNAWRPU)
 9005 RETURN
      END

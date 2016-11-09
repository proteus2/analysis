C   IMSL ROUTINE NAME   - NAWRPE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - WILSONS ANOVA (1, 2, OR 3 WAY DESIGNS) WITH
C                           EQUAL REPLICATION
C
C   USAGE               - CALL NAWRPE (Y,N,NREPS,IR,NF,NDF,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR CONTAINING THE RESPONSE VALUES.
C                           IF N(1)=1, THEN Y IS OF LENGTH
C                             N(2)*NREPS.
C                           IF N(1)=2, THEN Y IS OF LENGTH
C                             N(2)*N(3)*NREPS.
C                           IF N(1)=3, THEN Y IS OF LENGTH
C                             N(2)*N(3)*N(4)*NREPS.
C                           SEE THE DECSRIPTION OF N BELOW. SEE THE
C                           PROGRAMMING NOTES IN THE MANUAL DOCUMENT
C                           FOR THE ORDERING OF THE ELEMENTS IN Y.
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
C                         N(4) CONTAINS THE NUMBER OF LEVELS FOR FACTOR
C                           THREE. N(4) IS REQUIRED ONLY WHEN N(1) IS
C                           THREE AND MUST BE GREATER THAN OR EQUAL TO
C                           TWO.
C                NREPS  - INPUT VALUE CONTAINING THE NUMBER OF REPLI-
C                           CATIONS PER CELL. NREPS MUST BE GREATER
C                           THAN OR EQUAL TO TWO.
C                IR     - WORK VECTOR OF SAME LENGTH AS Y.
C                           THE ORIGINAL POSITION OF EACH ELEMENT IN THE
C                           VECTOR Y IS CONTAINED IN THE CORRESPONDING
C                           ELEMENT IN IR.
C                NF     - WORK VECTOR CONTAINING THE FREQUENCIES OF
C                           OBSERVATIONS BELOW THE MEDIAN. THE LENGTH
C                           OF NF IS DEPENDENT UPON N(1) AS FOLLOWS
C                           FOR N(1)=1, THE LENGTH OF NF IS N(2)
C                           FOR N(1)=2, THE LENGTH OF NF IS
C                             N(2)*N(3)+N(2)+N(3)
C                           FOR N(1)=3, THE LENGTH OF NF IS
C                             N(2)*N(3)*N(4)+N(2)+N(3)
C                             +N(4)+N(2)*N(3)+N(2)*N(4)+N(3)*N(4)
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
C                           STATISTIC FOR FACORS ONE, TWO, AND THREE.
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
C                           IER=131 INDICATES NREPS IS LESS THAN 2
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
      SUBROUTINE NAWRPE (Y,N,NREPS,IR,NF,NDF,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N(1),NREPS,IR(1),NF(1),NDF(1),IER
      REAL               Y(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L,LB,LC,LL,LR,NAA,NB,NBB,NC,NFF,NFTOT,
     1                   NR,NSUM,NTOT,MIDEL
      REAL               CHI,CHITA,EXPC,EXPCA,FNDF,XF,XNAA,XNBB,XNSUM,
     1                   XNTOT,XREP,YMED
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
      N(4) = 1
      NBB = 0
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
C                                  CHECK NREPS
      IF (NREPS .GE. 2) GO TO 35
      IER = 131
      GO TO 9000
   35 NSUM = NREPS * NTOT
      IER = 0
C                                  INITIALIZE IR
      DO 40 I=1,NSUM
         IR(I) = I
   40 CONTINUE
C                                  SORT DATA VECTOR
      CALL VSRTR (Y,NSUM,IR)
C                                  FIND MEDIAN
      MIDEL = NSUM/2 + 1
      YMED = Y(MIDEL)
C                                  INITIALIZE FREQUENCY ARRAY
      NFTOT = NTOT + LR + LC + LB + NBB*LC + NR*LB + LC*LB
      DO 45 I=1,NFTOT
         NF(I) = 0
   45 CONTINUE
C                                  COUNT Y BELOW MEDIAN
      LL = NTOT/NR
      NBB = 0
      MIDEL = MIDEL - 1
      DO 50 I=1,MIDEL
C                                  CHECK FOR TIE WITH MEDIAN
         IF (Y(I) .EQ. YMED) GO TO 55
         NBB = NBB + 1
         LC = IR(I)
C                                  FIND CELL NUMBER
         LB = LC/NREPS
         IF (LB*NREPS .NE. LC) LB = LB+1
C                                  INCREMENT COUNT FOR THAT CELL
         NF(LB) = NF(LB) + 1
C                                  FIND ROW NUMBER FOR THIS Y
         IF (NFF .EQ. 1) GO TO 50
         LR = LB/LL
         IF (LR*LL .NE. LB) LR=LR+1
         J = NTOT + LR
         NF(J) = NF(J) + 1
C                                  FIND COLUMN NUMBER FOR THIS Y
         J = LB/NB
         IF (J*NB .NE. LB) J=J+1
         LC = MOD(J,NC)
         IF (LC .EQ. 0) LC = NC
         J = NTOT + NR + LC
         NF(J) = NF(J) + 1
C                                  FIND BLOCK NUMBER FOR THIS Y
         IF (NFF .EQ. 2) GO TO 50
         L = MOD(LB,NB)
         IF (L .EQ. 0) L=NB
         J = NTOT + NR + NC + L
         NF(J) = NF(J) + 1
C                                  FIND R X C NUMBER FOR THIS Y
         J = NTOT + NR + NC + NB
         LB = J + (LR-1)*NC + LC
         NF(LB) = NF(LB) + 1
C                                  FIND R X B NUMBER FOR THIS Y
         J = J + NR*NC
         LB = J + (LR-1)*NB + L
         NF(LB) = NF(LB) + 1
C                                  FIND C X B NUMBER FOR THIS Y
         J = J + NR*NB
         LB = J + (LC-1)*NB + L
         NF(LB) = NF(LB) + 1
   50 CONTINUE
C                                  INITIALIZE CHI-SQ VALUES
   55 DO 60 I=1,7
         STAT(I) = 0.0
   60 CONTINUE
C                                  DEFINE DEGREES OF FREEDOM
      NDF(1) = NR - 1
      IF (NFF .EQ. 1) GO TO 65
      NDF(2) = NC - 1
      NDF(4) = NDF(1) * NDF(2)
      IF (NFF .EQ. 2) GO TO 65
      NDF(3) = NB - 1
      NDF(5) = NDF(1) * NDF(3)
      NDF(6) = NDF(2) * NDF(3)
      NDF(7) = NDF(1) * NDF(6)
   65 NAA = NSUM - NBB
      XNAA = NAA
      XNBB = NBB
      XNSUM = NSUM
      XREP = NREPS
      LL = 1
      LR = NTOT
      LB = 1
      LC = NTOT
C                                  CALCULATION AND TEST LOOP
      DO 155 I=1,8
         IF ((NFF .EQ. 2) .AND. (I .EQ. 4)) GO TO 155
         K = I - 1
         CHI = 0.D0
         GO TO (70,125,125,125,130,135,135,140), I
C                                  CHI-SQ FOR NAA .EQ. NBB
   70    XNTOT = LR
         IF (NAA .NE. NBB) GO TO 85
         EXPC = XNSUM * 0.5/XNTOT
         DO 75 J=LB,LC
            XF = NF(J)
            DIF = XF - EXPC
            CHI = CHI + DIF*DIF
   75    CONTINUE
         CHI = (CHI+CHI)/EXPC
   80    IF (I .GT. 1) GO TO 95
         CHIT = CHI
         IF (NFF .GT. 1) GO TO 155
         K = 1
         LL = 0
         GO TO 115
C                                  CHI-SQ FOR NAA .NE. NBB
   85    CHITA = 0.D0
         EXPC = XNBB/XNTOT
         EXPCA = XNAA/XNTOT
         DO 90 J=LB,LC
            XF = NF(J)
            DIF = XF - EXPC
            CHI = CHI + DIF*DIF
            DIF = XREP - XF - EXPCA
            CHITA = CHITA + DIF*DIF
   90    CONTINUE
         CHI = CHI/EXPC + CHITA/EXPCA
         GO TO 80
   95    GO TO (115,115,115,115,100,105,110,115), I
C                                  SUBTRACT MAIN EFFECT FROM INTERACTION
  100    CHI = CHI - STAT(1) - STAT(2)
         GO TO 115
  105    CHI = CHI - STAT(1) - STAT(3)
         GO TO 115
  110    CHI = CHI - STAT(2) - STAT(3)
  115    STAT(K) = CHI
C                                  CHI-SQ TEST MDCH INOVATION
         FNDF = NDF(K)
         CALL MDCH(CHI,FNDF,CHITA,J)
  120    STAT(K+7) = 1. - CHITA
         GO TO 150
C                                  R,C, AND B INITIALIZATION FOR CHI-SQ
  125    LB = LC + 1
         LC = LC + N(I)
         LR = N(I)
         XREP = NSUM/LR
         GO TO 70
C                                  2-FACTOR INTERACTION CHI-SQ
  130    CHINT = CHIT - STAT(1) - STAT(2) - STAT(3)
         IF (NFF .EQ. 3) GO TO 135
         CHI = CHINT
         LL = 0
         IF(CHI .LT. 0.0) GO TO 145
         GO TO 115
  135    LB = LC + 1
         LC = LC + NTOT/N(9-I)
         LR = LC - LB + 1
         XREP = NSUM/LR
         GO TO 70
C                                  R X C X B  INTERACTION CHI-SQ
  140    CHI = CHINT - STAT(4) - STAT(5) - STAT(6)
         IF(CHI .LT. 0.0) GO TO 145
         GO TO 115
C                                  NEGATIVE CHI-SQ BY SUBTRACTION
  145    STAT(K+7) = -1.0
         IER = 36
         GO TO 9000
  150    IF (LL .EQ. 0) GO TO 9005
  155 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HNAWRPE)
 9005 RETURN
      END

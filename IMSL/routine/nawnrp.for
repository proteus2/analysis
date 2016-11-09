C   IMSL ROUTINE NAME   - NAWNRP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - WILSONS ANOVA (2 OR 3 WAY DESIGNS)
C                           WITHOUT REPLICATES
C
C   USAGE               - CALL NAWNRP (Y,N,IR,NF,NDF,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR CONTAINING THE RESPONSES.
C                           IF N(1)=2, THEN Y IS OF LENGTH
C                             N(2)*N(3).
C                           IF N(1)=3, THEN Y IS OF LENGTH
C                             N(2)*N(3)*N(4).
C                           SEE THE DESCRIPTION OF N BELOW. SEE THE
C                           PROGRAMMING NOTES IN THE MANUAL DOCUMENT
C                           FOR THE ORDERING OF THE ELEMENTS IN Y.
C                         ON OUTPUT, Y WILL BE SORTED INTO ASCENDING
C                           ORDER.
C                N      - INPUT VECTOR OF LENGTH 4 INDICATING THE NUMBER
C                           OF FACTORS AND THE NUMBER OF LEVELS OF EACH
C                           FACTOR.
C                         N(1) CONTAINS THE NUMBER OF FACTORS. N MUST BE
C                           TWO OR THREE.
C                         N(2) CONTAINS THE NUMBER OF LEVELS FOR FACTOR
C                           ONE. N(2) MUST BE GREATER THAN OR EQUAL TO
C                           TWO.
C                         N(3) CONTAINS THE NUMBER OF LEVELS FOR FACTOR
C                           TWO. N(2) MUST BE GREATER THAN OR EQUAL TO
C                           TWO.
C                         N(4) CONTAINS THE NUMBER OF LEVELS FOR FACTOR
C                           THREE. N(3) IS REQUIRED ONLY WHEN N(1)=3
C                           AND MUST BE GREATER THAN OR EQUAL TO TWO.
C                IR     - WORKING VECTOR OF SAME LENGTH AS Y.
C                           THE ORIGINAL POSITION OF EACH ELEMENT IN THE
C                           VECTOR Y IS CONTAINED IN THE CORRESPONDING
C                           ELEMENT IN IR.
C                NF     - OUTPUT VECTOR OF LENGTH N(2)+N(3)+..+N(N(1)+1)
C                           CONTAINING COUNT OF OBSERVATIONS BELOW THE
C                           MEDIAN FOR EACH LEVEL OF EVERY FACTOR.
C                NDF    - OUTPUT VECTOR OF LENGTH 3 CONTAINING DEGREES
C                           OF FREEDOM FOR CHI-SQUARE TESTS.
C                         NDF(1) CONTAINS THE DEGREES OF FREEDOM FOR
C                           TEST OF FACTOR ONE.
C                         NDF(2) CONTAINS THE DEGREES OF FREEDOM FOR
C                           TEST OF FACTOR TWO.
C                         NDF(3) CONTAINS THE DEGREES OF FREEDOM FOR
C                           TEST OF FACTOR THREE. NDF(3) IS DEFINED ONLY
C                           WHEN N(1)=3.
C                STAT   - OUTPUT VECTOR OF LENGTH 6 CONTAINING CHI-
C                           SQUARE TEST STATISTIC VALUES AND CORRESPOND-
C                           ING PROBABILITIES OF OBSERVING THESE VALUES
C                           OR GREATER IF THE NULL HYPOTHESIS OF NO MAIN
C                           EFFECT OF EACH FACTOR IS TRUE.
C                         STAT(1) CONTAINS THE CHI-SQUARE STATISTIC FOR
C                           FACTOR ONE.
C                         STAT(2) CONTAINS THE CHI-SQUARE STATISTIC FOR
C                           FACTOR TWO.
C                         STAT(3) CONTAINS THE CHI-SQUARE STATISTIC FOR
C                           FACTOR THREE. STAT(3) IS DEFINED ONLY WHEN
C                           N(1)=3.
C                         STAT(4) CONTAINS THE PROBABILITY CORRESPONDING
C                           TO STAT(1).
C                         STAT(5) CONTAINS THE PROBABILITY CORRESPONDING
C                           TO STAT(2).
C                         STAT(6) CONTAINS THE PROBABILITY CORRESPONDING
C                           TO STAT(3). STAT(3) IS DEFINED ONLY WHEN
C                           N(1)=3.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES N(1) IS NOT 2 OR 3
C                           IER=130 INDICATES N(2),N(3), OR N(4) (IF
C                             N(1)=3) IS LESS THAN 2
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
C   REMARKS      IF THE DEGREES OF FREEDOM, NDF, ARE SMALL,
C                PROBABILITIES PRODUCED WILL BE VERY APPROXIMATE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NAWNRP (Y,N,IR,NF,NDF,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N(1),IR(1),NF(1),NDF(1),IER
      REAL               Y(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IER2,L,LL,NAA,NALL,NB,NBB,NC,NEL,NFF,NR,NTOT,
     1                   MIDEL
      REAL               CHI,CHIA,DF,EXPC,XNA,XNB,YMED,YTOT
C                                  FIRST EXECUTABLE STATEMENT
      NFF = N(1)
C                                  CHECK N(1)
      IF ((NFF .GE. 2) .AND. (NFF .LE. 3)) GO TO 5
      IER = 129
      GO TO 9000
C                                  CHECK N(2) AND N(3)
    5 NR = N(2)
      NC = N(3)
      IF ((NR .GE. 2) .AND. (NC .GE. 2)) GO TO 10
      IER = 130
      GO TO 9000
C                                  CHECK N(4)
   10 IF (NFF .EQ. 2) GO TO 20
      NB = N(4)
      IF (NB .GE. 2) GO TO 15
      IER = 130
      GO TO 9000
   15 NALL = NB
      GO TO 25
   20 NB = 1
      NALL = 0
      N(4) = 1
   25 IER = 0
      NALL = NR + NC + NALL
      NTOT = NR * NC * NB
C                                  INITIALIZE WORKING VECTOR
      DO 30 I=1,NTOT
         IR(I) = I
   30 CONTINUE
C                                  INVOKE VSRTR FOR SORTING
      CALL VSRTR (Y,NTOT,IR)
C                                  DETERMINE MEDIAN
      MIDEL = NTOT/2 + 1
      YMED = Y(MIDEL)
C                                  INITIALIZE FREQUENCY VECTOR
      DO 35 I=1,NALL
         NF(I) = 0
   35 CONTINUE
C                                  LOOP TO COUNT Y BELOW MEDIAN AND TO
C                                  ACCUMULATE FOR EACH LEVEL OF FACTORS
      NBB = 0
      LL = NTOT/NR
      MIDEL = MIDEL - 1
      DO 40 I=1,MIDEL
C                                  CHECK FOR TIE WITH MEDIAN
         IF (Y(I) .EQ. YMED) GO TO 45
         NBB = NBB + 1
         L = IR(I)
C                                  DETERMINE ROW LEVEL OF L-TH Y
         NEL = L/LL
         IF ((NEL*LL) .NE. L) NEL=NEL+1
         NF(NEL) = NF(NEL) + 1
C                                  DETERMINE COLUMN LEVEL OF L-TH Y
         NEL = L/NB
         IF ((NEL*NB) .NE. L) NEL=NEL+1
         NEL = MOD(NEL,NC)
         IF (NEL .EQ. 0) NEL = NC
         NF(NR+NEL) = NF(NR+NEL) + 1
C                                  DETERMINE BLOCK LEVEL OF L-TH Y
         IF (NB .EQ. 1) GO TO 40
         NEL = MOD(L,NB)
         IF (NEL .EQ. 0) NEL = NB
         NF(NR+NC+NEL) = NF(NR+NC+NEL) + 1
   40 CONTINUE
   45 NAA = NTOT - NBB
C                                  ASSIGN DEGREES OF FREEDOM
      NDF(1) = NR-1
      NDF(2) = NC-1
      IF (NFF .EQ. 3) NDF(3)=NB-1
      IF (NAA .NE. NBB) GO TO 60
C                                  LOOP TO DETERMINE CHI-SQ AND PROB
C                                  WHEN MEDIAN SPLITS DATA EXACTLY
      NEL = 0
      YTOT = NTOT
      DO 55 I=1,NFF
C                                  EXPECTED FREQUENCY FOR THIS FACTOR
         EXPC = YTOT * 0.5 / N(I+1)
         CHI = 0.0
         LL = NEL + 1
         NEL = LL + N(I+1) - 1
C                                  CHI-SQ ACCUMULATION
         DO 50 L=LL,NEL
            CHI = CHI + (NF(L)-EXPC)**2
   50    CONTINUE
         CHI =  2.0/EXPC * CHI
         STAT(I) = CHI
C                                  INVOKE MDCH FOR PROBABILITY
         DF = NDF(I)
         CALL MDCH(CHI,DF,CHIA,IER2)
         STAT(I+3) = 1.-CHIA
   55 CONTINUE
      GO TO 9005
C                                  LOOP FOR CHI-SQ AND PROBABILITY WHEN
C                                  SPLIT BY MEDIAN IS NOT EVEN
   60 XNA = NAA
      XNB = NBB
      NALL = 0
      NEL = NR * NB
      N(1) = NC
      DO 70 I=1,NFF
C                                  EXPECTED COUNTS ABOVE AND BELOW MED.
         NEL = NEL/N(I+1) * N(I)
         YTOT = XNA/N(I+1)
         EXPC = XNB/N(I+1)
         CHIA = 0.0
         CHI = 0.0
         LL = NALL + 1
         NALL = LL + N(I+1) - 1
C                                  CHI-SQ ACCUMULATION
         DO 65 L=LL,NALL
            NTOT = NEL - NF(L)
            CHIA = CHIA + (NTOT-YTOT)**2
            CHI = CHI + (NF(L)-EXPC)**2
   65    CONTINUE
         CHI = CHI/EXPC + CHIA/YTOT
         STAT(I) = CHI
C                                  INVOKE MDCH FOR PROBABILITY
         DF = NDF(I)
         CALL MDCH(CHI,DF,CHIA,IER2)
         STAT(I+3) = 1.-CHIA
   70 CONTINUE
      N(1) = NFF
      GO TO 9005
 9000 CONTINUE
       CALL UERTST(IER,6HNAWNRP)
 9005 RETURN
      END

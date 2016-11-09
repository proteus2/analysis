C   IMSL ROUTINE NAME   - RLEAP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - LEAPS AND BOUNDS ALGORITHM FOR DETERMINING
C                           A NUMBER OF BEST REGRESSION SUBSETS FROM
C                           A FULL REGRESSION MODEL
C
C   USAGE               - CALL RLEAP (RR,KZ,IJOB,IXS,STAT,IXV,NVAR,IXB,
C                           BEST,IB,WK,IW,IER)
C
C   ARGUMENTS    RR     - INPUT VECTOR OF LENGTH (2*KZ**3+4*KZ)/3.
C                           THE FIRST KZ*(KZ+1)/2 LOCATIONS OF RR
C                           CONTAIN THE CORRELATION MATRIX OR SUMS OF
C                           SQUARES AND CROSS PRODUCTS MATRIX FOR THE
C                           KZ VARIABLES IN SYMMETRIC STORAGE MODE.
C                           VARIABLE KZ IS THE DEPENDENT VARIABLE. THE
C                           REMAINDER OF THE VECTOR IS USED AS WORK
C                           STORAGE. SEE REMARKS.
C                KZ     - INPUT NUMBER OF VARIABLES (NUMBER OF
C                           INDEPENDENT VARIABLES PLUS ONE). KZ MUST
C                           BE GREATER THAN 3 AND LESS THAN OR EQUAL
C                           TO IJOB(1). SEE THE DESCRIPTION OF IJOB
C                           BELOW.
C                IJOB   - INPUT AND OUTPUT VECTOR OF LENGTH 4.
C                         IJOB(1) CONTAINS THE INPUT NUMBER OF DEGREES
C                           OF FREEDOM FOR RR (NUMBER OF DATA POINTS
C                           MINUS ONE). IJOB(1) MUST BE GREATER THAN OR
C                           EQUAL TO KZ.
C                         IJOB(2) CONTAINS THE INPUT MODEL SELECTION
C                           CRITERION.
C                           IF IJOB(2) EQUALS 1 OR A NEGATIVE INTEGER
C                           THEN THE R-SQUARED CRITERION IS USED. IF
C                           IJOB(2)=-K, WHERE K IS A POSITIVE INTEGER
C                           LESS THAN KZ, THE MAXIMUM NUMBER OF
C                           INDEPENDENT VARIABLES TO BE CONSIDERED
C                           IN ANY MODEL IS K.
C                           IF IJOB(2) EQUALS 2, THEN THE ADJUSTED
C                           R-SQUARED IS USED.
C                           IF IJOB(2) EQUALS 3, THEN THE MALLOWS
C                           CP STATISTIC IS USED.
C                           THE THREE CRITERIA ARE SIMPLE FUNCTIONS OF
C                           THE RESIDUAL SUM OF SQUARES.
C                         IJOB(3) CONTAINS, ON INPUT, THE MAXIMUM NUMBER
C                           OF BEST REGRESSIONS TO BE COMPUTED.
C                           IF IJOB(2) = 1, THE BEST IJOB(3) OF EACH
C                           SUBSET SIZE WILL BE RETURNED.
C                           OTHERWISE, THE IJOB(3) BEST OVERALL WILL BE
C                           RETURNED. IJOB(3) MUST BE GREATER THAN ZERO.
C                           ON OUTPUT, IJOB(3) CONTAINS THE TOTAL
C                           NUMBER OF BEST SUBSETS FOUND.
C                         IJOB(4) CONTAINS THE INPUT MAXIMUM NUMBER OF
C                           GOOD REGRESSIONS OF EACH SUBSET SIZE TO BE
C                           SAVED IN FINDING THE IJOB(3) BEST
C                           REGRESSIONS. IJOB(4) MUST BE GREATER THAN
C                           OR EQUAL TO IJOB(3). NORMALLY IJOB(4)
C                           SHOULD BE LESS THAN OR EQUAL TO 10. IT
C                           NEED NOT EVER BE SET LARGER THAN THE
C                           MAXIMUM NUMBER OF SUBSETS FOR ANY SUBSET
C                           SIZE. COMPUTING TIME REQUIRED IS INVERSELY
C                           RELATED TO IJOB(4).
C                IXS    - OUTPUT VECTOR OF LENGTH KZ CONTAINING THE
C                           LOCATION OF THE FIRST ELEMENT FOR EACH
C                           SUBSET SIZE IN STAT. SEE THE PROGRAMMING
C                           NOTES IN THE MANUAL DOCUMENT FOR FURTHER
C                           DETAILS.
C                STAT   - OUTPUT VECTOR OF LENGTH THE MAXIMUN OF KZ
C                           AND IJOB(4)*(KZ-1) CONTAINING THE CRITERION
C                           VALUES FOR EACH SUBSET CONSIDERED, IN
C                           INCREASING SUBSET SIZE ORDER.  WITHIN EACH
C                           SUBSET SIZE, RESULTS ARE RETURNED IN MONO-
C                           TONE ORDER ACCORDING TO CRITERION VALUE,
C                           WITH THE RESULTS FOR THE BETTER REGRESSIONS
C                           LISTED FIRST.  SEE THE PROGRAMMING NOTES IN
C                           THE MANUAL DOCUMENT FOR FURTHER DETAILS.
C                IXV    - OUTPUT VECTOR OF LENGTH KZ CONTAINING THE
C                           LOCATIONS OF THE FIRST ELEMENT FOR EACH
C                           SUBSET SIZE IN NVAR. SEE THE PROGRAMMING
C                           NOTES IN THE MANUAL DOCUMENT FOR FURTHER
C                           DETAILS. THE LAST LOCATION IS WORK STORAGE.
C                NVAR   - OUTPUT VECTOR OF LENGTH IJOB(4)*(KZ-1)*KZ/2
C                           CONTAINING THE VARIABLE NUMBERS FOR EACH
C                           SUBSET CONSIDERED AND IN THE SAME ORDER AS
C                           RETURNED IN STAT.
C                IXB    - OUTPUT VECTOR CONTAINING THE ROW NUMBER OF
C                           THE FIRST ROW FOR EACH SUBSET IN BEST.
C                           IF IJOB(2) IS EQUAL TO 2 OR 3, THEN THE
C                           LENGTH OF IXB IS IJOB(3)+1 (IJOB(3) AS SET
C                           ON INPUT). IF IJOB(2) EQUALS 1, THEN THE
C                           LENGTH OF IXB IS IJOB(3)*(KZ-1)+1 (IJOB(3)
C                           AS SET ON INPUT). SEE THE PROGRAMMING NOTES
C                           IN THE MANUAL DOCUMENT FOR FURTHER DETAILS.
C                BEST   - OUTPUT MATRIX CONTAINING THE RESULTS FOR BEST
C                           REGRESSIONS BEGINNING WITH ONE VARIABLE
C                           REGRESSIONS AND INCREASING THE SUBSET SIZE.
C                           COLUMNS 1,2,3, AND 4 CONTAIN VARIABLE
C                           NUMBER, REGRESSION COEFFICIENT, F VALUE, AND
C                           TAIL AREA OF THE F DISTRIBUTION
C                           RESPECTIVELY.
C                           IF IJOB(2) IS EQUAL TO 1, THEN BEST IS
C                           OF DIMENSION IJOB(3)*(KZ-1)*KZ/2 BY 4
C                           (IJOB(3) AS SET ON INPUT). IF IJOB(2) IS
C                           EQUAL TO 2 OR 3, THEN BEST IS OF DIMENSION
C                           IJOB(3)*(KZ-1) BY 4 (IJOB(3) AS SET ON
C                           INPUT). SEE THE PROGRAMMING NOTES IN THE
C                           MANUAL DOCUMENT FOR FURTHER DETAILS.
C                IB     - INPUT ROW DIMENSION OF THE MATRIX BEST
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                WK     - REAL WORK AREA OF LENGTH KZ*(2*IJOB(4)+6).
C                IW     - INTEGER WORK AREA OF LENGTH 3*KZ**2+9*KZ+12.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES KZ OR AN ELEMENT OF IJOB
C                             WAS SPECIFIED INCORRECTLY.
C                           IER=130 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL SUBROUTINE MDFD.
C                           IER=131 INDICATES THAT THE ERROR SUM OF
C                             SQUARES WAS NON-POSITIVE.
C                           IER=132 INDICATES DELETION OF VARIABLES
C                             (BECAUSE OF SINGULARITY OF RR) RESULTED IN
C                             FEWER THAN THREE INDEPENDENT VARIABLES
C                             REMAINING. THE FIRST KZ-1 ELEMENTS OF IW
C                             INDICATE, BY VARIABLE NUMBER, WHICH
C                             VARIABLES WERE DELETED.
C                         WARNING ERROR
C                           IER=37 INDICATES AT LEAST ONE VARIABLE WAS
C                             DELETED BECAUSE OF SINGULARITY OF RR, BUT
C                             THREE OR MORE WERE STILL AVAILABLE. IF
C                             ANY IW(1),...,IW(KZ-1) ARE NON-ZERO, THEN
C                             THE NON-ZERO ELEMENT GIVES THE VARIABLE
C                             NUMBER WHICH CAUSED DELETION FOR
C                             SINGULARITY REASONS. THOUGH A VARIABLE
C                             MAY BE DELETED AND THUS NOT APPEAR IN THE
C                             FULL MODEL, IT MAY ENTER IN SOME OF THE
C                             SMALLER SUBMODELS.
C
C   REQD. IMSL ROUTINES - MDFD,MERRC=ERFC,RLEAP1,RLEAP2,RLEAP3,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      INPUT VECTOR RR MAY BE COMPUTED FROM THE RAW DATA
C                MATRIX USING IMSL ROUTINES IN CHAPTER B. THE
C                RESULTANT RR WILL BE STORED PROPERLY FOR INPUT TO
C                RLEAP. IF SOME OTHER ROUTINE WAS USED TO PRODUCE RR
C                (NOT IN SYMMETRIC STORAGE MODE), IMSL ROUTINES IN
C                CHAPTER V FOR STORAGE MODE CONVERSION MAY BE USEFUL.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLEAP  (RR,KZ,IJOB,IXS,STAT,IXV,NVAR,IXB,BEST,IB,WK,
     *                   IW,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            KZ,IJOB(4),IXS(1),IXV(1),NVAR(1),IXB(1),IB,
     *                   IW(KZ,1),IER
      REAL               RR(1),STAT(1),BEST(IB,4),WK(KZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBIT,ICL,II,IIER,IJ,IK,ILI,ILN,IN,INC,INDEX1,
     *                   INDEX2,INDEX3,INDEX4,INDI,IP,IS,IWK,IXC,I1,I2,
     *                   J,JC,JER,K,KM,KO,KP,KX,KZR,KZ1,K5,L,LA,LB,LC,
     *                   LI,LL,LN,LOW,LS,M,MBST,MI,MMM,MMMP1,MN,MP,MT,
     *                   MV,N,NCOF,NDEF,NX
      REAL               F,P
      REAL               CAB,TOL,RS,BOUND,ZERO,T2,ONE,HUND,TWO,TWODF,
     *                   VAR,R2,SIGM,XNDEF,SS,SIG
      DATA               TOL /1.E-4/
      DATA               HUND /100./,ZERO /0./,TWO /2./,ONE /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  TEST INPUT
      NDEF = IJOB(1)
      IBIT = IJOB(2)
      MMM = KZ-1
      IF (IBIT.GT.0) GO TO 5
      MMM = -IBIT
      IBIT = 1
    5 MBST = IJOB(3)
      KO = IJOB(4)
      KX = KZ-1
      IWK = KZ
      IF (KX.GE.3 .AND. NDEF.GT.KX .AND. MBST.GT.0 .AND. MBST.LE.KO
     *.AND. IBIT.GE.1 .AND. IBIT.LE.3) GO TO 10
C                                  TERMINAL ERROR - INPUT PARAMETER IS
C                                  OUT OF BOUNDS
      IER = 129
      GO TO 9000
C                                  INITIALIZE
   10 ILI = 5+KZ
      ILN = ILI+KZ
      ICL = KO+6
      K5 = KO+5
      IXC = (KZ**3+KZ+KZ)/3
      IIER = 9+3*KX
      LS = (KZ*(KZ+1))/2
      SS = RR(LS)/HUND
      IF (IBIT.EQ.2) SS = SS/NDEF
      WK(1,1) = ZERO
      WK(1,2) = ZERO
      IJ = LS+1
      IK = 3
      IW(1,5) = 0
      DO 15 I=2,KZ
         IW(I,5) = IJ-IK+IW(I-1,5)
         IK = IK+I+1
   15 CONTINUE
      TWODF = (NDEF+NDEF)*RR(LS)
      LOW = KO-MBST+1
      I2 = LS*KX
      I2 = I2+I2
      T2 = TWO
      KZ1 = KZ+1
      IK = 0
      MMMP1 = MMM+1
      DO 25 LL=1,KZ
         L = KZ1-LL
         IW(L,2) = 1
         IW(1,ILI+L) = L
         WK(L,ICL) = -TWODF
         WK(L,3) = T2
         T2 = T2+T2
         IW(1,L+5) = L
         IK = IK+LL
         STAT(LL) = TOL*RR(IK)
         DO 20 M=1,KO
            WK(L,ICL+M) = ZERO
            WK(L,M+5) = TWODF
   20    CONTINUE
   25 CONTINUE
C                                  SAVE RR
      DO 30 I=1,LS
         RR(IXC+I) = RR(I)
   30 CONTINUE
C                                  INVERT MATRIX STEPWISE
      KZR = (KX*KZ)/2+IXC
      DO 40 N=1,KX
         DO 35 LA=N,KX
            LL = IW(1,ILI+LA)
            L = (LL*(LL-1))/2
            IF (RR(L+LL+IXC).LE.STAT(LL)) GO TO 35
            RS = RR(IXC+LS)-RR(KZR+LL)*RR(KZR+LL)/RR(IXC+L+LL)
            IF (RS.LT.WK(N,K5)) J = LA
            IF (RS.LT.WK(N,6)) CALL RLEAP3(RS,WK(1,2)+WK(LL,3),KO,WK(1,
     *      ICL+1),WK(1,6),N,IWK)
   35    CONTINUE
         IF (WK(N,K5).EQ.TWODF) GO TO 45
         M = IW(1,ILI+J)
         IW(1,ILI+J) = IW(1,ILI+N)
         IW(1,ILI+N) = M
         IW(1,ILN+N) = M
         WK(1,2) = WK(1,2)+WK(M,3)
         CALL RLEAP2(RR(IXC+1),KZ,M)
   40 CONTINUE
      N = KZ
   45 K = N-1
      KP = K+1
      IF (K.EQ.KX) GO TO 60
      IER = 37
      DO 50 I=1,KP
         IW(I,IIER) = 0
   50 CONTINUE
      DO 55 I=KP,KX
         IW(I,IIER) = IW(1,ILI+I)
   55 CONTINUE
      IF (K.LT.3) GO TO 200
   60 KM = K-1
      SIG = (RR(LS+IXC)+RR(LS+IXC))/(NDEF-K)
      WK(1,4) = RR(LS+IXC)
      WK(1,5) = RR(LS)
      IW(1,3) = K
      IW(1,4) = K
      IF (IBIT.EQ.1) GO TO 75
      DO 70 M=1,K
         SIGM = SIG*M
         XNDEF = ONE/(NDEF-M)
         DO 65 L=1,KO
            IF (IBIT.EQ.2) RS = WK(M,L+5)*XNDEF
            IF (IBIT.EQ.3) RS = WK(M,L+5)+SIGM
            IF (RS.LT.WK(KZ,6)) CALL RLEAP3(RS,WK(M,ICL+L),KO,WK(1,
     *      ICL+1),WK(1,6),KZ,IWK)
   65    CONTINUE
   70 CONTINUE
   75 MN = 2
      MV = -1
C                                  STAGE  LOOP
   80 IF (MN.EQ.1) GO TO 130
      IP = IW(MN,2)
      IW(MN,2) = IP+1
      MV = MV-IW(MN+1,2)+IP+2
      IW(MV,1) = IP
      MN = MN-1
      IN = IW(MN,2)
      JC = MV
      BOUND = WK(IP,4)
      WK(IP,4) = TWODF
C                                  FIND LEAP FROM BOUNDS
      DO 85 LB=IP,KM
         MT = MN+KM-LB
         IF (IBIT.EQ.1 .AND. WK(MT,LOW+5).GT.BOUND) GO TO 90
         IF (IBIT.EQ.2 .AND. WK(KZ,LOW+5).GT.BOUND/(NDEF-MT)) GO TO 90
         IF (IBIT.EQ.3 .AND. WK(KZ,LOW+5).GT.BOUND+SIG*MT) GO TO 90
   85 CONTINUE
      GO TO 80
   90 LC = KM+IP-LB
      IF (IP.EQ.1) LC = K
C                                  REGRESSIONS FROM INVERSE MATRIX
      DO 110 LB=IP,LC
         IS = LB+1
         CALL RLEAP1(IW(1,6),LB,LI,IW(1,1),MV,RS,BOUND,IW(1,ILI+1),JC,
     *   IW(1,5),RR(IXC+1),1,IW(1,3),IWK,KZ)
C                                  RE-ORDER VARIABLES
         M = LB
         IF (LB.GT.IW(IN,4)) GO TO 105
         LN = IW(IN,ILN+LB)
   95    IF (RS.LE.WK(M,4)) GO TO 100
         WK(M+1,4) = WK(M,4)
         IW(IP,ILI+M) = IW(IP,ILI+M-1)
         IW(IN,ILN+M) = IW(IN,ILN+M-1)
         M = M-1
         GO TO 95
  100    IW(IP,ILI+M) = LI
         IW(IN,ILN+M) = LN
  105    WK(M+1,4) = RS
         IW(IS,3) = LB
         IW(IS,4) = LB
  110 CONTINUE
      IF (LC.EQ.K) LC = KM
      MI = K-MV
      JC = MN
C                                  REGRESSIONS FROM PRODUCT MATRIX
      SIGM = SIG*MI
      XNDEF = ONE/(NDEF-MI)
      DO 125 LB=IP,LC
         IS = LB+1
         CALL RLEAP1(IW(1,6),LB,L,IW(1,2),MN,WK(IS,5),WK(IN,5),IW(1,
     *   ILN+1),JC,IW(1,5),RR,0,IW(1,4),IWK,KZ)
         INC = IW(IN,L+5)
         WK(IS,2) = WK(IP,2)-WK(INC,3)
         WK(IS,1) = WK(IN,1)+WK(INC,3)
         IF (WK(IS,4).GE.WK(MI,6)) GO TO 115
         CALL RLEAP3(WK(IS,4),WK(IS,2),KO,WK(1,ICL+1),WK(1,6),MI,IWK)
         IF (IBIT.EQ.1) GO TO 115
         IF (IBIT.EQ.2) RS = WK(IS,4)*XNDEF
         IF (IBIT.EQ.3) RS = WK(IS,4)+SIGM
         IF (RS.LT.WK(KZ,6)) CALL RLEAP3(RS,WK(IS,2),KO,WK(1,ICL+1),
     *   WK(1,6),KZ,IWK)
  115    IF (WK(IS,5).GE.WK(MN,6)) GO TO 120
         CALL RLEAP3(WK(IS,5),WK(IS,1),KO,WK(1,ICL+1),WK(1,6),MN,IWK)
         IF (IBIT.EQ.1) GO TO 120
         IF (IBIT.EQ.2) RS = WK(IS,5)/(NDEF-MN)
         IF (IBIT.EQ.3) RS = WK(IS,5)+MN*SIG
         IF (RS.LT.WK(KZ,6)) CALL RLEAP3(RS,WK(IS,1),KO,WK(1,ICL+1),
     *   WK(1,6),KZ,IWK)
  120    MN = MN+1
         IW(MN+1,2) = IW(MN,2)+1
         IN = IS
  125 CONTINUE
      IF (LC.EQ.KM) MN = MN-1
      GO TO 80
C                                  OUTPUT
  130 INDEX1 = 1
      INDEX2 = 1
      INDEX3 = 1
      INDEX4 = 1
      DO 175 M=1,MMM
         SIGM = SIG*M
         NX = NDEF-M
         XNDEF = ONE/NX
         IXS(M) = INDEX1
         IXV(M) = INDEX2
         DO 170 LA=1,KO
            NCOF = 1
            L = KO-LA+1
            IF (WK(M,L+5).EQ.TWODF) GO TO 175
            IF (IBIT.EQ.1) R2 = HUND-WK(M,L+5)/SS
            IF (IBIT.EQ.2) RS = WK(M,L+5)*XNDEF
            IF (IBIT.EQ.3) RS = WK(M,L+5)+SIGM
            IF (IBIT.EQ.1 .AND. LA.LE.MBST .OR. IBIT.GT.1 .AND.
     *      RS.LE.WK(KZ,LOW+5)) NCOF = 0
            IF (IBIT.EQ.2) R2 = HUND-RS/SS
            IF (IBIT.EQ.3) R2 = (RS+RS)/SIG-NDEF+ONE
C                                  DECODE LABELS
            CAB = WK(M,ICL+L)
            MP = 1
            DO 135 I=1,KX
               IF (CAB.LT.WK(I,3)) GO TO 135
               IW(MP,2) = I
               MP = MP+1
               CAB = CAB-WK(I,3)
  135       CONTINUE
            IF (NCOF.NE.0) GO TO 160
            STAT(INDEX1) = R2
            INDEX1 = INDEX1+1
            IXB(INDEX4) = INDEX3
            IJOB(3) = INDEX4
            INDEX4 = INDEX4+1
C                                  FORM SUBMATRIX
            IW(MP,2) = KZ
            II = 0
            DO 145 I=1,MP
               II = II+I
               IJ = II
               DO 140 J=I,MP
                  IK = MAX0(IW(I,2),IW(J,2))
                  I1 = MIN0(IW(I,2),IW(J,2))
                  INDI = (IK*(IK-1))/2+I1
                  RR(IXC+IJ) = RR(INDI)
                  IJ = IJ+J
  140          CONTINUE
  145       CONTINUE
C                                  INVERT THE SUBMATRIX
            DO 150 N=1,M
               CALL RLEAP2(RR(IXC+1),MP,N)
  150       CONTINUE
            IJ = (MP*(MP+1))/2
            VAR = RR(IJ+IXC)*XNDEF
            INDI = 0
            DO 155 I=1,M
               INDI = INDI+I
               IJ = (MP*(MP-1))/2+I
               BEST(INDEX3,2) = -RR(IXC+IJ)
               IF (RR(INDI+IXC)*VAR.GE.ZERO) GO TO 195
               BEST(INDEX3,3) = -BEST(INDEX3,2)*BEST(INDEX3,2)/
     *         (RR(INDI+IXC)*VAR)
               F = BEST(INDEX3,3)
               NVAR(INDEX2) = IW(I,2)
               INDEX2 = INDEX2+1
               BEST(INDEX3,1) = IW(I,2)
               CALL MDFD(F,1,NX,P,JER)
               IF (JER.NE.0) GO TO 190
               BEST(INDEX3,4) = 1.-P
               INDEX3 = INDEX3+1
               IXB(INDEX4) = INDEX3
  155       CONTINUE
            GO TO 170
  160       STAT(INDEX1) = R2
            INDEX1 = INDEX1+1
            DO 165 I=1,M
               NVAR(INDEX2) = IW(I,2)
               INDEX2 = INDEX2+1
  165       CONTINUE
  170    CONTINUE
  175 CONTINUE
C                                  ZERO OUT UNUSED ELEMENTS OF IXS AND I
      IF (K.EQ.KZ) GO TO 185
      K = K+1
      DO 180 I=K,KZ
         IXS(I) = 0.0
         IXV(I) = 0.0
  180 CONTINUE
  185 IF (IER-37) 9005, 205, 9005
  190 IER = 130
      GO TO 9000
  195 IER = 131
      GO TO 9000
  200 IER = 132
  205 DO 210 I=1,KX
         IW(I,1) = IW(I,IIER)
  210 CONTINUE
 9000 CONTINUE
      CALL UERTST(IER,6HRLEAP )
 9005 RETURN
      END

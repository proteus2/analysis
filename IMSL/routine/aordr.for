C   IMSL ROUTINE NAME   - AORDR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REORDERING OF THE DATA OBTAINED FROM ANY
C                           BALANCED COMPLETE EXPERIMENTAL DESIGN
C
C   USAGE               - CALL AORDR (NF,NL,IORD,YIN,JORD,YOUT,IWK,IER)
C
C   ARGUMENTS    NF     - INPUT NUMBER OF FACTORS (NUMBER OF DISTINCT
C                           SUBSCRIPTS IN THE MODEL), INCLUDING A
C                           REPLICATION FACTOR IF PRESENT. NF MUST BE
C                           GREATER THAN OR EQUAL TO TWO.
C                NL     - INPUT VECTOR OF LENGTH NF CONTAINING THE
C                           NUMBER OF LEVELS FOR EACH OF THE NF FACTORS.
C                           NL(I) MUST BE GREATER THAN OR EQUAL TO TWO
C                           FOR I=1,2,...,NF.
C                IORD   - INPUT VECTOR OF LENGTH NF INDICATING THE
C                           ORDERING OF THE DATA IN VECTOR YIN.
C                           IORD(I) = J IMPLIES THAT THE MODEL SUB-
C                           SCRIPT CORRESPONDING TO FACTOR I IS
C                           ALTERING J-TH MOST RAPIDLY. J = 1,2,...,NF
C                YIN    - INPUT VECTOR OF LENGTH NL(1)*NL(2)*...*NL(NF)
C                           CONTAINING THE DATA IN THE PATTERN IMPLIED
C                           BY VECTOR IORD.
C                JORD   - INPUT VECTOR OF LENGTH NF INDICATING THE
C                           DESIRED ORDERING OF THE DATA IN VECTOR YOUT
C                           (THE REORDERED YIN VECTOR). IORD(K) = L
C                           IMPLIES THAT THE MODEL SUBSCRIPT CORRESPOND-
C                           ING TO FACTOR K IS ALTERING L-TH MOST
C                           RAPIDLY. L = 1,2,...,NF
C                YOUT   - OUTPUT VECTOR OF LENGTH NL(1)*NL(2)*...*NL(NF)
C                           CONTAINING THE DATA IN THE PATTERN IMPLIED
C                           BY VECTOR JORD.
C                IWK    - WORK VECTOR OF LENGTH 4*NF.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NF WAS SPECIFIED LESS
C                             THAN TWO.
C                           IER=130 INDICATES THAT SOME NL(I) WAS
C                             SPECIFIED LESS THAN TWO FOR I=1,2,...,NF.
C                           IER=131 INDICATES THAT IORD OR JORD OR BOTH
C                             DID NOT CONTAIN THE INTEGERS 1,2,...,NF.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AORDR  (NF,NL,IORD,YIN,JORD,YOUT,IWK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NF,NL(NF),IORD(NF),JORD(NF),IWK(1),IER
      REAL               YIN(1),YOUT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ID,II,IJ,ITOT,J,JJ,KK,LJ,ND,NF2,NF3,NJ,NK,
     1                   NT,NT1
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NF.LT.2) GO TO 45
      DO 5 I=1,NF
         IF (NL(I).LT.2) GO TO 50
         IF (IORD(I).LT.1.OR.IORD(I).GT.NF) GO TO 55
         IF (JORD(I).LT.1.OR.JORD(I).GT.NF) GO TO 55
    5 CONTINUE
      JJ = NF-1
      DO 15 I=1,JJ
         J = I+1
         DO 10 IJ=J,NF
            IF (IORD(I).EQ.IORD(IJ)) GO TO 55
            IF (JORD(I).EQ.JORD(IJ)) GO TO 55
   10    CONTINUE
   15 CONTINUE
C                                  REORDER NL SUBJECT TO IORD AND JORD
C                                  AND STORE THE RESULTS BEGINNING AT
C                                  IWK(1) AND IWK(NF+1), RESPECTIVELY.
C                                  IORD IS REORDERED SUBJECT TO NL AND
C                                  STORED BEGINNING AT IWK(2*NF+1).
      NF2 = NF+NF
      NF3 = NF2+NF
      DO 20 I=1,NF
         II = IORD(I)
         IWK(II) = NL(I)
         JJ = JORD(I)+NF
         IWK(JJ) = NL(I)
         KK = IORD(I)+NF2
         IWK(KK) = I
   20 CONTINUE
C                                  COMPUTE THE TOTAL NUMBER OF DATA
C                                  POINTS.
      NT = 1
      DO 25 I=1,NF
         NT = NT*NL(I)
   25 CONTINUE
C                                  COMPUTE MODEL SUBSCRIPTS ASSOCIATED
C                                  WITH YIN(M), M=1,2,...,NT.
      NT1 = NT-1
      DO 40 I=2,NT1
         LJ = I
C                                  REORDER JORD SUBJECT TO NL AND STORE
C                                  THE RESULTS BEGINNING AT IWK(3*NF+1)
         DO 30 J=1,NF
            ID = IWK(NF2+J)
            ND = NF3+JORD(ID)
            IWK(ND) = MOD(LJ,IWK(J))
            IF (IWK(ND).EQ.0) IWK(ND) = IWK(J)
            LJ = (LJ-IWK(ND))/IWK(J)+1
   30    CONTINUE
C                                  COMPUTE THE CORRESPONDING SUBSCRIPT
C                                  IN YOUT.
         ITOT = IWK(NF3+1)
         KK = 1
         DO 35 J=2,NF
            NK = NF+J-1
            KK = KK*IWK(NK)
            NJ = NF3+J
            ITOT = ITOT+(IWK(NJ)-1)*KK
   35    CONTINUE
         IF (ITOT.LT.2.OR.ITOT.GT.NT1) GO TO 45
         YOUT(ITOT) = YIN(I)
   40 CONTINUE
      YOUT(1) = YIN(1)
      YOUT(NT) = YIN(NT)
      GO TO 9005
   45 IER = 129
      GO TO 9000
   50 IER = 130
      GO TO 9000
   55 IER = 131
 9000 CONTINUE
      CALL UERTST (IER,'AORDR ')
 9005 RETURN
      END

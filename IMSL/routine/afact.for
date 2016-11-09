C   IMSL ROUTINE NAME   - AFACT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - FULL FACTORIAL PLAN ANALYSIS - EASY
C                           TO USE VERSION
C
C   USAGE               - CALL AFACT (NF,NL,IJOB,Y,GMEAN,YMEANS,INDEX,
C                           STAT,IS,IER)
C
C   ARGUMENTS    NF     - INPUT NUMBER OF FACTORS (NUMBER OF DISTINCT
C                           SUBSCRIPTS IN THE MODEL), INCLUDING A
C                           FACTOR FOR REPLICATE OBSERVATIONS PER CELL,
C                           IF PRESENT.
C                NL     - INPUT VECTOR OF LENGTH NF CONTAINING THE
C                           NUMBER OF LEVELS FOR EACH OF THE NF FACTORS.
C                IJOB   - INPUT OPTION PARAMETER
C                           IF IJOB IS ZERO, ONLY ONE OBSERVATION PER
C                             CELL IS AVAILABLE.
C                           IF IJOB IS NONZERO, THE LAST FACTOR
C                             (SUBSCRIPT) IS A REPLICATION EFFECT.
C                Y      - INPUT VECTOR OF LENGTH 2**NF+MS CONTAINING
C                           THE MT OBSERVATIONS BEGINNING IN LOCATION
C                           2**NF+1, WHERE MS IS (NL(1)+1)*(NL(2)+1)*...
C                           *(NL(NF)+1) AND MT IS NL(1)*NL(2)*...
C                           *NL(NF).  THE REMAINING LOCATIONS ARE WORK
C                           STORAGE.  INPUT Y IS DESTROYED.  PRIOR TO
C                           CALLING AFACT, IMSL ROUTINE AORDR MAY BE
C                           USED TO REORDER Y AS DESIRED.
C                GMEAN  - OUTPUT GRAND MEAN
C                YMEANS - OUTPUT VECTOR OF LENGTH MS-MT-1 IF IJOB IS
C                           ZERO, OR MS/(NL(NF)+1)-1 IF IJOB IS
C                           NONZERO, CONTAINING THE FULL SET OF MEANS.
C                INDEX  - OUTPUT VECTOR OF LENGTH 2**(NF+4) INDICATING
C                           THE LOCATION OF THE FIRST MEAN FOR EACH SET
C                           OF EFFECT MEANS IN YMEANS.  IF IJOB IS ZERO,
C                           THE FIRST (2**NF)-2 LOCATIONS ARE DEFINED.
C                           IF IJOB IS NONZERO THE FIRST (2**(NF-1))-1
C                           LOCATIONS ARE DEFINED.  THE REMAINING
C                           LOCATIONS ARE WORK STORAGE.
C                STAT   - OUTPUT MATRIX OF DIMENSION 2**NF BY 4.  IF
C                           IJOB IS ZERO, ONLY THE FIRST 2**NF-1 ROWS
C                           ARE DEFINED.  IF IJOB IS NONZERO, ONLY THE
C                           FIRST 2**(NF-1) ROWS ARE DEFINED.  IN EITHER
C                           CASE COLUMN 4 IS WORK STORAGE.  COLUMNS 1,
C                           2, AND 3 CONTAIN THE DEGREES OF FREEDOM,
C                           SUMS OF SQUARES, AND MEAN SQUARES,
C                           RESPECTIVELY.
C                IS     - INPUT ROW DIMENSION OF THE MATRIX STAT EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT A TERMINAL ERROR
C                             OCCURRED IN IMSL ROUTINE AFACN
C                           IER=130 INDICATES IS WAS SPECIFIED LESS
C                             THAN REQUIRED FOR STAT
C
C   REQD. IMSL ROUTINES - SINGLE/AFACN,UERTST,UGETIO
C                       - DOUBLE/AFACN,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AFACT  (NF,NL,IJOB,Y,GMEAN,YMEANS,INDEX,STAT,IS,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NF,NL(NF),IJOB,INDEX(1),IS,IER
      REAL               Y(1),GMEAN,YMEANS(1),STAT(IS,4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            KF,KM,KF1,KM1,KM11,NFP1,I,KM2,KS,J,K,ID,JJ,KV,
     1                   KF2,JK,JD,KK,KK1,LL,ISUM,IIF,LF,MT,MS,IT,IJ,
     2                   K1,KD,ISW
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      KF = 2**NF
      IF (IJOB.EQ.0) GO TO 5
      IF (IS.GE.KF/2) GO TO 10
      GO TO 110
    5 IF (IS.LT.KF) GO TO 110
C                                  FILL IN INDEX(II) WITH INFO FOR
C                                  AFACN
   10 KM = KF+KF
      KF1 = KF-1
      KM1 = KM-1
      KM11 = KM-1
      NFP1 = NF+1
      INDEX(KM1) = KF
      DO 20 I=2,NF
         KM1 = KM1-1
         KM2 = KM1
         INDEX(KM1) = 2**(NFP1-I)
         KS = KM11-KM1
         DO 15 J=1,KS
            KM1 = KM1-1
            INDEX(KM1) = INDEX(KM2)+INDEX(KM-J)
   15    CONTINUE
   20 CONTINUE
      IF (IJOB.EQ.0) GO TO 30
      DO 25 I=1,KF1,2
         K = KF+I
         INDEX(K) = -INDEX(K)
   25 CONTINUE
C                                  FILL IN INDEX(I) WITH FIRST MEAN
C                                  LOCATIONS FOR AFACT
   30 INDEX(1) = NL(1)
      INDEX(2) = NL(2)
      INDEX(3) = NL(1)*NL(2)
      IF (NF.EQ.2) GO TO 45
      ID = 3
      KM1 = 2
      DO 40 I=3,NF
         ID = ID+1
         INDEX(ID) = NL(I)
         KM1 = KM1+KM1
         JJ = KM1-1
         DO 35 J=1,JJ
            ID = ID+1
            INDEX(ID) = NL(I)*INDEX(J)
   35    CONTINUE
   40 CONTINUE
C                                  FILL IN INDEX(III) WITH FIRST MEAN
C                                  LOCATIONS FOR AFACN
   45 KV = KF/2
      KF2 = KF-2
      DO 60 I=1,KF1
         KK = 0
         JK = 0
         K = KF+I
         JD = INDEX(K)
   50    JK = JK+1
         KK = KK+2
         KK1 = KK
         IF (IJOB.EQ.0) GO TO 55
         IF (JK.GE.KV) KK1 = -KK1
   55    IF (JD.NE.KK1) GO TO 50
         LL = KM+I
         INDEX(LL) = INDEX(JK)
   60 CONTINUE
C                                  BACK TO INDEX(I)
      ISUM = 1
      DO 65 I=1,KF1
         K = INDEX(I)
         INDEX(I) = ISUM
         ISUM = ISUM+K
   65 CONTINUE
C                                  AND TO INDEX (III) AGAIN
      IIF = KM+KF1
      KM2 = KM+1
      LF = KM2
      ISUM = 1
      DO 70 I=LF,IIF
         K = INDEX(I)
         INDEX(I) = ISUM
         ISUM = ISUM+K
   70 CONTINUE
      MT = 1
      MS = 1
      DO 75 I=1,NF
         MT = MT*NL(I)
         MS = MS*(NL(I)+1)
   75 CONTINUE
      ISW = 0
      CALL AFACN (ISW,NF,NL,Y,STAT(1,4),INDEX,IER)
      IF (IER.NE.0) GO TO 105
      IF (IJOB.NE.0) KF2 = KV
      IT = 0
      DO 90 I=1,KF2
         IJ = INDEX(I+1)-INDEX(I)
         K1 = INDEX(I)-1
         IT = IT+2
         DO 85 J=1,KF1
            IF (IT.NE.INDEX(KF+J)) GO TO 85
            KK = INDEX(KM+J)+KF1
            DO 80 K=1,IJ
               YMEANS(K1+K) = Y(KK+K)
   80       CONTINUE
            GO TO 90
   85    CONTINUE
   90 CONTINUE
      GMEAN = Y(KF+MS)
      KD = KM+1
      ISW = 1
      CALL AFACN (ISW,NF,NL,Y,STAT(1,4),INDEX(KD),IER)
      IF (IER.NE.0) GO TO 105
      IT = 0
      DO 100 I=1,KF2
         IT = IT+2
         DO 95 J=1,KF1
            IF (IT.NE.INDEX(KF+J)) GO TO 95
            STAT(I,1) = INDEX(KM+J)
            STAT(I,2) = STAT(J,4)
            STAT(I,3) = STAT(I,2)/STAT(I,1)
            GO TO 100
   95    CONTINUE
  100 CONTINUE
      IF (IJOB.EQ.0) KF2 = KF1
      STAT(KF2,1) = INDEX(KM+KF)
      STAT(KF2,2) = STAT(KF,4)
      STAT(KF2,3) = STAT(KF2,2)/STAT(KF2,1)
      GO TO 9005
  105 IER = 129
      GO TO 9000
  110 IER = 130
 9000 CONTINUE
      CALL UERTST (IER,'AFACT ')
 9005 RETURN
      END

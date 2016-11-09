C   IMSL ROUTINE NAME   - FTAUTO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - MEAN, VARIANCE, AUTOCOVARIANCES,
C                           AUTOCORRELATIONS, AND PARTIAL
C                           AUTOCORRELATIONS FOR A STATIONARY
C                           TIME SERIES.
C
C   USAGE               - CALL FTAUTO (W,LW,K,L,ISW,AMEAN,VAR,ACV,AC,
C                           PACV,WKAREA)
C
C   ARGUMENTS    W      - INPUT VECTOR OF LENGTH LW CONTAINING THE TIME
C                           SERIES.
C                LW     - INPUT LENGTH OF W.
C                K      - INPUT NUMBER OF AUTOCOVARIANCES AND
C                           AUTOCORRELATIONS TO BE COMPUTED.
C                L      - INPUT NUMBER OF PARTIAL AUTOCORRELATIONS TO
C                           BE COMPUTED.  L MUST BE LESS THAN OR EQUAL
C                           TO K.
C                ISW    - INPUT CONTROL PARAMETER FOR DETERMINING TASK
C                           TO BE PERFORMED. IF
C                             ISW = 1  FIND MEAN AND VARIANCE.
C                             ISW = 2  FIND AUTOCOVARIANCE.
C                             ISW = 3  FIND MEAN, VARIANCE, AND
C                                        AUTOCOVARIANCES.
C                             ISW = 4  FIND AUTOCOVARIANCES AND
C                                        AUTOCORRELATIONS.
C                             ISW = 5  FIND MEAN, VARIANCE, AUTO-
C                                        COVARIANCES, AND AUTOCORRELAT-
C                                        IONS.
C                             ISW = 6  FIND AUTOCOVARIANCES, AUTOCORREL-
C                                        ATIONS, AND PARTIAL AUTOCO-
C                                        RRELATIONS.
C                             ISW = 7  FIND MEAN, VARIANCE, AUTOCOVAR-
C                                        IANCES, AUTOCORRELATIONS, AND
C                                        PARTIAL AUTOCORRELATIONS.
C                AMEAN  - OUTPUT FOR ISW = 1,3,5, AND 7.
C                           INPUT FOR ISW = 2,4, AND 6. MEAN VALUE OF
C                           THE TIME SERIES W.
C                VAR    - OUTPUT FOR ISW = 1,3,5, AND 7.
C                           INPUT FOR ISW = 2,4, AND 6. VARIANCE OF
C                           TIME SERIES W.
C                ACV    - VECTOR OF LENGTH K.
C                           OUTPUT FOR ISW = 2,3,4,5,6 AND 7.
C                           AUTOCOVARIANCES FOR TIME SERIES W.  ACV(I)
C                           CORRESPONDS TO A TIME LAG OF I TIME UNITS.
C                AC     - VECTOR OF LENGTH K.
C                           OUTPUT FOR ISW = 4,5,6, AND 7.
C                           AUTOCORRELATIONS FOR TIME SERIES W.  AC(I)
C                           CORRESPONDS TO A TIME LAG OF I TIME UNITS.
C                PACV   - VECTOR OF LENGTH L.
C                           OUTPUT FOR ISW = 6 AND 7.
C                           PARTIAL AUTOCORRELATIONS OF TIME SERIES W.
C                           PACV(1) = AC(1).
C                WKAREA - WORK AREA VECTOR OF LENGTH L.
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IN SOME ROUTINES OF THIS CHAPTER ACV(I) CORRESPONDS
C                TO A TIME LAG OF (I-1) TIME UNITS RATHER THAN I TIME
C                UNITS. THUS, IN THESE SUBROUTINES, ACV(1) IS THE SAME
C                AS THE VARIANCE VAR. IN THE CALLING PROGRAM TO FTAUTO,
C                IF THE USER WISHES THE VARIANCE TO BE THE FIRST ENTRY
C                IN HIS AUTOCOVARIANCE ARRAY THE FOLLOWING CALL CAN BE
C                MADE
C                  CALL FTAUTO(W,LW,K,L,ISW,AMEAN,ACV(1),ACV(2),AC,PACV,
C                 *            WKAREA)
C                THE USER SHOULD NOTE THAT IN THIS CASE, ACV MUST BE
C                DIMENSIONED K+1 IN THE MAIN PROGRAM.
C            2.  IF THE TIME SERIES W IS CONSTANT, THEN ANY OF ACV, AC,
C                AND PACV THAT ARE OUTPUT, ACCORDING TO THE ISW SETTING,
C                ARE SET TO ZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTAUTO  (W,LW,K,L,ISW,AMEAN,VAR,ACV,AC,PACV,WKAREA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LW,K,L,ISW
      REAL               W(LW),ACV(1),AC(1),PACV(1),WKAREA(1),AMEAN,VAR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,KEND,J1,K0,J2,J1MK,IFLAG,IM
      REAL               TEMP2,ZERO
      DOUBLE PRECISION   TEMP,TEMP1,ONE,DZERO
      DATA               ZERO/0.0/,ONE/1.0D0/,DZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (ISW.EQ.1) GO TO 40
      IFLAG = 0
      IF (K.LE.0) GO TO 15
      DO 5 I=1,K
    5 ACV(I) = ZERO
      IF (ISW.LT.4) GO TO 25
      DO 10 I=1,K
   10 AC(I) = ZERO
   15 IF (L.LE.0) GO TO 25
      DO 20 I=1,L
   20 PACV(I) = ZERO
   25 IF ((ISW/2)*2.NE.ISW) VAR = ZERO
      DO 30 I=2,LW
         IF (W(I).NE.W(I-1)) GO TO 35
   30 CONTINUE
      IFLAG = 1
      AMEAN = W(1)
   35 IF (IFLAG.EQ.1) GO TO 9005
   40 IM = (ISW/2)*2-ISW
C                                  COMPUTE THE MEAN
      IF (IM.EQ.0) GO TO 55
      TEMP = DZERO
      DO 45 I=1,LW
         TEMP = TEMP+DBLE(W(I))
   45 CONTINUE
      AMEAN = TEMP/LW
C                                  COMPUTE THE VARIANCE
      TEMP = DZERO
      DO 50 I=1,LW
         TEMP = TEMP+(DBLE(W(I)-AMEAN)*DBLE(W(I)-AMEAN))
   50 CONTINUE
      VAR = TEMP/LW
      IF (ISW.EQ.1) GO TO 9005
C                                  COMPUTE AUTOCOVARIANCES
   55 DO 65 J=1,K
         KEND = LW-J
         IF (KEND .LE. 0) GO TO 65
         TEMP = DZERO
         DO 60 I=1,KEND
            TEMP = TEMP+(DBLE(W(I)-AMEAN)*DBLE(W(I+J)-AMEAN))
   60    CONTINUE
         ACV(J) = TEMP/LW
   65 CONTINUE
      IF (ISW.LT.4) GO TO 9005
C                                  COMPUTE AUTOCORRELATIONS
      DO 70 J=1,K
         AC(J) = ACV(J)/VAR
   70 CONTINUE
      IF (ISW.LT.6) GO TO 9005
C                                  COMPUTE PARTIAL AUTOCOVARIANCE
      PACV(1) = AC(1)
      DO 90 J=2,L
         J1 = J-1
         WKAREA(J1) = PACV(J1)
         J2 = (J1)/2
         IF (J.EQ.2) GO TO 80
         DO 75 K0=1,J2
            J1MK = J1-K0
            TEMP2 = WKAREA(K0)-PACV(J1)*WKAREA(J1MK)
            WKAREA(J1MK) = WKAREA(J1MK)-PACV(J1)*WKAREA(K0)
            WKAREA(K0) = TEMP2
   75    CONTINUE
   80    CONTINUE
         TEMP = DZERO
         TEMP1 = DZERO
         DO 85 I=1,J1
            TEMP = TEMP+(DBLE(AC(J-I))*DBLE(WKAREA(I)))
            TEMP1 = TEMP1+(DBLE(AC(I))*DBLE(WKAREA(I)))
   85    CONTINUE
         PACV(J) = (AC(J)-TEMP)/(ONE-TEMP1)
   90 CONTINUE
 9005 RETURN
      END

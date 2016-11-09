C   IMSL ROUTINE NAME   - FTWENM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - MULTICHANNEL WIENER FORECAST
C
C   USAGE               - CALL FTWENM (AC,NI,M,LF,IA,IB,CC,IC,ID,EPS,TR,
C                           IF,IG,F,ERP,LP,WK,IER)
C
C   ARGUMENTS    AC     - INPUT NI BY NI BY LF ARRAY CONTAINING
C                           MULTICHANNEL AUTOCOVARIANCES IN MULTIPLEXED
C                           MODE.
C                NI     - INPUT NUMBER OF INPUT CHANNELS.
C                           NI MUST BE GREATER THAN ZERO.
C                M      - INPUT NUMBER OF OUTPUT CHANNELS.
C                           M MUST BE GREATER THAN ZERO AND LESS
C                           THAN OR EQUAL TO NI.
C                LF     - INPUT MAXIMUM ALLOWABLE LENGTH OF THE
C                           FILTER.  LF MUST BE GREATER THAN ZERO.
C                IA     - INPUT FIRST DIMENSION OF AC EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                IB     - INPUT SECOND DIMENSION OF AC EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                CC     - INPUT M BY NI BY LF ARRAY CONTAINING
C                           MULTICHANNEL CROSS-COVARIANCES IN
C                           MULTIPLEXED MODE.
C                IC     - INPUT FIRST DIMENSION OF CC EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                ID     - INPUT SECOND DIMENSION OF CC EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                EPS    - INPUT CONSTANT IN THE INCLUSIVE INTERVAL (0,1)
C                           FOR CONTROLLING THE FILTER LENGTH.  SEE
C                           REMARKS.
C                TR     - INPUT TRACE OF THE ZERO TIME LAG
C                           AUTOCOVARIANCE MATRIX OF THE DESIRED OUTPUT.
C                IF     - INPUT FIRST DIMENSION OF F EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                IG     - INPUT SECOND DIMENSION OF F EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                F      - OUTPUT M BY NI BY LF ARRAY CONTAINING
C                           THE MULTICHANNEL FILTER COEFFICIENTS IN
C                           MULTIPLEXED MODE.
C                ERP    - OUTPUT VECTOR OF LENGTH LF+M.
C                           THE FIRST LP LOCATIONS CONTAIN THE
C                           NORMALIZED MEAN SQUARE ERROR FOR FILTERS OF
C                           LENGTH 1,2,...,LP.
C                           THE LAST M LOCATIONS ARE USED AS
C                           WORK STORAGE.
C                LP     - OUTPUT CONSTANT GIVING THE LENGTH OF THE
C                           FILTER.
C                WK     - WORK AREA OF LENGTH NI*NI*(2*LF+12).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THE FILTER LENGTH WAS
C                             SPECIFIED LESS THAN ONE.
C                           IER=130 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL ROUTINE LINV3F.  AN INPUT MATRIX TO
C                             THE ROUTINE APPEARS TO BE SINGULAR.
C                           IER=131 INDICATES THE NUMBER OF OUTPUT
C                             CHANNELS IS LESS THAN ONE OR THAT THE
C                             NUMBER OF OUTPUT CHANNELS EXCEEDS THE
C                             NUMBER OF INPUT CHANNELS.
C
C   REQD. IMSL ROUTINES - SINGLE/LINV3F,LUDATN,LUELMN,UERTST,UGETIO,
C                           VMULFF
C                       - DOUBLE/LINV3F,LUDATN,LUELMN,UERTST,UGETIO,
C                           VMULFF,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE LENGTH OF THE FILTER IS CONTROLLED BY THE TWO
C                PARAMETERS LF AND EPS. RECURSION TO A LONGER FILTER
C                STOPS WHENEVER EITHER OF THE FOLLOWING TWO CONDITIONS
C                OCCUR;
C                  A. THE NORMALIZED MEAN SQUARE ERROR FALLS BELOW EPS
C                  B. THE FILTER REACHES THE MAXIMUM ALLOWABLE LENGTH,
C                     LF
C            2.  IMSL SUBROUTINE FTCROS MAY BE USED TO OBTAIN THE
C                INPUT PARAMETERS AC, CC, AND TR.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTWENM (AC,NI,M,LF,IA,IB,CC,IC,ID,EPS,TR,IF,IG,F,ERP,
     1                   LP,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NI,M,LF,IA,IB,IC,ID,IF,IG,LP,IER
      REAL               AC(IA,IB,1),CC(IC,ID,1),F(IF,IG,1),
     1                   WK(NI,NI,1),ERP(1),EPS,TR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IE,IH,IHP,II,INO,IP,IQ,IR,IRP,IS,IT,J,JJ,K,
     1                   KK,KT,L,LF1,LT
      REAL               D1,D2,ZERO,ONE
      DOUBLE PRECISION   TEMP,GM
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF(LF .GT. 0) GO TO 5
      IER = 129
      GO TO 9000
    5 IF(M .GT. 0 .AND. NI .GE. M) GO TO 10
      IER = 131
      GO TO 9000
   10 LT = LF+11
C                                  INITIALIZE WORK AREAS
      DO 20 I=1,NI
         DO 15 J=1,NI
            WK(I,J,11) = ZERO
            WK(I,J,LT) = ZERO
            WK(I,J,8) = AC(I,J,1)
   15    WK(I,J,7) = AC(I,J,1)
         WK(I,I,11) = ONE
         WK(I,I,LT) = ONE
   20 CONTINUE
      IH = 2*LF+11
      IHP = IH+1
      D1 = ZERO
      D2 = ZERO
C                                  INVERT AC(I,J,1) INTO WK(I,J,1)
      DO 25 I=1,NI
         DO 25 J=1,NI
   25 WK(J,I,1) = AC(J,I,1)
      CALL LINV3F(WK,WK(1,1,2),1,NI,NI,D1,D2,WK(1,1,3),IER)
      IF(IER .GT. 127) GO TO 9000
      DO 30 J=1,NI
         DO 30 I=1,NI
   30 WK(I,J,2) = WK(I,J,1)
      CALL VMULFF(CC,WK,M,NI,NI,IC,NI,F,IF,IER)
      GM = TR
      DO 40 II=1,M
         TEMP = 0.0D0
         DO 35 JJ=1,NI
         TEMP = TEMP+DBLE(CC(II,JJ,1))*DBLE(F(II,JJ,1))
   35 CONTINUE
         ERP(II+LF) = TEMP
   40 CONTINUE
      TEMP = 0.0D0
      DO 45 II=1,M
         TEMP = TEMP+DBLE(ERP(II+LF))
   45 CONTINUE
      ERP(1) = (GM-TEMP)/GM
      LP = 1
      IF(ERP(1) .LE. EPS .OR. LF .LE. 1) GO TO 9005
      LF1 = LF-1
C                                  COMPUTE THE FILTER COEFFICIENTS.
      DO 180 I=1,LF1
         DO 50 JJ=1,NI
            DO 50 II=1,NI
   50    WK(II,JJ,3) = ZERO
         DO 60 K=1,I
            CALL VMULFF(WK(1,1,K+10),AC(1,1,I-K+2),NI,NI,NI,NI,IA,
     *                  WK(1,1,10),NI,IER)
            DO 55 JJ=1,NI
               DO 55 II=1,NI
   55       WK(II,JJ,3) = WK(II,JJ,3)+WK(II,JJ,10)
   60    CONTINUE
         DO 65 JJ=1,NI
            DO 65 II=1,NI
   65    WK(JJ,II,4) = WK(II,JJ,3)
         CALL VMULFF(WK(1,1,3),WK,NI,NI,NI,NI,NI,WK(1,1,5),NI,IER)
         CALL VMULFF(WK(1,1,4),WK(1,1,2),NI,NI,NI,NI,NI,WK(1,1,6),NI,
     *               IER)
         IQ = I+11
         IR = I+LT
         DO 75 JJ=1,NI
            DO 70 II=1,NI
               WK(II,JJ,IR) = -WK(II,JJ,6)
   70       WK(II,JJ,IQ) = -WK(II,JJ,5)
   75    CONTINUE
         IF(I .EQ. 1) GO TO 95
         IE = I-1
         DO 90 K=1,IE
            IS = IR-K
            CALL VMULFF(WK(1,1,5),WK(1,1,IS),NI,NI,NI,NI,NI,WK(1,1,9),
     *                  NI,IER)
            IT = K+11
            CALL VMULFF(WK(1,1,6),WK(1,1,IT),NI,NI,NI,NI,NI,WK(1,1,10),
     *                  NI,IER)
            DO 85 JJ=1,NI
               DO 80 II=1,NI
                  WK(II,JJ,IT) = WK(II,JJ,IT) - WK(II,JJ,9)
   80          WK(II,JJ,IS) = WK(II,JJ,IS) - WK(II,JJ,10)
   85       CONTINUE
   90    CONTINUE
   95    CALL VMULFF(WK(1,1,5),WK(1,1,4),NI,NI,NI,NI,NI,WK(1,1,9),
     *               NI,IER)
         CALL VMULFF(WK(1,1,6),WK(1,1,3),NI,NI,NI,NI,NI,WK(1,1,10),
     *               NI,IER)
         DO 105 JJ=1,NI
            DO 100 II=1,NI
               WK(II,JJ,8) = WK(II,JJ,8) - WK(II,JJ,9)
  100       WK(II,JJ,7) = WK(II,JJ,7) - WK(II,JJ,10)
  105    CONTINUE
         INO = 2
         IF(I .EQ. LF1) INO = 1
         DO 115 K=1,INO
            DO 110 JJ=1,NI
               DO 110 II=1,NI
  110       WK(II,JJ,K) = WK(II,JJ,K+6)
            CALL LINV3F(WK(1,1,K),WK(1,1,5),1,NI,NI,D1,D2,WK(1,1,6),IER)
            IF(IER .GT. 127) GO TO 9000
  115    CONTINUE
         IP = I+1
         KT = I+2
         DO 120 JJ=1,NI
            DO 120 II=1,M
  120    WK(II,JJ,IH) = ZERO
         DO 130 K=1,I
            CALL VMULFF(F(1,1,K),AC(1,1,KT-K),M,NI,NI,IF,IA,WK(1,1,10),
     *                  NI,IER)
            DO 125 JJ=1,NI
               DO 125 II=1,M
  125       WK(II,JJ,IH) = WK(II,JJ,IH)+WK(II,JJ,10)
  130    CONTINUE
         DO 135 JJ=1,NI
            DO 135 II=1,M
  135    WK(II,JJ,IH) = WK(II,JJ,IH)-CC(II,JJ,IP)
         CALL VMULFF(WK(1,1,IH),WK,M,NI,NI,NI,NI,WK(1,1,IHP),NI,IER)
         DO 140 JJ=1,NI
            DO 140 II=1,M
  140    F(II,JJ,IP) = -WK(II,JJ,IHP)
         IRP = IR+1
         DO 155 L=1,I
            CALL VMULFF(WK(1,1,IHP),WK(1,1,IRP-L),M,NI,NI,NI,NI,
     *                  WK(1,1,10),NI,IER)
            DO 150 JJ=1,NI
               DO 145 II=1,M
  145          F(II,JJ,L) = F(II,JJ,L)-WK(II,JJ,10)
  150       CONTINUE
  155    CONTINUE
C                                  COMPUTE THE MEAN SQUARE ERROR
         DO 170 II=1,M
            TEMP = 0.0D0
            DO 165 KK=1,IP
               DO 160 JJ=1,NI
  160          TEMP = TEMP+DBLE(CC(II,JJ,KK))*DBLE(F(II,JJ,KK))
  165       CONTINUE
            ERP(II+LF) = TEMP
  170    CONTINUE
         TEMP = 0.0D0
         DO 175 II=1,M
  175    TEMP = TEMP+DBLE(ERP(II+LF))
         ERP(I+1) = (GM-TEMP)/GM
         LP = LP+1
         IF(ERP(I+1) .LE. EPS) GO TO 9005
  180 CONTINUE
      IF(IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'FTWENM')
 9005 CONTINUE
      RETURN
      END

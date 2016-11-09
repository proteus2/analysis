C   IMSL ROUTINE NAME   - AGXPM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EXPECTED MEAN SQUARES FOR BALANCED COMPLETE
C                           DESIGN MODELS
C
C   USAGE               - CALL AGXPM (IOPT,NF,M,IA,NL,INL,CMS,IORD,IEMS,
C                           STAT,IS,ERTM,IE,IER)
C
C   ARGUMENTS    IOPT   - INPUT OPTION SWITCH.
C                         IF IOPT IS EQUAL TO ZERO, THE ROUTINE COMPUTES
C                           FOR EACH MODEL TERM, THE EXPECTED MEAN
C                           SQUARE, DEGREES OF FREEDOM, TEST TERM, AND
C                           TEST TERM DEGREES OF FREEDOM.
C                         IF IOPT IS NOT EQUAL TO ZERO, THE ROUTINE
C                           COMPUTES FOR EACH MODEL TERM, THE EXPECTED
C                           MEAN SQUARE, DEGREES OF FREEDOM, TEST TERM,
C                           TEST TERM DEGREES OF FREEDOM, VARIANCE
C                           COMPONENT ESTIMATE, AND F-VALUE.
C                NF     - INPUT NUMBER OF FACTORS (NUMBER OF DISTINCT
C                           SUBSCRIPTS IN THE MODEL). NF MUST BE
C                           GREATER THAN OR EQUAL TO TWO.
C                M      - INPUT NUMBER OF TERMS SPECIFIED IN THE STATIS-
C                           TICAL MODEL CONTAINED IN VECTOR IA.
C                IA     - INPUT/OUTPUT VECTOR.
C                         ON INPUT, IA CONTAINS THE STATISTICAL MODEL.
C                           IA IS OF LENGTH EQUAL TO TWICE THE
C                           NUMBER OF LETTERS AND SYMBOLS IN THE
C                           ALGEBRAIC MODEL SPECIFICATION. IA MUST HAVE
C                           ONE CHARACTER LEFT-JUSTIFIED PER WORD.
C                           THE LAST TERM IN THE MODEL IS RESERVED FOR
C                           ERROR. IT IS REPRESENTED BY THE LETTER E AND
C                           MUST HAVE A FULL SET OF SUBSCRIPTS. THE
C                           REMAINDER OF IA IS USED AS WORK STORAGE.
C                           FOR EXAMPLE, THE FIRST HALF (33 LOCATIONS)
C                           OF IA WOULD CONTAIN
C                           S(I)T(IJ)P(K)TP(IJK)SP(IK)E(IJKL)
C                           FOR THE MODEL
C                           Y(IJKL)=M+S(I)+T(IJ)+P(K)+TP(IJK)+SP(IK)+
C                                   E(IJKL).
C                         ON OUTPUT, THE ELEMENTS OF IA WILL BE
C                           REORDERED AS INDICATED BY VECTOR IORD
C                NL     - INPUT/OUTPUT MATRIX OF DIMENSION (NF+1) BY
C                           (2M+5).
C                         ON INPUT, THE FIRST NF LOCATIONS
C                           OF COLUMN 1 MUST CONTAIN THE NUMBER OF
C                           LEVELS FOR EACH OF THE NF FACTORS.
C                           NL(I,1) MUST BE GREATER THAN OR EQUAL TO
C                           TWO FOR I=1,2,...,NF. THE FIRST
C                           NF LOCATIONS OF COLUMN 2 MUST CONTAIN
C                           THE FIXED OR RANDOM FACTOR INDICATOR. IF
C                           NL(I,2) = 0, THEN THE I-TH MAIN EFFECT IS
C                           ASSUMED TO BE FIXED FOR I=1,2,...,NF.
C                           IF NL(I,2) = 1, THEN THE THE I-TH MAIN
C                           EFFECT IS ASSUMED TO BE RANDOM FOR
C                           I=1,2,...,NF. THE REMAINDER OF NL IS
C                           USED AS WORK STORAGE.
C                         ON OUTPUT, COLUMNS 1 AND 2 MAY BE REARRANGED.
C                           REARRANGEMENT IS INDICATED BY VECTOR IORD.
C                INL    - INPUT ROW DIMENSION OF THE MATRIX NL EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           CALLING PROGRAM.
C                CMS    - INPUT/OUTPUT VECTOR OF LENGTH M USED ONLY
C                           WHEN IOPT IS NOT EQUAL TO ZERO.
C                         ON INPUT, CMS CONTAINS THE CALCULATED MEAN
C                           SQUARE FOR EACH TERM.
C                         ON OUTPUT, THE ELEMENTS OF CMS WILL BE
C                           REORDERED AS INDICATED BY VECTOR IORD.
C                IORD   - OUTPUT VECTOR OF LENGTH M, INDICATING THE
C                           REORDERING OF THE TERMS IN IA. THE OUTPUT
C                           PARAMETERS ARE CALCULATED USING THE
C                           REORDERED MODEL. IORD(I) = J IMPLIES THE
C                           J-TH TERM IN THE INPUT MODEL IS THE I-TH
C                           TERM IN THE OUTPUT MODEL.
C                IEMS   - OUTPUT VECTOR OF LENGTH M(M+1)/2 CONTAINING
C                           THE VARIANCE COMPONENT COEFFICIENTS FOR
C                           EACH MODEL TERM. SEE REMARKS.
C                STAT   - OUTPUT MATRIX.
C                         STAT IS OF DIMENSION M BY (M+3) FOR
C                           IOPT EQUAL TO ZERO.
C                         STAT IS OF DIMENSION M BY (M+5) FOR
C                           IOPT NOT EQUAL TO ZERO.
C                         STAT CONTAINS THE DEGREES OF FREEDOM FOR THE
C                           I-TH TERM AND ITS ERROR TERM IN STAT(I,1)
C                           AND STAT(I,2) RESPECTIVELY.
C                           STAT(I,3) INDICATES WHETHER OR NOT THE I-TH
C                           ERROR TERM IS A COMPOSITE. IF STAT(I,3)=1.0,
C                           THEN THE ERROR TERM IS A COMPOSITE OF TWO
C                           OR MORE MEAN SQUARES. IF STAT(I,3)=0.0,
C                           THEN THE ERROR TERM IS NOT A COMPOSITE
C                           OF TWO OR MORE MEAN SQUARES.
C                           FOR IOPT NOT EQUAL TO ZERO,
C                           STAT(I,4) CONTAINS THE VARIANCE COMPONENT
C                           ESTIMATE WHILE STAT(I,5) CONTAINS THE
C                           F-VALUE FOR THE I-TH TERM FOR I=1,2,...,M.
C                           STAT(M,2), STAT(M,3) AND STAT(M,5) ARE
C                           UNDEFINED.
C                           FOR IOPT EQUAL TO ZERO, IF STAT(I,3) EQUALS
C                           1.0, STAT(I,2) IS UNDEFINED.
C                           THE REMAINDER OF STAT IS USED AS WORK
C                           STORAGE.
C                IS     - INPUT ROW DIMENSION OF THE MATRIX STAT
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                ERTM   - OUTPUT MATRIX OF DIMENSION M BY (STAT(1,3)+
C                           STAT(2,3)+...+STAT(M-1,3)) CONTAINING THE
C                           MEAN SQUARE COEFFICIENTS INDICATING THE
C                           COMPOSITION OF THE COMPOSITE ERROR TERMS.
C                           SEE REMARKS.
C                IE     - INPUT ROW DIMENSION OF THE MATRIX ERTM
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NF WAS SPECIFIED LESS
C                             THAN 2.
C                           IER=130 INDICATES THAT SOME NL(I,1), I = 1,
C                             2,...,NF WAS SPECIFIED LESS THAN 2.
C                           IER=131 INDICATES THAT THE STATISTICAL MODEL
C                             CONTAINED IN VECTOR IA WAS SPECIFIED
C                             INCORRECTLY. EITHER THE SAME MAIN EFFECT
C                             (SINGLE EFFECT SYMBOL) APPEARS MORE THAN
C                             ONCE IN THE MODEL OR THE SAME SUBSCRIPT
C                             IS ASSOCIATED WITH MORE THAN ONE MAIN
C                             EFFECT.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO,VSRTR
C                       - DOUBLE/UERTST,UGETIO,VSRTRD,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  PARAMETERS NL, CMS, IEMS, STAT, AND ERTM WILL BE
C                ORDERED IN A MANNER REFLECTING THE ORDERING OF THE
C                TERMS IN IA. THE ORDERING OF TERMS IN IA, WHICH MAY
C                BE CHANGED INSIDE OF THE ROUTINE FOR COMPUTATIONAL
C                CONVENIENCE, IS INDICATED BY VECTOR IORD. THE
C                ORDERING OF THE ELEMENTS IN IEMS, STAT, AND ERTM
C                IS ILLUSTRATED IN THE REFERENCE MANUAL EXAMPLE.
C            2.  THE MEAN SQUARES MUST BE KNOWN TO COMPUTE THE DEGREES
C                OF FREEDOM FOR THE COMPOSITE ERROR TERMS. THEREFORE,
C                FOR IOPT EQUAL TO ZERO, STAT(I,2) IS UNDEFINED WHEN
C                THE ERROR TERM FOR TERM I IS COMPOSITE.
C            3.  CALLING AGXPM WITH VARIABLE IOPT EQUAL TO ZERO IS
C                DESIGNED FOR USE BEFORE THE DATA IS AVAILABLE AND
C                RETURNS THE SAME INFORMATION AS CALLING AGXPM WITH
C                IOPT NOT EQUAL TO ZERO, WITH THE EXCEPTION OF THE
C                VARIANCE COMPONENT ESTIMATES AND F-VALUES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AGXPM  (IOPT,NF,M,IA,NL,INL,CMS,IORD,IEMS,STAT,IS,ERTM,
     1                  IE,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER           NF,M,INL,IS,IE,IER,NL(INL,1),IA(1),IORD(1),
     1                  IEMS(1),IOPT
      REAL              CMS(1),STAT(IS,1),ERTM(IE,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICNT,II,IRP,ISTAT,ISW,I1,J,J1,K,K1,L,LP,MI,
     1                   M1,M25,NF1,NS,NT
      REAL               ZERO,ONE
      DOUBLE PRECISION   SUM,SUM1
      DATA               ZERO/0.0/,ONE/1.0/
      DATA               LP/1H(/,IRP/1H)/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
C                                  TERMINAL ERROR - NF IS LESS THAN 2
      IF (NF .LT. 2) GO TO 9000
C                                  TERMINAL ERROR - SOME NL(I,1) WAS
C                                  SPECIFIED INCORRECTLY
      IER = 130
      M25 = M+M+5
      DO 10 I = 1,NF
         IF (NL(I,1)  .LT. 2) GO TO 9000
         DO 10 J = 6,M25
            NL(I,J) = 0
   10 CONTINUE
      IER = 131
      ISW = 0
C                                  ANALYZE THE MODEL
      ICNT = 0
      I = 1
      K = 1
      NF1 = NF+1
   15 IF (IA(I) .EQ. LP) GO TO 30
      IF (IA(I) .NE. IRP) GO TO 20
C                                  STORE ASSOCIATED SUBSCRIPTS
      IF (ICNT .NE. 0) GO TO 35
      NL(K,4) = IA(I-1)
      K = K+1
      IF (K .GT. NF) GO TO 40
      ISW = 0
      I = I+1
      GO TO 15
   20 IF (ISW .EQ. 0) ICNT = ICNT+1
   25 I = I+1
      GO TO 15
   30 IF (ICNT .NE. 1) GO TO 25
C                                  STORE EFFECT SYMBOL
      NL(K,3) = IA(I-1)
      ISW = 1
   35 ICNT = 0
      GO TO 25
C                                  ASSIGN POWER OF TWO TO EACH
C                                  SUBSCRIPT
   40 NL(1,5) = 1
      I1 = 1
      J = NL(1,4)
      ICNT = 1
      NT = NL(1,1)
      DO 50 I = 2,NF
         NT = NT*NL(I,1)
         NL(I,5) = ICNT+ICNT
         ICNT = NL(I,5)
C                                  TERMINAL ERROR - THE MODEL CONTAINED
C                                  IN VECTOR IA WAS SPECIFIED
C                                  INCORRECTLY
         DO 45 K = I,NF
            IF ((NL(I1,4) .EQ. NL(K,4)) .OR. (NL(I1,3) .EQ. NL(K,3)))
     *      GO TO 9000
   45    CONTINUE
         I1 = I
   50 CONTINUE
      IER = 0
C                                  SUM NUMBERS ASSOCIATED WITH
C                                  SUBSCRIPTS FOR EACH INDIVIDUAL TERM
      M1 = M-1
      K = 1
      I = 0
      NL(NF1,1) = 0
      IF (IOPT .EQ. 1) STAT(K,2) = CMS(K)
      ICNT = 1
   55 I = I+1
      IF (IA(I) .NE. LP) GO TO 55
      I = I+1
      IORD(ICNT) = ICNT
      ICNT = ICNT+1
      STAT(K,1) = ZERO
   60 DO 65 J = 1,NF
         IF (IA(I) .NE. NL(J,4)) GO TO 65
C                                  SUM THE SUBSCRIPTS
         STAT(K,1) = STAT(K,1)+NL(J,5)
         I = I+1
         IF (IA(I) .NE. IRP) GO TO 60
         K = K+1
         IF (K .GT. M) GO TO 70
C                                  STORE INDEX FOR THE BEGINNING OF
C                                  EACH TERM
         NL(NF1,K) = I
         IF (IOPT .EQ. 1) STAT(K,2) = CMS(K)
         GO TO 55
   65 CONTINUE
   70 J = NL(NF1,M)
      ICNT = J+NF+3
      DO 75 I = 1,J
         IA(ICNT+I) = IA(I)
   75 CONTINUE
C                                  SORT THE MODEL SO THAT NUMBERS
C                                  ASSIGNED TO SUBSCRIPTS ARE IN
C                                  ASCENDING ORDER
      CALL VSRTR (STAT(1,1),M1,IORD)
      I = 1
C                                  REARRANGE THE TERMS
      DO 85 J = 1,M1
         I1 = IORD(J)
         IF (IOPT .EQ. 1) CMS(J) = STAT(I1,2)
         K = NL(NF1,I1)+ICNT+1
         K1 = NL(NF1,I1+1)+ICNT
         DO 80 L = K,K1
            IA(I) = IA(L)
            I = I+1
   80    CONTINUE
   85 CONTINUE
C                                  RECORD EFFECT SYMBOLS AND SUBSCRIPTS
C                                  FOR EACH TERM
      I = 1
      MI = 6
  100 DO 105 K = 1,NF
         IF (IA(I) .NE. NL(K,3)) GO TO 105
         NL(K,MI) = 1
         GO TO 110
  105 CONTINUE
  110 I = I+1
      IF (IA(I) .NE. LP) GO TO 100
      MI = MI+1
      I = I+1
  115 DO 120 K = 1,NF
         IF (IA(I) .NE. NL(K,4)) GO TO 120
         NL(K,MI) = 1
         GO TO 125
  120 CONTINUE
  125 I = I+1
      IF (IA(I) .NE. IRP) GO TO 115
      MI = MI+1
      I = I+1
      IF (MI .LT. M25) GO TO 100
C                                  OBTAIN THE EMS STRUCTURE
      ISTAT = 5
      IF (IOPT .NE. 1) ISTAT = 3
      J = ISTAT+1
      DO 130 I = 1,M
         STAT(I,J) = ONE
         STAT(M,J) = ONE
         J = J+1
  130 CONTINUE
      I1 = 5
      J1 = 1
  135 J = I1+1
      MI = J1
C                                  LOGICAL .OR. OF SYMBOL AND RANDOM
C                                  FACTORS
      DO 140 K = 1,NF
         NL(K,5) = NL(K,J)
         IF (NL(K,2) .EQ. 1) NL(K,5) = 1
  140 CONTINUE
      J1 = J1+1
      I1 = I1+2
      II = I1
      DO 165 L = J1,M1
         II = II+2
C                                  COMPARE SUBSCRIPTS
         DO 145 ICNT = 1,NF
            IF (NL(ICNT,I1) .EQ. 0) GO TO 145
            IF (NL(ICNT,II) .EQ. 0) GO TO 155
  145    CONTINUE
         DO 150 ICNT = 1,NF
            IF (NL(ICNT,II-1) .EQ. 0) GO TO 150
            IF (NL(ICNT,5) .EQ. 0) GO TO 155
  150    CONTINUE
         STAT(L,MI+ISTAT) = ONE
         GO TO 160
  155    STAT(L,MI+ISTAT) = ZERO
  160    STAT(MI,L+ISTAT) = ZERO
  165 CONTINUE
      IF (J1 .LT. M1) GO TO 135
      STAT(M1,M+ISTAT) = -ONE
C                                  FIND MEAN SQUARE COEFFICIENTS FOR
C                                  INDICATING STRUCTURE OF VARIANCE
C                                  COMPONENT ESTIMATES
      K = M
      MI = M+ISTAT
      DO 185 I = 2,M1
         II = K-I
         ICNT = II+ISTAT
         I1 = II+1
         DO 170 J = I1,M
            STAT(II,J+ISTAT)= -STAT(J,ICNT)
  170    CONTINUE
         DO 180 J = I1,M1
            IF (STAT(J,ICNT) .EQ. ZERO) GO TO 180
            J1 = J+1+ISTAT
            DO 175 L = J1,MI
               STAT(II,L) = STAT(II,L)-STAT(J,L)
  175       CONTINUE
  180    CONTINUE
  185 CONTINUE
      K = 1
      NS = 0
      DO 205 I = 7,M25,2
         I1 = 1
         ICNT = 1
         DO 200 J = 1,NF
            IF (NL(J,I) .EQ. 0) GO TO 195
            IF (NL(J,I-1) .EQ. 1) GO TO 190
            I1 = I1*NL(J,1)
            GO TO 200
  190       I1 = I1*(NL(J,1)-1)
            GO TO 200
  195       ICNT = ICNT*NL(J,1)
  200    CONTINUE
         STAT(K,3) = ICNT
         STAT(K,1) = I1
         IF (K .NE. M) NS = NS+I1
         K = K+1
  205 CONTINUE
      STAT(M,1) = NT-NS-1
      IF (IOPT .EQ. 0) GO TO 220
C                                  FIND VARIANCE COMPONENT ESTIMATE
      DO 215 I = 1,M
         SUM = 0.D0
         DO 210 J = I,M
            SUM = SUM+DBLE(CMS(J))*DBLE(STAT(I,J+ISTAT))
  210    CONTINUE
         STAT(I,4) = SUM/STAT(I,3)
  215 CONTINUE
C                                  FIND MEAN SQUARE COEFFICIENTS FOR
C                                  DETERMINING THE STRUCTURE OF THE
C                                  COMPOSITE ERROR TERMS
  220 STAT(2,ISTAT+2) = STAT(2,ISTAT+2)*STAT(2,ISTAT+1)
      J = 3
  225 I1 = J-1
      J1 = J+ISTAT
      DO 235 I = 1,I1
         SUM = 0.D0
         II = I+1
         ICNT = I+ISTAT
         DO 230 K = II,J
            SUM = SUM+STAT(K,ICNT)*STAT(K,J1)
  230    CONTINUE
         STAT(1,ICNT) = SUM
  235 CONTINUE
      K = 2
      DO 240 I = 1,I1
         STAT(K,J1) = STAT(1,I+ISTAT)
         K = K+1
  240 CONTINUE
      J = J+1
      IF (J .LE. M) GO TO 225
      DO 245 J = 1,M1
         J1 = J+1
         I1 = J+ISTAT
         DO 243 I = J1,M
            STAT(I,I1) = STAT(I,I1)*STAT(I,3)
  243    CONTINUE
  245 CONTINUE
C                                  OBTAIN THE VARIANCE COMPONENT
C                                  COEFFICIENTS (EMS STRUCTURE) FOR
C                                  EACH MODEL TERM
      K = 1
      ICNT = M+1
      DO 255 I = 1,M1
         II = M-I
         DO 250 J = 1,II
            J1 = ICNT-J
            IEMS(K) = STAT(J1,I+ISTAT)
            K = K+1
  250    CONTINUE
         IEMS(K) = STAT(I,3)
         K = K+1
  255 CONTINUE
      IEMS(K) = 1
      K = 1
      DO 280 I = 2,M
         I1 = I-1
         ICNT = 0
         SUM1 = 0.0D0
         SUM = 0.D0
         II = I+ISTAT
         DO 265 J = I,M
            IF (STAT(I,II) .EQ. ZERO) GO TO 260
            ICNT = ICNT+1
C                                  FIND ERROR TERM DEGREES OF FREEDOM
C                                  FOR NON-COMPOSITE TERMS
            STAT(I1,2) = STAT(J,1)
            IF (IOPT .NE. 1) GO TO 260
            SUM1 = SUM1+DBLE(CMS(J))*DBLE(CMS(J))/DBLE(STAT(J,1))
            SUM = SUM+DBLE(STAT(I,II))*DBLE(CMS(J))
  260       II = II+1
  265    CONTINUE
         STAT(I1,3) = ZERO
C                                  COMPUTE F-VALUE FOR EACH TERM
         IF (IOPT .EQ. 1) STAT(I1,5) = CMS(I1)/SUM
         IF (ICNT .EQ. 1) GO TO 280
C                                  COMPOSITE TERM
         STAT(I1,3) = ONE
C                                  FIND MEAN SQUARE COEFFICIENTS FOR
C                                  THE COMPOSITE ERROR TERMS
         DO 270 J = 1,I1
            ERTM(J,K) = ZERO
  270    CONTINUE
         DO 275 J = I,M
            ERTM(J,K) = STAT(I,J+ISTAT)
  275    CONTINUE
C                                  FIND ERROR TERM DEGREES OF FREEDOM
C                                  FOR COMPOSITE TERMS
         STAT(I1,2) = ZERO
         IF (IOPT .EQ. 1) STAT(I1,2) = SUM*SUM/SUM1
         K = K+1
  280 CONTINUE
      IF (IOPT .EQ. 1) STAT(M,5) = ZERO
      STAT(M,3) = ZERO
      STAT(M,2) = ZERO
      DO 285 I = 1,NF
         NL(I,6) = NL(I,1)
         NL(I,7) = NL(I,2)
  285 CONTINUE
      K = 1
      DO 290 I = 1,M
         IF (IORD(I) .GT. NF) GO TO 290
         NL(NF1,K) = IORD(I)
         K = K+1
  290 CONTINUE
      DO 295 I = 1,NF
         J = NL(NF1,I)
         NL(I,1) = NL(J,6)
         NL(I,2) = NL(J,7)
  295 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'AGXPM ')
 9005 RETURN
      END

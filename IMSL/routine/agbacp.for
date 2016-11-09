C   IMSL ROUTINE NAME   - AGBACP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ANALYSIS OF BALANCED COMPLETE EXPERIMENTAL
C                           DESIGN STRUCTURE DATA
C
C   USAGE               - CALL AGBACP (NF,NL,IA,Y,IWK,LST,LOC,SS,NDF,
C                           IER)
C
C   ARGUMENTS    NF     - INPUT NUMBER OF FACTORS (NUMBER OF DISTINCT
C                           SUBSCRIPTS IN THE MODEL).
C                NL     - INPUT VECTOR OF LENGTH NF CONTAINING THE
C                           NUMBER OF LEVELS FOR EACH OF THE NF FACTORS.
C                IA     - INPUT VECTOR CONTAINING THE STATISTICAL MODEL.
C                           IA IS OF LENGTH EQUAL TO THE NUMBER
C                           OF LETTERS AND SYMBOLS IN THE ALGEBRAIC
C                           MODEL SPECIFICATION. IA MUST HAVE ONE
C                           CHARACTER LEFT-JUSTIFIED PER WORD.
C                           THE LAST TERM IN THE MODEL IS RESERVED FOR
C                           ERROR. IT IS REPRESENTED BY THE LETTER E AND
C                           MUST HAVE A FULL SET OF SUBSCRIPTS.
C                Y      - VECTOR OF LENGTH 2**NF+(NL(1)+1)*(NL(2)+1)*
C                           ...*(NL(NF)+1). THE FIRST 2*NF LOCATIONS
C                           OF Y ARE WORK STORAGE.
C                         ON INPUT, Y CONTAINS THE RESPONSES IN
C                           LOCATIONS 2**NF+1, 2**NF+2,...,
C                           2**NF+NL(1)*NL(2)*...*NL(NF).
C                           THE INPUT RESPONSES ARE DESTROYED ON OUTPUT.
C                           THE REMAINING COMPONENTS OF Y ARE UNDEFINED.
C                         ON OUTPUT, Y CONTAINS THE FULL SET OF MEANS
C                           (I.E. AS IF THE MODEL WERE A COMPLETE
C                           FACTORIAL). THIS FULL SET OF MEANS
C                           IS CONTAINED IN
C                           LOCATIONS 2**NF+NL(1)*NL(2)*...*NL(NF)+1,
C                           2**NF+NL(1)*NL(2)*...*NL(NF)+2, ... ,
C                           2**NF+(NL(1)+1)*(NL(2)+1)*...*(NL(NF)+1).
C                           THE LOCATION OF THE FIRST MEAN IN Y FOR
C                           EACH MODEL TERM SPECIFIED IS GIVEN IN
C                           VECTOR LOC.
C                IWK    - WORK AREA OF LENGTH 9*NF.
C                LST    - WORK AREA OF LENGTH 2**NF.
C                LOC    - OUTPUT VECTOR OF LENGTH M. M IS EQUAL TO THE
C                           NUMBER OF MODEL TERMS SPECIFIED IN VECTOR
C                           IA. LOC CONTAINS THE LOCATION OF THE FIRST
C                           MEAN IN THE GROUP OF MEANS FOR EACH TERM IN
C                           THE MODEL, INCLUDING THE GRAND MEAN.
C                SS     - OUTPUT VECTOR OF LENGTH M+1 (M IS DESCRIBED
C                           ABOVE). SS CONTAINS THE SUM OF SQUARES
C                           FOR EACH MODEL TERM. THE CORRECTED TOTAL
C                           SUM OF SQUARES IS IN SS(M+1).
C                NDF    - OUTPUT VECTOR OF LENGTH M+1 (M IS DECSRIBED
C                           ABOVE) CONTAINING THE DEGREES OF FREEDOM
C                           ASSOCIATED WITH EACH COMPONENT OF SS.
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NF WAS SPECIFIED
C                             LESS THAN 2
C                           IER=130 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL SUBROUTINE AFACN
C                           IER=131 INDICATES THAT THE MODEL CONTAINED
C                             IN VECTOR IA WAS SPECIFIED INCORRECTLY.
C                             EITHER THE SAME MAIN EFFECT (SINGLE
C                             EFFECT SYMBOL) APPEARS MORE THAN ONCE
C                             IN THE MODEL OR THE SAME SUBSCRIPT IS
C                             ASSOCIATED WITH MORE THAN ONE MAIN
C                             EFFECT.
C
C   REQD. IMSL ROUTINES - SINGLE/AFACN,UERTST,UGETIO
C                       - DOUBLE/AFACN,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE ORDER OF THE ELEMENTS OF LOC, SS, AND NDF
C                CORRESPONDS TO THE ORDER OF THE TERMS IN IA.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AGBACP (NF,NL,IA,Y,IWK,LST,LOC,SS,NDF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NF,NL(NF),IA(1),IWK(NF,9),LST(1),LOC(1),NDF(1),
     1                   IER
      REAL               Y(1),SS(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LP,IRP,IE,K,I,NF1,J,II,NF2,I1,I2,KX,NA,ICNT,M,
     1                   ISUM,NA2,KK,JJ,ILM,KKK,ISW,NDFT,NF21,IOPT,NALL,
     2                   NT,L,MM,K1
      REAL               X,XX
      DOUBLE PRECISION   SUM,SUM1,XY
      DATA               LP/'('/,IRP/')'/,IE/'E'/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
C                                  TERMINAL ERROR - NF WAS SPECIFIED
C                                    LESS THAN 2
      IF (NF .LT. 2) GO TO 9000
      ISW = 0
      ICNT = 0
      I = 1
      K = 1
      NDFT = 1
      NF1 = NF-1
      NF21 = 2**NF1
      NF2 = NF21+NF21
    5 IF (IA(I) .EQ. LP) GO TO 20
      IF (IA(I) .NE. IRP) GO TO 10
C                                  STORE ASSOCIATED SUBSCRIPTS
      IF (ICNT .NE. 0) GO TO 25
      IWK(K,2) = IA(I-1)
      NDFT = NDFT*NL(K)
      K = K+1
      IF (K .GT. NF) GO TO 30
      ISW = 0
      I = I+1
      GO TO 5
   10 IF (ISW .EQ. 0) ICNT = ICNT+1
   15 I = I+1
      GO TO 5
   20 IF (ICNT .NE. 1) GO TO 15
C                                  STORE THE EFFECT SYMBOL
      IWK(K,1) = IA(I-1)
      ISW = 1
   25 ICNT = 0
      GO TO 15
C                                  ASSIGN POWERS OF TWO TO EACH
C                                  SUBSCRIPT
   30 IWK(NF,3) = 1
      DO 40 I = 1,NF1
C                                  TERMINAL ERROR - THE MODEL WAS
C                                  SPECIFIED INCORRECTLY
         II = I+1
         DO 35 J = II,NF
            IF ((IWK(J,2) .NE. IWK(I,2)) .OR. (IWK(J,1) .NE. IWK(I,1)))
     *      GO TO 35
            IER = 131
            GO TO 9000
   35    CONTINUE
         K = NF-I
         ICNT = IWK(K+1,3)
         IWK(K,3) = ICNT+ICNT
   40 CONTINUE
C                                  COMPUTE COMPLETE SET OF MEANS
      IOPT = 0
      CALL AFACN (IOPT,NF,NL,Y,SS,NDF,IER)
      IF (IER .NE. 0) GO TO 9000
      X = Y(1)
      Y(1) = NF2+1
      DO 45 I = 2,NF2
         XX = Y(I)
         Y(I) = Y(I-1)+X
         X = XX
   45 CONTINUE
C                                  COMPUTE TOTAL SUM OF SQUARES
      SUM1 = 0.0D0
      I1 = Y(1)
      I2 = Y(2)-1
      J = Y(NF2)
      DO 50 II = I1,I2
         Y(II) = Y(II)-Y(J)
         SUM1 = SUM1+DBLE(Y(II))*DBLE(Y(II))
   50 CONTINUE
      NDFT = NDFT-1
      NALL = 0
      XY = 0.D0
      KX = 1
      I = 1
      ISW = 0
   55 ILM = 1
      NDF(KX) = 1
      DO 60 K = 1,NF
         IWK(K,5) = 0
         IWK(K,4) = 0
         IWK(K,6) = 0
   60 CONTINUE
      NA = 1
      ISUM = 0
   65 IF (ISW .NE. 0) GO TO 75
      IF (IA(I) .EQ. IE) GO TO 220
      DO 70 K = 1,NF
         IF (IWK(K,1) .NE. IA(I)) GO TO 70
         IWK(K,4) = 1
         GO TO 95
   70 CONTINUE
      GO TO 95
   75 DO 85 K = 1,NF
         IF (IWK(K,2) .NE. IA(I)) GO TO 85
         IWK(K,5) = 1
         ILM = ILM*NL(K)
         IF (IWK(K,4) .NE. 1) GO TO 80
         IWK(NA,6) = IWK(K,3)
         NDF(KX) = NDF(KX)*(NL(K)-1)
         NA = NA+1
         GO TO 90
   80    ISUM = ISUM+IWK(K,3)
         NDF(KX) = NDF(KX)*NL(K)
         GO TO 90
   85 CONTINUE
C                                  OBTAIN STATISTICS ON SUBSCRIPTS
   90 NT = NT+1
   95 I = I+1
      IF (IA(I) .NE. LP) GO TO 100
      ISW = 1
      NT = 0
      GO TO 95
  100 IF (IA(I) .NE. IRP) GO TO 65
      NA = NA-1
C                                  FIND THE APPROPRIATE STORAGE AREA
      ICNT = 2**NA
      NA2 = ICNT
      J = ICNT+1
      LST(ICNT) = ISUM
      ISUM = ICNT
      M = ICNT
      ICNT = ICNT-1
      I1 = NA+1
      DO 110 K = 1,NA
         KK = IWK(I1-K,6)
         L = ISUM+M
         DO 105 II = M,ISUM
            MM = L-II
            LST(ICNT) = LST(MM)+KK
            I2 = ICNT
            ICNT = ICNT-1
  105    CONTINUE
         M = I2
  110 CONTINUE
      DO 115 K = 1,NA2
         LST(K) = NF2-LST(K)
  115 CONTINUE
C                                  COMPUTE INDICES OF PRIMARY STORAGE
C                                  ARRAY
      J = LST(1)
      LOC(KX) = Y(J)
      KK = Y(J)
      KKK = Y(1)
      I1 = KK+ILM-1
      K = Y(1)
      DO 120 II = KK,I1
         Y(K) = Y(II)
         K = K+1
  120 CONTINUE
      I1 = K-1
C                                  COMPUTE DEVIATES
      DO 205 JJ = 2,NA2
         KK = KKK
         ILM = 0
C                                  FIND INDICES APPEARING IN SECONDARY
C                                  ARRAY
         K = NF21
         J = NF2-LST(JJ)
         DO 125 II = 1,NF
            IWK(II,9) = J/K
            ILM = ILM+IWK(II,9)
            IF (IWK(II,9) .EQ. 1) J = J-K
            K = K/2
  125    CONTINUE
C                                  DETERMINE SIGN
         M = NT+ILM
         IF (M .EQ. 2*(M/2)) GO TO 130
         M = -1
         GO TO 135
  130    M = 1
C                                  FIND APPROPRIATE SUBSCRIPT WITHIN
C                                  SECONDARY ARRAY
  135    DO 150 II = 1,NF
            IWK(II,7) = 0
            IF (IWK(II,5) .NE. 1) GO TO 150
            IF (IWK(II,9) .EQ. 0) GO TO 150
            K = II+1
            ISUM = 1
            IF (K .GT. NF) GO TO 145
            DO 140 J = K,NF
               IF (IWK(J,9) .EQ. 1) ISUM = ISUM*NL(J)
  140       CONTINUE
  145       IWK(II,7) = ISUM
  150    CONTINUE
C                                  COMPUTE INDICES OF SECONDARY ARRAY
         K1 = LST(JJ)
         K1 = Y(K1)
         DO 155 II = 1,NF
            IWK(II,8) = IWK(II,5)
            IF (IWK(II,8) .EQ. 1) ICNT = II
  155    CONTINUE
         GO TO 195
  160    K = 0
         NA = ICNT+1
         DO 165 II = 1,ICNT
            J = NA-II
            IF (IWK(J,5) .EQ. 0) GO TO 165
            IF (IWK(J,8) .NE. NL(J)) GO TO 170
            K = J
  165    CONTINUE
  170    IF (K .NE. 0) GO TO 175
         IWK(ICNT,8) = IWK(ICNT,8)+1
         GO TO 195
  175    ILM = K-1
  180    IF (ILM .EQ. 0) GO TO 195
         IF (IWK(ILM,5) .EQ. 1) GO TO 185
         ILM = ILM-1
         GO TO 180
  185    IWK(ILM,8) = IWK(ILM,8)+1
         ILM = ILM+1
         DO 190 II = ILM,NF
            IF (IWK(II,5) .EQ. 1) IWK(II,8) = 1
  190    CONTINUE
  195    L = 0
         DO 200 II = 1,NF
            L = L+(IWK(II,8)-1)*IWK(II,7)
  200    CONTINUE
         X = M
         Y(KK) = Y(KK)+Y(K1+L)*X
         KK = KK+1
         IF (KK .LE. I1) GO TO 160
  205 CONTINUE
C                                  FIND SUM OF SQUARES
      SUM = 0.0D0
      DO 210 II = KKK,I1
         SUM = SUM+DBLE(Y(II))*DBLE(Y(II))
  210 CONTINUE
      DO 215 II = 1,NF
         IF (IWK(II,5) .NE. 1) SUM = SUM*NL(II)
  215 CONTINUE
      XY = XY+SUM
      NALL = NALL+NDF(KX)
      SS(KX) = SUM
      KX = KX+1
      I = I+1
      ISW = 0
      GO TO 55
  220 SS(KX) = SUM1-XY
      NDF(KX) = NDFT-NALL
      KX = KX+1
      SS(KX) = SUM1
      NDF(KX) = NDFT
      LOC(KX-1) = Y(NF2)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'AGBACP')
 9005 RETURN
      END

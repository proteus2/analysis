C   IMSL ROUTINE NAME   - AFACN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - FULL FACTORIAL PLAN ANALYSIS
C
C   USAGE               - CALL AFACN (IOPT,NF,NL,Y,SS,NDF,IER)
C
C   ARGUMENTS    IOPT   - INPUT OPTION SWITCH.
C                         IF IOPT IS EQUAL TO ZERO, THE ROUTINE COMPUTES
C                           A COMPLETE SET OF MEANS FOR ALL EFFECTS IN
C                           THE FULL FACTORIAL PLAN.
C                         IF IOPT IS NOT EQUAL TO ZERO, THE ROUTINE
C                           COMPUTES THE SUMS OF SQUARES AND DEGREES OF
C                           FREEDOM FOR ALL EFFECTS IN THE FULL
C                           FACTORIAL PLAN.
C                NF     - INPUT NUMBER OF FACTORS (NUMBER OF DISTINCT
C                           SUBSCRIPTS IN THE MODEL) INCLUDING
C                           REPLICATION
C                NL     - INPUT VECTOR OF LENGTH NF CONTAINING THE
C                           NUMBER OF LEVELS FOR EACH OF THE NF FACTORS
C                Y      - VECTOR OF LENGTH 2**NF+(NL(1)+1)*(NL(2)+1)*
C                           ...*(NL(NF)+1). THE FIRST 2**NF
C                           LOCATIONS ARE WORK STORAGE.
C                         ON INPUT, Y CONTAINS THE RESPONSES IN
C                           LOCATIONS 2**NF+1,2**NF+2,...,2**NF+NL(1)*
C                           NL(2)*...*NL(NF). THE REMAINING COMPONENTS
C                           ARE UNDEFINED.
C                         ON OUTPUT, FOR IOPT EQUAL TO ZERO, THE INPUT
C                           RESPONSES ARE UNCHANGED. THE FULL SET OF
C                           MEANS IS CONTAINED IN LOCATIONS 2**NF+NL(1)*
C                           NL(2)*...*NL(NF)+1,2**NF+NL(1)*NL(2)*...*
C                           NL(NF)+2,...,2**NF+(NL(1)+1)*(NL(2)+1)*...*
C                           (NL(NF)+1).
C                           FOR IOPT NOT EQUAL TO ZERO, THE DEVIATES
C                           REQUIRED FOR COMPUTING SUMS OF SQUARES ARE
C                           CONTAINED IN LOCATIONS 2**NF+1,2**NF+2,...,
C                           2**NF+(NL(1)+1)*(NL(2)+1)*...*(NL(NF)+1).
C                SS     - OUTPUT VECTOR OF LENGTH 2**NF, CONTAINING
C                           THE SUMS OF SQUARES FOR ALL MAIN EFFECTS,
C                           INTERACTIONS, AND THE CORRECTED TOTAL SUM
C                           OF SQUARES. USED ONLY WHEN IOPT IS NOT EQUAL
C                           TO ZERO.
C                NDF    - OUTPUT VECTOR OF LENGTH 2**NF CONTAINING THE
C                           DEGREES OF FREEDOM CORRESPONDING TO
C                           COMPONENTS OF THE SS VECTOR. USED ONLY WHEN
C                           IOPT IS NOT EQUAL TO ZERO.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NF WAS SPECIFIED LESS
C                             THAN 2.
C                           IER=130 INDICATES THAT SOME NL(I) WAS
C                             SPECIFIED LESS THAN 2 FOR I=1,2,...,NF
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      TYPICAL USAGE OF AFACN WOULD HAVE THE USER CALLING
C                THE ROUTINE WITH IOPT EQUAL TO ZERO TO COMPUTE THE
C                COMPLETE SET OF MEANS FOLLOWED BY ANOTHER CALL TO
C                THE ROUTINE WITH IOPT NOT EQUAL TO ZERO TO COMPUTE
C                THE SUMS OF SQUARES AND DEGREES OF FREEDOM.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AFACN  (IOPT,NF,NL,Y,SS,NDF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NF,NL(1),NDF(1),IER
      REAL               Y(1),SS(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,LL,NF2,MM,NN,LT,KK,II,L,J,M,I2,NPM,N,INC,
     1                   J2,LD,I1,K,L1
      REAL               X,ZERO,ONE
      DOUBLE PRECISION   TEMP,DZERO
      DATA               DZERO/0.0D0/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IOPT .EQ. 0) GO TO 10
      NF2 = 2**NF
      DO 5 I = 1,NF2
         SS(I) = ZERO
    5 CONTINUE
   10 IF (NF .GE. 2) GO TO 15
C                                  TERMINAL ERROR - NUMBER OF FACTORS
C                                  IS LESS THAN 2
      IER = 129
      GO TO 9000
   15 DO 20 I = 1,NF
         IF (NL(I) .GE. 2) GO TO 20
C                                  TERMINAL ERROR - SOME NL(I) IS LESS
C                                  THAN 2
         IER = 130
         GO TO 9000
   20 CONTINUE
C                                  COMPUTE PARAMETERS FOR ARRAY
C                                  MANIPULATION AND STORE IN Y(1),Y(2),
C                                  ...,Y(2**NF)
      LL = 2**NF
      NF2 = LL+1
      Y(LL) = ONE
      IF (IOPT .NE. 0) NDF(LL) = 1
      MM = LL
      NN = LL
      LL = LL-1
      LT = NF+1
      DO 30 I = 1,NF
         KK = NL(LT-I)
         IF (IOPT .NE. 0) II = NL(LT-I)-1
         L = MM+NN
         DO 25 J = MM,NN
            M = L-J
            IF (IOPT .NE. 0) NDF(LL) = NDF(M)*II
            Y(LL) = Y(M)*KK
            I2 = LL
            LL = LL-1
   25    CONTINUE
         MM = I2
   30 CONTINUE
C                                  INITIALATION TO FIND COMPLETE SET OF
C                                  MEANS
      LL = 1
      MM = 1
      NN = 1
      LT = NF2
   35 L1 = NF2
      KK = LL
C                                  FIND NUMBER OF ELEMENTS IN THE MEAN
      NPM = NL(NF+1-NN)
   40 LT = LT+Y(MM)
C                                  FIND NUMBER OF MEANS FOR EACH FACTOR
      N = Y(MM+1)
C                                  FIND NUMBER OF ENTRIES BETWEEN
C                                  ELEMENTS OF THE MEAN
      INC = Y(NF2-KK)
      M = 1
      DO 60 I = 1,N,INC
         J2 = I+INC-1
         DO 55 J = I,J2
            L = M
            LD = M
            I1 = LT+J-1
C                                  COMPUTE THE MEAN AND STORE IN Y(I1)
            TEMP = DZERO
            DO 45 K = 1,NPM
               I2 = L1+L-1
               TEMP = TEMP+Y(I2)
               L = L+INC
   45       CONTINUE
            X = NPM
            Y(I1) = TEMP/X
            IF (IOPT .EQ. 0) GO TO 55
C                                  COMPUTE DEVIATES
            DO 50 K = 1,NPM
               I2 = L1+LD-1
               Y(I2) = Y(I2)-Y(I1)
               LD = LD+INC
   50       CONTINUE
   55    M = M+1
         M = L-INC+1
   60 CONTINUE
      IF (KK .EQ. 1) GO TO 65
      KK = KK-1
      MM = MM+1
      L1 = L1+Y(LL-KK)
      GO TO 40
   65 IF (NN .EQ. NF) GO TO 70
      LL = LL+LL
      NN = NN+1
      MM = MM+1
      GO TO 35
   70 IF (IOPT .EQ. 0) GO TO 9005
C                                  COMPUTE SUMS OF SQUARES
      L1 = NF2
      LL = L1-1
      KK = NF2-2
      DO 80 M = 1,KK
         LT = L1+Y(M)-1
         TEMP = DZERO
         DO 75 II = L1,LT
            TEMP = TEMP+DBLE(Y(II))*DBLE(Y(II))
   75    CONTINUE
         SS(M) = TEMP*Y(LL)
         LL = LL-1
         L1 = LT+1
   80 CONTINUE
      NF2 = NF2-1
      NDF(NF2) = 0
      DO 85 I = 1,KK
         SS(NF2) = SS(NF2)+SS(I)
         NDF(NF2) = NDF(NF2)+NDF(I)
   85 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'AFACN ')
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - ACRDAN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ANALYSIS OF ONE-WAY CLASSIFICATION DESIGN DATA
C
C   USAGE               - CALL ACRDAN (Y,NT,N,TM,WTV,S,GM,NDF,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N(1)+N(2)+...+N(NT)
C                           CONTAINING THE RESPONSES FOR EACH OF THE NT
C                           TREATMENTS.
C                NT     - INPUT NUMBER OF TREATMENTS.
C                N      - INPUT VECTOR OF LENGTH NT CONTAINING THE
C                           NUMBER OF RESPONSES FOR EACH OF THE NT
C                           TREATMENTS. (INPUT)
C                TM     - OUTPUT VECTOR OF LENGTH NT CONTAINING THE
C                           TREATMENT MEANS FOR EACH OF THE NT
C                           TREATMENTS.
C                WTV    - OUTPUT VECTOR OF LENGTH NT CONTAINING
C                           ESTIMATES OF THE ERROR VARIANCE WITHIN EACH
C                           TREATMENT.
C                S      - OUTPUT VECTOR OF OF LENGTH 3.
C                           S(1) CONTAINS THE TREATMENT SUM OF SQUARES
C                           S(2) CONTAINS THE ERROR SUM OF SQUARES
C                           S(3) CONTAINS THE CORRECTED TOTAL SUM OF
C                             SQUARES
C                GM     - OUTPUT GRAND MEAN OF THE RESPONSES.
C                NDF    - OUTPUT VECTOR OF LENGTH 3.
C                           NDF(1) CONTAINS THE TREATMENT DEGREES OF
C                             FREEDOM
C                           NDF(2) CONTAINS THE ERROR DEGREES OF FREEDOM
C                           NDF(3) CONTAINS THE CORRECTED TOTAL DEGREES
C                             OF FREEDOM
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NT IS LESS THAN 2.
C                           IER=130 INDICATES THAT SOME N(I),
C                             I=1,2,...,NT, IS LESS THAN 1.
C                           IER=131 INDICATES THAT ALL THE ELEMENTS IN
C                             THE VECTOR N ARE EQUAL TO 1, SO THERE
C                             ARE NO DEGREES OF FREEDOM FOR ERROR.
C                         WARNING ERROR
C                           IER=36 INDICATES THAT ONE OR MORE (BUT NOT
C                             ALL) OF THE N(I) EQUAL 1. CORRESPONDING
C                             ELEMENTS IN VECTOR WTV ARE SET TO ZERO.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE INPUT VECTOR Y HAS A SPECIFIC ARRANGEMENT OF
C                RESPONSES. THE FIRST N(1) COMPONENTS CONTAIN THE
C                RESPONSES FOR TREATMENT 1. THE NEXT N(2) COMPONENTS
C                CONTAIN THE RESPONSES FOR TREATMENT 2. THIS PATTERN
C                CONTINUES FOR THE NT TREATMENTS ENDING WITH THE
C                LAST N(NT) COMPONENTS CONTAINING THE RESPONSES FOR
C                TREATMENT NT.
C            2.  THE ORDER OF VALUES IN OUTPUT VECTORS TM AND WTV
C                CORRESPONDS TO THE ORDER OF TREATMENTS REPRESENTED
C                IN THE VECTOR Y. THAT IS, ELEMENTS TM(I) AND WTV(I)
C                CONTAIN THE TREATMENT MEAN AND ESTIMATE OF ERROR
C                VARIANCE, RESPECTIVELY, FOR TREATMENT I, I=1,2,...,NT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ACRDAN (Y,NT,N,TM,WTV,S,GM,NDF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NT,N(NT),NDF(3),IER
      REAL               Y(1),TM(NT),WTV(NT),S(3),GM
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NN,I,K,KK,J,NI
      REAL               ZERO
      DOUBLE PRECISION   TEMP,Z,SS,XX,TMI
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF(NT .GE. 2) GO TO 5
C                                  TERMINAL ERROR - NUMBER OF
C                                  TREATMENTS IS LESS THAN 2
      IER=129
      GO TO 9000
    5 NN=0
      DO 10 I = 1,NT
C                                  TERMINAL ERROR - SOME N(I) IS LESS
C                                  THAN 1
         IF (N(I) .LT. 1) GO TO 15
         NN = NN+N(I)
   10 CONTINUE
C                                  TERMINAL ERROR - ALL ELEMENTS IN
C                                    VECTOR N ARE EQUAL 1
      IF (NN .NE. NT) GO TO 20
      IER=131
      GO TO 9000
   15 IER=130
      GO TO 9000
C                                  COMPUTE TREATMENT MEANS
   20 Z=0.D0
      K=1
      KK=0
      DO 30 I = 1,NT
         TEMP=0.D0
         KK=N(I)+KK
         DO 25 J = K,KK
            TEMP = TEMP+Y(J)
   25    CONTINUE
         TM(I)=TEMP/N(I)
         K=KK+1
C                                  FIND GRAND MEAN
         Z = Z+TEMP
   30 CONTINUE
      GM=Z/NN
      SS=0.D0
      TEMP=0.D0
      KK=0
      K=1
C                                  FIND WITHIN TREATMENT ERROR
C                                  VARIANCE
      DO 50 I = 1,NT
         TMI=TM(I)
         NI=N(I)
         XX=TMI-GM
         SS=SS+NI*XX*XX
         Z=0.D0
         KK=KK+NI
         DO 35 J = K,KK
            XX=Y(J)-GM
            TEMP=TEMP+XX*XX
            XX=Y(J)-TMI
            Z = Z+XX*XX
   35    CONTINUE
         IF (NI .EQ. 1) GO TO 40
         WTV(I)=Z/(NI-1)
         GO TO 45
   40    WTV(I) = ZERO
         IER = 36
   45    K = KK+1
   50 CONTINUE
C                                  FIND SUMS OF SQUARES
      S(1)=SS
      S(2)=TEMP-SS
      S(3)=TEMP
C                                  FIND DEGREES OF FREEDOM
      NDF(1)=NT-1
      NDF(3)=NN-1
      NDF(2)=NDF(3)-NDF(1)
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'ACRDAN')
 9005 RETURN
      END

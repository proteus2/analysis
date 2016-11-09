C   IMSL ROUTINE NAME   - ANCOV1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COVARIANCE ANALYSIS FOR ONE-WAY
C                           CLASSIFICATION DESIGN DATA
C
C   USAGE               - CALL ANCOV1 (XY,NOP,IXY,XYM,IXYM,TM,SXY,VARB,
C                           VART,SS,NDF,WK,IER)
C
C   ARGUMENTS    XY     - INPUT (NOP(4)+NOP(5)+...+NOP(NOP(1)+3)) BY
C                           (NOP(2)+1) MATRIX CONTAINING THE DATA FOR
C                           EACH COVARIATE AND THE RESPONSE VARIABLE
C                           AT EACH CLASSIFICATION LEVEL (TREATMENT).
C                           DATA FOR EACH CLASSIFICATION LEVEL MUST
C                           APPEAR IN CONTIGUOUS ROWS OF XY WITH THE
C                           RESPONSE VARIABLE SETTINGS APPEARING IN
C                           THE LAST COLUMN.
C                NOP    - INPUT VECTOR OF LENGTH NOP(1)+3.
C                           NOP(1) CONTAINS THE NUMBER OF TREATMENTS.
C                             NOP(1) MUST BE GREATER THAN OR EQUAL TO
C                             TWO.
C                           NOP(2) CONTAINS THE NUMBER OF COVARIATES.
C                             NOP(2) MUST BE GREATER THAN OR EQUAL TO
C                             ONE.
C                           NOP(3) CONTAINS AN OPTION FOR TESTING
C                             THE HOMOGENEITY OF THE REGRESSION
C                             RELATIONSHIPS WITHIN TREATMENTS.
C                             NOP(3)=0 IMPLIES TEST IS NOT DESIRED.
C                             NOP(3) NOT EQUAL TO ZERO IMPLIES TEST
C                             IS DESIRED.
C                           NOP(4),NOP(5),...,NOP(NOP(1)+3) CONTAINS
C                             THE NUMBER OF OBSERVATIONS FOR EACH OF
C                             THE NOP(1) TREATMENTS.
C                             EACH OF NOP(4),NOP(5),...,NOP(NOP(1)+3)
C                             MUST BE GREATER THAN OR EQUAL TO ONE AND
C                             THEIR SUM MUST BE GREATER THAN OR EQUAL
C                             TO NOP(1)+NOP(2)+1. IN ADDITION, WHEN
C                             NOP(3) IS NOT ZERO, EACH OF NOP(4),NOP(5),
C                             ...,NOP(NOP(1)+3) MUST BE GREATER THAN OR
C                             EQUAL TO NOP(2)+1 AND THEIR SUM MUST BE
C                             GREATER THAN OR EQUAL TO
C                             NOP(1)*(NOP(2)+1)+1.
C                IXY    - INPUT ROW DIMENSION OF THE MATRIX XY
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                XYM    - OUTPUT (NOP(1)+1) BY (NOP(2)+1) MATRIX
C                           CONTAINING THE COVARIATE AND RESPONSE
C                           VARIABLE MEANS FOR EACH TREATMENT. THE
C                           LAST ROW CONTAINS THE OVERALL MEANS
C                           FOR EACH VARIABLE.
C                IXYM   - INPUT ROW DIMENSION OF THE MATRIX XYM
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           IN THE CALLING PROGRAM.
C                TM     - OUTPUT VECTOR OF LENGTH NOP(1)+NOP(2)
C                           CONTAINING THE ADJUSTED TREATMENT MEANS
C                           AND THE REGRESSION COEFFICIENT
C                           ESTIMATES.
C                SXY    - OUTPUT 3 BY (NOP(2)+1)*(NOP(2)+2)/2
C                           MATRIX CONTAINING THREE
C                           (NOP(2)+1) BY (NOP(2)+1) MATRICES EACH
C                           STORED IN SYMMETRIC STORAGE MODE. MATRIX 1,
C                           STORED IN THE FIRST ROW OF SXY,
C                           CONTAINS THE BETWEEN-TREATMENT SUMS OF
C                           SQUARES AND CROSS PRODUCTS.
C                           MATRIX 2, STORED IN THE SECOND ROW OF SXY,
C                           CONTAINS THE POOLED WITHIN-
C                           TREATMENT SUMS OF SQUARES AND CROSS
C                           PRODUCTS.
C                           MATRIX 3, STORED IN THE THIRD ROW OF SXY,
C                           CONTAINS THE CORRECTED TOTAL SUMS
C                           OF SQUARES AND CROSS PRODUCTS.
C                VARB   - OUTPUT VECTOR OF LENGTH NOP(2)*(NOP(2)+1)/2
C                           CONTAINING THE NOP(2) BY NOP(2)
C                           INVERSE OF THE POOLED WITHIN-TREATMENT
C                           MATRIX OF SUMS OF SQUARES AND CROSS
C                           PRODUCTS FOR THE COVARIATES, STORED IN
C                           SYMMETRIC STORAGE MODE.
C                VART   - OUTPUT VECTOR OF LENGTH NOP(1)*(NOP(1)+1)/2
C                           CONTAINING THE NOP(1) BY NOP(1)
C                           VARIANCE-COVARIANCE MATRIX OF THE
C                           ADJUSTED TREATMENT MEANS, DIVIDED BY
C                           THE ERROR MEAN SQUARE. THE MATRIX IS
C                           STORED IN SYMMETRIC STORAGE MODE.
C                SS     - OUTPUT VECTOR OF LENGTH 6 CONTAINING SUMS
C                           OF SQUARES.
C                           SS(1) CONTAINS THE SUM OF SQUARES FOR
C                             ADJUSTED TREATMENTS.
C                           SS(2) CONTAINS THE SUM OF SQUARES FOR
C                             REGRESSION.
C                           SS(3) CONTAINS THE ERROR SUM OF SQUARES.
C                           SS(4) CONTAINS THE SUM OF SQUARES FOR
C                             REGRESSION HOMOGENEITY.
C                           SS(5) CONTAINS THE ERROR SUM OF SQUARES FOR
C                             TESTING REGRESSION HOMOGENEITY.
C                           SS(6) CONTAINS THE CORRECTED TOTAL SUM OF
C                             SQUARES. IF NOP(3) =0, SS(4) AND SS(5) ARE
C                             NOT DEFINED AND ARE SET TO ZERO IN THE
C                             PROGRAM.
C                NDF    - OUTPUT VECTOR OF LENGTH 6 CONTAINING
C                           DEGREES OF FREEDOM CORRESPONDING TO
C                           COMPONENTS OF THE SS VECTOR. IF NOP(3) = 0,
C                           NDF(4) AND NDF(5) ARE NOT DEFINED AND ARE
C                           SET TO ZERO IN THE PROGRAM.
C                WK     - WORK AREA OF LENGTH
C                           (NOP(2)+1)*(NOP(2)+2)/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 IMPLIES THE NUMBER OF TREATMENTS
C                             (NOP(1)) WAS SPECIFIED LESS THAN 2.
C                           IER=130 IMPLIES THE NUMBER OF COVARIATES
C                             (NOP(2)) WAS SPECIFIED LESS THAN 1.
C                           IER=131 IMPLIES THAT ONE OR MORE OF NOP(4),
C                             NOP(5),...,NOP(NOP(1)+3) ARE LESS THAN
C                             ONE, OR THAT THEIR SUM IS LESS THAN
C                             NOP(1)+NOP(2)+1.
C                           IER=132 IMPLIES WHEN NOP(3) IS NOT ZERO,
C                             THAT ONE OR MORE OF NOP(4),NOP(5),...,
C                             NOP(NOP(1)+3) ARE LESS THAN NOP(2)+1,
C                             OR THAT THEIR SUM IS LESS THAN
C                             NOP(1)*(NOP(2)+1)+1.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE ESTIMATED VARIANCE-COVARIANCE MATRIX OF THE
C                REGRESSION COEFFICIENT ESTIMATES IS EASILY AVAILABLE
C                FROM THE OUTPUT VECTOR VARB. MULTIPLICATION OF EACH
C                ELEMENT OF VARB BY SS(3)/NDF(3) IS REQUIRED.
C            2.  WHEN HIGH CORRELATIONS EXIST AMONG SOME OF THE
C                COVARIATES, AN OVERFLOW CAN OCCUR. DELETION OF ONE
C                OR MORE COVARIATES CAN BE DONE TO ALLEVIATE THE
C                PROBLEM.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ANCOV1 (XY,NOP,IXY,XYM,IXYM,TM,SXY,VARB,VART,SS,
     1                   NDF,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NOP(1),IXY,IXYM,NDF(6),IER
      REAL               XY(IXY,1),XYM(IXYM,1),TM(1),SXY(3,1),
     1                   VARB(1),VART(1),SS(6),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBACK,IH,II,IJ,IL,INX,INY,ISUM,IX,J,
     1                   JJ,JX,J1,K,KK,KX,L,LL,L1,M,MD,M1,NA,NABP1,
     2                   NAPI,NAP1,NAP3,NATNB,NB,NBP1,NN,NSUM,ITEMP
      REAL               ZERO,X
      DOUBLE PRECISION   SUM,GSUM,FJJ,SNM
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
C                                  CHECK FOR TERMINAL ERRORS
      IF (NOP(1) .GE. 2) GO TO 5
      IER=129
      GO TO 9000
    5 IF (NOP(2) .GE. 1) GO TO 10
      IER=130
      GO TO 9000
   10 NA=NOP(1)
      NB=NOP(2)
      NAP1=NA+1
      NBP1=NB+1
      NABP1=NA+NB+1
      NAP3=NA+3
      NSUM=0
      DO 15 I=4,NAP3
         IF (NOP(I) .LT. 1) GO TO 20
         NSUM=NSUM+NOP(I)
   15 CONTINUE
      IF (NSUM .GE. NABP1) GO TO 25
   20 IER=131
      GO TO 9000
   25 IF (NOP(3) .EQ. 0) GO TO 40
      NATNB=NA*NBP1+1
      ISUM=0
      DO 30 I=4,NAP3
         IF (NOP(I) .LT. NBP1) GO TO 35
         ISUM=ISUM+NOP(I)
   30 CONTINUE
      IF (ISUM .GE. NATNB) GO TO 40
   35 IER=132
      GO TO 9000
C                                  COMPUTE COVARIATE AND RESPONSE
C                                    VARIABLE MEANS FOR EACH TREATMENT.
   40 DO 55 K=1,NBP1
         J1=1
         JJ=0
         DO 50 I=1,NA
            FJJ=NOP(I+3)
            JJ=NOP(I+3)+JJ
            SUM=0.0D0
            DO 45 J=J1,JJ
               SUM=SUM+DBLE(XY(J,K))
   45       CONTINUE
            XYM(I,K)=SUM/FJJ
            J1=JJ+1
   50    CONTINUE
   55 CONTINUE
C                                  COMPUTE GRAND MEAN FOR EACH VARIABLE
      DO 65 K=1,NBP1
         GSUM=0.0D0
         SUM=0.0D0
         DO 60 I=1,NA
            SNM=NOP(I+3)
            SUM=SUM+SNM
            GSUM=GSUM+SNM*DBLE(XYM(I,K))
   60    CONTINUE
         XYM(NAP1,K)=GSUM/SUM
   65 CONTINUE
C                                  COMPUTE THE WITHIN-TREATMENT SUMS
C                                    OF SQUARES AND CROSS PRODUCTS
C                                    FOR EACH TREATMENT AND COMPUTE
C                                    SS(5)
      SS(4)=ZERO
      SS(5)=ZERO
      NN=(NBP1*(NBP1+1))/2
      DO 70 I=1,NN
         SXY(2,I)=ZERO
   70 CONTINUE
      LL=0
      L1=1
      K=1
   75 LL=NOP(3+K)+LL
      DO 85 J=1,NBP1
         INY=(J*(J-1))/2
         DO 85 I=1,J
            INX=INY+I
            SUM=0.0D0
            DO 80 L=L1,LL
               SUM=SUM+DBLE(XY(L,I)-XYM(K,I))*DBLE(XY(L,J)-XYM(K,J))
   80       CONTINUE
            WK(INX)=SUM
            SXY(2,INX)=SXY(2,INX)+WK(INX)
   85 CONTINUE
      IF (NOP(3) .EQ. 0) GO TO 130
C                                  PERFORM NOP(2) JORDAN REDUCTIONS
C                                    ON EACH WITHIN-TREATMENT SUMS
C                                    OF SQUARES AND CROSS PRODUCTS
C                                    MATRIX
      ASSIGN 125 TO IBACK
   90 DO 120 M=1,NB
         M1=(M*(M-1))/2
         MD=M1+M
         X=1.0/WK(MD)
         SXY(3,MD)=X
         DO 110 I=1,NBP1
            IL=1
            IF (I .NE. M) IL=-1
            II=(I*(I-1))/2
            ITEMP = I
            DO 105 J=1,ITEMP
               IH=1
               JJ=(J*(J-1))/2
               IF (ITEMP.EQ.J .AND.ITEMP.EQ.M) GO TO 110
               IX=II+J
               IF (ITEMP.NE.M .AND. J.NE.M) GO TO 95
               SXY(3,IX)=WK(IX)*X*IL
               GO TO 105
   95          JX=II+M
               KX=JJ+M
               IF (M.GT.ITEMP) JX=M1+ITEMP
               IF (M .GT. J) KX=M1+J
               IF (M.EQ.1) GO TO 100
               IF (ITEMP .GE. J .AND. M .GT. ITEMP) IH=-1
  100          SXY(3,IX)=WK(IX)-(WK(JX)*IH*WK(KX)*X)
  105       CONTINUE
  110    CONTINUE
         DO 115 I=1,NN
            WK(I)=SXY(3,I)
  115    CONTINUE
  120 CONTINUE
      GO TO IBACK, (125,160,180)
  125 SS(5)=SS(5)+SXY(3,NN)
  130 L1=LL+1
      IF (K.EQ.NA) GO TO 135
      K=K+1
      GO TO 75
C                                  COMPUTE MATRIX 1 OF SXY
  135 DO 150 J=1,NBP1
         INY=(J*(J-1))/2
         DO 145 I=1,J
            SUM=0.0D0
            INX=INY+I
            DO 140 K=1,NA
               FJJ=NOP(K+3)
               SUM=SUM+FJJ*DBLE(XYM(K,I)-XYM(NAP1,I))
     1           *DBLE(XYM(K,J)-XYM(NAP1,J))
  140       CONTINUE
            SXY(1,INX)=SUM
  145    CONTINUE
  150 CONTINUE
C                                  PERFORM NOP(2) JORDAN REDUCTIONS
C                                    ON MATRIX 2 OF SXY AND COMPUTE
C                                    TM(NOP(1)+1),...,TM(NOP(1)+NOP(2))
      ASSIGN 160 TO IBACK
      DO 155 I=1,NN
         WK(I)=SXY(2,I)
  155 CONTINUE
      GO TO 90
  160 INY=(NB*(NB+1))/2
      SS(2)=ZERO
      DO 165 I=1,NB
         NAPI=NA+I
         INX=INY+I
         TM(NAPI)=-SXY(3,INX)
         SS(2)=SS(2)+TM(NAPI)*SXY(2,INX)
  165 CONTINUE
C                                  COMPUTE VARB VECTOR
      SS(3)=SXY(3,NN)
      DO 170 I=1,INY
         VARB(I)=SXY(3,I)
  170 CONTINUE
C                                  COMPUTE SS(1) AND PERFORM NOP(2)
C                                    JORDAN REDUCTIONS ON MATRIX 3
C                                    OF SXY
      ASSIGN 180 TO IBACK
      DO 175 I=1,NN
         WK(I)=SXY(1,I)+SXY(2,I)
  175 CONTINUE
      GO TO 90
  180 SS(1)=SXY(3,NN)-SS(3)
C                                  COMPUTE SS(4) IF NOP(3) IS NOT
C                                    ZERO
      IF (NOP(3) .EQ. 0) GO TO 185
      SS(4)=SS(3)-SS(5)
C
C
C                                  COMPUTE MATRIX 3 OF SXY
  185 DO 190 I=1,NN
         SXY(3,I)=SXY(1,I)+SXY(2,I)
  190 CONTINUE
C                                  COMPUTE TM(I),I=1,2,...,NOP(1)
      SS(6)=SXY(3,NN)
      DO 200 I=1,NA
         TM(I)=ZERO
         DO 195 J=1,NB
            TM(I)=TM(I)+TM(NA+J)*DBLE(XYM(I,J)-XYM(NAP1,J))
  195    CONTINUE
         TM(I)=XYM(I,NBP1)-TM(I)
  200 CONTINUE
C                                  COMPUTE VART VECTOR
      DO 210 J=1,NA
         JJ=(J*(J-1))/2
         DO 210 I=1,J
            IJ=JJ+I
            SUM=0.0D0
            SNM=NOP(I+3)
            DO 205 K=1,NB
               KK=(K*(K-1))/2
               DO 205 L=1,NB
                  LL=(L*(L-1))/2
                  INX=KK+L
                  IF (L.GT.K) INX=LL+K
                  SUM=SUM+DBLE(XYM(I,K)-XYM(NAP1,K))*DBLE((XYM(J,L)-
     1              XYM(NAP1,L))*VARB(INX))
  205       CONTINUE
            VART(IJ)=SUM
            IF (I.EQ.J) VART(IJ)=SUM+DBLE(1.0)/SNM
  210 CONTINUE
C
C                                  COMPUTE NDF VECTOR
  215 NDF(1)=NA-1
      NDF(2)=NB
      NN=0
      DO 220 I=1,NA
         NN=NN+NOP(I+3)
  220 CONTINUE
      NDF(3)=NN-NA-NB
      NDF(6)=NN-1
      IF (NOP(3).NE.0) GO TO 225
      NDF(4)=0
      NDF(5)=0
      GO TO 9005
  225 NDF(5)=NN-NA*NBP1
      NDF(4)=NDF(3)-NDF(5)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'ANCOV1')
 9005 RETURN
      END

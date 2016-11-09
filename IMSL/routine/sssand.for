C   IMSL ROUTINE NAME   - SSSAND
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SIMPLE RANDOM SAMPLING WITH CONTINUOUS
C                           DATA - INFERENCES REGARDING THE POPULATION
C                           MEAN AND TOTAL
C
C   USAGE               - CALL SSSAND (Y,NBR,ALPHA,TEMP,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT SUBVECTOR OF LENGTH NBR(2) OF THE VECTOR
C                           (CALL IT YY) CONTAINING THE ENTIRE RANDOM
C                           SAMPLE. THE LAST SUBVECTOR OF YY MAY HAVE
C                           FEWER THAN NBR(2) ELEMENTS.
C                NBR    - INPUT VECTOR OF LENGTH 8. NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF OBSERVATIONS IN YY.
C                           I=2, NUMBER OF OBSERVATIONS IN EACH SUB-
C                             VECTOR Y, NOT INCLUDING THE LAST SUB-
C                             VECTOR, WHERE THE NUMBER MAY BE LESS THAN
C                             OR EQUAL TO NBR(2). HOWEVER NBR(2)
C                             SHOULD BE THE SAME FOR ALL CALLS.
C                           I=3, THE NUMBER OF THE SUBVECTOR STORED
C                             IN Y. SEE REMARKS.
C                           I=4, THE TEMPORARY MEAN INDICATOR. IF
C                             NBR(4) = 0, THE USER SUPPLIES THE TEMPOR-
C                             ARY MEAN IN TEMP. OTHERWISE, THE FIRST
C                             ELEMENT OF YY (OR FIRST ELEMENT OF Y WHEN
C                             NBR(3) = 1) IS UTILIZED.
C                           I=5, SUBPOPULATION INDICATOR.
C                             IF NBR(5) = 0, THE INPUT DATA IS A SAMPLE
C                               FROM A POPULATION.
C                             IF NBR(5) IS NEGATIVE THE SAMPLE IS FROM
C                               A SUBPOPULATION OF UNKNOWN SIZE.
C                             IF NBR(5) IS POSITIVE THE SAMPLE IS FROM
C                               A SUBPOPULATION OF KNOWN SIZE.
C                           I=6, SIZE OF THE SAMPLED POPULATION. NOT
C                             REQUIRED IF NBR(5) IS POSITIVE.
C                           I=7, SIZE OF THE SAMPLED SUBPOPULATION.
C                             REQUIRED ONLY WHEN NBR(5) IS POSITIVE.
C                           I=8, SIZE OF THE SAMPLE FROM THE POPULATION
C                             FOR WHICH NBR(1) WERE TAKEN TO CONSTITUTE
C                             THE SAMPLE FROM THE SUBPOPULATION OF IN-
C                             TEREST. REQUIRED ONLY WHEN NBR(5) IS
C                             NEGATIVE.
C                ALPHA  - INPUT VALUE IN THE EXCLUSIVE INTERVAL (0,1)
C                           USED FOR COMPUTING 100(1-ALPHA) PERCENT CON-
C                           FIDENCE INTERVALS FOR THE MEAN AND TOTAL
C                           PARAMETERS. 0.05 IS A COMMON CHOICE.
C                TEMP   - INPUT TEMPORARY MEAN. REQUIRED ONLY IF
C                           NBR(4) = 0.
C                STAT   - OUTPUT VECTOR OF LENGTH 9. STAT(I) CONTAINS,
C                         WHEN
C                           I=1, ESTIMATE OF THE MEAN.
C                           I=2, ESTIMATE OF THE TOTAL.
C                           I=3, WITHIN SAMPLE VARIANCE ESTIMATE.
C                           I=4, VARIANCE ESTIMATE OF THE MEAN ESTIMATE.
C                           I=5, VARIANCE ESTIMATE OF THE TOTAL ESTIMATE
C                           I=6, LOWER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=7, UPPER CONFIDENCE LIMIT FOR THE MEAN.
C                           I=8, LOWER CONFIDENCE LIMIT FOR THE TOTAL.
C                           I=9, UPPER CONFIDENCE LIMIT FOR THE TOTAL.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES ONE OR MORE OF NBR(1),
C                             NBR(6) (WHEN NBR(5) IS NONPOSITIVE),
C                             NBR(7) (WHEN NBR(5) IS POSITIVE), AND
C                             NBR(8) (WHEN NBR(5) IS NEGATIVE), ARE
C                             LESS THAN 2.
C                           IER = 130 INDICATES NBR(3) IS LESS THAN ONE
C                             OR THAT NBR(2)*(NBR(3)-1) EXCEEDS NBR(1).
C                           IER = 131 INDICATES THAT NBR(2) IS LESS
C                             THAN 1 OR THAT ALPHA IS NOT IN THE EX-
C                             CLUSIVE INTERVAL (0,1).
C
C   REQD. IMSL ROUTINES - MDNRIS,MDSTI,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      BETWEEN THE FIRST CALL AND THE LAST CALL (M-TH CALL) TO
C                SSSAND ONLY NBR(3) MAY BE MODIFIED AND IT SHOULD FOLLOW
C                THE PATTERN 1,2,...,M. THOUGH THIS PATTERN IS THE OB-
C                VIOUS ONE TO FOLLOW, IT IS NOT NECESSARY IN ITS EN-
C                TIRETY. FOR CALLS 2,3,...,M-1, NBR(3) MAY TAKE ANY
C                VALUE IN THE SET (2,3,...,M-1). ON THE FIRST CALL
C                NBR(3) MUST EQUAL 1, AND ON THE LAST CALL NBR(3) MUST
C                EQUAL M.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSSAND (Y,NBR,ALPHA,TEMP,STAT,IER)
C
      DIMENSION          Y(1),NBR(1),STAT(1)
      DOUBLE PRECISION   WK,WKA
      EQUIVALENCE        (F,SS),(T1,SSQ),(T2,FBR)
      DATA               ZERO,ONE,TWO/0.0,1.0,2.0/
C                                  INITIALIZE IER
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IF (NBR(1) .LT. 2 .OR. (NBR(5) .LE. 0 .AND. NBR(6) .LT. 2))
     1   GO TO 5
      IF (NBR(5) .GT. 0 .AND. NBR(7) .LT. 2) GO TO 5
      IF (.NOT. (NBR(5) .LT. 0 .AND. NBR(8) .LT. 2)) GO TO 10
C                                  TERMINAL ERROR - INVALID NBR VALUES
    5 IER = 129
      GO TO 9000
   10 IF (NBR(3) .LT. 1) GO TO 15
      IF (NBR(2)*(NBR(3)-1) .LE. NBR(1)) GO TO 20
   15 IER = 130
      GO TO 9000
   20 IF (NBR(2) .LT. 1) GO TO 25
      IF (ALPHA .GT. ZERO .AND. ALPHA .LT. ONE) GO TO 30
   25 IER = 131
      GO TO 9000
   30 MM = (NBR(1)+NBR(2)-1)/NBR(2)
      NN=NBR(2)
      IF (NBR(3) .EQ. MM) NN = NBR(1)-(MM-1)*NN
      IF (NBR(3) .NE. 1) GO TO 40
C                                  MOVE THE FIRST ELEMENT OF Y INTO TEMP
      IF (NBR(4) .NE. 0) TEMP = Y(1)
C                                  INITIALIZE STAT VECTOR
      DO 35 I=1,9
         STAT(I) = ZERO
   35 CONTINUE
   40 WK = STAT(1)
C                                  COMPUTE ESTIMATE OF THE MEAN
      DO 45 J=1,NN
         WK = WK +DBLE(Y(J))
   45 CONTINUE
      STAT(1) = WK
      WKA = STAT(3)
C                                  COMPUTE APPROXIMATION TO WITHIN
C                                    SAMPLE CORRECTED SUM OF SQUARES
      DO 50 I=1,NN
         WK = DBLE(Y(I))-DBLE(TEMP)
         WKA = WKA+WK**2
   50 CONTINUE
      STAT(3) = WKA
C                                  RETURN - IF THIS IS NOT THE FINAL CAL
      IF (NBR(3) .NE. MM) GO TO 9005
      FBR1 = NBR(1)
      RFBR1 = ONE/FBR1
      YBAR = STAT(1)*RFBR1
      WK = YBAR-TEMP
      SS = (STAT(3)-DBLE(FBR1)*WK*WK)/(FBR1-ONE)
C                                  CHECK SUBPOPULATION INDICATOR(NBR(5))
C                                    TO COMPUTE THE MEAN ESTIMATE AND TH
C                                    TOTAL ESTIMATE
      FBR = NBR(6)
      IF (NBR(5) .GT. 0) FBR = NBR(7)
      RFBR = ONE/FBR
      IF (NBR(5) .GE. 0) GO TO 55
      FBR8 = NBR(8)
      RFBR8 = ONE/FBR8
      YBARP = STAT(1)*RFBR8
      WK = YBARP-TEMP
      SSQ = STAT(3)-FBR1*WK**2-TWO*WK*(DBLE(STAT(1))*(ONE-FBR1*RFBR8))
      SSQ = (SSQ+DBLE(FBR8-FBR1)*DBLE(YBARP)*DBLE(YBARP))/(FBR8-ONE)
      STAT(2) = FBR*RFBR8*STAT(1)
      FBR8 = ONE-FBR8*RFBR
      STAT(4) = (SS*RFBR1)*FBR8
      STAT(5) = ((FBR**2)*RFBR8)*SSQ*FBR8
      GO TO 60
   55 STAT(2) = FBR*RFBR1*STAT(1)
      STAT(4) = (SS*RFBR1)*(ONE-FBR1*RFBR)
      STAT(5) = (FBR**2)*STAT(4)
   60 STAT(1) = YBAR
      STAT(3) =SS
      F = NBR(1)-1
C                                  CALL ROUTINE MDSTI
      CALL MDSTI (ALPHA,F,X,JER)
      T1 = SQRT(STAT(4))*X
      T2 = SQRT(STAT(5))*X
      STAT(6) = STAT(1)-T1
      STAT(7) = STAT(1)+T1
      STAT(8) = STAT(2)-T2
      STAT(9) = STAT(2)+T2
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HSSSAND)
 9005 RETURN
      END

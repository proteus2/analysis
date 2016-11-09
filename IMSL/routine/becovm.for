C   IMSL ROUTINE NAME   - BECOVM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MEANS AND VARIANCE-COVARIANCE MATRIX
C
C   USAGE               - CALL BECOVM (X,IX,NBR,TEMP,XM,VCV,IER)
C
C   ARGUMENTS    X      - ON INPUT, X IS AN NBR(3) BY NBR(1) SUBMATRIX
C                           OF THE MATRIX (CALL IT XX) OF DATA FOR
C                           WHICH MEANS, VARIANCES AND COVARIANCES,
C                           OR CORRECTED SUMS OF SQUARES AND
C                           CROSS-PRODUCTS ARE DESIRED. THE LAST
C                           SUBMATRIX IN XX MAY HAVE FEWER THAN NBR(3)
C                           ROWS. SEE EXAMPLE.
C                         ON OUTPUT, THE ROWS OF X HAVE BEEN ADJUSTED
C                           BY THE TEMPORARY MEANS.
C                IX     - INPUT, ROW DIMENSION OF X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NBR    - INPUT VECTOR OF LENGTH 6. NBR(I) CONTAINS,
C                         WHEN
C                           I=1, NUMBER OF VARIABLES.
C                           I=2, NUMBER OF OBSERVATIONS PER VARIABLE
C                             IN XX.
C                           I=3, NUMBER OF OBSERVATIONS PER VARIABLE IN
C                             EACH SUBMATRIX X, NOT INCLUDING THE LAST
C                             SUBMATRIX WHERE THE NUMBER MAY BE LESS
C                             THAN OR EQUAL TO NBR(3). HOWEVER, NBR(3)
C                             SHOULD BE THE SAME FOR ALL CALLS.
C                           I=4, THE NUMBER OF THE SUBMATRIX STORED IN
C                             X. SEE EXAMPLE.
C                           I=5, THE TEMPORARY MEAN INDICATOR. IF
C                             NBR(5) = 0, THE USER SUPPLIES TEMPORARY
C                             MEANS IN TEMP. OTHERWISE, THE FIRST ROW
C                             OF XX (OR FIRST ROW OF X WHEN NBR(4) = 1)
C                             IS UTILIZED.
C                           I=6, THE VCV OPTION. IF NBR(6) = 0, VCV
C                             CONTAINS THE VARIANCE-COVARIANCE MATRIX.
C                             OTHERWISE, VCV CONTAINS THE CORRECTED
C                             SUMS OF SQUARES AND CROSS-PRODUCTS
C                             MATRIX.
C                TEMP   - INPUT VECTOR OF LENGTH NBR(1). IF NBR(5) = 0,
C                           TEMP MUST CONTAIN THE TEMPORARY MEANS WHEN
C                           NBR(4) = 1. OTHERWISE, TEMP IS WORK STORAGE.
C                XM     - OUTPUT VECTOR OF LENGTH NBR(1) CONTAINING THE
C                           VARIABLE MEANS.
C                VCV    - OUTPUT NBR(1) BY NBR(1) MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE REQUIRING
C                           (NBR(1)*(NBR(1)+1))/2 STORAGE LOCATIONS.
C                           VCV CONTAINS THE VARIANCE-COVARIANCE MATRIX
C                           OR THE CORRECTED SUMS OF SQUARES AND
C                           CROSS-PRODUCTS MATRIX, AS CONTROLLED BY THE
C                           VCV OPTION, NBR(6).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NBR(4) IS LESS THAN
C                             1 OR NBR(3)*(NBR(4)-1) EXCEEDS NBR(2).
C                           IER=130 INDICATES THAT NBR(1) IS LESS THAN
C                             1 OR NBR(2) IS LESS THAN 2 OR NBR(3)
C                             EXCEEDS NBR(2).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BECOVM (X,IX,NBR,TEMP,XM,VCV,IER)
C
      DIMENSION          X(IX,1),NBR(1),TEMP(1),XM(1),VCV(1)
      DOUBLE PRECISION   WK
      DATA               ZERO,ONE/0.0,1.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
      IF (NBR(4) .LT. 1) GO TO 10
    5 IF (NBR(3)*(NBR(4)-1) .LE. NBR(2)) GO TO 15
C                                  TERMINAL ERROR - INVALID NBR VALUES
   10 IER = 129
      GO TO 9000
   15 IF (NBR(1) .GE. 1 .AND. NBR(2) .GE. 2 .AND. NBR(3) .LE. NBR(2))
     1   GO TO 20
      IER = 130
      GO TO 9000
   20 MM = (NBR(2)+NBR(3)-1)/NBR(3)
      N1 = NBR(1)
      NN1 = (N1*(N1+1))/2
      IF (NBR(4) .NE. 1) GO TO 40
C                                  INITIALIZE XM, VCV VECTORS
      DO 25 I=1,N1
         XM(I) = ZERO
   25 CONTINUE
      DO 30 I=1,NN1
         VCV(I) = ZERO
   30 CONTINUE
      IF (NBR(5) .EQ. 0) GO TO 40
C                                  MOVE THE FIRST ROW OF X INTO TEMP
      DO 35 J=1,N1
         TEMP(J) = X(1,J)
   35 CONTINUE
   40 NN = NBR(3)
   45 IF (NBR(4) .EQ. MM) NN = NBR(2)-(MM-1)*NBR(3)
C                                  COMPUTE VARIABLE MEANS
      DO 55 J=1,N1
         WK = XM(J)
         DO 50 I=1,NN
            WK = WK+DBLE(X(I,J))
   50    CONTINUE
         XM(J) = WK
   55 CONTINUE
C                                  REPLACE THE ELEMENTS OF SUBMATRIX
C                                    X WITH THEIR DEVIATIONS FROM THE
C                                    TEMPORARY MEANS
      DO 60 I=1,NN
         DO 60 J=1,N1
            X(I,J) = X(I,J)-TEMP(J)
   60 CONTINUE
C                                  COMPUTE THE APPROXIMATE SUMS OF
C                                    SQUARES AND CROSS-PRODUCTS MATRIX
      INX = 0
      DO 70 L=1,N1
         DO 70 J=1,L
            INX = INX+1
            WK = VCV(INX)
            DO 65 I=1,NN
               WK = WK+DBLE(X(I,J))*DBLE(X(I,L))
   65       CONTINUE
            VCV(INX) = WK
   70 CONTINUE
      IF (NBR(4) .NE. MM) GO TO 9005
      FNBR2 = NBR(2)
      RNBR2 = ONE/FNBR2
      RNBR3 = ONE/(FNBR2-ONE)
C                                  ADJUST THE APPROXIMATE SUMS OF
C                                    SQUARES AND CROSS-PRODUCTS MATRIX
C                                    FOR THE TEMPORARY MEAN EFFECT
C                                    AND COMPUTE THE VARIANCE-
C                                    COVARIANCE MATRIX ON OPTION
      INX = 0
      DO 75 L=1,N1
         XM(L) = XM(L)*RNBR2
         DO 75 J=1,L
            INX = INX+1
            WK = DBLE(FNBR2)*DBLE(XM(J)-TEMP(J))*DBLE(XM(L)-TEMP(L))
            VCV(INX) = VCV(INX)-WK
            IF (NBR(6) .EQ. 0) VCV(INX) = VCV(INX)*RNBR3
   75 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BECOVM')
 9005 RETURN
      END

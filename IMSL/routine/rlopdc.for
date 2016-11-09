C   IMSL ROUTINE NAME   - RLOPDC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - RESPONSE PREDICTION USING AN ORTHOGONAL
C                           POLYNOMIAL REGRESSION MODEL
C
C   USAGE               - CALL RLOPDC (X,N,A,B,C,ID,IOPT,P,YHAT,IER)
C
C   ARGUMENTS    X      - INPUT AND OUTPUT VECTOR OF LENGTH N.
C                         ON INPUT, VECTOR OF N INDEPENDENT VARIABLE
C                           SETTINGS.
C                         ON OUTPUT, X IS SCALED BY THE LINEAR
C                           TRANSFORMATION SX=C(ID+2)*X+C(ID+3) IF
C                           SCALING WAS REQUESTED. OTHERWISE, X IS
C                           UNCHANGED.
C                N      - NUMBER OF INDEPENDENT VARIABLE SETTINGS.
C                A      - INPUT VECTOR OF LENGTH ID CONTAINING CONSTANTS
C                           USED IN GENERATING THE ORTHOGONAL
C                           POLYNOMIALS.
C                B      - INPUT VECTOR OF LENGTH ID CONTAINING
C                           ADDITIONAL CONSTANTS USED IN GENERATING THE
C                           ORTHOGONAL POLYNOMIALS.
C                C      - INPUT VECTOR OF LENGTH ID+3 CONTAINING THE
C                           INTERCEPT AND THE (ORTHOGONAL POLYNOMIAL)
C                           COEFFICIENTS OF THE FITTED POLYNOMIAL MODEL
C                           IN THE FIRST ID+1 LOCATIONS. THE LAST TWO
C                           LOCATIONS CONTAIN THE SCALING CONSTANTS USED
C                           IN THE LINEAR TRANSFORMATION
C                           SX = C(ID+2)*X+C(ID+3).
C                ID     - INPUT DEGREE OF THE POLYNOMIAL MODEL USED FOR
C                           PREDICTION. THE ROUTINE HAS NOT BEEN TESTED
C                           FOR DEGREES GREATER THAN 10.
C                IOPT   - INPUT SCALING OPTION. SEE REMARKS.
C                           IF IOPT =1, NO SCALING WILL BE PERFORMED.
C                           OTHERWISE, SCALING WILL BE PERFORMED.
C                P      - WORK VECTOR OF LENGTH 2*N LOCATIONS.
C                           P MUST BE TYPED DOUBLE PRECISION.
C                YHAT   - OUTPUT VECTOR OF LENGTH N CONTAINING
C                           THE PREDICTED RESPONSES FOR THE
C                           INPUT X VALUES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT ID IS GREATER THAN 10.
C                             RLOPDC HAS NOT BEEN TESTED FOR ACCURACY
C                             ABOVE DEGREE 10.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE INPUT X VECTOR WILL BE TRANSFORMED, IF NECESSARY,
C                TO THE SAME SCALE AS THE DATA USED IN FITTING THE
C                MODEL. IF THE MODEL WAS FITTED USING IMSL ROUTINES
C                RLFOTH OR RLFOTW AND THE PREDICTED RESPONSES ARE
C                DESIRED FOR THE INDEPENDENT VARIABLE SETTINGS USED IN
C                THE FITTING PROCESS, THEN THE SCALED X VALUES ARE
C                AVAILABLE FROM RLFOTH OR RLFOTW AND NO SCALING IS
C                NEEDED IN RLOPDC. HOWEVER, TO PREDICT FOR VALUES OF THE
C                UNSCALED X VARIABLE, SCALING IS REQUIRED AND THE NEEDED
C                SCALING CONSTANTS ARE AVAILABLE FROM RLFOTH OR RLFOTW
C                IN LOCATIONS C(ID+2) AND C(ID+3).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLOPDC (X,N,A,B,C,ID,IOPT,P,YHAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,ID,IOPT,IER
      REAL               X(1),A(1),B(1),C(1),YHAT(1)
      DOUBLE PRECISION   P(N,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IXC,J,LAST,NEXT
      DOUBLE PRECISION   Q,W,Z,Y
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
C                                  WARNING ERROR
      IF (ID .GT. 10) IER=33
      IF (IOPT .EQ. 1) GO TO 10
C                                  SCALE X
      Q=C(ID+2)
      Y=C(ID+3)
      DO 5 I=1,N
    5    X(I)=Q*X(I)+Y
   10 LAST=1
      NEXT=2
      Z=A(1)
      Q=C(1)
      Y=C(2)
C                                  COMPUTE THE ORTHOGONAL POLYNOMIALS
C                                  AND THE PREDICTED RESPONSES
      DO 15 I=1,N
         P(I,1)=1.D0
         W=X(I)-Z
         YHAT(I)=Q+Y*W
         P(I,2)=W
   15 CONTINUE
      IF (ID .LT. 2) GO TO 9005
      DO 25 I=2,ID
         Z=A(I)
         Q=B(I)
         Y=C(I+1)
         DO 20 J=1,N
            W=(X(J)-Z)*P(J,NEXT)-Q*P(J,LAST)
            YHAT(J)=YHAT(J)+Y*W
            P(J,LAST)=W
   20    CONTINUE
         IXC=NEXT
         NEXT=LAST
         LAST=IXC
   25 CONTINUE
      IF (IER) 9000,9005,9000
 9000 CONTINUE
      CALL UERTST (IER,6HRLOPDC)
 9005 RETURN
      END

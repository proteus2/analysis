C   IMSL ROUTINE NAME   - RLPOL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - GENERATE ORTHOGONAL POLYNOMIALS WITH THE
C                           ASSOCIATED CONSTANTS AA AND BB
C
C   USAGE               - CALL RLPOL (X,N,ID,SM,SA,AA,BB,P,IP,IER)
C
C   ARGUMENTS    X      - DATA VECTOR OF LENGTH N. (INPUT)
C                         ON OUTPUT, X IS SCALED TO THE INTERVAL (-2,2)
C                           BY THE LINEAR TRANSFORMATION XS=SM*X+SA.
C                N      - NUMBER OF DATA POINTS. (INPUT)
C                ID     - DEGREE OF THE HIGHEST DEGREE ORTHOGONAL
C                           POLYNOMIAL TO BE GENERATED. (INPUT)
C                SM     - MULTIPLICATIVE CONSTANT USED FOR SCALING
C                           X TO (-2,2). (OUTPUT)
C                SA     - ADDITIVE CONSTANT USED FOR SCALING X TO (-2,2)
C                           (OUTPUT)
C                AA     - VECTOR OF LENGTH ID CONTAINING THE
C                           CONSTANTS USED IN GENERATING THE ORTHOGONAL
C                           POLYNOMIALS. (OUTPUT)
C                BB     - VECTOR OF LENGTH ID CONTAINING
C                           ADDITIONAL CONSTANTS USED IN GENERATING THE
C                           ORTHOGONAL POLYNOMIALS. (OUTPUT)
C                P      - N BY ID MATRIX CONTAINING THE ID ORTHOGONAL
C                           POLYNOMIALS. (OUTPUT)
C                IP     - ROW DIMENSION OF MATRIX P EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE X VECTOR IS
C                             CONSTANT.
C                           IER = 130 INDICATES THAT THE NUMBER OF
C                             POINTS IS INSUFFICIENT FOR THE DEGREE
C                             SPECIFIED (THAT IS, N IS LESS THAN OR
C                             EQUAL TO ID) OR THAT ID WAS SPECIFIED
C                             LESS THAN 1.
C                         WARNING ERROR
C                           IER = 35 INDICATES THAT ID WAS SPECIFIED
C                             GREATER THAN 10. THE ACCURACY OF THE
C                             RESULTS MAY BE QUESTIONABLE.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE NUMBER OF DISTINCT VALUES IN THE X VECTOR
C                MUST BE GREATER THAN OR EQUAL TO ID+1.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLPOL  (X,N,ID,SM,SA,AA,BB,P,IP,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,ID,IP,IER
      REAL               X(N),P(IP,ID),AA(ID),BB(ID),SM,SA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,KK
      REAL               XN,X1,ZERO
      DOUBLE PRECISION   A,B,DZERO,FOUR,HALFM,SX,T(2),TEMP,TEMP1,TEMP2
      DATA               ZERO/0.0/
      DATA               DZERO/0.0D0/,FOUR/4.0D0/,HALFM/-0.5D0/
C                                  FIRST EXECUTABLE STATEMENT
      BB(1) = ZERO
      IER = 0
      IF (ID .LT. N .AND. ID .GE. 1) GO TO 10
C                                  TERMINAL ERROR - DEGREE IS LESS THAN
C                                  1 OR THERE ARE NOT ENOUGH DATA POINTS
C                                  FOR DESIRED DEGREE
      IER=130
      GO TO 9000
C                                  WARNING ID IS GREATER THAN 10
   10 IF (ID .GT. 10) IER=35
      XN = X(1)
      X1 = X(1)
C                                  FIND MAXIMUM AND MINIMUM X
      DO 15 I = 2,N
         IF (X(I) .GT. XN) XN=X(I)
         IF (X(I) .LT. X1) X1=X(I)
   15 CONTINUE
C                                  TERMINAL ERROR X VECTOR IS CONSTANT
      IF (X1 .NE. XN) GO TO 20
      IER=129
      GO TO 9000
C                                  SCALE X TO (-2,2)
   20 A=FOUR/(XN-X1)
      B=HALFM*(X1+XN)*A
      SM=A
      SA=B
      TEMP=DZERO
      DO 25 I=1,N
         SX=A*X(I)+B
         TEMP=TEMP+SX
         X(I)=SX
   25 CONTINUE
      A=TEMP/N
      AA(1) = A
      DO 30 I=1,N
         P(I,1)=X(I)-A
   30 CONTINUE
      IF (ID .LT. 2) GO TO 9005
      K=1
      TEMP2=N
      DO 70 J=2,ID
C                                  COMPUTE ALPHA AND BETA
         TEMP=DZERO
         DO 35 I=1,N
C                                  T CONTAINS P(I,K)**2
C                                  TOT CONTAINS THE SUMMATION OF
C                                  X(I)*P(I,K)**2
            T(1)=DBLE(P(I,K))
            TEMP=TEMP+DBLE(X(I))*T(1)*T(1)
   35    CONTINUE
         TEMP1=DZERO
         DO 40 I=1,N
            T(1)=DBLE(P(I,K))
            TEMP1=TEMP1+T(1)*T(1)
   40    CONTINUE
         A=TEMP/TEMP1
         B=TEMP1/TEMP2
         AA(J) = A
         BB(J) = B
   45    TEMP2=TEMP1
         IF (J .EQ. 2) GO TO 55
C                                  COMPUTE NEXT ORTHOGONAL POLYNOMIAL
         DO 50 I=1,N
            P(I,J)=(X(I)-A)*P(I,K)-B*P(I,KK)
   50    CONTINUE
         GO TO 65
   55    DO 60 I=1,N
            P(I,J)=(X(I)-A)*P(I,K)-B
   60    CONTINUE
   65    KK=K
         K=J
   70 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HRLPOL )
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - RLFOTH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - FIT A UNIVARIATE CURVILINEAR REGRESSION MODEL
C                           USING ORTHOGONAL POLYNOMIALS
C
C   USAGE               - CALL RLFOTH (X,Y,N,RSQ,MD,ID,P,C,S,A,B,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           INDEPENDENT VARIABLE SETTINGS. THE X VECTOR
C                           MUST NOT BE CONSTANT.
C                         ON OUTPUT, THE ELEMENTS OF X ARE SCALED TO
C                           THE INTERVAL (-2,2), INCLUSIVELY.
C                Y      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           OBSERVED RESPONSES. THE Y VECTOR MUST NOT
C                           BE CONSTANT.
C                N      - INPUT NUMBER OF OBSERVATIONS ON X AND Y.
C                RSQ    - INPUT VALUE IN THE INTERVAL (0,100) USED TO
C                           CONTROL THE DEGREE OF THE FITTED MODEL.
C                           FITTING STOPS WHEN THE COEFFICIENT OF
C                           DETERMINATION EXCEEDS RSQ, UNLESS DEGREE MD
C                           HAS ALREADY BEEN REACHED.
C                MD     - INPUT MAXIMUM DEGREE ALLOWED FOR THE FITTED
C                           MODEL. MD MUST BE GREATER THAN OR
C                           EQUAL TO 1 AND LESS THAN N.
C                ID     - OUTPUT DEGREE OF THE MODEL.
C                P      - WORK VECTOR OF LENGTH 2*N. P MUST BE TYPED
C                           DOUBLE PRECISION IN THE CALLING PROGRAM.
C                C      - OUTPUT VECTOR OF LENGTH MD+3 WHOSE FIRST ID+1
C                           COMPONENTS CONTAIN THE REGRESSION
C                           COEFFICIENTS OF THE FITTED MODEL IN
C                           ASCENDING DEGREE ORDER. THE NEXT TWO
C                           COMPONENTS CONTAIN THE SCALING CONSTANTS D
C                           AND E USED IN THE TRANSFORMATION SX=D*X+E.
C                S      - OUTPUT VECTOR OF LENGTH MD+3 WHOSE FIRST TWO
C                           COMPONENTS CONTAIN THE UNCORRECTED TOTAL
C                           SUM OF SQUARES AND THE SUM OF SQUARES
C                           ATTRIBUTABLE TO THE MEAN, RESPECTIVELY. THE
C                           NEXT ID COMPONENTS CONTAIN IN ASCENDING
C                           DEGREE ORDER, THE SUMS OF SQUARES
C                           ATTRIBUTABLE TO ORTHOGONAL POLYNOMIALS OF
C                           EACH DEGREE LESS THAN OR EQUAL TO ID. THE
C                           LAST COMPONENT CONTAINS THE ERROR SUM OF
C                           SQUARES.
C                A      - OUTPUT VECTOR OF LENGTH MD WHOSE FIRST ID
C                           COMPONENTS CONTAIN CONSTANTS USED IN
C                           GENERATING THE ORTHOGONAL POLYNOMIALS.
C                B      - OUTPUT VECTOR OF LENGTH MD WHOSE FIRST ID
C                           COMPONENTS CONTAIN ADDITIONAL CONSTANTS
C                           USED IN GENERATING THE ORTHOGONAL
C                           POLYNOMIALS.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT THE DEGREE OF THE
C                             MODEL WAS SPECIFIED LESS THAN 1 OR
C                             THAT NOT ENOUGH DATA POINTS WERE GIVEN
C                             FOR THE DEGREE REQUESTED.
C                           IER=130 INDICATES THAT THE X VECTOR OR THE
C                             Y VECTOR WAS CONSTANT.
C                           IER=131 INDICATES THAT RSQ WAS SPECIFIED
C                             LESS THAN OR EQUAL TO ZERO.
C                         WARNING WITH FIX ERROR
C                           IER=68 INDICATES THAT RSQ WAS SPECIFIED
C                             INCORRECTLY. IF RSQ WAS SPECIFIED
C                             GREATER THAN 100., RSQ IS SET TO 100.
C                         WARNING ERROR
C                           IER=37 INDICATES THAT THE MAXIMUM DEGREE
C                             SPECIFIED FOR THE FITTED POLYNOMIAL
C                             WAS GREATER THAN 10. RLFOTH HAS NOT BEEN
C                             TESTED FOR ACCURACY ABOVE DEGREE 10.
C                           IER=38 INDICATES THAT RSQ WAS 100.,
C                             BUT THAT A PERFECT FIT WAS OBTAINED
C                             WITH POLYNOMIAL OF LOWER DEGREE THAN MD.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF IT IS DESIRED TO FIT A MODEL OF A SPECIFIED
C                DEGREE, THEN SET RSQ EQUAL TO 100.0. AND MD EQUAL TO
C                THE DESIRED DEGREE. ALTERNATIVELY, THE DEGREE IS
C                CONTROLLED BY RSQ, WITH RSQ IN THE INTERVAL (1,100),
C                IF MD IS SUFFICIENTLY LARGE.
C            2.  IN ORDER TO FIT A MODEL OF DEGREE ID, THE NUMBER OF
C                DISTINCT SETTINGS OF THE INDEPENDENT VARIABLE MUST
C                BE GREATER THAN OR EQUAL TO ID+1. RLFOTH DOES NOT
C                CHECK FOR FAILURE OF THIS REQUIREMENT. IT DOES
C                CHECK TO SEE THAT N IS NOT LESS THAN ID+1.
C            3.  VECTOR C CONTAINS THE FITTED COEFFICIENTS. HOWEVER,
C                THESE COEFFICIENTS ARE NOT THOSE FOR THE VARIOUS
C                POWERS OF THE X VARIABLE. THEY ARE COEFFICIENTS OF
C                MORE COMPLICATED FUNCTIONS (THE ORTHOGONAL
C                POLYNOMIALS) OF X. THE USER MAY DECODE VECTOR C
C                BY CALLING IMSL ROUTINE RLDOPM. THE RESULTANT
C                COEFFICIENTS WILL BE COEFFICIENTS OF POWERS OF X.
C                HOWEVER, TO COMPUTE THE PREDICTED RESPONSE FOR A
C                GIVEN VALUE OF X WITHOUT DECODING, THE USER MAY
C                CALL IMSL ROUTINE RLOPDC.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLFOTH (X,Y,N,RSQ,MD,ID,P,C,S,A,B,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MD,ID,IER
      REAL               X(N),Y(N),RSQ,C(1),S(1),A(1),B(1)
      DOUBLE PRECISION   P(N,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,ISW,I1,I2,J,K,MD1
      REAL               R,ST,ZERO,HUND,P01,XX1,XXN
      DOUBLE PRECISION   SY,SY2,X1,XN,D,SN,YY,YN,SX
      DATA               ZERO/0.0/,HUND/100.0/,P01/.01/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF (MD .LT. N .AND. MD .GE. 1) GO TO 5
C                                  TERMINAL ERROR - DEGREE IS LESS THAN
C                                  1 OR THERE ARE NOT ENOUGH DATA POINTS
C                                  TO FIT DESIRED DEGREE
      IER=129
      GO TO 9000
C                                  RSQ WAS SPECIFIED INCORRECTLY
    5 IF (RSQ .GT. ZERO) GO TO 10
      IER=131
      GO TO 9000
   10 IF (RSQ .LE. HUND) GO TO 15
      IER=68
      RSQ=HUND
C                                  WARNING - MD IS OUTSIDE THE RANGE
C                                  TESTED
   15 IF (MD .GT. 10) IER=37
C                                  FIND MAXIMUM AND MINIMUM X
      XXN = X(1)
      XX1 = X(1)
      ISW = 1
      DO 20 I = 2,N
         IF (Y(I-1) .EQ. Y(I)) ISW = ISW+1
         IF (X(I) .GT. XXN) XXN=X(I)
         IF (X(I) .LT. XX1) XX1=X(I)
   20 CONTINUE
      X1 = XX1
      XN = XXN
C                                  TERMINAL ERROR - X VECTOR IS CONSTANT
      IF (X1 .NE. XN .AND. ISW .NE. N) GO TO 25
      IER=130
      GO TO 9000
C                                  SCALE X TO (-2,2)
   25 YY=(X1+XN)*.5D0
      YN=4.D0/(XN-X1)
      I1=1
      I2=2
      D=0.D0
      DO 30 I = 1,N
         P(I,I1)=1.D0
         SX=YN*(X(I)-YY)
         D=D+SX
         X(I) = SX
   30 CONTINUE
      D=D/N
      A(1)=D
      B(1)=ZERO
      SY=0.D0
      DO 35 I = 1,N
         SY=SY+Y(I)
         P(I,I2) = X(I)-D
   35 CONTINUE
      XN=0.D0
      DO 40 I = 1,N
         XN=XN+P(I,I2)*P(I,I2)
   40 CONTINUE
      SY2=0.D0
      DO 45 I = 1,N
         SY2=SY2+DBLE(Y(I))*DBLE(Y(I))
   45 CONTINUE
      S(1)=SY2
      R=RSQ*P01
      ST=ZERO
      S(2)=SY*SY/N
      SY2=SY2-S(2)
C                                  FIND FIRST COEFFICIENT
      D=N
      C(1)=SY/N
      MD1=MD+1
      ID=1
      SY=0.D0
      DO 50 I = 1,N
         SY=SY+P(I,I2)*Y(I)
   50 CONTINUE
      C(2)=SY/XN
      IF (MD1 .LT. 2) GO TO 80
      DO 75 I = 2,MD1
         J=I+1
         S(J)=SY*SY/XN
         ST=ST+S(J)
         SY2=SY2-S(J)
         IF (SY2 .LE. 0.D0) SY2=0.D0
         IF (ST.GT.R*(S(1)-S(2))) GO TO 80
         IF (I .EQ. MD1) GO TO 80
         ID=I
C                                  COMPUTE CONSTANTS
         SN = 0.D0
         DO 55 K = 1,N
            SN=SN+X(K)*P(K,I2)*P(K,I2)
   55    CONTINUE
         SN=SN/XN
         X1=XN/D
C                                  FIND NEXT ORTHOGONAL POLYNOMIAL
         DO 60 K = 1,N
            P(K,I1) = (X(K)-SN)*P(K,I2)-X1*P(K,I1)
   60    CONTINUE
         A(I)=SN
         B(I)=X1
         D=XN
         SY=0.D0
         DO 65 K = 1,N
            SY=SY+P(K,I1)*Y(K)
   65    CONTINUE
         XN=0.D0
         DO 70 K=1,N
            XN=XN+P(K,I1)*P(K,I1)
   70    CONTINUE
C                                  COMPUTE NEW COEFFICIENTS
         C(J)=SY/XN
         II=I1
         I1=I2
         I2=II
   75 CONTINUE
   80 S(MD+3)=SY2
      C(ID+2)=YN
      C(ID+3)=-YN*YY
      IF (RSQ.EQ.100..AND.ID.LT.MD) IER=38
      IF (IER) 9000,9005,9000
 9000 CONTINUE
      CALL UERTST (IER,6HRLFOTH)
 9005 RETURN
      END

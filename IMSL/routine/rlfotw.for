C   IMSL ROUTINE NAME   - RLFOTW
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - FIT A UNIVARIATE CURVILINEAR REGRESSION MODEL
C                           USING ORTHOGONAL POLYNOMIALS WITH WEIGHTING
C
C   USAGE               - CALL RLFOTW (X,Y,N,RSQ,MD,W,ID,P,C,S,A,B,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           INDEPENDENT VARIABLE SETTINGS. THE X VECTOR
C                           MUST NOT BE CONSTANT.
C                         ON OUTPUT, THE ELEMENTS OF X ARE SCALED TO
C                           THE INTERVAL (-2,2).
C                Y      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           OBSERVED RESPONSES.
C                N      - INPUT NUMBER OF OBSERVATIONS ON X AND Y.
C                RSQ    - INPUT VALUE IN THE INTERVAL (0,100) USED TO
C                           CONTROL THE DEGREE OF THE FITTED MODEL.
C                           FITTING STOPS WHEN THE COEFFICIENT OF
C                           DETERMINATION EXCEEDS RSQ, UNLESS DEGREE MD
C                           HAS ALREADY BEEN REACHED.
C                MD     - INPUT MAXIMUM DEGREE ALLOWED FOR THE FITTED
C                           MODEL. MD MUST BE GREATER THAN OR
C                           EQUAL TO 1 AND LESS THAN N.
C                W      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           WEIGHTS OF THE DATA POINTS.
C                         ON OUTPUT, W CONTAINS THE SQUARE ROOTS OF
C                           THE WEIGHTS.
C                ID     - OUTPUT DEGREE OF THE MODEL.
C                P      - WORK VECTOR OF LENGTH 4*N. P MUST BE TYPED
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
C                           IER=130 INDICATES THAT THE X VECTOR WAS
C                             CONSTANT.
C                           IER=131 INDICATES THAT RSQ WAS SPECIFIED
C                             LESS THAN OR EQUAL TO ZERO.
C                         WARNING WITH FIX ERROR
C                           IER=68 INDICATES THAT RSQ WAS SPECIFIED
C                             INCORRECTLY. IF RSQ WAS SPECIFIED
C                             GREATER THAN 100., RSQ IS SET TO 100.
C                         WARNING ERROR
C                           IER=36 INDICATES THAT THE Y VECTOR WAS
C                             CONSTANT.
C                           IER=37 INDICATES THAT THE DEGREE
C                             SPECIFIED WAS GREATER THAN 10. RLFOTW
C                             HAS NOT BEEN TESTED FOR ACCURACY ABOVE
C                             DEGREE 10.
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
C                BE GREATER THAN OR EQUAL TO ID+1. RLFOTW DOES NOT
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
      SUBROUTINE RLFOTW (X,Y,N,RSQ,MD,W,ID,P,C,S,A,B,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MD,ID,IER
      REAL               X(1),Y(1),RSQ,C(1),S(1),A(1),B(1),W(1)
      DOUBLE PRECISION   P(N,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,ISW,I1,I2,I3,I4,J,K,MD1
      REAL               R,ST,SU,ZERO,XX1,XXN,P01,HUND
      DOUBLE PRECISION   SY,SY2,X1,XN,D,SN,YY,YN,SX,SJ
      DATA               ZERO /0.0/,HUND /100./,P01 /.01/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (MD.LT.N .AND. MD.GE.1) GO TO 5
C                                  TERMINAL ERROR - DEGREE IS LESS THAN
C                                  1 OR THERE ARE NOT ENOUGH DATA POINTS
C                                  TO FIT DESIRED DEGREE
      IER = 129
      GO TO 9000
C                                  RSQ WAS SPECIFIED INCORRECTLY
    5 IF (RSQ.GT.ZERO) GO TO 10
      IER = 131
      GO TO 9000
   10 IF (RSQ.LE.HUND) GO TO 15
      IER = 68
      RSQ = HUND
C                                  WARNING - MD IS OUTSIDE THE RANGE
C                                  TESTED
   15 IF (MD.GT.10) IER = 37
C                                  FIND MAXIMUM AND MINIMUM X
      XXN = X(1)
      XX1 = X(1)
      ISW = 1
      DO 20 I=2,N
         IF (Y(I-1).EQ.Y(I)) ISW = ISW+1
         IF (X(I).GT.XXN) XXN = X(I)
         IF (X(I).LT.XX1) XX1 = X(I)
   20 CONTINUE
      X1 = XX1
      XN = XXN
C                                  WARNING ERROR - Y IS A CONSTANT
      IF (ISW.EQ.N) IER = 36
C                                  TERMINAL ERROR - X IS A CONSTANT
      IF (X1.NE.XN) GO TO 25
      IER = 130
      GO TO 9000
C                                  SCALE X TO (-2,2)
   25 YY = .5D0*(X1+XN)
      YN = 4.D0/(XN-X1)
      I1 = 1
      I2 = 2
      I3 = 3
      I4 = 4
      SY = 0.D0
      DO 30 I=1,N
         SY = SY+DBLE(Y(I))*DBLE(W(I))
   30 CONTINUE
      SY2 = 0.D0
      DO 35 I=1,N
         SY2 = SY2+DBLE(Y(I))*DBLE(Y(I))*DBLE(W(I))
   35 CONTINUE
      XN = 0.D0
      DO 40 I=1,N
         XN = XN+W(I)
         P(I,I1) = 0.D0
         P(I,I3) = 0.D0
         P(I,I2) = 1.D0
         SX = YN*(X(I)-YY)
         X(I) = SX
   40 CONTINUE
      DO 45 I=1,N
         W(I) = SQRT(W(I))
         P(I,I4) = W(I)
   45 CONTINUE
      X1 = 0.D0
      B(1) = ZERO
      S(1) = SY2
      R = RSQ*P01
      ST = ZERO
C                                  FIND FIRST COEFFICIENT
      D = XN
      C(1) = SY/XN
      MD1 = MD+1
      ID = 0
      DO 80 I=1,MD1
         J = I+1
         SJ = SY*SY/XN
         S(J) = SJ
         SY2 = SY2-SJ
         IF (SY2.LE.0.D0) SY2 = 0.D0
         IF (J.GT.2) GO TO 50
         SU = SY2
         IF (ISW.EQ.N .OR. SU.LE.ZERO) GO TO 85
         GO TO 55
   50    ST = ST+S(J)
         IF (ST.GT.SU*R) GO TO 85
   55    IF (I.EQ.MD1) GO TO 85
         ID = I
C                                  COMPUTE CONSTANTS
         SN = 0.D0
         DO 60 K=1,N
            SN = SN+X(K)*P(K,I4)*P(K,I4)
   60    CONTINUE
         SN = SN/XN
         IF (I.NE.1) X1 = XN/D
C                                  FIND NEXT ORTHOGONAL POLYNOMIAL
         DO 65 K=1,N
            P(K,I1) = (X(K)-SN)*P(K,I2)-X1*P(K,I1)
            P(K,I3) = P(K,I1)*W(K)
   65    CONTINUE
         A(I) = SN
         B(I) = X1
         D = XN
         SY = 0.D0
         DO 70 K=1,N
            SY = SY+P(K,I3)*Y(K)*W(K)
   70    CONTINUE
         XN = 0.D0
         DO 75 K=1,N
            XN = XN+P(K,I3)*P(K,I3)
   75    CONTINUE
C                                  COMPUTE NEW COEFFICIENTS
         C(J) = SY/XN
         II = I1
         I1 = I2
         I2 = II
         II = I3
         I3 = I4
         I4 = II
   80 CONTINUE
   85 S(MD+3) = SY2
      C(ID+2) = YN
      C(ID+3) = -YN*YY
      IF (RSQ.EQ.100..AND.ID.LT.MD) IER=38
      IF (IER) 9000, 9005, 9000
 9000 CONTINUE
      CALL UERTST(IER,6HRLFOTW)
 9005 RETURN
      END

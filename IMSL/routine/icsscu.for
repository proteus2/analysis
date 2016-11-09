C   IMSL ROUTINE NAME   - ICSSCU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CUBIC SPLINE DATA SMOOTHER
C
C   USAGE               - CALL ICSSCU (X,F,DF,NX,SM,Y,C,IC,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),F(I)) I=1,...,
C                           NX. (INPUT) X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                F      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS.
C                           (INPUT)
C                DF     - VECTOR OF LENGTH NX. (INPUT)
C                           DF(I) IS THE RELATIVE WEIGHT OF DATA
C                           POINT I (SEE PARAMETER SM BELOW).
C                NX     - NUMBER OF ELEMENTS IN X, F, DF, AND Y. (INPUT)
C                           NX MUST BE .GE. 2.
C                SM     - A NON-NEGATIVE NUMBER WHICH CONTROLS THE
C                           EXTENT OF SMOOTHING. (INPUT) THE SPLINE
C                           FUNCTION S IS DETERMINED SUCH THAT THE
C                           SUM FROM 1 TO NX OF
C                           ((S(X(I))-F(I))/DF(I))**2
C                           IS LESS THAN OR EQUAL TO SM,
C                           WHERE EQUALITY HOLDS UNLESS S DESCRIBES
C                           A STRAIGHT LINE.
C                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y IS A VECTOR
C                           OF LENGTH NX. C IS AN NX-1 BY 3 MATRIX.
C                           THE VALUE OF THE SPLINE APPROXIMATION
C                           AT T IS
C                           S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I) .LE. T .LT. X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                WK     - WORK AREA VECTOR OF LENGTH 7*NX+14.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN NX-1
C                           IER = 130, NX IS LESS THAN 2
C                           IER = 131, INPUT ABSCISSAE ARE NOT ORDERED
C                             SO THAT X(1) .LT. X(2) ... .LT. X(NX)
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE ROUTINE PRODUCES A NATURAL CUBIC SPLINE. HENCE,
C                THE SECOND DERIVATIVE OF THE SPLINE FUNCTION S AT
C                X(1) AND X(NX) IS ZERO.
C            2.  FOR EACH SET OF DATA POINTS THERE EXISTS A MAXIMUM
C                VALUE FOR THE SMOOTHING PARAMETER. LET US CALL THIS
C                                *
C                MAXIMUM VALUE SM . IT IS DEFINED BY THE FOLLOWING
C                FORMULA;
C                    *
C                  SM  = THE SUM FROM I EQUAL 1 TO NX OF
C                                           2
C                        ((Y(I)-F(I))/DF(I))
C                WHERE Y IS THE SET OF FUNCTION VALUES DEFINING THE
C                STRAIGHT LINE WHICH BEST APPROXIMATES THE DATA IN
C                THE LEAST SQUARES SENSE (WITH WEIGHTS DF).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSSCU (X,F,DF,NX,SM,Y,C,IC,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IER
      REAL               X(1),F(1),DF(1),Y(1),C(IC,1),WK(7,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,M2,NP1,NP3
      REAL               E,FF,F2,G,H,HMG,ONE,ONEDH,ONED3,P,SM,
     1                   TWOD3,ZERO
      DATA               TWOD3/.6666667/,
     1                   ONED3/.3333333/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CHECK ERROR CONDITIONS
      IF (IC .LT. NX-1) GO TO 65
      IF (NX .LT. 2) GO TO 70
C                                  SET UP WORKING AREAS
      M2 = NX+2
      NP1 = NX+1
      WK(1,1) = ZERO
      WK(1,2) = ZERO
      WK(2,NP1) = ZERO
      WK(3,M2) = ZERO
      WK(3,NP1) = ZERO
      WK(6,1) = ZERO
      WK(6,2) = ZERO
      WK(6,M2) = ZERO
      WK(6,NP1) = ZERO
      P = ZERO
      H = X(2)-X(1)
      IF (H .LE. ZERO) GO TO 75
      F2 = -SM
      FF = (F(2)-F(1))/H
      IF (NX .LT. 3) GO TO 30
      DO 5 I=3,NX
         G = H
         H = X(I)-X(I-1)
         IF (H .LE. ZERO) GO TO 75
         ONEDH = ONE/H
         E = FF
         FF = (F(I)-F(I-1))*ONEDH
         Y(I) = FF-E
         WK(4,I) = (G+H)*TWOD3
         WK(5,I) = H*ONED3
         WK(3,I) = DF(I-2)/G
         WK(1,I) = DF(I)*ONEDH
         WK(2,I) = -DF(I-1)/G-DF(I-1)*ONEDH
    5 CONTINUE
      DO 10 I=3,NX
         C(I-1,1) = WK(1,I)*WK(1,I)+WK(2,I)*WK(2,I)+WK(3,I)*WK(3,I)
         C(I-1,2) = WK(1,I)*WK(2,I+1)+WK(2,I)*WK(3,I+1)
         C(I-1,3) = WK(1,I)*WK(3,I+2)
   10 CONTINUE
C                                  NEXT ITERATION
   15 IF (NX .LT. 3) GO TO 30
      DO 20 I=3,NX
         WK(2,I-1) = FF*WK(1,I-1)
         WK(3,I-2) = G*WK(1,I-2)
         WK(1,I) = ONE/(P*C(I-1,1)+WK(4,I)-FF*WK(2,I-1)-G*WK(3,I-2))
         WK(6,I) = Y(I)-WK(2,I-1)*WK(6,I-1)-WK(3,I-2)*WK(6,I-2)
         FF = P*C(I-1,2)+WK(5,I)-H*WK(2,I-1)
         G = H
         H = C(I-1,3)*P
   20 CONTINUE
      NP3 = NX+3
      DO 25 I=3,NX
         J = NP3-I
         WK(6,J) = WK(1,J)*WK(6,J)-WK(2,J)*WK(6,J+1)-WK(3,J)*WK(6,J+2)
   25 CONTINUE
   30 E = ZERO
      H = ZERO
C                                  COMPUTE U AND ACCUMULATE E
      DO 35 I=2,NX
         G = H
         H = (WK(6,I+1)-WK(6,I))/(X(I)-X(I-1))
         HMG = H-G
         WK(7,I) = HMG*DF(I-1)*DF(I-1)
         E = E+WK(7,I)*HMG
   35 CONTINUE
      G = -H*DF(NX)*DF(NX)
      WK(7,NP1) = G
      E = E-G*H
      G = F2
      F2 = E*P*P
      IF (F2 .GE. SM .OR. F2 .LE. G) GO TO 50
      FF = ZERO
      H = (WK(7,3)-WK(7,2))/(X(2)-X(1))
      IF (NX .LT. 3) GO TO 45
      DO 40 I=3,NX
         G = H
         H = (WK(7,I+1)-WK(7,I))/(X(I)-X(I-1))
         G = H-G-WK(2,I-1)*WK(1,I-1)-WK(3,I-2)*WK(1,I-2)
         FF = FF+G*WK(1,I)*G
         WK(1,I) = G
   40 CONTINUE
   45 H = E-P*FF
      IF (H .LE. ZERO) GO TO 50
C                                  UPDATE THE LAGRANGE MULTIPLIER P
C                                     FOR THE NEXT ITERATION
      P = P+(SM-F2)/((SQRT(SM/E)+P)*H)
      GO TO 15
C                                  IF E LESS THAN OR EQUAL TO S,
C                                  COMPUTE THE COEFFICIENTS AND RETURN.
   50 NP1 = NX-1
      DO 55 I=1,NP1
         Y(I) = F(I)-P*WK(7,I+1)
         C(I,2) = WK(6,I+1)
         WK(1,I) = Y(I)
   55 CONTINUE
      WK(1,NX) = F(NX)-P*WK(7,NX+1)
      Y(NX) = WK(1,NX)
      DO 60 I=2,NX
         H = X(I)-X(I-1)
         C(I-1,3) = (WK(6,I+1)-C(I-1,2))/(H+H+H)
         C(I-1,1) = (WK(1,I)-Y(I-1))/H-(H*C(I-1,3)+C(I-1,2))*H
   60 CONTINUE
      GO TO 9005
   65 IER = 129
      GO TO 9000
   70 IER = 130
      GO TO 9000
   75 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,6HICSSCU)
 9005 RETURN
      END

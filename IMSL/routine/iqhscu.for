C   IMSL ROUTINE NAME   - IQHSCU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ONE-DIMENSIONAL QUASI-CUBIC HERMITE
C                           INTERPOLATION
C
C   USAGE               - CALL IQHSCU (X,Y,NX,C,IC,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),Y(I)) I=1,...,
C                           NX. (INPUT) X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS.
C                           (INPUT)
C                NX     - NUMBER OF ELEMENTS IN X AND Y. (INPUT) NX
C                           MUST BE .GE. 4.
C                C      - SPLINE COEFFICIENTS. (OUTPUT) C IS AN NX BY
C                           3 MATRIX. THE VALUE OF THE SPLINE
C                           APPROXIMATION AT T IS
C                           S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I) .LE. T .LT. X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN NX.
C                           IER = 130, NX IS LESS THAN 4.
C                           IER = 131, INPUT ABSCISSAE ARE NOT ORDERED
C                             SO THAT X(1) .LT. X(2) ... .LT. X(NX).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE COEFFICIENTS OF THE INTERPOLATING POLYNOMIAL
C                ARE CONTAINED IN THE FIRST NX-1 ROWS OF THE MATRIX C.
C                THE ONE REMAINING ROW OF THE MATRIX C IS USED AS A
C                WORK STORAGE AREA.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IQHSCU (X,Y,NX,C,IC,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IER
      REAL               X(NX),Y(NX),C(IC,3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,N0
      REAL               T1,T2,B,RM1,RM2,RM3,RM4,HALF,ONE,ZERO
      DATA               HALF/0.50/,ONE/1.0/,ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IC .GE. NX) GO TO 5
C                                  TERMINAL - IC LESS THAN NX
      IER = 129
      GO TO 9000
    5 IF (NX .GT. 3) GO TO 10
C                                  TERMINAL - TOO FEW DATA POINTS
      IER = 130
      GO TO 9000
   10 DO 15 I=2,NX
         IF (X(I) .GT. X(I-1)) GO TO 15
C                                  TERMINAL - X NOT MONOTONE INCREASING
         IER = 131
         GO TO 9000
   15 CONTINUE
      IER = 0
      RM3 = (Y(2)-Y(1))/(X(2)-X(1))
      T1 = RM3-(Y(2)-Y(3))/(X(2)-X(3))
      RM2 = RM3+T1
      RM1 = RM2+T1
C                                  NOW GET THE SLOPES
      N0 = NX-2
      DO 40 I=1,NX
         IF (I .GT. N0) GO TO 20
         RM4 = (Y(I+2)-Y(I+1))/(X(I+2)-X(I+1))
         GO TO 25
   20    RM4 = RM3-RM2+RM3
   25    T1 = ABS(RM4-RM3)
         T2 = ABS(RM2-RM1)
         B = T1+T2
         IF (B .NE. ZERO) GO TO 30
C                                  IF DENOMINATOR IS ZERO, GET AVERAGE
         C(I,1) = HALF*(RM2+RM3)
         GO TO 35
   30    C(I,1) = (T1*RM2+T2*RM3)/B
   35    RM1 = RM2
         RM2 = RM3
         RM3 = RM4
   40 CONTINUE
      N0 = NX-1
C                                  COMPUTE THE COEFFICIENTS FOR THE
C                                    NX-1 INTERVALS
      DO 45 I=1,N0
         T1 = ONE/(X(I+1)-X(I))
         T2 = (Y(I+1)-Y(I))*T1
         B = (C(I,1)+C(I+1,1)-T2-T2)*T1
         C(I,3) = B*T1
         C(I,2) = -B+(T2-C(I,1))*T1
   45 CONTINUE
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST(IER,6HIQHSCU)
 9005 RETURN
      END

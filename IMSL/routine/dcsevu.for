C   IMSL ROUTINE NAME   - DCSEVU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - CUBIC SPLINE FIRST AND SECOND DERIVATIVE
C                           EVALUATOR
C
C   USAGE               - CALL DCSEVU (X,Y,NX,C,IC,U,DS,M1,DDS,M2,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),Y(I)) I=1,...,
C                           NX. (INPUT) X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS.
C                           (INPUT)
C                NX     - NUMBER OF ELEMENTS IN X AND Y. (INPUT)
C                           NX MUST BE .GE. 2.
C                C      - SPLINE COEFFICIENTS. (INPUT) C IS AN NX-1 BY
C                           3 MATRIX.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                U      - VECTOR OF LENGTH MAX(M1,M2) CONTAINING
C                           THE ABSCISSAE OF THE POINTS AT WHICH THE
C                           FIRST DERIVATIVE AND/OR THE SECOND
C                           DERIVATIVE OF THE CUBIC SPLINE IS TO BE
C                           EVALUATED. (INPUT)
C                DS     - VECTOR OF LENGTH M1. (OUTPUT)
C                           THE VALUE OF THE FIRST DERIVATIVE OF
C                           THE SPLINE APPROXIMATION AT U(I) IS
C                           DS(I) = (3.0*C(J,3)*D+2.0*C(J,2))*D+C(J,1)
C                           WHERE X(J) .LE. U(I) .LT. X(J+1) AND
C                           D = U(I)-X(J), FOR I=1,...,M1.
C                M1     - NUMBER OF ELEMENTS IN DS. (INPUT)
C                DDS    - VECTOR OF LENGTH M2. (OUTPUT)
C                           THE VALUE OF THE SECOND DERIVATIVE OF
C                           THE SPLINE APPROXIMATION AT U(I) IS
C                           DDS(I) = 6.0*C(J,3)*D+2.0*C(J,2)
C                           WHERE X(J) .LE. U(I) .LT. X(J+1) AND
C                           D = U(I)-X(J), FOR I=1,...,M2.
C                M2     - NUMBER OF ELEMENTS IN DDS. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, U(I) IS LESS THAN X(1).
C                           IER = 34, U(I) IS GREATER THAN X(NX).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE CONDITIONS REQUIRED FOR INPUT ARGUMENTS X AND NX
C                ARE NOT CHECKED IN THE ROUTINE. IF ANY OF THOSE
C                CONDITIONS ARE NOT MET THE ROUTINE WILL NOT PERFORM
C                CORRECTLY.
C            2.  DCSEVU ALLOWS THE USER THE OPTION OF EVALUATING FIRST
C                DERIVATIVES AND/OR SECOND DERIVATIVES. IF SECOND
C                DERIVATIVES ARE NOT DESIRED, THE INPUT PARAMETER
C                M2 SHOULD BE SET TO ZERO AND THE DIMENSION OF DDS
C                SHOULD BE SET TO ONE. IF FIRST DERIVATIVES ARE
C                NOT DESIRED, THE INPUT PARAMETER M1 SHOULD BE SET
C                TO ZERO AND THE DIMENSION OF DS SHOULD BE SET TO ONE.
C            3.  THE ORDINATE Y(NX) IS NOT USED BY THE ROUTINE. FOR
C                U(K) GREATER THAN X(NX-1), THE VALUE OF THE FIRST
C                DERIVATIVE OF THE SPLINE, DS(K), IS GIVEN BY
C                  DS(K) = (3.0*C(NX-1,3)*D+2.0*C(NX-1,2))*D+C(NX-1,1)
C                AND THE VALUE OF THE SECOND DERIVATIVE OF THE SPLINE,
C                DDS(K), IS GIVEN BY
C                  DDS(K) = 6.0*C(NX-1,3)*D+2.0*C(NX-1,2)
C                WHERE D = U(K)-X(NX-1).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DCSEVU (X,Y,NX,C,IC,U,DS,M1,DDS,M2,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,M1,M2,IER
      REAL               X(NX),Y(NX),C(IC,3),U(1),DS(1),DDS(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JER,K,KER,MM,NXM1
      REAL               D,DD,SPP,THREE,ZERO
      DATA               I/1/,ZERO/0.0/,THREE/3.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  INITIALIZE ERROR PARAMETERS
      JER = 0
      KER = 0
      IF (M1 .LE. 0 .AND. M2 .LE. 0) GO TO 9005
      NXM1 = NX-1
      IF (I .GT. NXM1) I = 1
C                                  EVALUATE FIRST DERIVATIVE OF SPLINE
C                                  AT M1 POINTS AND EVALUATE SECOND
C                                  DERIVATIVE OF SPLINE AT M2 POINTS
      MM = MAX0(M1,M2)
      DO 40 K=1,MM
C                                  FIND THE PROPER INTERVAL
         D = U(K)-X(I)
         IF (D) 5,25,15
    5    IF (I .EQ. 1) GO TO 30
         I = I-1
         D = U(K)-X(I)
         IF (D) 5,25,20
   10    I = I+1
         D = DD
   15    IF (I .GE. NX) GO TO 35
         DD = U(K)-X(I+1)
         IF (DD .GE. ZERO) GO TO 10
         IF (D .EQ. ZERO) GO TO 25
C                                  PERFORM FIRST AND SECOND
C                                  DERIVATIVE EVALUATIONS
   20    SPP = THREE*C(I,3)*D+C(I,2)
         IF (K .LE. M1) DS(K) = (SPP+C(I,2))*D+C(I,1)
         IF (K .LE. M2) DDS(K) = SPP+SPP
         GO TO 40
   25    IF (K .LE. M1) DS(K) = C(I,1)
         IF (K .LE. M2) DDS(K) = C(I,2)+C(I,2)
         GO TO 40
C                                  WARNING - U(I) .LT. X(1)
   30    JER = 33
         GO TO 20
C                                  IF U(I) .GT. X(NX) - WARNING
   35    IF (DD .GT. ZERO) KER = 34
         D = U(K)-X(NXM1)
         I = NXM1
         GO TO 20
   40 CONTINUE
      IER = MAX0(JER,KER)
 9000 CONTINUE
      IF (JER .GT. 0) CALL UERTST(JER,6HDCSEVU)
      IF (KER .GT. 0) CALL UERTST(KER,6HDCSEVU)
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - ICSEVU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EVALUATION OF A CUBIC SPLINE
C
C   USAGE               - CALL ICSEVU(X,Y,NX,C,IC,U,S,M,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),Y(I)) I=1,...,
C                           NX (INPUT). X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS
C                           (INPUT).
C                NX     - NUMBER OF ELEMENTS IN X AND Y (INPUT).
C                           NX MUST BE .GE. 2.
C                C      - SPLINE COEFFICIENTS (INPUT). C IS AN NX-1 BY
C                           3 MATRIX.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM (INPUT).
C                           IC MUST BE .GE. NX-1.
C                U      - VECTOR OF LENGTH M CONTAINING THE ABSCISSAE
C                           OF THE M POINTS AT WHICH THE CUBIC SPLINE
C                           IS TO BE EVALUATED (INPUT).
C                S      - VECTOR OF LENGTH M (OUTPUT).
C                           THE VALUE OF THE SPLINE APPROXIMATION AT
C                           U(I) IS
C                           S(I) = ((C(J,3)*D+C(J,2))*D+C(J,1))*D+Y(J)
C                           WHERE X(J) .LE. U(I) .LT. X(J+1) AND
C                           D = U(I)-X(J).
C                M      - NUMBER OF ELEMENTS IN U AND S (INPUT).
C                IER    - ERROR PARAMETER (OUTPUT).
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
C   REMARKS  1.  THE ROUTINE ASSUMES THAT THE ABSCISSAE OF THE NX
C                DATA POINTS ARE ORDERED SUCH THAT X(I) IS LESS THAN
C                X(I+1) FOR I=1,...,NX-1. NO CHECK OF THIS CONDITION
C                IS MADE IN THE ROUTINE. UNORDERED ABSCISSAE WILL CAUSE
C                THE ALGORITHM TO PRODUCE INCORRECT RESULTS.
C            2.  THE ROUTINE GENERATES TWO WARNING ERRORS. ONE ERROR
C                OCCURS IF U(I) IS LESS THAN X(1), FOR SOME I IN THE
C                THE INTERVAL (1,M) INCLUSIVELY. THE OTHER ERROR OCCURS
C                IF U(I) IS GREATER THAN X(NX), FOR SOME I IN THE
C                INTERVAL (1,M) INCLUSIVELY.
C            3.  THE ORDINATE Y(NX) IS NOT USED BY THE ROUTINE. FOR
C                U(K) .GT. X(NX-1), THE VALUE OF THE SPLINE, S(K), IS
C                GIVEN BY
C                 S(K)=((C(NX-1,3)*D+C(NX-1,2))*D+C(NX-1,1))*D+Y(NX-1)
C                WHERE D=U(K)-X(NX-1).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSEVU  (X,Y,NX,C,IC,U,S,M,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,M,IER
      REAL               X(NX),Y(NX),C(IC,3),U(M),S(M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JER,KER,NXM1,K
      REAL               D,DD,ZERO
      DATA               I/1/,ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      JER = 0
      KER = 0
      IF (M .LE. 0) GO TO 9005
      NXM1 = NX-1
      IF (I .GT. NXM1) I = 1
C                                  EVALUATE SPLINE AT M POINTS
      DO 40 K=1,M
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
C                                  PERFORM EVALUATION
   20    S(K) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
         GO TO 40
   25    S(K) = Y(I)
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
      IF (JER .GT. 0) CALL UERTST(JER,'ICSEVU')
      IF (KER .GT. 0) CALL UERTST(KER,'ICSEVU')
 9005 RETURN
      END

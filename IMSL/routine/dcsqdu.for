C   IMSL ROUTINE NAME   - DCSQDU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - CUBIC SPLINE QUADRATURE
C
C   USAGE               - CALL DCSQDU (X,Y,NX,C,IC,A,B,Q,IER)
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
C                           3 MATRIX. THE VALUE OF THE SPLINE
C                           APPROXIMATION AT A POINT Z IS GIVEN BY
C                           S(Z) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I) .LE. Z .LT. X(I+1) AND
C                           D = Z-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                A,B    - LIMITS OF INTEGRATION. (INPUT)
C                Q      - INTEGRAL FROM A TO B. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, A AND/OR B IS LESS THAN X(1).
C                           IER = 34, A AND/OR B IS GREATER THAN X(NX).
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
C            2.  WHEN THE LIMITS OF INTEGRATION ARE OUTSIDE OF THE
C                INTERVAL (X(1),X(NX)), THE INTEGRATION IS PERFORMED
C                ON THE EXTENDED SPLINE FUNCTION. THE SPLINE IS
C                EXTENDED ON THE LEFT BY THE FOLLOWING FORMULA-
C                  S(Z) = ((C(1,3)*D+C(1,2))*D+C(1,1))*D+Y(1)
C                WHERE Z IS LESS THAN X(1) AND D = Z-X(1).
C                THE SPLINE IS EXTENDED ON THE RIGHT BY THE
C                FOLLOWING FORMULA-
C                 S(U) = ((C(NX-1,3)*D+C(NX-1,2))*D+C(NX-1,1))*D+Y(NX-1)
C                WHERE U IS GREATER THAN X(NX) AND D = U-X(NX-1).
C            3.  THE ORDINATE Y(NX) IS NOT USED BY THE ROUTINE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DCSQDU (X,Y,NX,C,IC,A,B,Q,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IER
      REAL               X(NX),Y(NX),C(IC,3),A,B,Q
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IA,IB,IBM1,IPT,IV,JER,KER,NXM1
      REAL               D,DA,DB,DD,DX,FOURTH,HALF,QA,QAB,QB,
     1                   SIXTH,THIRD,V,ZERO
      DATA               SIXTH/.1666667/,
     1                   THIRD/.3333333/
      DATA               ZERO/0.0/,FOURTH/.25/,HALF/.5/
C                                  FIRST EXECUTABLE STATEMENT
      JER = 0
      KER = 0
C                                  FIND THE INTERVAL FOR A AND B
      NXM1 = NX-1
      IPT = 1
      IA = 1
      V = AMIN1(A,B)
    5 D = V-X(IA)
      DO 10 I=IA,NXM1
         IV = I
         DD = V-X(I+1)
         IF (DD .LT. ZERO) GO TO 15
         IF (I .LT. NXM1) D = DD
   10 CONTINUE
      IV = NXM1
C                                  IF V .GT. X(NX) - WARNING
      IF (DD .GT. ZERO) KER = 34
   15 CONTINUE
C                                  CHECK FOR V .LT. X(1)
      IF (D .LT. ZERO) JER = 33
      IF (IPT .EQ. 2) GO TO 20
      IPT = 2
      IA = IV
      DA = D
      V = AMAX1(A,B)
      GO TO 5
   20 IB = IV
      DB = D
C                                  INTEGRATE FROM X(IA) TO MIN(A,B)
C
      QA = (((FOURTH*C(IA,3)*DA+THIRD*C(IA,2))*DA+HALF*C(IA,1))*DA
     *        +Y(IA))*DA
      QAB = ZERO
      IBM1 = IB-1
      IF (IBM1 .LT. IA) GO TO 30
C                                  INTEGRATE FROM X(IA) TO X(IB)
      DO 25 I=IA,IBM1
         DX = X(I+1)-X(I)
         QAB = QAB+HALF*DX*(Y(I+1)+Y(I)-(C(I+1,2)+C(I,2))*DX*DX*SIXTH)
   25 CONTINUE
C                                  INTEGRATE FROM X(IB) TO MAX(A,B)
C
   30 QB = (((FOURTH*C(IB,3)*DB+THIRD*C(IB,2))*DB+HALF*C(IB,1))*DB
     *        +Y(IB))*DB
C                                  Q = INTEGRAL FROM A TO B
      Q = QB+QAB-QA
C                                  IF B .LT. A, CHANGE SIGN OF Q
      IF (B .LT. A) Q = -Q
      IER = MAX0(JER,KER)
 9000 CONTINUE
      IF (JER .GT. 0) CALL UERTST(JER,6HDCSQDU)
      IF (KER .GT. 0) CALL UERTST(KER,6HDCSQDU)
 9005 RETURN
      END

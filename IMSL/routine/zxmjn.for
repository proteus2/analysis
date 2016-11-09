C   IMSL ROUTINE NAME   - ZXMJN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES ZXMIN AND
C                           ZXMWD
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXMJN(A,N,Z,SIG,W,IR,MK,EPS)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IR,MK
      REAL               A(1),Z(N),SIG,W(N),EPS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,JJ,IJ,JP1,I,II,MM
      REAL               ZERO,ONE,FOUR,TI,V,TIM,AL,R,B,GM,Y,RINF,SQRINF
      DATA               ZERO/0.0/,ONE/1.0/,FOUR/4.0/
      DATA               RINF/Z7FFFFFFF/
C                                  UPDATE FACTORS GIVEN IN A
C                                    SIG*Z*Z-TRANSPOSE IS ADDED
C                                  FIRST EXECUTABLE STATEMENT
      SQRINF = SQRT(RINF)
      IF (N.GT.1) GO TO 5
C                                  N .EQ. 1
      A(1) = A(1) + SIG*Z(1)*Z(1)
      IR = 1
      IF (A(1).GT.ZERO) GO TO 9005
      A(1) = ZERO
      IR = 0
      GO TO 9005
C                                  N .GT. 1
    5 IF (SIG.GT.ZERO) GO TO 65
      IF (SIG.EQ.ZERO .OR. IR.EQ.0) GO TO 9005
      TI = ONE/SIG
      JJ = 0
      IF (MK.EQ.0) GO TO 15
C                                  L*W = Z ON INPUT
      DO 10 J = 1, N
        JJ = JJ + J
        IF (A(JJ).NE.ZERO) TI = TI + (W(J)*W(J))/A(JJ)
   10 CONTINUE
      GO TO 40
C                                  SOLVE L*W = Z
   15 DO 20 J = 1, N
        W(J) = Z(J)
   20 CONTINUE
      DO 35 J = 1, N
        JJ = JJ + J
        V = W(J)
        IF (A(JJ).GT.ZERO) GO TO 25
        W(J) = ZERO
        GO TO 35
   25   TI = TI + (V*V)/A(JJ)
        IF (J.EQ.N) GO TO 35
        IJ = JJ
        JP1 = J + 1
        DO 30 I = JP1, N
          IJ = IJ + I - 1
          W(I) = W(I) - V*A(IJ)
   30   CONTINUE
   35 CONTINUE
C                                  SET TI, TIM AND W
   40 IF (IR.LE.0) GO TO 45
      IF (TI.GT.ZERO) GO TO 50
      IF (MK-1) 65, 65, 55
   45 TI = ZERO
      IR = -IR - 1
      GO TO 55
   50 TI = EPS/SIG
      IF (EPS.EQ.ZERO) IR = IR - 1
   55 TIM = TI
      II = JJ
      I = N
      DO 60 J = 1, N
        IF (A(II).NE.ZERO) TIM = TI - (W(I)*W(I))/A(II)
        W(I) = TI
        TI = TIM
        II = II - I
        I = I - 1
   60 CONTINUE
      MM = 1
      GO TO 70
   65 MM = 0
      TIM = ONE/SIG
   70 JJ = 0
C                                  UPDATE A
      DO 120 J = 1, N
        JJ = JJ + J
        IJ = JJ
        JP1 = J + 1
C                                  UPDATE A(J,J)
        V = Z(J)
        IF (A(JJ).GT.ZERO) GO TO 95
C                                  A(J,J) .EQ. ZERO
        IF (IR.GT.0 .OR. SIG.LT.ZERO .OR. V.EQ.ZERO) GO TO 90
        IR = 1 - IR
        IF (V.GE.SQRINF) GO TO 75
        A(JJ) = (V*V)/TIM
        GO TO 80
   75   A(JJ) = RINF/TIM
   80   IF (J.EQ.N) GO TO 9005
        DO 85 I = JP1, N
          IJ = IJ + I - 1
          A(IJ) = Z(I)/V
   85   CONTINUE
        GO TO 9005
   90   TI = TIM
        GO TO 120
C                                  A(J,J) .GT. ZERO
   95   AL = V/A(JJ)
        TI = W(J)
        IF (MM.EQ.0) TI = TIM + V*AL
        R = TI/TIM
        A(JJ) = R*A(JJ)
        IF (R.EQ.ZERO) GO TO 125
        IF (J.EQ.N) GO TO 125
C                                  UPDATE REMAINDER OF COLUMN J
        B = AL/TI
        IF (R.GT.FOUR) GO TO 105
        DO 100 I = JP1, N
          IJ = IJ + I - 1
          Z(I) = Z(I) - V*A(IJ)
          A(IJ) = A(IJ) + B*Z(I)
  100   CONTINUE
        GO TO 115
  105   GM = TIM/TI
        DO 110 I = JP1, N
          IJ = IJ + I - 1
          Y = A(IJ)
          A(IJ) = B*Z(I) + Y*GM
          Z(I) = Z(I) - V*Y
  110   CONTINUE
  115   TIM = TI
  120 CONTINUE
  125 IF (IR.LT.0) IR = -IR
 9005 CONTINUE
      RETURN
      END

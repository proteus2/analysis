C   IMSL ROUTINE NAME   - EREDU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE EIGZS
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EREDU (N,A,B,AL,BL,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      REAL               A(1),B(1),AL(1),BL(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,II,IJ,IK,I,J1,JK,J,KJ,K,IK0,JK0
      REAL               X,Y
C                                  FIRST EXECUTABLE STATEMENT
C                                  FORM L IN THE ARRAY BL
      DO 25 I=1,N
         I1 = I-1
         DO 20 J=I,N
            IJ = J*(J-1)/2+I
            X = B(IJ)
            IF (I.EQ.1) GO TO 10
            IK0 = I*(I-1)/2
            JK0 = J*(J-1)/2
            DO 5 K=1,I1
               IK = IK0+K
               JK = JK0+K
               X = X-BL(IK)*BL(JK)
    5       CONTINUE
   10       IF (J.NE.I) GO TO 15
C                                  TEST B FOR POSITIVE DEFINITENESS
            IER = 129
            IF (X.LE.0.0) GO TO 9000
            Y = SQRT(X)
            II = I*(I+1)/2
            BL(II) = Y
            GO TO 20
   15       BL(IJ) = X/Y
   20    CONTINUE
   25 CONTINUE
C                                  FORM THE TRANSPOSE OF THE UPPER
C                                    TRIANGLE OF INV(L)*A IN AL
      DO 45 I=1,N
         I1 = I-1
         II = I*(I+1)/2
         Y = BL(II)
         DO 40 J=I,N
            IJ = J*(J-1)/2+I
            X = A(IJ)
            IF (I.EQ.1) GO TO 35
            IK0 = I*(I-1)/2
            JK0 = J*(J-1)/2
            DO 30 K=1,I1
               IK = IK0+K
               JK = JK0+K
               X = X-BL(IK)*AL(JK)
   30       CONTINUE
   35       AL(IJ) = X/Y
   40    CONTINUE
   45 CONTINUE
C                                  PRE MULTIPLY BY INV(L)
      DO 75 J=1,N
         J1 = J-1
         DO 70 I=J,N
            IJ = I*(I-1)/2+J
            X = AL(IJ)
            IF (I.EQ.J) GO TO 55
            I1 = I-1
            IK0 = I*(I-1)/2
            DO 50 K=J,I1
               KJ = K*(K-1)/2+J
               IK = IK0+K
               X = X-AL(KJ)*BL(IK)
   50       CONTINUE
   55       IF (J.EQ.1) GO TO 65
            IK0 = I*(I-1)/2
            JK0 = J*(J-1)/2
            DO 60 K=1,J1
               JK = JK0+K
               IK = IK0+K
               X = X-AL(JK)*BL(IK)
   60       CONTINUE
   65       II = I*(I+1)/2
            AL(IJ) = X/BL(II)
   70    CONTINUE
   75 CONTINUE
      IER = 0
 9000 RETURN
      END

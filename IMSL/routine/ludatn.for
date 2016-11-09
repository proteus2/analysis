C   IMSL ROUTINE NAME   - LUDATN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE LEQT2F
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUDATN (A,IA,N,LU,ILU,IDGT,D1,D2,APVT,EQUIL,WA,IER)
C
      DIMENSION          A(IA,1),LU(ILU,1),APVT(1),EQUIL(1)
      REAL               LU
      DATA               ZERO,ONE,FOUR,SIXTN,SIXTH/0.0,1.,4.,16.,.0625/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZATION
      IER = 0
      RN = N
      WREL = ZERO
      D1 = ONE
      D2 = ZERO
      BIGA = ZERO
      DO 10 I=1,N
         BIG = ZERO
         DO 5 J=1,N
            P = A(I,J)
            LU(I,J) = P
            P = ABS(P)
            IF (P .GT. BIG) BIG = P
    5    CONTINUE
         IF (BIG .GT. BIGA) BIGA = BIG
         IF (BIG .EQ. ZERO) GO TO 110
         EQUIL(I) = ONE/BIG
   10 CONTINUE
      DO 105 J=1,N
         JM1 = J-1
         IF (JM1 .LT. 1) GO TO 40
C                                  COMPUTE U(I,J), I=1,...,J-1
         DO 35 I=1,JM1
            SUM = LU(I,J)
            IM1 = I-1
            IF (IDGT .EQ. 0) GO TO 25
C                                  WITH ACCURACY TEST
            AI = ABS(SUM)
            WI = ZERO
            IF (IM1 .LT. 1) GO TO 20
            DO 15 K=1,IM1
               T = LU(I,K)*LU(K,J)
               SUM = SUM-T
               WI = WI+ABS(T)
   15       CONTINUE
            LU(I,J) = SUM
   20       WI = WI+ABS(SUM)
            IF (AI .EQ. ZERO) AI = BIGA
            TEST = WI/AI
            IF (TEST .GT. WREL) WREL = TEST
            GO TO 35
C                                  WITHOUT ACCURACY
   25       IF (IM1 .LT. 1) GO TO 35
            DO 30 K=1,IM1
               SUM = SUM-LU(I,K)*LU(K,J)
   30       CONTINUE
            LU(I,J) = SUM
   35    CONTINUE
   40    P = ZERO
C                                  COMPUTE U(J,J) AND L(I,J), I=J+1,...,
         DO 70 I=J,N
            SUM = LU(I,J)
            IF (IDGT .EQ. 0) GO TO 55
C                                  WITH ACCURACY TEST
            AI = ABS(SUM)
            WI = ZERO
            IF (JM1 .LT. 1) GO TO 50
            DO 45 K=1,JM1
               T = LU(I,K)*LU(K,J)
               SUM = SUM-T
               WI = WI+ABS(T)
   45       CONTINUE
            LU(I,J) = SUM
   50       WI = WI+ABS(SUM)
            IF (AI .EQ. ZERO) AI = BIGA
            TEST = WI/AI
            IF (TEST .GT. WREL) WREL = TEST
            GO TO 65
C                                  WITHOUT ACCURACY TEST
   55       IF (JM1 .LT. 1) GO TO 65
            DO 60 K=1,JM1
               SUM = SUM-LU(I,K)*LU(K,J)
   60       CONTINUE
            LU(I,J) = SUM
   65       Q = EQUIL(I)*ABS(SUM)
            IF (P .GE. Q) GO TO 70
            P = Q
            IMAX = I
   70    CONTINUE
C                                  TEST FOR ALGORITHMIC SINGULARITY
         IF (RN+P .EQ. RN) GO TO 110
         IF (J .EQ. IMAX) GO TO 80
C                                  INTERCHANGE ROWS J AND IMAX
         D1 = -D1
         DO 75 K=1,N
            P = LU(IMAX,K)
            LU(IMAX,K) = LU(J,K)
            LU(J,K) = P
   75    CONTINUE
         EQUIL(IMAX) = EQUIL(J)
   80    APVT(J) = IMAX
         D1 = D1*LU(J,J)
   85    IF (ABS(D1) .LE. ONE) GO TO 90
         D1 = D1*SIXTH
         D2 = D2+FOUR
         GO TO 85
   90    IF (ABS(D1) .GE. SIXTH) GO TO 95
         D1 = D1*SIXTN
         D2 = D2-FOUR
         GO TO 90
   95    CONTINUE
         JP1 = J+1
         IF (JP1 .GT. N) GO TO 105
C                                  DIVIDE BY PIVOT ELEMENT U(J,J)
         P = LU(J,J)
         DO 100 I=JP1,N
            LU(I,J) = LU(I,J)/P
  100    CONTINUE
  105 CONTINUE
C                                  PERFORM ACCURACY TEST
      IF (IDGT .EQ. 0) GO TO 9005
      P = 3*N+3
      WA = P*WREL
      IF (WA+10.0**(-IDGT) .NE. WA) GO TO 9005
      IER = 34
      GO TO 9000
C                                  ALGORITHMIC SINGULARITY
  110 IER = 129
      D1 = ZERO
      D2 = ZERO
 9000 CONTINUE
C                                  PRINT ERROR
      CALL UERTST(IER,6HLUDATN)
 9005 RETURN
      END

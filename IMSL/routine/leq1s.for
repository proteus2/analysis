C   IMSL ROUTINE NAME   - LEQ1S
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - INDEFINITE
C                           MATRIX - SYMMETRIC STORAGE MODE - SPACE
C                           ECONOMIZER SOLUTION
C
C   USAGE               - CALL LEQ1S (A,N,B,M,IB,IJOB,ICHNG,DET,IER)
C
C   ARGUMENTS    A      - INPUT/OUTPUT VECTOR OF DIMENSION N*(N+1)/2.
C                           SEE PARAMETER IJOB.
C                N      - ORDER OF MATRIX A AND THE NUMBER OF ROWS IN
C                           B. (INPUT)
C                B      - INPUT/OUTPUT MATRIX OF DIMENSION N BY M.
C                           ON INPUT, B CONTAINS THE M RIGHT-HAND SIDES
C                           OF THE EQUATION AX = B. ON OUTPUT, THE
C                           SOLUTION MATRIX X REPLACES B. IF IJOB = 1,
C                           B IS NOT USED.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES WHEN
C                           I = 0, FACTOR THE MATRIX A AND SOLVE THE
C                             EQUATION AX = B. ON INPUT, A CONTAINS THE
C                             COEFFICIENT MATRIX OF THE EQUATION AX = B,
C                             WHERE A IS ASSUMED TO BE AN N BY N
C                             SYMMETRIC MATRIX. A IS STORED IN SYMMETRIC
C                             STORAGE MODE AND THEREFORE HAS DIMENSION
C                             N*(N+1)/2. ON OUTPUT, A IS REPLACED BY
C                             ITS FACTORIZED FORM.
C                           I = 1, FACTOR THE MATRIX A. A CONTAINS THE
C                             SAME INPUT/OUTPUT INFORMATION AS IF
C                             IJOB = 0. B IS NOT USED.
C                           I = 2, SOLVE THE EQUATION AX = B. THIS
C                             OPTION IMPLIES THAT LEQ1S HAS ALREADY
C                             BEEN CALLED USING IJOB = 0 OR 1 SO THAT
C                             THE MATRIX A HAS ALREADY BEEN FACTORED.
C                             IN THIS CASE, OUTPUT MATRIX A MUST HAVE
C                             BEEN SAVED FOR REUSE IN THE CALL TO
C                             LEQ1S.
C                ICHNG  - WORK AREA OF LENGTH 2N.
C                DET    - WORK AREA OF LENGTH N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQ1S  (A,N,B,M,IB,IJOB,ICHNG,DET,IER)
C
      DIMENSION          A(1),B(IB,1),ICHNG(1),DET(N)
      DATA               ZERO/0./,ONE/1./
      DATA               ALPHA/.6403882/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
C                                  DECOMPOSE A INTO THE PRODUCT
C                                  M*D*M-TRANSPOSE WHERE M IS UNIT LOWER
C                                  TRIANGULAR AND D IS BLOCK DIAGONAL
C                                  WITH BLOCKS OF ORDER 1 OR 2.
C                                  M(I+1,I) .EQ. 0 WHEN D(I+1,I) .NE. 0
C                                  M AND D ARE WRITTEN OVER A.
      N0 = N-1
      N1 = N+1
      IF (IJOB.EQ.2) GO TO 125
C                                  CALCULATE EQUILIBRATION FACTORS
      RN = 16.0*N
      LL = 1
      DO 10 I=1,N
         DET(I) = ZERO
         IC = 1
         L = LL
         DO 5 J=1,N
            TEMP = ABS(A(L))
            IF (TEMP.GT.DET(I)) DET(I) = TEMP
            IF (J.GE.I) IC = J
            L = L+IC
    5    CONTINUE
         DET(I) = DET(I)*RN
         LL = LL+I
   10 CONTINUE
C                                  CALCULATE THE MAXIMUM DIAGONAL OF A
      IDET = 1
      I = 1
      L = 0
   15 LL = L
      IX = L+I
      IXI = IX+I
      IX1 = IXI+1
      X0 = ZERO
      DO 20 J=I,N
         L = L+J
         X1 = ABS(A(L))
         IF (X1.LE.X0) GO TO 20
         X0 = X1
         K = J
   20 CONTINUE
      X1 = X0
      IR = K
      IS = K
      ICHNG(I) = K
      I1 = I+1
      IF (I.GT.N0) GO TO 75
C                                  CALCULATE MAXIMUM OFF DIAGONAL
C                                  ELEMENT IN SUBMATRIX
      J = 1
      LX1 = IXI
      DO 35 L=I1,N
         LX = LX1
         DO 30 KK=1,J
            TEMP = ABS(A(LX))
            IF (TEMP.LE.X0) GO TO 25
            X0 = TEMP
            IR = L
            IS = KK
   25       LX = LX+1
   30    CONTINUE
         LX1 = LX1+L
         J = J+1
   35 CONTINUE
      IF (IS.NE.K) IS = IS+I-1
      INX = LL
      II = I
      KK = K
      ASSIGN 75 TO IBACK
      IF (X1.GE.X0*ALPHA) GO TO 40
      ASSIGN 70 TO IBACK
      ICHNG(I) = IS
      KK = IS
C                                  INTERCHANGE ROWS II AND KK
   40 IF (KK.EQ.II) GO TO IBACK, (75,70,95)
      TEMP = DET(II)
      DET(II) = DET(KK)
      DET(KK) = TEMP
      K1 = KK+1
      JK = (K1*KK)/2
      KX = JK-KK
      IF (K1.GT.N) GO TO 50
      JI = JK+II
      JK = JK+KK
      DO 45 J=K1,N
         TEMP = A(JK)
         A(JK) = A(JI)
         A(JI) = TEMP
         JK = JK+J
         JI = JI+J
   45 CONTINUE
   50 K0 = KK-1
      II1 = II+1
      IF (II1.GT.K0) GO TO 60
      JK = KX+II1
      JI = INX+II+II
      DO 55 J=II1,K0
         TEMP = A(JI)
         A(JI) = A(JK)
         A(JK) = TEMP
         JK = JK+1
         JI = JI+J
   55 CONTINUE
   60 INXJ = INX+II
      KNXJ = KX+KK
      TEMP = A(INXJ)
      A(INXJ) = A(KNXJ)
      A(KNXJ) = TEMP
      IF (II.EQ.1) GO TO IBACK, (75,70,95)
      IM1 = II-1
      DO 65 J=1,IM1
         INXJ = INX+J
         KNXJ = KX+J
         TEMP = A(INXJ)
         A(INXJ) = A(KNXJ)
         A(KNXJ) = TEMP
   65 CONTINUE
      GO TO IBACK, (75,70,95)
   70 ICHNG(I1) = IR
      ASSIGN 95 TO IBACK
      KK = IR
      II = I1
      INX = INX+I
      GO TO 40
C                                  USE A 1 BY 1 PIVOT
   75 TEMP = DET(I)+A(IX)
      IF (TEMP.EQ.DET(I)) GO TO 200
      IF (I1.GT.N) GO TO 90
      JX = IX
      L = JX+I
      TEMP = ONE/A(IX)
      DO 85 J=I1,N
         INXJ = JX+I
         SAVE = A(INXJ)
         A(INXJ) = SAVE*TEMP
C                                  A(J,I) HAS BEEN SET EQUAL TO THE
C                                  MULTIPLIER
         KX = L
         DO 80 K=I1,J
            INXJ = JX+K
            A(INXJ) = A(INXJ)-SAVE*A(KX)
            KX = KX+K
   80    CONTINUE
         JX = JX+J
   85 CONTINUE
C                                  ICHNG(N+I) = 1 IF A 1 BY 1 PIVOT WAS
C                                  USED AT ROW I
   90 ICHNG(N+I) = 1
      I = I1
      L = IX
      GO TO 120
C                                  WE USE A 2 BY 2 PIVOT
   95 TEMP = A(IX)*A(IX1)
      SAVE = ABS(TEMP)+A(IXI)*A(IXI)
      DET(IDET) = TEMP-A(IXI)*A(IXI)
      IF (SAVE.EQ.ZERO) SAVE = ONE
      SAVE = N*SAVE
      TEMP = SAVE+ABS(DET(IDET))
      IF (TEMP.EQ.SAVE) GO TO 200
      I2 = I+2
      IF (I2.GT.N) GO TO 115
      JX = IX+I1
      L = JX
      TEMP1 = ONE/DET(IDET)
      IDET = IDET+1
      J0 = I1
      II = L+I
      DO 110 J=I2,N
         JXI = JX+I
         JXI1 = JX+I1
         SAVE = A(JXI)
         TEMP = A(JXI1)
         IF (I2.GT.J0) GO TO 105
         KX = II
C                                  SET A(J,K) TO A NEW VALUE
         DO 100 K=I2,J0
            JXK = JX+K
            A(JXK) = A(JXK)-A(KX)*SAVE-A(KX+1)*TEMP
            KX = KX+K
  100    CONTINUE
C                                  SET A(J,I) AND A(J,I+1) TO THE
C                                  APPROPRIATE MULTIPLIER
  105    A(JXI) = (A(IX1)*SAVE-A(IXI)*TEMP)*TEMP1
         A(JXI1) = (A(IX)*TEMP-A(IXI)*SAVE)*TEMP1
         JX = JX+J
         A(JX) = A(JX)-A(JXI)*SAVE-A(JXI1)*TEMP
         J0 = J
  110 CONTINUE
C                                  IF A 2 BY 2 PIVOT IS USED AT ROW I,
C                                  THEN ICHNG(N+I+1) = 0
  115 ICHNG(N+I) = 2
      ICHNG(N+I1) = 0
      L = IX+I1
      I = I2
  120 IF (I.LE.N) GO TO 15
C                                  SOLVE M*D*MT*X = B WHERE MT =
C                                  M-TRANSPOSE
  125 IF (IJOB.EQ.1) GO TO 9005
C                                  SOLVE MC = B AND STORE C IN B
      DO 195 JC=1,M
         IF (N .EQ. 1) GO TO 145
         DO 130 I=1,N
            SAVE = B(I,JC)
            J = ICHNG(I)
            B(I,JC) = B(J,JC)
            B(J,JC) = SAVE
  130    CONTINUE
         L = 1
         DO 143 I=1,N0
            JX = L+I
            I1 = I+1
            L = L+I1
            SAVE = B(I,JC)
            IF (SAVE .EQ. ZERO) GO TO 143
            IF (ICHNG(N+I).NE.2) GO TO 135
            JX = JX+I1
            I1 = I1+1
            IF (I1 .GT. N) GO TO 143
  135       DO 140 J=I1,N
               B(J,JC) = B(J,JC)-A(JX)*SAVE
               JX = JX+J
  140       CONTINUE
  143    CONTINUE
C                                  SOLVE DY=C AND STORE Y IN B
  145    L = 0
         IDET = 1
         DO 165 I=1,N
            L = L+I
            IF (ICHNG(N+I)-1) 150,155,160
  150       B(I,JC) = (B(I,JC)*A(L-I)-SAVE*A(L-1))/DET(IDET)
            IDET = IDET+1
            GO TO 165
  155       B(I,JC) = B(I,JC)/A(L)
            GO TO 165
  160       SAVE = B(I,JC)
            B(I,JC) = (SAVE*A(L+I+1)-B(I+1,JC)*A(L+I))/DET(IDET)
  165    CONTINUE
C                                  SOLVE MT*X = Y AND STORE X IN B
         NN = N0
         IF (ICHNG(N+N).EQ.0) NN = NN-1
         IF (NN.LT.1) GO TO 185
         NN1 = NN+1
         L = (NN*NN1)/2+NN
         DO 180 K=1,NN
            I = NN1-K
            I1 = I+1
            JX = L
            L = L-I1
            SAVE = B(I,JC)
            IF (ICHNG(N+I).NE.2) GO TO 170
            JX = JX+I1
            I1 = I1+1
  170       DO 175 J=I1,N
               SAVE = SAVE-A(JX)*B(J,JC)
               JX = JX+J
  175       CONTINUE
            B(I,JC) = SAVE
  180    CONTINUE
  185    DO 190 K=1,N
            I = N1-K
            SAVE = B(I,JC)
            J = ICHNG(I)
            B(I,JC) = B(J,JC)
            B(J,JC) = SAVE
  190    CONTINUE
  195 CONTINUE
      GO TO 9005
  200 IER = 129
      CALL UERTST (IER,6HLEQ1S )
 9005 RETURN
      END

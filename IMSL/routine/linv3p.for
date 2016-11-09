C   IMSL ROUTINE NAME   - LINV3P
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - IN PLACE INVERSE, EQUATION SOLUTION, POSITIVE
C                           DEFINITE MATRIX - SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL LINV3P (A,B,IJOB,N,IER)
C
C   ARGUMENTS    A      - INPUT/OUTPUT VECTOR OF LENGTH N(N+1)/2.
C                         ON INPUT, A CONTAINS THE N BY N POSITIVE
C                           DEFINITE SYMMETRIC MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE.
C                         ON OUTPUT, SEE PARAMETER IJOB.
C                B      - INPUT/OUTPUT VECTOR OF LENGTH N WHEN IJOB =
C                           2 OR 3.
C                         ON INPUT, B CONTAINS THE RIGHT HAND SIDE OF
C                           THE EQUATION AX=B.
C                         ON OUTPUT, THE SOLUTION X REPLACES B WHEN
C                           IJOB=2 OR 3. OTHERWISE, B IS NOT USED.
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES:
C                           I = 1, INVERT MATRIX A. A IS REPLACED BY ITS
C                             INVERSE STORED IN SYMMETRIC STORAGE MODE.
C                           I = 2, SOLVE THE EQUATION AX = B. A IS
C                             REPLACED BY THE DECOMPOSED MATRIX L SUCH
C                             THAT A = L*L-TRANSPOSE. L IS STORED IN
C                             SYMMETRIC STORAGE MODE. THE DIAGONAL OF L
C                             CONTAINS THE RECIPROCALS OF THE ACTUAL
C                             DIAGONAL ELEMENTS OF L.
C                           I = 3, SOLVE AX = B AND INVERT MATRIX A.
C                             A IS REPLACED BY ITS INVERSE AND THE
C                             SOLUTION X REPLACES B.
C                N      - ORDER OF A. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IJOB WAS LESS THAN
C                             1 OR GREATER THAN 3.
C                           IER = 130 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C                             (SEE THE CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LINV3P (A,B,IJOB,N,IER)
C
      REAL               A(1),B(N),ONE,Q,RN,S,SIXTN,T,X,ZERO
      DATA               ZERO/0.0/,ONE/1.0/,SIXTN/16.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IJOB.LT.1.OR.IJOB.GT.3) GO TO 115
      ISW = 0
C                                  CHOLESKY DECOMPOSITION OF A
C                                     A = L*L-TRANSPOSE
      RN = ONE/(N*SIXTN)
      IP = 1
      IER = 0
      DO 30 I=1,N
         IQ = IP
         IR = 1
         DO 25 J=1,I
            X = A(IP)
            IF (J.EQ.1) GO TO 10
            DO 5 K=IQ,IP1
               X = X-A(K)*A(IR)
               IR = IR+1
    5       CONTINUE
   10       IF (I.NE.J) GO TO 15
            Q = A(IP)+X*RN
            IF (Q.LE.A(IP)) GO TO 120
            A(IP) = ONE/SQRT(X)
            GO TO 20
   15       A(IP) = X*A(IR)
   20       IP1 = IP
            IP = IP+1
            IR = IR+1
   25    CONTINUE
   30 CONTINUE
      IF (ISW.EQ.0.AND.IJOB.NE.1) GO TO 75
C                                  FORM A-INVERSE
   35 NM1 = N-1
      IF (N.EQ.1) GO TO 55
C                                  COMPUTE L-INVERSE FIRST
      II = 1
      DO 50 I=1,NM1
         IP1 = I+1
         JM1 = I
         JJ = II
         DO 45 J=IP1,N
            S = ZERO
            LI = II
            JI = JJ+I
            JL = JI
            DO 40 L=I,JM1
               S = S+A(LI)*A(JL)
               JL = JL+1
               LI = LI+L
   40       CONTINUE
            JJ = JJ+J
            A(JI) = -A(JJ)*S
            JM1 = J
   45    CONTINUE
         II = II+IP1
   50 CONTINUE
   55 II = 0
C                                  NOW FORM A-INVERSE
      DO 70 I=1,N
         JJ = II
         DO 65 J=I,N
            S = ZERO
            JI = JJ+I
            LI = JI
            JJ = JJ+J
            LJ = JJ
            DO 60 L=J,N
               S = S+A(LI)*A(LJ)
               LI = LI+L
               LJ = LJ+L
   60       CONTINUE
            A(JI) = S
   65    CONTINUE
         II = II+I
   70 CONTINUE
      IF (IJOB.EQ.1.OR.ISW.EQ.1) GO TO 9005
C                                  SOLVE AX = B
   75 ISW = 1
      IP = 1
      IW = 0
C                                  SOLUTION OF LY = B
      IM1 = 0
      DO 95 I=1,N
         T = B(I)
         IF (IW.EQ.0) GO TO 85
         IP = IP+IW-1
         DO 80 K=IW,IM1
            T = T-A(IP)*B(K)
            IP = IP+1
   80    CONTINUE
         GO TO 90
   85    IF (T.NE.ZERO) IW = I
         IP = IP+IM1
   90    B(I) = T*A(IP)
         IP = IP+1
         IM1 = I
   95 CONTINUE
C                                  SOLUTION OF UX = Y
C                                  WHERE U = L-TRANSPOSE
      N1 = N+1
      DO 110 I=1,N
         II = N1-I
         IP = IP-1
         IS = IP
         IQ = II+1
         T = B(II)
         IF (N.LT.IQ) GO TO 105
         KK = N
         DO 100 K=IQ,N
            T = T-A(IS)*B(KK)
            KK = KK-1
            IS = IS-KK
  100    CONTINUE
  105    B(II) = T*A(IS)
  110 CONTINUE
      IF (IJOB.EQ.3) GO TO 35
      GO TO 9005
C                                  IJOB OUT OF RANGE
  115 IER = 129
      GO TO 9000
C                                  A IS ALGORITHMICALLY NOT
C                                     POSITIVE DEFINITE
  120 IER = 130
 9000 CONTINUE
      CALL UERTST (IER,6HLINV3P)
 9005 RETURN
      END

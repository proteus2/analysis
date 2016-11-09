      SUBROUTINE LUELMP (A,B,N,X)
      DIMENSION      A(1),B(1),X(1)
      DATA           ZERO/0.0/
C                                   first executable statement
C                                   solution of lY = B
      IP = 1
      IW = 0
      DO 15 I =1,N
         T = B(I)
         IM1 = I - 1
         IF (IW .EQ. 0) GO TO 9
         IP = IP + IW - 1
         DO 5 K = IW,IM1
            T = T - A(IP) * X(K)
            IP = IP + 1
5        CONTINUE
         GO TO 10
9        IF (T .NE. ZERO) IW = I
         IP = IP + IM1
10       X(I) = T * A(IP)
         IP = IP + 1
15    CONTINUE
C                                               solution of UX = Y
      N1 = N + 1
      DO 30 I = 1,N
         II = N1 - I
         IP = IP - 1
         IS = IP
         IQ = II + 1
         T = X(II)
         IF (N .LT. IQ) GO TO 25
         KK = N
         DO 20 K = IQ,N
            T = T -A(IS) * X(KK)
            KK = KK - 1
            IS = IS - KK
20       CONTINUE
25       X(II) = T * A(IS)
30    CONTINUE
      RETURN
      END



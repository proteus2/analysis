C   IMSL ROUTINE NAME   - IBCDCU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE IBCCCU
C
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IBCDCU (TAU,GTAU,N,M,W,VS,IC1,IC2,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IC1,IC2,IER
      REAL               TAU(N),GTAU(IC1,1),W(N,2),VS(IC2,2,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JJ,JM1,JP1,J,K,LIM,LL,LP1,NM1
      REAL               AA,BB,C1,C2,CC,DD,DTAU,G,H,RATIO,U,XILIM
C                                  FIRST EXECUTABLE STATEMENT
      LIM = N-3
      NM1 = N-1
      LP1 = LIM+1
      IER = 132
      W(2,1) = TAU(3)-TAU(1)
      IF (W(2,1).LE.0.0) RETURN
      DO 5 K=1,M
         VS(K,1,1) = GTAU(1,K)
    5 CONTINUE
      XILIM = TAU(1)
      IF (LIM.LT.2) GO TO 20
      XILIM = TAU(N-2)
      DO 15 I=2,LIM
         J = I+1
         W(J,1) = TAU(I+2)-TAU(J)
         IF (W(J,1).LE.0.0) RETURN
         DO 10 K=1,M
   10    VS(K,1,I) = GTAU(J,K)
   15 CONTINUE
   20 W(LP1,1) = TAU(N)-XILIM
      IF (W(LP1,1).LE.0.0) RETURN
      DO 25 K=1,M
   25 VS(K,1,LP1) = GTAU(N,K)
      DO 35 I=2,LP1
         DO 30 K=1,M
   30    VS(K,2,I) = (VS(K,1,I)-VS(K,1,I-1))/W(I,1)
   35 CONTINUE
      DTAU = TAU(2)-TAU(1)
      RATIO = DTAU/W(2,1)
      W(1,2) = (RATIO-1.)**2
      W(1,1) = RATIO*(RATIO-1.)
      C1 = RATIO*(2.*RATIO-3.)
      DO 40 K=1,M
   40 VS(K,2,1) = (GTAU(2,K)-GTAU(1,K))/DTAU+VS(K,2,2)*C1
      IF (LIM.LT.2) GO TO 55
      DO 50 I=2,LIM
         J = I+1
         JJ = I-1
         G = -W(J,1)/W(JJ,2)
         C1 = 3.*W(I,1)
         C2 = 3.*W(J,1)
         DO 45 K=1,M
   45    VS(K,2,I) = G*VS(K,2,JJ)+C1*VS(K,2,J)+C2*VS(K,2,I)
         W(I,2) = G*W(JJ,1)+2.*(W(I,1)+W(J,1))
   50 CONTINUE
   55 DTAU = TAU(N-1)-XILIM
      RATIO = DTAU/W(LP1,1)
      G = -(RATIO-1.)**2/W(LIM,2)
      W(LP1,2) = RATIO*(RATIO-1.)
      C1 = RATIO*(2.*RATIO-3.)
      DO 60 K=1,M
   60 VS(K,2,LP1) = (GTAU(N-1,K)-VS(K,1,LIM))/DTAU+VS(K,2,LP1)*C1
      W(LP1,2) = G*W(LIM,1)+W(LP1,2)
      DO 65 K=1,M
   65 VS(K,2,LP1) = (G*VS(K,2,LIM)+VS(K,2,LP1))/W(LP1,2)
      J = LIM
   70 DO 75 K=1,M
   75 VS(K,2,J) = (VS(K,2,J)-W(J,1)*VS(K,2,J+1))/W(J,2)
      J = J-1
      IF (J.GT.0) GO TO 70
      DO 95 K=1,M
         DO 85 JJ=1,N
            J = N+1-JJ
            JM1 = J-1
            IF (J.EQ.N) JM1 = J-2
            IF (J.EQ.1) JM1 = J
            DO 80 LL=1,2
               VS(K,LL,J) = VS(K,LL,JM1)
   80       CONTINUE
   85    CONTINUE
         DO 90 J=2,NM1,LIM
            JM1 = J-1
            JP1 = J+1
            IF (JM1.EQ.2) JM1 = 1
            IF (JP1.EQ.NM1) JP1 = N
            H = TAU(JP1)-TAU(JM1)
            U = TAU(J)-TAU(JM1)
            AA = VS(K,1,JM1)
            BB = VS(K,2,JM1)
            CC = (3.*(VS(K,1,JP1)-VS(K,1,JM1))/H-(VS(K,2,JP1)+
     *      2.*VS(K,2,JM1)))/H
            DD = (2.*(VS(K,1,JM1)-VS(K,1,JP1))/H+(VS(K,2,JP1)+
     *      VS(K,2,JM1)))/H**2
            VS(K,1,J) = AA+U*(BB+U*(CC+DD*U))
            VS(K,2,J) = BB+U*(2.*CC+3.*DD*U)
   90    CONTINUE
   95 CONTINUE
      IER = 0
      RETURN
      END

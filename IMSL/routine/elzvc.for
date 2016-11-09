C   IMSL ROUTINE NAME   - ELZVC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGZC
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ELZVC  (N,A,IA,B,IB,Z,IZ,IJOB,EIGA,EIGB,INFER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IB,IZ,IJOB,INFER,IER
      COMPLEX            A(IA,N),B(IB,N),EIGA(N),EIGB(N),Z(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NN,IM1,I,ITS,NNP2,LB,L,LM1,NL,MB,M,MN1,
     *                   LOR1,NNORN,NM1,JP1,J,K,L1
      REAL               EPSA,EPSB,SS,R,ANORM,BNORM,ANI,BNI,C,
     *                   D0,D1,D2,E0,E1,ONE,PT8,TWO,ZRO
      COMPLEX            ANNM1,ALFM,BETM,D,SL,DEN,NUM,ANM1M1,S,W,Y,
     *                   ZC,ZERO
      DATA               PT8/.8/,TWO/2.0/,ZRO/0.0/,ONE/1.0/,
     *                   ZERO/(0.0,0.0)/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      INFER = 0
      NN = N
C                                  COMPUTE THE MACHINE PRECISION TIMES
C                                    THE NORM OF A AND B
      ANORM = ZRO
      BNORM = ZRO
      IM1 = 1
      DO 15 I=1,N
         ANI = ZRO
         IF (I.EQ.1) GO TO 5
         Y = A(I,IM1)
         ANI = ANI+ABS(REAL(Y))+ABS(AIMAG(Y))
    5    BNI = ZRO
         DO 10 J=I,N
            ANI = ANI+ABS(REAL(A(I,J)))+ABS(AIMAG(A(I,J)))
            BNI = BNI+ABS(REAL(B(I,J)))+ABS(AIMAG(B(I,J)))
   10    CONTINUE
         IF (ANI.GT.ANORM) ANORM = ANI
         IF (BNI.GT.BNORM) BNORM = BNI
         IM1 = I
   15 CONTINUE
      IF (ANORM.EQ.ZRO) ANORM = ONE
      IF (BNORM.EQ.ZRO) BNORM = ONE
      EPSB = BNORM
      EPSA = ANORM
   20 EPSA = EPSA/TWO
      EPSB = EPSB/TWO
      C = ANORM+EPSA
      IF (C.GT.ANORM) GO TO 20
      IF (N.LE.1) GO TO 160
   25 ITS = 0
      NM1 = NN-1
C                                  CHECK FOR NEGLIGIBLE
C                                    SUBDIAGONAL ELEMENTS
   30 D2 = ABS(REAL(A(NN,NN)))+ABS(AIMAG(A(NN,NN)))
      NNP2 = NN+2
      DO 35 LB=2,NN
         L = NNP2-LB
         LM1 = L-1
         SS = D2
         Y = A(LM1,LM1)
         D2 = ABS(REAL(Y))+ABS(AIMAG(Y))
         SS = SS+D2
         Y = A(L,LM1)
         R = SS+ABS(REAL(Y))+ABS(AIMAG(Y))
         IF (R.EQ.SS) GO TO 40
   35 CONTINUE
      L = 1
   40 IF (L.EQ.NN) GO TO 160
      IF (ITS.LT.30) GO TO 45
      IF(INFER.EQ.0) INFER = NN
      IER = 129
      GO TO 9000
   45 IF (ITS.EQ.10.OR.ITS.EQ.20) GO TO 55
C                                  COMPUTE SHIFT AS EIGENVALUE
C                                    OF LOWER 2 BY 2
      ANNM1 = A(NN,NM1)
      ANM1M1 = A(NM1,NM1)
      S = A(NN,NN)*B(NM1,NM1)-ANNM1*B(NM1,NN)
      W = ANNM1*B(NN,NN)*(A(NM1,NN)*B(NM1,NM1)-B(NM1,NN)*ANM1M1)
      Y = (ANM1M1*B(NN,NN)-S)/TWO
      ZC = CSQRT(Y*Y+W)
      IF (REAL(ZC).EQ.ZRO.AND.AIMAG(ZC).EQ.ZRO) GO TO 50
      D0 = REAL(Y/ZC)
      IF (D0.LT.ZRO) ZC = -ZC
   50 DEN = (Y+ZC)*B(NM1,NM1)*B(NN,NN)
      IF (REAL(DEN).EQ.ZRO.AND.AIMAG(DEN).EQ.ZRO) DEN =
     *CMPLX(EPSA,ZRO)
      NUM = (Y+ZC)*S-W
      GO TO 60
C                                  AD-HOC SHIFT
   55 Y = A(NM1,NN-2)
      NUM = CMPLX(ABS(REAL(ANNM1))+ABS(AIMAG(ANNM1)),ABS(REAL(Y))+
     *ABS(AIMAG(Y)))
      DEN = CMPLX(ONE,ZRO)
C                                  CHECK FOR 2 CONSECUTIVE SMALL
C                                    SUBDIAGONAL ELEMENTS
   60 IF (NN.EQ.L+1) GO TO 70
      D2 = ABS(REAL(A(NM1,NM1)))+ABS(AIMAG(A(NM1,NM1)))
      E1 = ABS(REAL(ANNM1))+ABS(AIMAG(ANNM1))
      D1 = ABS(REAL(A(NN,NN)))+ABS(AIMAG(A(NN,NN)))
      NL = NN-(L+1)
      DO 65 MB=1,NL
         M = NN-MB
         MN1 = M-1
         E0 = E1
         Y = A(M,MN1)
         E1 = ABS(REAL(Y))+ABS(AIMAG(Y))
         D0 = D1
         D1 = D2
         Y = A(MN1,MN1)
         D2 = ABS(REAL(Y))+ABS(AIMAG(Y))
         Y = A(M,M)*DEN-B(M,M)*NUM
         D0 = (D0+D1+D2)*(ABS(REAL(Y))+ABS(AIMAG(Y)))
         E0 = E0*E1*(ABS(REAL(DEN))+ABS(AIMAG(DEN)))+D0
         IF (E0.EQ.D0) GO TO 75
   65 CONTINUE
   70 M = L
   75 CONTINUE
      ITS = ITS+1
      W = A(M,M)*DEN-B(M,M)*NUM
      ZC = A(M+1,M)*DEN
      D1 = ABS(REAL(ZC))+ABS(AIMAG(ZC))
      D2 = ABS(REAL(W))+ABS(AIMAG(W))
C                                  FIND L AND M AND SET A=LAM AND B=LBM
      LOR1 = L
      NNORN = NN
      IF (IJOB.NE.1) GO TO 80
      LOR1 = 1
      NNORN = N
   80 IM1 = M
      DO 155 I=M,NM1
         J = I+1
         JP1 = J+1
C                                  FIND ROW TRANSFORMATIONS TO RESTORE
C                                    A TO UPPER HESSENBERG FORM.
C                                    APPLY TRANSFORMATIONS TO A AND B
         IF (I.EQ.M) GO TO 85
         W = A(I,IM1)
         ZC = A(J,IM1)
         D1 = ABS(REAL(ZC))+ABS(AIMAG(ZC))
         D2 = ABS(REAL(W))+ABS(AIMAG(W))
         IF (D1.EQ.ZRO) GO TO 30
   85    IF (D2.GT.D1) GO TO 95
C                                  MUST INTERCHANGE ROWS
         DO 90 K=I,NNORN
            Y = A(I,K)
            A(I,K) = A(J,K)
            A(J,K) = Y
            Y = B(I,K)
            B(I,K) = B(J,K)
            B(J,K) = Y
   90    CONTINUE
         IF (I.GT.M) A(I,IM1) = A(J,IM1)
         IF (D2.EQ.ZRO) GO TO 110
C                                  THE SCALING OF W AND Z IS DESIGNED
C                                    TO AVOID A DIVISION BY ZERO
C                                    WHEN THE DENOMINATOR IS SMALL
         Y = CMPLX(REAL(W)/D1,AIMAG(W)/D1)/CMPLX(REAL(ZC)/D1,
     *   AIMAG(ZC)/D1)
         GO TO 100
   95    Y = CMPLX(REAL(ZC)/D2,AIMAG(ZC)/D2)/CMPLX(REAL(W)/D2,
     *   AIMAG(W)/D2)
  100    DO 105 K=I,NNORN
            A(J,K) = A(J,K)-Y*A(I,K)
            B(J,K) = B(J,K)-Y*B(I,K)
  105    CONTINUE
  110    IF (I.GT.M) A(J,IM1) = ZERO
C                                  PERFORM TRANSFORMATIONS FROM RIGHT
C                                    TO RESTORE B TO TRIANGULAR FORM
C                                  APPLY TRANSFORMATIONS TO A
         ZC = B(J,I)
         IM1 = I
         W = B(J,J)
         D2 = ABS(REAL(W))+ABS(AIMAG(W))
         D1 = ABS(REAL(ZC))+ABS(AIMAG(ZC))
         IF (D1.EQ.ZRO) GO TO 30
         IF (D2.GT.D1) GO TO 135
C                                  MUST INTERCHANGE COLUMNS
         DO 115 K=LOR1,J
            Y = A(K,J)
            A(K,J) = A(K,I)
            A(K,I) = Y
            Y = B(K,J)
            B(K,J) = B(K,I)
            B(K,I) = Y
  115    CONTINUE
         IF (I.EQ.NM1) GO TO 120
         Y = A(JP1,J)
         A(JP1,J) = A(JP1,I)
         A(JP1,I) = Y
  120    IF (IJOB.NE.1) GO TO 130
         DO 125 K=1,N
            Y = Z(K,J)
            Z(K,J) = Z(K,I)
            Z(K,I) = Y
  125    CONTINUE
  130    B(J,I) = ZERO
         IF (D2.EQ.ZRO) GO TO 155
         ZC = CMPLX(REAL(W)/D1,AIMAG(W)/D1)/CMPLX(REAL(ZC)/D1,
     *   AIMAG(ZC)/D1)
         GO TO 140
  135    ZC = CMPLX(REAL(ZC)/D2,AIMAG(ZC)/D2)/CMPLX(REAL(W)/D2,
     *   AIMAG(W)/D2)
  140    DO 145 K=LOR1,J
            A(K,I) = A(K,I)-ZC*A(K,J)
            B(K,I) = B(K,I)-ZC*B(K,J)
  145    CONTINUE
         B(J,I) = ZERO
         IF (I.LT.NM1) A(JP1,I) = A(JP1,I)-ZC*A(JP1,J)
         IF (IJOB.NE.1) GO TO 155
         DO 150 K=1,N
            Z(K,I) = Z(K,I)-ZC*Z(K,J)
  150    CONTINUE
  155 CONTINUE
      GO TO 30
  160 CONTINUE
      EIGA(NN) = A(NN,NN)
      EIGB(NN) = B(NN,NN)
      IF (NN.EQ.1) GO TO 165
      NN = NM1
      IF (NN.GT.1) GO TO 25
      GO TO 160
C                                  FIND EIGENVECTORS USING B FOR
C                                    INTERMEDIATE STORAGE
  165 IF (IJOB.NE.1) GO TO 9000
      M = N
  170 CONTINUE
      ALFM = A(M,M)
      BETM = B(M,M)
      B(M,M) = CMPLX(ONE,ZRO)
      L = M-1
      IF (L.EQ.0) GO TO 185
  175 CONTINUE
      L1 = L+1
      SL = ZERO
      DO 180 J=L1,M
         SL = SL+(BETM*A(L,J)-ALFM*B(L,J))*B(J,M)
  180 CONTINUE
      Y = BETM*A(L,L)-ALFM*B(L,L)
      IF (REAL(Y).EQ.ZRO.AND.AIMAG(Y).EQ.ZRO) Y =
     *CMPLX((EPSA+EPSB)/TWO,ZRO)
      B(L,M) = -SL/Y
      L = L-1
  185 IF (L.GT.0) GO TO 175
      M = M-1
      IF (M.GT.0) GO TO 170
C                                  TRANSFORM TO ORIGINAL COORDINATE
C                                    SYSTEM
      M = N
  190 CONTINUE
      DO 200 I=1,N
         S = ZERO
         DO 195 J=1,M
            S = S+Z(I,J)*B(J,M)
  195    CONTINUE
         Z(I,M) = S
  200 CONTINUE
      M = M-1
      IF (M.GT.0) GO TO 190
C                                  NORMALIZE SO THAT LARGEST
C                                    COMPONENT = 1.0
      M = N
  205 CONTINUE
      SS = ZRO
      DO 210 I=1,N
         R = ABS(REAL(Z(I,M)))+ABS(AIMAG(Z(I,M)))
         IF (R.LT.SS) GO TO 210
         SS = R
         D = Z(I,M)
  210 CONTINUE
      IF (SS.EQ.ZRO) GO TO 220
      DO 215 I=1,N
         Z(I,M) = Z(I,M)/D
  215 CONTINUE
  220 M = M-1
      IF (M.GT.0) GO TO 205
      GO TO 9005
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HELZVC )
 9005 RETURN
      END

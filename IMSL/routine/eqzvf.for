C   IMSL ROUTINE NAME   - EQZVF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGZF
C
C   REQD. IMSL ROUTINES - VHSH2C,VHSH2R
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQZVF  (A,IA,B,IB,N,EPSA,EPSB,ALFR,ALFI,BETA,Z,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,N,IZ
      REAL               A(IA,N),B(IB,N),ALFR(N),ALFI(N),BETA(N),Z(IZ,1)
      REAL               EPSA,EPSB
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            M,L,I,J,L1,K,MR,MI,IRET
      REAL               E,CZ,TI,V1,A21,B11,SQR,SZR,A11R,A21R,
     1                   AN,EI,TR,V2,A22,B12,SSI,A12I,A22I,
     2                   C,R,BN,ER,U1,A11,BDI,B22,SSR,A12R,
     3                   A22R,D,T,CQ,U2,A12,BDR,SQI,SZI,A11I,A21I,
     4                   DI,SLI,TLK,TKLR,TLLR,DR,SLR,TLL,ALMI,
     5                   TKKI,TLKI,S,SK,SKI,TKK,ALMR,TKKR,TLKR,SL,SKR,
     6                   TKL,ALFM,TKLI,TLLI,XR,XI,YR,YI,ZR,ZI,
     7                   H,BETM,ZERO,ONE,F,HALF
      LOGICAL            WANTX,FLIP
      DATA               ZERO/0.0/,HALF/0.5/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      WANTX = .FALSE.
      IF (IZ.GE.N) WANTX = .TRUE.
C                                  FIND EIGENVALUES OF QUASI-TRIANGULAR
C                                     MATRICES
      M = N
    5 CONTINUE
      IF (M.EQ.1) GO TO 10
      IF (A(M,M-1).NE.ZERO) GO TO 15
C                                  ONE-BY-ONE SUBMATRIX, ONE REAL ROOT
   10 ALFR(M) = A(M,M)
      IF(B(M,M).LT.ZERO) ALFR(M) = -ALFR(M)
      BETA(M) = ABS(B(M,M))
      ALFI(M) = ZERO
      M = M-1
      GO TO 70
C                                  TWO-BY-TWO SUBMATRIX
   15 L = M-1
      IF (ABS(B(L,L)).GT.EPSB) GO TO 20
      B(L,L) = ZERO
      CALL VHSH2R (A(L,L),A(M,L),U1,U2,V1,V2)
      GO TO 50
   20 IF (ABS(B(M,M)).GT.EPSB) GO TO 25
      B(M,M) = ZERO
      CALL VHSH2R (A(M,M),A(M,L),U1,U2,V1,V2)
      BN = ZERO
      GO TO 30
   25 AN = ABS(A(L,L))+ABS(A(L,M))+ABS(A(M,L))+ABS(A(M,M))
      BN = ABS(B(L,L))+ABS(B(L,M))+ABS(B(M,M))
      F = ONE/AN
      A11 = A(L,L)*F
      A12 = A(L,M)*F
      A21 = A(M,L)*F
      A22 = A(M,M)*F
      F = ONE/BN
      B11 = B(L,L)*F
      B12 = B(L,M)*F
      B22 = B(M,M)*F
      E = A11/B11
      C = ((A22-E*B22)/B22-(A21*B12)/(B11*B22))*HALF
      D = C*C+(A21*(A12-E*B12))/(B11*B22)
      IF (D.LT.ZERO) GO TO 65
C                                  TWO REAL ROOTS
C                                  ZERO BOTH A(M,L) AND B(M,L)
      IF (C.GE.ZERO) E = E+(C+SQRT(D))
      IF (C.LT.ZERO) E = E+(C-SQRT(D))
      A11 = A11-E*B11
      A12 = A12-E*B12
      A22 = A22-E*B22
      FLIP = (ABS(A11)+ABS(A12)).GE.(ABS(A21)+ABS(A22))
      IF (FLIP) CALL VHSH2R (A12,A11,U1,U2,V1,V2)
      IF (.NOT.FLIP) CALL VHSH2R (A22,A21,U1,U2,V1,V2)
   30 IF (U1.NE.ONE) GO TO 45
      DO 35 I=1,M
         T = A(I,M)+U2*A(I,L)
         A(I,M) = A(I,M)+V1*T
         A(I,L) = A(I,L)+V2*T
         T = B(I,M)+U2*B(I,L)
         B(I,M) = B(I,M)+V1*T
         B(I,L) = B(I,L)+V2*T
   35 CONTINUE
      IF (.NOT.WANTX) GO TO 45
      DO 40 I=1,N
         T = Z(I,M)+U2*Z(I,L)
         Z(I,M) = Z(I,M)+V1*T
         Z(I,L) = Z(I,L)+V2*T
   40 CONTINUE
   45 IF (BN.EQ.ZERO) GO TO 60
      FLIP = AN.GE.ABS(E)*BN
      IF (FLIP) CALL VHSH2R (B(L,L),B(M,L),U1,U2,V1,V2)
      IF (.NOT.FLIP) CALL VHSH2R (A(L,L),A(M,L),U1,U2,V1,V2)
   50 IF (U1.NE.ONE) GO TO 60
      DO 55 J=L,N
         T = A(L,J)+U2*A(M,J)
         A(L,J) = A(L,J)+V1*T
         A(M,J) = A(M,J)+V2*T
         T = B(L,J)+U2*B(M,J)
         B(L,J) = B(L,J)+V1*T
         B(M,J) = B(M,J)+V2*T
   55 CONTINUE
   60 A(M,L) = ZERO
      B(M,L) = ZERO
      ALFR(L) = A(L,L)
      ALFR(M) = A(M,M)
      IF(B(L,L).LT.ZERO) ALFR(L) = -ALFR(L)
      IF(B(M,M).LT.ZERO) ALFR(M) = -ALFR(M)
      BETA(L) = ABS(B(L,L))
      BETA(M) = ABS(B(M,M))
      ALFI(M) = ZERO
      ALFI(L) = ZERO
      M = M-2
      GO TO 70
C                                  TWO COMPLEX ROOTS
   65 ER = E+C
      EI = SQRT(-D)
      A11R = A11-ER*B11
      A11I = EI*B11
      A12R = A12-ER*B12
      A12I = EI*B12
      A21R = A21
      A21I = ZERO
      A22R = A22-ER*B22
      A22I = EI*B22
      FLIP = (ABS(A11R)+ABS(A11I)+ABS(A12R)+ABS(A12I))
     1.GE.(ABS(A21R)+ABS(A22R)+ABS(A22I))
      IF (FLIP) CALL VHSH2C (A12R,A12I,-A11R,-A11I,CZ,SZR,SZI)
      IF (.NOT.FLIP) CALL VHSH2C (A22R,A22I,-A21R,-A21I,CZ,SZR,SZI)
      FLIP = AN.GE.(ABS(ER)+ABS(EI))*BN
      IF (FLIP) CALL VHSH2C (CZ*B11+SZR*B12,SZI*B12,SZR*B22,SZI*B22,CQ,S
     1QR,SQI)
      IF (.NOT.FLIP) CALL VHSH2C (CZ*A11+SZR*A12,SZI*A12,CZ*A21+SZR*A22,
     1SZI*A22,CQ,SQR,SQI)
      SSR = SQR*SZR+SQI*SZI
      SSI = SQR*SZI-SQI*SZR
      TR = CQ*CZ*A11+CQ*SZR*A12+SQR*CZ*A21+SSR*A22
      TI = CQ*SZI*A12-SQI*CZ*A21+SSI*A22
      BDR = CQ*CZ*B11+CQ*SZR*B12+SSR*B22
      BDI = CQ*SZI*B12+SSI*B22
      R = SQRT(BDR*BDR+BDI*BDI)
      BETA(L) = BN*R
      ALFR(L) = AN*(TR*BDR+TI*BDI)/R
      ALFI(L) = AN*(TR*BDI-TI*BDR)/R
      TR = SSR*A11-SQR*CZ*A12-CQ*SZR*A21+CQ*CZ*A22
      TI = -SSI*A11-SQI*CZ*A12+CQ*SZI*A21
      BDR = SSR*B11-SQR*CZ*B12+CQ*CZ*B22
      BDI = -SSI*B11-SQI*CZ*B12
      R = SQRT(BDR*BDR+BDI*BDI)
      BETA(M) = BN*R
      ALFR(M) = AN*(TR*BDR+TI*BDI)/R
      ALFI(M) = AN*(TR*BDI-TI*BDR)/R
      M = M-2
   70 IF (M.GT.0) GO TO 5
      IF (.NOT.WANTX) GO TO 240
C                                  FIND EIGENVECTORS OF QUASI-TRIANGULAR
C                                     MATRICES BY BACKSUBSTITUTION
C                                  USE B FOR INTERMEDIATE STORAGE,
C                                     M-TH VECTOR IN B(.,M)
      M = N
   75 CONTINUE
      IF (ALFI(M).NE.ZERO) GO TO 110
C                                  REAL VECTOR
      ALFM = ALFR(M)
      BETM = BETA(M)
      B(M,M) = ONE
      L = M-1
      IF (L.EQ.0) GO TO 105
   80 CONTINUE
      L1 = L+1
      SL = ZERO
      DO 85 J=L1,M
         SL = SL+(BETM*A(L,J)-ALFM*B(L,J))*B(J,M)
   85 CONTINUE
      IF (L.EQ.1) GO TO 90
      IF (BETM*A(L,L-1).NE.ZERO) GO TO 95
C                                  REAL 1-BY-1 BLOCK
   90 D = BETM*A(L,L)-ALFM*B(L,L)
      IF (D.EQ.ZERO) D = (EPSA+EPSB)*HALF
      B(L,M) = -SL/D
      L = L-1
      GO TO 105
C                                  REAL 2-BY-2 BLOCK
   95 K = L-1
      SK = ZERO
      DO 100 J=L1,M
         SK = SK+(BETM*A(K,J)-ALFM*B(K,J))*B(J,M)
  100 CONTINUE
      TKK = BETM*A(K,K)-ALFM*B(K,K)
      TKL = BETM*A(K,L)-ALFM*B(K,L)
      TLK = BETM*A(L,K)
      TLL = BETM*A(L,L)-ALFM*B(L,L)
      D = TKK*TLL-TKL*TLK
      IF (D.EQ.ZERO) D = (EPSA+EPSB)*HALF
      B(L,M) = (TLK*SK-TKK*SL)/D
      FLIP = ABS(TKK).GE.ABS(TLK)
      IF (FLIP) B(K,M) = -(SK+TKL*B(L,M))/TKK
      IF (.NOT.FLIP) B(K,M) = -(SL+TLL*B(L,M))/TLK
      L = L-2
  105 IF (L.GT.0) GO TO 80
      M = M-1
      GO TO 165
C                                  COMPLEX VECTOR
  110 ALMR = ALFR(M-1)
      ALMI = ALFI(M-1)
      BETM = BETA(M-1)
      MR = M-1
      MI = M
C                                  LET M-TH COMPONENT = (0.0,-1.0) SO
C                                     THAT B IS TRIANGULAR
C                                  (M-1)ST = -(BETM*A(M,M)-ALFM*B(M,M))*
C                                     *(M-TH)/(BETM*A(M,M-1))
      B(M-1,MR) = ALMI*B(M,M)/(BETM*A(M,M-1))
      B(M-1,MI) = (BETM*A(M,M)-ALMR*B(M,M))/(BETM*A(M,M-1))
      B(M,MR) = ZERO
      B(M,MI) = -ONE
      L = M-2
      IF (L.EQ.0) GO TO 160
  115 CONTINUE
      L1 = L+1
      SLR = ZERO
      SLI = ZERO
      DO 120 J=L1,M
         TR = BETM*A(L,J)-ALMR*B(L,J)
         TI = -ALMI*B(L,J)
         SLR = SLR+TR*B(J,MR)-TI*B(J,MI)
         SLI = SLI+TR*B(J,MI)+TI*B(J,MR)
  120 CONTINUE
      IF (L.EQ.1) GO TO 125
      IF (BETM*A(L,L-1).NE.ZERO) GO TO 135
C                                  COMPLEX 1-BY-1 BLOCK
  125 DR = BETM*A(L,L)-ALMR*B(L,L)
      DI = -ALMI*B(L,L)
      IRET = 1
      XR = -SLR
      XI = -SLI
      YR = DR
      YI = DI
      GO TO 225
  130 B(L,MR) = ZR
      B(L,MI) = ZI
      L = L-1
      GO TO 160
C                                  COMPLEX 2-BY-2 BLOCK
  135 K = L-1
      SKR = ZERO
      SKI = ZERO
      DO 140 J=L1,M
         TR = BETM*A(K,J)-ALMR*B(K,J)
         TI = -ALMI*B(K,J)
         SKR = SKR+TR*B(J,MR)-TI*B(J,MI)
         SKI = SKI+TR*B(J,MI)+TI*B(J,MR)
  140 CONTINUE
      TKKR = BETM*A(K,K)-ALMR*B(K,K)
      TKKI = -ALMI*B(K,K)
      TKLR = BETM*A(K,L)-ALMR*B(K,L)
      TKLI = -ALMI*B(K,L)
      TLKR = BETM*A(L,K)
      TLKI = ZERO
      TLLR = BETM*A(L,L)-ALMR*B(L,L)
      TLLI = -ALMI*B(L,L)
      DR = TKKR*TLLR-TKKI*TLLI-TKLR*TLKR
      DI = TKKR*TLLI+TKKI*TLLR-TKLI*TLKR
      IF (DR.EQ.ZERO.AND.DI.EQ.ZERO) DR = (EPSA+EPSB)*HALF
      IRET = 2
      XR = TLKR*SKR-TKKR*SLR+TKKI*SLI
      XI = TLKR*SKI-TKKR*SLI-TKKI*SLR
      YR = DR
      YI = DI
      GO TO 225
  145 B(L,MR) = ZR
      B(L,MI) = ZI
      FLIP = (ABS(TKKR)+ABS(TKKI)).GE.ABS(TLKR)
      IRET = 3
      IF (FLIP) GO TO 150
      XR = -SLR-TLLR*B(L,MR)+TLLI*B(L,MI)
      XI = -SLI-TLLR*B(L,MI)-TLLI*B(L,MR)
      YR = TLKR
      YI = TLKI
      GO TO 225
  150 XR = -SKR-TKLR*B(L,MR)+TKLI*B(L,MI)
      XI = -SKI-TKLR*B(L,MI)-TKLI*B(L,MR)
      YR = TKKR
      YI = TKKI
      GO TO 225
  155 B(K,MR) = ZR
      B(K,MI) = ZI
      L = L-2
  160 IF (L.GT.0) GO TO 115
      M = M-2
  165 IF (M.GT.0) GO TO 75
C                                  TRANSFORM TO ORIGINAL COORDINATE
C                                     SYSTEM
      M = N
  170 CONTINUE
      DO 180 I=1,N
         S = ZERO
         DO 175 J=1,M
            S = S+Z(I,J)*B(J,M)
  175    CONTINUE
         Z(I,M) = S
  180 CONTINUE
      M = M-1
      IF (M.GT.0) GO TO 170
C                                  NORMALIZE SO THAT LARGEST
C                                     COMPONENT = 1.
      M = N
  185 CONTINUE
      S = ZERO
      IF (ALFI(M).NE.ZERO) GO TO 200
      DO 190 I=1,N
         R = ABS(Z(I,M))
         IF (R.LT.S) GO TO 190
         S = R
         D = Z(I,M)
  190 CONTINUE
      DO 195 I=1,N
         Z(I,M) = Z(I,M)/D
  195 CONTINUE
      M = M-1
      GO TO 220
  200 DO 205 I=1,N
         R = ABS(Z(I,M-1))+ABS(Z(I,M))
         IF (R .EQ. ZERO) GO TO 205
         R = R*SQRT((Z(I,M-1)/R)**2+(Z(I,M)/R)**2)
         IF (R.LT.S) GO TO 205
         S = R
         DR = Z(I,M-1)
         DI = Z(I,M)
  205 CONTINUE
      IRET = 4
      I = 0
  210 I = I+1
      XR = Z(I,M-1)
      XI = Z(I,M)
      YR = DR
      YI = DI
      GO TO 225
  215 Z(I,M-1) = ZR
      Z(I,M) = ZI
      IF (I.LT.N) GO TO 210
      M = M-2
  220 IF (M.GT.0) GO TO 185
      GO TO 240
C                                  AVOID DIVISION BY ZERO
  225 IF (ABS(YR)+ABS(YI) .EQ. 0.0) YR = (EPSA+EPSB)/2.0
      IF (ABS(YR).LT.ABS(YI)) GO TO 230
      H = YI/YR
      F = YR+H*YI
      ZR = (XR+H*XI)/F
      ZI = (XI-H*XR)/F
      GO TO 235
  230 H = YR/YI
      F = YI+H*YR
      ZR = (H*XR+XI)/F
      ZI = (H*XI-XR)/F
  235 GO TO (130,145,155,215), IRET
  240 RETURN
      END

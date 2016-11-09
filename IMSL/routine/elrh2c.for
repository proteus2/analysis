C   IMSL ROUTINE NAME   - ELRH2C
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCC
C
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ELRH2C (HR,HI,K,L,N,IH,WR,WI,ZR,ZI,ID,INFER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,N,IH,INFER,IER
      REAL               HR(IH,1),HI(IH,1),WR(1),WI(1),ZR(IH,1),
     *                   ZI(IH,1),ID(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,IEND,II,IP1,IM1,M,NN,ITS,NNM1,NNM2,NPL,KK,
     *                   MM1,NNMJ,NM,MM,MMM1,MP1,NM1
      REAL               T1(2),T2(2),T3(2),ZERO,ONE,TWO,EPS,TR,TI,SR,SI,
     *                   XR,XI,YR,YI,ZZR,ZZI,FNORM
      COMPLEX            X,Y,Z
      EQUIVALENCE        (X,T1(1)),(T1(1),XR),(T1(2),XI),(Y,T2(1)),
     *                   (T2(1),YR),(T2(2),YI),(Z,T3(1)),(T3(1),ZZR),
     *                   (T3(2),ZZI)
      DATA               ZERO,ONE,TWO /0.0,1.0,2.0/
      DATA               EPS/Z3C100000/
C                                  INITIALIZE IER
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      INFER = 0
      TR = ZERO
      TI = ZERO
      DO 10 I=1,N
         DO 5 J=1,N
            ZR(I,J) = ZERO
            ZI(I,J) = ZERO
    5    CONTINUE
         ZR(I,I) = ONE
   10 CONTINUE
C                                  FORM THE MATRIX OF ACCUMULATED
C                                    TRANSFORMATIONS FROM THE INFOR-
C                                    MATION LEFT BY ROUTINE EHESSC
      IEND = L-K-1
      IF (IEND.LE.0) GO TO 30
C                                  FOR I=L-1 STEP -1 UNTIL K+1 DO
      DO 25 II=1,IEND
         I = L-II
         IP1 = I+1
         IM1 = I-1
         DO 15 M=IP1,L
            ZR(M,I) = HR(M,IM1)
            ZI(M,I) = HI(M,IM1)
   15    CONTINUE
         J = ID(I)
         IF (I.EQ.J) GO TO 25
         DO 20 M=I,L
            ZR(I,M) = ZR(J,M)
            ZI(I,M) = ZI(J,M)
            ZR(J,M) = ZERO
            ZI(J,M) = ZERO
   20    CONTINUE
         ZR(J,I) = ONE
   25 CONTINUE
   30 DO 35 I=1,N
         IF (I.GE.K .AND. I.LE.L) GO TO 35
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
   35 CONTINUE
      NN = L
C                                  SEARCH FOR NEXT EIGENVALUE
   40 IF (NN.LT.K) GO TO 155
      ITS = 0
      NNM1 = NN-1
      NNM2 = NN-2
      IF (NN.EQ.K) GO TO 55
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                    ELEMENT
   45 NPL = NN+K
C                                  FOR M=NN STEP -1 UNTIL K+1 DO
      DO 50 KK=K,NNM1
         M = NPL-KK
         MM1 = M-1
         IF (ABS(HR(M,MM1))+ABS(HI(M,MM1)).LE.EPS*(ABS(HR(MM1,MM1))
     *   +ABS(HI(MM1,MM1))+ABS(HR(M,M))+ABS(HI(M,M)))) GO TO 60
   50 CONTINUE
   55 M = K
   60 IF (M.EQ.NN) GO TO 150
      IF (ITS.EQ.30) GO TO 215
C                                  FORM SHIFT
      IF (ITS.EQ.10 .OR. ITS.EQ.20) GO TO 65
      SR = HR(NN,NN)
      SI = HI(NN,NN)
      XR = HR(NNM1,NN)*HR(NN,NNM1)-HI(NNM1,NN)*HI(NN,NNM1)
      XI = HR(NNM1,NN)*HI(NN,NNM1)+HI(NNM1,NN)*HR(NN,NNM1)
      IF (XR.EQ.ZERO .AND. XI.EQ.ZERO) GO TO 70
      YR = (HR(NNM1,NNM1)-SR)/TWO
      YI = (HI(NNM1,NNM1)-SI)/TWO
      Z = CSQRT(CMPLX(YR**2-YI**2+XR,TWO*YR*YI+XI))
      IF (YR*ZZR+YI*ZZI.LT.ZERO) Z = -Z
      X = X/(Y+Z)
      SR = SR-XR
      SI = SI-XI
      GO TO 70
   65 SR = ABS(HR(NN,NNM1))+ABS(HR(NNM1,NNM2))
      SI = ABS(HI(NN,NNM1))+ABS(HI(NNM1,NNM2))
   70 DO 75 I=K,NN
         HR(I,I) = HR(I,I)-SR
         HI(I,I) = HI(I,I)-SI
   75 CONTINUE
      TR = TR+SR
      TI = TI+SI
      ITS = ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                    SUB-DIAGONAL ELEMENTS
      XR = ABS(HR(NNM1,NNM1))+ABS(HI(NNM1,NNM1))
      YR = ABS(HR(NN,NNM1))+ABS(HI(NN,NNM1))
      ZZR = ABS(HR(NN,NN))+ABS(HI(NN,NN))
      NNMJ = NNM1-M
      IF (NNMJ.EQ.0) GO TO 85
C                                  FOR MM=NN-1 STEP -1 UNTIL M+1 DO
      DO 80 NM=1,NNMJ
         MM = NN-NM
         MMM1 = MM-1
         YI = YR
         YR = ABS(HR(MM,MMM1))+ABS(HI(MM,MMM1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(MMM1,MMM1))+ABS(HI(MMM1,MMM1))
         IF (YR.LE.(EPS*ZZR/YI)*(ZZR+XR+XI)) GO TO 90
   80 CONTINUE
   85 MM = M
C                                  TRIANGULAR DECOMPOSITION H=L*R
   90 MP1 = MM+1
      DO 115 I=MP1,NN
         IM1 = I-1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR)+ABS(XI).GE.ABS(YR)+ABS(YI)) GO TO 100
C                                  INTERCHANGE ROWS OF HR AND HI
         DO 95 J=IM1,N
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
   95    CONTINUE
         Z = X/Y
         WR(I) = ONE
         GO TO 105
  100    Z = Y/X
         WR(I) = -ONE
  105    HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
         DO 110 J=I,N
            HR(I,J) = HR(I,J)-ZZR*HR(IM1,J)+ZZI*HI(IM1,J)
            HI(I,J) = HI(I,J)-ZZR*HI(IM1,J)-ZZI*HR(IM1,J)
  110    CONTINUE
  115 CONTINUE
C                                  COMPOSITION R*L=H
      DO 145 J=MP1,NN
         JM1 = J-1
         XR = HR(J,JM1)
         XI = HI(J,JM1)
         HR(J,JM1) = ZERO
         HI(J,JM1) = ZERO
C                                  INTERCHANGE COLUMNS OF HR, HI,
C                                    ZR, AND ZI IF NECESSARY
         IF (WR(J).LE.ZERO) GO TO 130
         DO 120 I=1,J
            ZZR = HR(I,JM1)
            HR(I,JM1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,JM1)
            HI(I,JM1) = HI(I,J)
            HI(I,J) = ZZI
  120    CONTINUE
         DO 125 I=K,L
            ZZR = ZR(I,JM1)
            ZR(I,JM1) = ZR(I,J)
            ZR(I,J) = ZZR
            ZZI = ZI(I,JM1)
            ZI(I,JM1) = ZI(I,J)
            ZI(I,J) = ZZI
  125    CONTINUE
C                                  END INTERCHANGE COLUMNS
  130    DO 135 I=1,J
            HR(I,JM1) = HR(I,JM1)+XR*HR(I,J)-XI*HI(I,J)
            HI(I,JM1) = HI(I,JM1)+XR*HI(I,J)+XI*HR(I,J)
  135    CONTINUE
         DO 140 I=K,L
            ZR(I,JM1) = ZR(I,JM1)+XR*ZR(I,J)-XI*ZI(I,J)
            ZI(I,JM1) = ZI(I,JM1)+XR*ZI(I,J)+XI*ZR(I,J)
  140    CONTINUE
C                                  END ACCUMULATE TRANSFORMATIONS
  145 CONTINUE
      GO TO 45
C                                  A ROOT FOUND
  150 WR(NN) = HR(NN,NN)+TR
      WI(NN) = HI(NN,NN)+TI
      NN = NNM1
      GO TO 40
C                                  ALL ROOTS FOUND. BACKSUBSTITUTE TO
C                                    FIND VECTORS OF UPPER TRIANGULAR
C                                    FORM
  155 IF (N.EQ.1) GO TO 9005
      FNORM = ZERO
      DO 165 I=1,N
         FNORM = FNORM+ABS(WR(I))+ABS(WI(I))
         IF (I.EQ.N) GO TO 165
         IP1 = I+1
         DO 160 J=IP1,N
            FNORM = FNORM+ABS(HR(I,J))+ABS(HI(I,J))
  160    CONTINUE
  165 CONTINUE
      IF (FNORM.EQ.ZERO) GO TO 9005
C                                  FOR NN=N STEP -1 UNTIL 2 DO
      DO 185 NM=2,N
         NN = N+2-NM
         XR = WR(NN)
         XI = WI(NN)
         NNM1 = NN-1
C                                  FOR I=NN-1 STEP -1 UNTIL 1 DO
         DO 180 II=1,NNM1
            I = NN-II
            ZZR = HR(I,NN)
            ZZI = HI(I,NN)
            IF (I.EQ.NNM1) GO TO 175
            IP1 = I+1
            DO 170 J=IP1,NNM1
               ZZR = ZZR+HR(I,J)*HR(J,NN)-HI(I,J)*HI(J,NN)
               ZZI = ZZI+HR(I,J)*HI(J,NN)+HI(I,J)*HR(J,NN)
  170       CONTINUE
  175       YR = XR-WR(I)
            YI = XI-WI(I)
            IF (YR.EQ.ZERO .AND. YI.EQ.ZERO) YR = EPS*FNORM
            Z = Z/Y
            HR(I,NN) = T3(1)
            HI(I,NN) = T3(2)
  180    CONTINUE
  185 CONTINUE
C                                  END BACKSUBSTITUTION
      NM1 = N-1
C                                  VECTOR OF ISOLATED ROOTS
      DO 195 I=1,NM1
         IF (I.GE.K .AND. I.LE.L) GO TO 195
         IP1 = I+1
         DO 190 J=IP1,N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  190    CONTINUE
  195 CONTINUE
      IF (L.EQ.0) GO TO 9005
C                                  MULTIPLY OF ORIGINAL FULL MATRIX
      NPL = N+K
C                                  FOR J=N STEP -1 UNTIL K+1 DO
      DO 210 JJ=K,NM1
         J = NPL-JJ
         DO 205 I=K,L
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            MM = J-1
            IF (L.LT.J) MM = L
            DO 200 M=K,MM
               ZZR = ZZR+ZR(I,M)*HR(M,J)-ZI(I,M)*HI(M,J)
               ZZI = ZZI+ZR(I,M)*HI(M,J)+ZI(I,M)*HR(M,J)
  200       CONTINUE
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  205    CONTINUE
  210 CONTINUE
      GO TO 9005
C                                  SET ERROR - NO CONVERGENCE TO AN
C                                    EIGENVALUE AFTER 30 ITERATIONS
  215 IER = 129
      INFER = NN
 9000 CONTINUE
      CALL UERTST(IER,6HELRH2C)
 9005 RETURN
      END

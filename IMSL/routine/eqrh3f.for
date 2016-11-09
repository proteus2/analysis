C   IMSL ROUTINE NAME   - EQRH3F
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQRH3F (H,N,IH,K,L,WR,WI,Z,IZ,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IH,K,L,IZ,IER
      REAL               H(IH,N),WR(N),WI(N),Z(IZ,N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEN,ITS,IENM2,NPL,LL,LB,NAML,MM,M,MP2,KA,NA,
     *                   J,JJ
      REAL               T3(2),RDELP,P4,P5,P7,ZERO,ONE,T,X,Y,W,S,ZZ,R,P,
     *                   Q,RNORM,RA,SA,SCALE,VR,VI
      COMPLEX            Z3
      LOGICAL            NOTLAS
      EQUIVALENCE        (Z3,T3(1))
      DATA               RDELP/Z3C100000/
      DATA               P4 /0.4375/,P5 /0.5/,P7 /0.75/,ZERO /0.0/,ONE
     *                   /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  STORE ROOTS ISOLATED BY EBALAF
      RNORM = 0.0
      KA = 1
      DO 10 I=1,N
         DO 5 J=KA,N
    5    RNORM = RNORM+ABS(H(I,J))
         KA = I
         IF (I.GE.K .AND. I.LE.L) GO TO 10
         WR(I) = H(I,I)
         WI(I) = ZERO
   10 CONTINUE
      IEN = L
      T = ZERO
C                                  SEARCH FOR NEXT EIGENVALUES
   15 IF (IEN.LT.K) GO TO 145
      ITS = 0
      NA = IEN-1
      IENM2 = NA-1
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                  ELEMENT
   20 NPL = IEN+K
      DO 25 LL=K,IEN
         LB = NPL-LL
         IF (LB.EQ.K) GO TO 30
         S = ABS(H(LB-1,LB-1))+ABS(H(LB,LB))
         IF (S.EQ.0.0) S = RNORM
         IF (ABS(H(LB,LB-1)).LE.RDELP*S) GO TO 30
   25 CONTINUE
C
   30 X = H(IEN,IEN)
      IF (LB.EQ.IEN) GO TO 110
      Y = H(NA,NA)
      W = H(IEN,NA)*H(NA,IEN)
      IF (LB.EQ.NA) GO TO 115
      IF (ITS.EQ.30) GO TO 250
C                                  FORM SHIFT
      IF (ITS.NE.10 .AND. ITS.NE.20) GO TO 40
      T = T+X
      DO 35 I=K,IEN
         H(I,I) = H(I,I)-X
   35 CONTINUE
      S = ABS(H(IEN,NA))+ABS(H(NA,IENM2))
      X = P7*S
      Y = X
      W = -P4*S*S
   40 ITS = ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                  SUB-DIAGONAL ELEMENTS
      NAML = IENM2+LB
      DO 45 MM=LB,IENM2
         M = NAML-MM
         ZZ = H(M,M)
         R = X-ZZ
         S = Y-ZZ
         P = (R*S-W)/H(M+1,M)+H(M,M+1)
         Q = H(M+1,M+1)-ZZ-R-S
         R = H(M+2,M+1)
         S = ABS(P)+ABS(Q)+ABS(R)
         P = P/S
         Q = Q/S
         R = R/S
         IF (M.EQ.LB) GO TO 50
         IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)).LE.RDELP*ABS(P)*(ABS(H(M-1,
     *   M-1))+ABS(ZZ)+ABS(H(M+1,M+1)))) GO TO 50
   45 CONTINUE
   50 MP2 = M+2
      DO 55 I=MP2,IEN
         H(I,I-2) = ZERO
         IF (I.EQ.MP2) GO TO 55
         H(I,I-3) = ZERO
   55 CONTINUE
C                                  DOUBLE QR STEP INVOLVING ROWS
C                                  L TO EN AND COLUMNS M TO EN
      DO 105 KA=M,NA
         NOTLAS = KA.NE.NA
         IF (KA.EQ.M) GO TO 60
         P = H(KA,KA-1)
         Q = H(KA+1,KA-1)
         R = ZERO
         IF (NOTLAS) R = H(KA+2,KA-1)
         X = ABS(P)+ABS(Q)+ABS(R)
         IF (X.EQ.ZERO) GO TO 105
         P = P/X
         Q = Q/X
         R = R/X
   60    CONTINUE
         S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (KA.EQ.M) GO TO 65
         H(KA,KA-1) = -S*X
         GO TO 70
   65    IF (LB.NE.M) H(KA,KA-1) = -H(KA,KA-1)
   70    P = P+S
         X = P/S
         Y = Q/S
         ZZ = R/S
         Q = Q/P
         R = R/P
C                                  ROW MODIFICATION
         DO 80 J=KA,N
            P = H(KA,J)+Q*H(KA+1,J)
            IF (.NOT.NOTLAS) GO TO 75
            P = P+R*H(KA+2,J)
            H(KA+2,J) = H(KA+2,J)-P*ZZ
   75       H(KA+1,J) = H(KA+1,J)-P*Y
            H(KA,J) = H(KA,J)-P*X
   80    CONTINUE
         J = MIN0(IEN,KA+3)
C                                  COLUMN MODIFICATION
         DO 90 I=1,J
            P = X*H(I,KA)+Y*H(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 85
            P = P+ZZ*H(I,KA+2)
            H(I,KA+2) = H(I,KA+2)-P*R
   85       H(I,KA+1) = H(I,KA+1)-P*Q
            H(I,KA) = H(I,KA)-P
   90    CONTINUE
         IF (IZ.LT.N) GO TO 105
C                                  ACCUMULATE TRANSFORMATIONS
         DO 100 I=K,L
            P = X*Z(I,KA)+Y*Z(I,KA+1)
            IF (.NOT.NOTLAS) GO TO 95
            P = P+ZZ*Z(I,KA+2)
            Z(I,KA+2) = Z(I,KA+2)-P*R
   95       Z(I,KA+1) = Z(I,KA+1)-P*Q
            Z(I,KA) = Z(I,KA)-P
  100    CONTINUE
  105 CONTINUE
      GO TO 20
C                                  ONE ROOT FOUND
  110 H(IEN,IEN) = X+T
      WR(IEN) = H(IEN,IEN)
      WI(IEN) = ZERO
      IEN = NA
      GO TO 15
C                                  TWO ROOTS FOUND
  115 P = (Y-X)*P5
      Q = P*P+W
      ZZ = SQRT(ABS(Q))
      H(IEN,IEN) = X+T
      X = H(IEN,IEN)
      H(NA,NA) = Y+T
      IF (Q.LT.ZERO) GO TO 135
C                                  REAL PAIR
      ZZ = P+SIGN(ZZ,P)
      WR(NA) = X+ZZ
      WR(IEN) = WR(NA)
      IF (ZZ.NE.ZERO) WR(IEN) = X-W/ZZ
      WI(NA) = ZERO
      WI(IEN) = ZERO
      X = H(IEN,NA)
C                                  EMPLOY SCALE FACTOR IN CASE X AND
C                                  ZZ ARE VERY SMALL
      SCALE = ABS(X) + ABS(ZZ)
      R = SCALE * SQRT( (X/SCALE)**2 + (ZZ/SCALE)**2 )
      P = X/R
      Q = ZZ/R
C                                  ROW MODIFICATION
      DO 120 J=NA,N
         ZZ = H(NA,J)
         H(NA,J) = Q*ZZ+P*H(IEN,J)
         H(IEN,J) = Q*H(IEN,J)-P*ZZ
  120 CONTINUE
C                                  COLUMN MODIFICATION
      DO 125 I=1,IEN
         ZZ = H(I,NA)
         H(I,NA) = Q*ZZ+P*H(I,IEN)
         H(I,IEN) = Q*H(I,IEN)-P*ZZ
  125 CONTINUE
      IF (IZ.LT.N) GO TO 140
C                                  ACCUMULATE TRANSFORMATIONS
      DO 130 I=K,L
         ZZ = Z(I,NA)
         Z(I,NA) = Q*ZZ+P*Z(I,IEN)
         Z(I,IEN) = Q*Z(I,IEN)-P*ZZ
  130 CONTINUE
      GO TO 140
C                                  COMPLEX PAIR
  135 WR(NA) = X+P
      WR(IEN) = X+P
      WI(NA) = ZZ
      WI(IEN) = -ZZ
  140 IEN = IENM2
      GO TO 15
C                                  ALL ROOTS FOUND, NOW
C                                  BACKSUBSTITUTE
  145 IF (IZ.LT.N) GO TO 9005
      IF (RNORM.EQ.ZERO) GO TO 9005
      DO 220 NN=1,N
         IEN = N+1-NN
         P = WR(IEN)
         Q = WI(IEN)
         NA = IEN-1
         IF (Q.GT.ZERO) GO TO 220
         IF (Q.LT.ZERO) GO TO 180
C                                  REAL VECTOR
         M = IEN
         H(IEN,IEN) = ONE
         IF (NA.EQ.0) GO TO 220
         DO 175 II=1,NA
            I = IEN-II
            W = H(I,I)-P
            R = H(I,IEN)
            IF (M.GT.NA) GO TO 155
            DO 150 J=M,NA
               R = R+H(I,J)*H(J,IEN)
  150       CONTINUE
  155       IF (WI(I).GE.ZERO) GO TO 160
            ZZ = W
            S = R
            GO TO 175
  160       M = I
            IF (WI(I).NE.ZERO) GO TO 165
            T = W
            IF (W.EQ.ZERO) T = RDELP*RNORM
            H(I,IEN) = -R/T
            GO TO 175
C                                  SOLVE REAL EQUATIONS
  165       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
            T = (X*S-ZZ*R)/Q
            H(I,IEN) = T
            IF (ABS(X).LE.ABS(ZZ)) GO TO 170
            H(I+1,IEN) = (-R-W*T)/X
            GO TO 175
  170       H(I+1,IEN) = (-S-Y*T)/ZZ
  175    CONTINUE
C                                  END REAL VECTOR
         GO TO 220
C                                  LAST VECTOR COMPONENT CHOSEN
C                                    IMAGINARY SO THAT EIGENVECTOR
C                                    MATRIX IS TRIANGULAR
  180    M = NA
C                                  COMPLEX VECTOR
         IF (ABS(H(IEN,NA)).LE.ABS(H(NA,IEN))) GO TO 185
         H(NA,NA) = Q/H(IEN,NA)
         H(NA,IEN) = -(H(IEN,IEN)-P)/H(IEN,NA)
         GO TO 190
  185    CONTINUE
         Z3 = CMPLX(ZERO,-H(NA,IEN))/CMPLX(H(NA,NA)-P,Q)
         H(NA,NA) = T3(1)
         H(NA,IEN) = T3(2)
  190    H(IEN,NA) = ZERO
         H(IEN,IEN) = ONE
         IENM2 = NA-1
         IF (IENM2.EQ.0) GO TO 220
         DO 215 II=1,IENM2
            I = NA-II
            W = H(I,I)-P
            RA = ZERO
            SA = H(I,IEN)
            DO 195 J=M,NA
               RA = RA+H(I,J)*H(J,NA)
               SA = SA+H(I,J)*H(J,IEN)
  195       CONTINUE
            IF (WI(I).GE.ZERO) GO TO 200
            ZZ = W
            R = RA
            S = SA
            GO TO 215
  200       M = I
            IF (WI(I).NE.ZERO) GO TO 205
            Z3 = CMPLX(-RA,-SA)/CMPLX(W,Q)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            GO TO 215
C                                  SOLVE COMPLEX EQUATIONS
  205       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
            VI = (WR(I)-P)*Q
            VI = VI+VI
            IF (VR.EQ.ZERO .AND. VI.EQ.ZERO) VR = RDELP*RNORM*(ABS(W)
     *      +ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
            Z3 = CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/CMPLX(VR,VI)
            H(I,NA) = T3(1)
            H(I,IEN) = T3(2)
            IF (ABS(X).LE.ABS(ZZ)+ABS(Q)) GO TO 210
            H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,IEN))/X
            H(I+1,IEN) = (-SA-W*H(I,IEN)-Q*H(I,NA))/X
            GO TO 215
  210       CONTINUE
            Z3 = CMPLX(-R-Y*H(I,NA),-S-Y*H(I,IEN))/CMPLX(ZZ,Q)
            H(I+1,NA) = T3(1)
            H(I+1,IEN) = T3(2)
  215    CONTINUE
C                                  END COMPLEX VECTOR
  220 CONTINUE
C                                  END BACKSUBSTITUTION
C                                  VECTORS OF ISOLATED ROOTS
      DO 230 I=1,N
         IF (I.GE.K .AND. I.LE.L) GO TO 230
         DO 225 J=I,N
            Z(I,J) = H(I,J)
  225    CONTINUE
  230 CONTINUE
      IF (L.EQ.0) GO TO 9005
C                                  MULTIPLY BY TRANSFORMATION MATRIX
      DO 245 JJ=K,N
         J = N+K-JJ
         M = MIN0(J,L)
         DO 240 I=K,L
            ZZ = ZERO
            DO 235 KA=K,M
               ZZ = ZZ+Z(I,KA)*H(KA,J)
  235       CONTINUE
            Z(I,J) = ZZ
  240    CONTINUE
  245 CONTINUE
      GO TO 9005
C                                  NO CONVERGENCE AFTER 30 ITERATIONS
C                                  SET ERROR INDICATOR  TO THE INDEX
C                                  OF THE CURRENT EIGENVALUE
  250 IER = 128+IEN
      DO 255 I=1,IEN
         WR(I) = ZERO
         WI(I) = ZERO
  255 CONTINUE
      IF (IZ.LT.N) GO TO 9000
      DO 265 I=1,N
         DO 260 J=1,N
            Z(I,J) = ZERO
  260    CONTINUE
  265 CONTINUE
 9000 CONTINUE
      CALL UERTST (IER,6HEQRH3F)
 9005 RETURN
      END

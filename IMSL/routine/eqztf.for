
 
C   IMSL ROUTINE NAME   - EQZTF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGZF
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,VHSH2R,VHSH3R
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQZTF  (A,IA,B,IB,N,EPSA,EPSB,Z,IZ,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,N,IZ,IER
      REAL               A(IA,N),B(IB,N),Z(IZ,1),EPSA,EPSB
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,M,ITER,LB,L,L1,M1,LOR1,MORN,K,K1,K2,
     *                   K3,KM1
      REAL               T,V1,A10,A21,A34,B11,B34,OLD1,BNORM,
     1                   U1,V2,A11,A22,A43,B12,B44,OLD2,CONST,
     2                   U2,V3,A12,A30,A44,B22,EPS,ZERO,ONE,
     3                   U3,ANI,A20,A33,BNI,B33,ANORM
      LOGICAL            WANTX,MID
      DATA               EPS/Z3C100000/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      WANTX = .FALSE.
      IF (IZ.GE.N) WANTX = .TRUE.
C                                  INITIALIZE ITER. COMPUTE EPSA,EPSB
      ANORM = ZERO
      BNORM = ZERO
      DO 10 I=1,N
         ANI = ZERO
         IF (I.NE.1) ANI = ABS(A(I,I-1))
         BNI = ZERO
         DO 5 J=I,N
            ANI = ANI+ABS(A(I,J))
            BNI = BNI+ABS(B(I,J))
    5    CONTINUE
         IF (ANI.GT.ANORM) ANORM = ANI
         IF (BNI.GT.BNORM) BNORM = BNI
   10 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = EPS
      IF (BNORM.EQ.ZERO) BNORM = EPS
      EPSA = EPS*ANORM
      EPSB = EPS*BNORM
C                                  REDUCE A TO QUASI-TRIANGULAR,
C                                     KEEP B TRIANGULAR
      M = N
      ITER = 0
   15 IF (M.LE.2) GO TO 110
C                                  CHECK FOR CONVERGENCE OR REDUCIBILITY
      DO 20 LB=1,M
         L = M+1-LB
         IF (L.EQ.1) GO TO 30
         IF (ABS(A(L,L-1)).LE.EPSA) GO TO 25
   20 CONTINUE
   25 A(L,L-1) = ZERO
      IF (L.LT.M-1) GO TO 30
      M = L-1
      ITER = 0
      GO TO 15
C                                  CHECK FOR SMALL TOP OF B
   30 IF (ABS(B(L,L)).GT.EPSB) GO TO 45
      B(L,L) = ZERO
      L1 = L+1
      CALL VHSH2R (A(L,L),A(L1,L),U1,U2,V1,V2)
      IF (U1.NE.ONE) GO TO 40
      DO 35 J=L,N
         T = A(L,J)+U2*A(L1,J)
         A(L,J) = A(L,J)+T*V1
         A(L1,J) = A(L1,J)+T*V2
         T = B(L,J)+U2*B(L1,J)
         B(L,J) = B(L,J)+T*V1
         B(L1,J) = B(L1,J)+T*V2
   35 CONTINUE
   40 L = L1
      GO TO 25
C                                  BEGIN ONE QZ STEP. ITERATION STRATEGY
   45 M1 = M-1
      L1 = L+1
      CONST = 0.75
      ITER = ITER+1
      IF (ITER.EQ.1) GO TO 50
      IF (ABS(A(M,M-1)).LT.CONST*OLD1) GO TO 50
      IF (ABS(A(M-1,M-2)).LT.CONST*OLD2) GO TO 50
      IF (ITER.EQ.10) GO TO 55
      IF (ITER.GT.30) GO TO 105
C                                  ZEROTH COLUMN OF A
   50 B11 = B(L,L)
      B22 = B(L1,L1)
      IF (ABS(B22).LT.EPSB) B22 = EPSB
      B33 = B(M1,M1)
      IF (ABS(B33).LT.EPSB) B33 = EPSB
      B44 = B(M,M)
      IF (ABS(B44).LT.EPSB) B44 = EPSB
      A11 = A(L,L)/B11
      A12 = A(L,L1)/B22
      A21 = A(L1,L)/B11
      A22 = A(L1,L1)/B22
      A33 = A(M1,M1)/B33
      A34 = A(M1,M)/B44
      A43 = A(M,M1)/B33
      A44 = A(M,M)/B44
      B12 = B(L,L1)/B22
      B34 = B(M1,M)/B44
      A10 = ((A33-A11)*(A44-A11)-A34*A43+A43*B34*A11)/A21+A12-A11*B12
      A20 = (A22-A11-A21*B12)-(A33-A11)-(A44-A11)+A43*B34
      A30 = A(L+2,L1)/B22
      GO TO 60
C
C                                  AD HOC SHIFT
   55 A10 = ZERO
      A20 = ONE
      A30 = 1.1605
   60 OLD1 = ABS(A(M,M-1))
      OLD2 = ABS(A(M-1,M-2))
      IF (.NOT.WANTX) LOR1 = L
      IF (WANTX) LOR1 = 1
      IF (.NOT.WANTX) MORN = M
      IF (WANTX) MORN = N
C                                  BEGIN MAIN LOOP
      DO 100 K=L,M1
         MID = K.NE.M1
         K1 = K+1
         K2 = K+2
         K3 = K+3
         IF (K3.GT.M) K3 = M
         KM1 = K-1
         IF (KM1.LT.L) KM1 = L
C                                  ZERO A(K+1,K-1) AND A(K+2,K-1)
         IF (K.EQ.L) CALL VHSH3R (A10,A20,A30,U1,U2,U3,V1,V2,V3)
         IF (K.GT.L.AND.K.LT.M1) CALL VHSH3R (A(K,KM1),A(K1,KM1),A(K2,KM
     1   1),U1,U2,U3,V1,V2,V3)
         IF (K.EQ.M1) CALL VHSH2R (A(K,KM1),A(K1,KM1),U1,U2,V1,V2)
         IF (U1.NE.ONE) GO TO 70
         DO 65 J=KM1,MORN
            T = A(K,J)+U2*A(K1,J)
            IF (MID) T = T+U3*A(K2,J)
            A(K,J) = A(K,J)+T*V1
            A(K1,J) = A(K1,J)+T*V2
            IF (MID) A(K2,J) = A(K2,J)+T*V3
            T = B(K,J)+U2*B(K1,J)
            IF (MID) T = T+U3*B(K2,J)
            B(K,J) = B(K,J)+T*V1
            B(K1,J) = B(K1,J)+T*V2
            IF (MID) B(K2,J) = B(K2,J)+T*V3
   65    CONTINUE
         IF (K.EQ.L) GO TO 70
         A(K1,K-1) = ZERO
         IF (MID) A(K2,K-1) = ZERO
C                                  ZERO B(K+2,K+1) AND B(K+2,K)
   70    IF (K.EQ.M1) GO TO 85
         CALL VHSH3R (B(K2,K2),B(K2,K1),B(K2,K),U1,U2,U3,V1,V2,V3)
         IF (U1.NE.ONE) GO TO 85
         DO 75 I=LOR1,K3
            T = A(I,K2)+U2*A(I,K1)+U3*A(I,K)
            A(I,K2) = A(I,K2)+T*V1
            A(I,K1) = A(I,K1)+T*V2
            A(I,K) = A(I,K)+T*V3
            T = B(I,K2)+U2*B(I,K1)+U3*B(I,K)
            B(I,K2) = B(I,K2)+T*V1
            B(I,K1) = B(I,K1)+T*V2
            B(I,K) = B(I,K)+T*V3
   75    CONTINUE
         B(K2,K) = ZERO
         B(K2,K1) = ZERO
         IF (.NOT.WANTX) GO TO 85
         DO 80 I=1,N
            T = Z(I,K2)+U2*Z(I,K1)+U3*Z(I,K)
            Z(I,K2) = Z(I,K2)+T*V1
            Z(I,K1) = Z(I,K1)+T*V2
            Z(I,K) = Z(I,K)+T*V3
   80    CONTINUE
C                                  ZERO B(K+1,K)
   85    CALL VHSH2R (B(K1,K1),B(K1,K),U1,U2,V1,V2)
         IF (U1.NE.ONE) GO TO 100
         DO 90 I=LOR1,K3
            T = A(I,K1)+U2*A(I,K)
            A(I,K1) = A(I,K1)+T*V1
            A(I,K) = A(I,K)+T*V2
            T = B(I,K1)+U2*B(I,K)
            B(I,K1) = B(I,K1)+T*V1
            B(I,K) = B(I,K)+T*V2
   90    CONTINUE
         B(K1,K) = ZERO
         IF (.NOT.WANTX) GO TO 100
         DO 95 I=1,N
            T = Z(I,K1)+U2*Z(I,K)
            Z(I,K1) = Z(I,K1)+T*V1
            Z(I,K) = Z(I,K)+T*V2
   95    CONTINUE
C                                  END MAIN LOOP
  100 CONTINUE
C                                  END ONE QZ STEP
      GO TO 15
  105 IER = 128+M
      CALL UERTST (IER,6HEQZTF )
  110 RETURN
      END
 
R; T=0.06/0.50 23:29:09
TEP
      GO TO 15
  105 IER = 128+M
      CALL UERTST (IER,6HEQZTF )
  110 R
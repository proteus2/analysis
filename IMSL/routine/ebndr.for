C   IMSL ROUTINE NAME   - EBNDR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE EIGBS
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBNDR  (A,IA,N,NC,IOPT,D,E,E2,Z,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,N,NC,IOPT,IZ
      REAL               A(IA,1),D(N),E(N),E2(N),Z(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,K,L,R,I1,I2,J1,J2,KR,MR,M1,N2,R1,UGL,MAXL,
     1                   MAXR,MB
      REAL               G,U,B1,B2,C2,F1,F2,S2,DMIN,DMINRT
C                                  FIRST EXECUTABLE STATEMENT
      DMIN = 2.0**(-64)
      DMINRT = 2.0**(-32)
C                                  INITIALIZE DIAGONAL SCALING MATRIX
      DO 5 J=1,N
    5 D(J) = 1.0
      IF (IOPT.NE.1) GO TO 20
      DO 15 J=1,N
         DO 10 K=1,N
   10    Z(J,K) = 0.0
         Z(J,J) = 1.0
   15 CONTINUE
   20 MB = NC+1
      M1 = NC
      IF (M1-1) 150,120,25
   25 N2 = N-2
      DO 115 K=1,N2
         MAXR = MIN0(M1,N-K)
C                                  FOR R=MAXR STEP -1 UNTIL 2 DO --
         DO 80 R1=2,MAXR
            R = MAXR+2-R1
            KR = K+R
            MR = MB-R
            G = A(KR,MR)
            A(KR-1,1) = A(KR-1,MR+1)
            UGL = K
            DO 75 J=KR,N,M1
               J1 = J-1
               J2 = J1-1
               IF (G.EQ.0.0) GO TO 80
               B1 = A(J1,1)/G
               B2 = B1*D(J1)/D(J)
               S2 = 1.0/(1.0+B1*B2)
               IF (S2.GE.0.5) GO TO 50
               B1 = G/A(J1,1)
               B2 = B1*D(J)/D(J1)
               C2 = 1.0-S2
               D(J1) = C2*D(J1)
               D(J) = C2*D(J)
               F1 = 2.0*A(J,M1)
               F2 = B1*A(J1,MB)
               A(J,M1) = -B2*(B1*A(J,M1)-A(J,MB))-F2+A(J,M1)
               A(J1,MB) = B2*(B2*A(J,MB)+F1)+A(J1,MB)
               A(J,MB) = B1*(F2-F1)+A(J,MB)
               DO 30 L=UGL,J2
                  I2 = MB-J+L
                  U = A(J1,I2+1)+B2*A(J,I2)
                  A(J,I2) = -B1*A(J1,I2+1)+A(J,I2)
                  A(J1,I2+1) = U
   30          CONTINUE
               UGL = J
               A(J1,1) = A(J1,1)+B2*G
               IF (J.EQ.N) GO TO 40
               MAXL = MIN0(M1,N-J1)
               DO 35 L=2,MAXL
                  I1 = J1+L
                  I2 = MB-L
                  U = A(I1,I2)+B2*A(I1,I2+1)
                  A(I1,I2+1) = -B1*A(I1,I2)+A(I1,I2+1)
                  A(I1,I2) = U
   35          CONTINUE
               I1 = J+M1
               IF (I1.GT.N) GO TO 40
               G = B2*A(I1,1)
   40          IF (IOPT.NE.1) GO TO 75
               DO 45 L=1,N
                  U = Z(L,J1)+B2*Z(L,J)
                  Z(L,J) = -B1*Z(L,J1)+Z(L,J)
                  Z(L,J1) = U
   45          CONTINUE
               GO TO 75
   50          U = D(J1)
               D(J1) = S2*D(J)
               D(J) = S2*U
               F1 = 2.0*A(J,M1)
               F2 = B1*A(J,MB)
               U = B1*(F2-F1)+A(J1,MB)
               A(J,M1) = B2*(B1*A(J,M1)-A(J1,MB))+F2-A(J,M1)
               A(J1,MB) = B2*(B2*A(J1,MB)+F1)+A(J,MB)
               A(J,MB) = U
               DO 55 L=UGL,J2
                  I2 = MB-J+L
                  U = B2*A(J1,I2+1)+A(J,I2)
                  A(J,I2) = -A(J1,I2+1)+B1*A(J,I2)
                  A(J1,I2+1) = U
   55          CONTINUE
               UGL = J
               A(J1,1) = B2*A(J1,1)+G
               IF (J.EQ.N) GO TO 65
               MAXL = MIN0(M1,N-J1)
               DO 60 L=2,MAXL
                  I1 = J1+L
                  I2 = MB-L
                  U = B2*A(I1,I2)+A(I1,I2+1)
                  A(I1,I2+1) = -A(I1,I2)+B1*A(I1,I2+1)
                  A(I1,I2) = U
   60          CONTINUE
               I1 = J+M1
               IF (I1.GT.N) GO TO 65
               G = A(I1,1)
               A(I1,1) = B1*A(I1,1)
   65          IF (IOPT.NE.1) GO TO 75
               DO 70 L=1,N
                  U = B2*Z(L,J1)+Z(L,J)
                  Z(L,J) = -Z(L,J1)+B1*Z(L,J)
                  Z(L,J1) = U
   70          CONTINUE
   75       CONTINUE
   80    CONTINUE
         IF (MOD(K,64).NE.0) GO TO 115
C                                  RESCALE TO AVOID UNDERFLOW OR
C                                    OVERFLOW
         DO 110 J=K,N
            IF (D(J).GE.DMIN) GO TO 110
            MAXL = MAX0(1,MB+1-J)
            DO 85 L=MAXL,M1
   85       A(J,L) = DMINRT*A(J,L)
            IF (J.EQ.N) GO TO 95
            MAXL = MIN0(M1,N-J)
            DO 90 L=1,MAXL
               I1 = J+L
               I2 = MB-L
               A(I1,I2) = DMINRT*A(I1,I2)
   90       CONTINUE
   95       IF (IOPT.NE.1) GO TO 105
            DO 100 L=1,N
  100       Z(L,J) = DMINRT*Z(L,J)
  105       A(J,MB) = DMIN*A(J,MB)
            D(J) = D(J)/DMIN
  110    CONTINUE
  115 CONTINUE
C                                  FORM SQUARE ROOT OF SCALING MATRIX
  120 DO 125 J=2,N
  125 E(J) = SQRT(D(J))
      IF (IOPT.NE.1) GO TO 140
      DO 135 J=1,N
         DO 130 K=2,N
  130    Z(J,K) = E(K)*Z(J,K)
  135 CONTINUE
  140 U = 1.0
      DO 145 J=2,N
         A(J,M1) = U*E(J)*A(J,M1)
         U = E(J)
         E2(J) = A(J,M1)**2
         A(J,MB) = D(J)*A(J,MB)
         D(J) = A(J,MB)
         E(J) = A(J,M1)
  145 CONTINUE
      D(1) = A(1,MB)
      E(1) = 0.0
      E2(1) = 0.0
      GO TO 160
  150 DO 155 J=1,N
         D(J) = A(J,MB)
         E(J) = 0.0
         E2(J) = 0.0
  155 CONTINUE
  160 RETURN
      END


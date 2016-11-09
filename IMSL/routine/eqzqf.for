C   IMSL ROUTINE NAME   - EQZQF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGZF
C
C   REQD. IMSL ROUTINES - VHSH2R
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EQZQF  (A,IA,B,IB,N,Z,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,IZ
      REAL               A(IA,N),B(IB,N),Z(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NM1,L
      REAL               S,U2,T,V1,V2,RHO,R,U1,SD,ZERO,ONE
      LOGICAL            WANTX
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      WANTX = .FALSE.
      IF (IZ.LT.N) GO TO 15
      WANTX = .TRUE.
C                                  INITIALIZE Z, USED TO SAVE
C                                     TRANSFORMATIONS
      DO 10 I=1,N
         DO 5 J=1,N
            Z(I,J) = ZERO
    5    CONTINUE
         Z(I,I) = ONE
   10 CONTINUE
C                                  REDUCE B TO UPPER TRIANGULAR FORM
   15 NM1 = N-1
      IF (N.LE.1) GO TO 110
      DO 65 L=1,NM1
         L1 = L+1
         S = ZERO
         DO 20 I=L1,N
            S = S+ABS(B(I,L))
   20    CONTINUE
         IF (S.EQ.ZERO) GO TO 65
         S = S+ABS(B(L,L))
         R = ZERO
         SD = ONE/S
         DO 25 I=L,N
            B(I,L) = B(I,L)*SD
            R = R+B(I,L)**2
   25    CONTINUE
         R = SQRT(R)
         IF (B(L,L).LT.ZERO) R = -R
         B(L,L) = B(L,L)+R
         RHO = R*B(L,L)
         SD = ONE/RHO
         DO 40 J=L1,N
            T = ZERO
            DO 30 I=L,N
               T = T+B(I,L)*B(I,J)
   30       CONTINUE
            T = -T*SD
            DO 35 I=L,N
               B(I,J) = B(I,J)+T*B(I,L)
   35       CONTINUE
   40    CONTINUE
         DO 55 J=1,N
            T = ZERO
            DO 45 I=L,N
               T = T+B(I,L)*A(I,J)
   45       CONTINUE
            T = -T*SD
            DO 50 I=L,N
               A(I,J) = A(I,J)+T*B(I,L)
   50       CONTINUE
   55    CONTINUE
         B(L,L) = -S*R
         DO 60 I=L1,N
            B(I,L) = ZERO
   60    CONTINUE
   65 CONTINUE
C                                  REDUCE A TO UPPER HESSENBERG,
C                                  KEEP B TRIANGULAR
      IF (N.LE.2) GO TO 110
      NM2 = N-2
      DO 105 K=1,NM2
         K1 = K+1
         NK1 = N-K1
         DO 100 LB=1,NK1
            L = N-LB
            L1 = L+1
            CALL VHSH2R (A(L,K),A(L1,K),U1,U2,V1,V2)
            IF (U1.NE.ONE) GO TO 80
            DO 70 J=K,N
               T = A(L,J)+U2*A(L1,J)
               A(L,J) = A(L,J)+T*V1
               A(L1,J) = A(L1,J)+T*V2
   70       CONTINUE
            A(L1,K) = ZERO
            DO 75 J=L,N
               T = B(L,J)+U2*B(L1,J)
               B(L,J) = B(L,J)+T*V1
               B(L1,J) = B(L1,J)+T*V2
   75       CONTINUE
   80       CALL VHSH2R (B(L1,L1),B(L1,L),U1,U2,V1,V2)
            IF (U1.NE.ONE) GO TO 100
            DO 85 I=1,L1
               T = B(I,L1)+U2*B(I,L)
               B(I,L1) = B(I,L1)+T*V1
               B(I,L) = B(I,L)+T*V2
   85       CONTINUE
            B(L1,L) = ZERO
            DO 90 I=1,N
               T = A(I,L1)+U2*A(I,L)
               A(I,L1) = A(I,L1)+T*V1
               A(I,L) = A(I,L)+T*V2
   90       CONTINUE
            IF (.NOT.WANTX) GO TO 100
            DO 95 I=1,N
               T = Z(I,L1)+U2*Z(I,L)
               Z(I,L1) = Z(I,L1)+T*V1
               Z(I,L) = Z(I,L)+T*V2
   95       CONTINUE
  100    CONTINUE
  105 CONTINUE
  110 RETURN
      END

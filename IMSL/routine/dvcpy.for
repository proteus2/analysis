C   IMSL ROUTINE NAME   - DVCPY
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DVCPR
C
C
C   REQD. IMSL ROUTINES - NONE
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPY (A,C,DEL,Y,M,N,P,R,IR,IC,U,MTNMAX,MMAX,NMAX,X)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER            M,N,P,R,MTNMAX,MMAX,NMAX,IR(NMAX,MMAX),IC(NMAX,
     *                   MMAX)
      REAL               A(MTNMAX,MMAX),C(MTNMAX,MMAX),DEL(MMAX,MTNMAX),
     *                   Y(MTNMAX),U(MTNMAX),X(MTNMAX)
C
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            I0J,I0,I11,I1,I2,II,IXN,IXP,IX,I,J1,JM,J,K1,K,
     *                   M1,MN,N1M,N1,N2,P1,P2
      REAL               TE
C                                  FIRST EXECUTABLE STATEMENT
      M1 = M-1
      N2 = N+1
      N1 = N-1
      P1 = P-1
      P2 = P+1
C                                  ROW INTERCHANGES ON THE RIGHT HAND
C                                    SIDE
      IF (P.EQ.0) GO TO 10
      DO 5 I=1,P
    5 X(I) = Y(I)
   10 MN = M*N1
      DO 15 I=1,M
   15 X(MN+I) = Y(MN+I)
      DO 25 I=1,N1
         IX = (I-1)*M+P
         DO 20 J=1,M
            IXP = IX+IR(I,J)
            X(IX+J) = Y(IXP)
   20    CONTINUE
   25 CONTINUE
C                                  SOLVE L * Y = X
      DO 40 I=2,N
         I1 = M*(I-2)+1
         I2 = I1+P1
         I11 = I1-1
         IF (P.EQ.0) GO TO 45
         DO 35 J=I1,I2
            JM = J+M
            TE = X(JM)
            DO 30 K=1,M
   30       TE = TE-C(J,K)*X(I11+K)
            X(JM) = TE
   35    CONTINUE
   40 CONTINUE
C
   45 N1M = N1*M
      IXN = N1M+P
      IF (R.EQ.0) GO TO 60
      DO 55 II=1,R
         I1 = IXN+II
         TE = X(I1)
         DO 50 J=1,N1M
   50    TE = TE-DEL(II,J)*X(J)
         X(I1) = TE
   55 CONTINUE
C                                  SOLVE U * Z = Y
   60 DO 100 I=1,N
         II = N2-I
         I0 = (II-1)*M
         I1 = I0+2
         I2 = II*M
         IF (I.EQ.1) GO TO 75
         DO 70 J=P2,M
            I0J = I0+J
            TE = X(I0J)
            DO 65 K=1,M
   65       TE = TE-C(I0J,K)*X(I2+K)
            X(I0J) = TE
   70    CONTINUE
   75    DO 85 J=I1,I2
            K1 = J-I1+1
            TE = X(J)
            DO 80 K=1,K1
   80       TE = TE-A(J,K)*X(I0+K)
            X(J) = TE
   85    CONTINUE
         X(I2) = X(I2)/A(I2,M)
         DO 95 J=1,M1
            J1 = I2-J
            TE = X(J1)
            K1 = M-J+1
            DO 90 K=K1,M
   90       TE = TE-A(J1,K)*X(I0+K)
            X(J1) = TE/A(J1,J1-I0)
   95    CONTINUE
  100 CONTINUE
C                                  INTERCHANGES IN X
      DO 110 I=1,N
         IX = (I-1)*M
         DO 105 J=1,M
            IXP = IX+IC(I,J)
            U(IXP) = X(IX+J)
  105    CONTINUE
  110 CONTINUE
C
      RETURN
      END

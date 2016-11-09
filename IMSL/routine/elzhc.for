C   IMSL ROUTINE NAME   - ELZHC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGZC
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ELZHC  (N,A,IA,B,IB,IJOB,Z,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IB,IJOB,IZ
      COMPLEX            A(IA,N),B(IB,N),Z(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NM1,NP1,I,IP1,K,II,J,NM2,JM2,JP1,IM1,IMJ
      REAL               C,D,ZRO
      COMPLEX            Y,W,ZC,ONE,ZERO
      DATA               ZERO/(0.0,0.0)/,ONE/(1.0,0.0)/
      DATA               ZRO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF(N.EQ.1) GO TO 135
      NM1 = N-1
      NP1 = N+1
C                                  REDUCE B TO TRIANGULAR FORM USING
C                                  ELEMENTARY TRANSFORMATIONS
      DO 40 I=1,NM1
         D = ZRO
         IP1 = I+1
         DO 5 K=IP1,N
            Y = B(K,I)
            C = ABS(REAL(Y))+ABS(AIMAG(Y))
            IF (C.LE.D) GO TO 5
            D = C
            II = K
    5    CONTINUE
         IF (D.EQ.ZRO) GO TO 40
         Y = B(I,I)
         IF (D.LE.ABS(REAL(Y))+ABS(AIMAG(Y))) GO TO 20
C                                  MUST INTERCHANGE
         DO 10 J=1,N
            Y = A(I,J)
            A(I,J) = A(II,J)
            A(II,J) = Y
   10    CONTINUE
         DO 15 J=I,N
            Y = B(I,J)
            B(I,J) = B(II,J)
            B(II,J) = Y
   15    CONTINUE
   20    DO 35 J=IP1,N
            Y = B(J,I)/B(I,I)
            IF (REAL(Y).EQ.ZRO.AND.AIMAG(Y).EQ.ZRO) GO TO 35
            DO 25 K=1,N
               A(J,K) = A(J,K)-Y*A(I,K)
   25       CONTINUE
            DO 30 K=IP1,N
               B(J,K) = B(J,K)-Y*B(I,K)
   30       CONTINUE
   35    CONTINUE
         B(IP1,I) = ZERO
   40 CONTINUE
C                                  INITIALIZE Z
      IF (IJOB.NE.1) GO TO 55
      DO 50 I=1,N
         DO 45 J=1,N
            Z(I,J) = ZERO
   45    CONTINUE
         Z(I,I) = ONE
   50 CONTINUE
C                                  REDUCE A TO UPPER HESSENBERG FORM
   55 NM2 = N-2
      IF (NM2.LT.1) GO TO 135
      DO 130 J=1,NM2
         JM2 = NM1-J
         JP1 = J+1
         DO 125 II=1,JM2
            I = NP1-II
            IM1 = I-1
            IMJ = I-J
            W = A(I,J)
            ZC = A(IM1,J)
            IF (ABS(REAL(W))+ABS(AIMAG(W)).LE.ABS(REAL(ZC))+
     *      ABS(AIMAG(ZC))) GO TO 70
C                                  MUST INTERCHANGE ROWS
            DO 60 K=J,N
               Y = A(I,K)
               A(I,K) = A(IM1,K)
               A(IM1,K) = Y
   60       CONTINUE
            DO 65 K=IM1,N
               Y = B(I,K)
               B(I,K) = B(IM1,K)
               B(IM1,K) = Y
   65       CONTINUE
   70       ZC = A(I,J)
            IF (REAL(ZC).EQ.ZRO.AND.AIMAG(ZC).EQ.ZRO) GO TO 85
            Y = ZC/A(IM1,J)
            DO 75 K=JP1,N
               A(I,K) = A(I,K)-Y*A(IM1,K)
   75       CONTINUE
            DO 80 K=IM1,N
               B(I,K) = B(I,K)-Y*B(IM1,K)
   80       CONTINUE
C                                  TRANSFORMATION FROM THE RIGHT
   85       W = B(I,IM1)
            ZC = B(I,I)
            IF (ABS(REAL(W))+ABS(AIMAG(W)).LE.ABS(REAL(ZC))+
     *      ABS(AIMAG(ZC))) GO TO 105
C                                  MUST INTERCHANGE COLUMNS
            DO 90 K=1,I
               Y = B(K,I)
               B(K,I) = B(K,IM1)
               B(K,IM1) = Y
   90       CONTINUE
            DO 95 K=1,N
               Y = A(K,I)
               A(K,I) = A(K,IM1)
               A(K,IM1) = Y
   95       CONTINUE
            IF (IJOB.NE.1) GO TO 105
            DO 100 K=IMJ,N
               Y = Z(K,I)
               Z(K,I) = Z(K,IM1)
               Z(K,IM1) = Y
  100       CONTINUE
  105       ZC = B(I,IM1)
            IF (REAL(ZC).EQ.ZRO.AND.AIMAG(ZC).EQ.ZRO) GO TO 125
            Y = ZC/B(I,I)
            DO 110 K=1,IM1
               B(K,IM1) = B(K,IM1)-Y*B(K,I)
  110       CONTINUE
            B(I,IM1) = ZERO
            DO 115 K=1,N
               A(K,IM1) = A(K,IM1)-Y*A(K,I)
  115       CONTINUE
            IF (IJOB.NE.1) GO TO 125
            DO 120 K=IMJ,N
               Z(K,IM1) = Z(K,IM1)-Y*Z(K,I)
  120       CONTINUE
  125    CONTINUE
         A(JP1+1,J) = ZERO
  130 CONTINUE
  135 RETURN
      END

C   IMSL ROUTINE NAME   - DVCPT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DVCPR
C
C   REQD. IMSL ROUTINES - DVCPU
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPT (K,P,Q,N,M,A,X,Y,S,IERROR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,P,Q,N,M,IERROR
      REAL               A(1),X(1),Y(1),S(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            I1,IAD,II1,II,IM,ITE,IT,I,J,KK1,KK,KMID1,
     *                   KMIDP1,KMID,L,NF
      REAL               ACUM,C(50)
C                                  ERROR EXIT
C                                  FIRST EXECUTABLE STATEMENT
      IF (K.GT.(N+1-Q)/P .OR. P.LT.1 .OR. K.LT.1) GO TO 65
      IF (Q.EQ.0) GO TO 10
      DO 5 I=1,Q
    5 A(I) = 0.
   10 KK1 = Q+P*K
      KK = KK1-1
      KMID = KK1/2
      IERROR = 0
      KMID1 = KMID-1
C                                  UNSYMMETRIC APPROXIMATION LEFT
C                                    BOUNDARY
      IF (KMID1.LT.1) GO TO 30
      DO 25 I=1,KMID1
         ITE = I
         CALL DVCPU(ITE,KK1,ITE,C,A,X,.5*(X(I+1)+X(I)))
         IM = (I-1)*M
         DO 20 L=1,M
            ACUM = 0.
            DO 15 J=1,KK1
   15       ACUM = ACUM+C(J)*Y((J-1)*M+L)
            S(IM+L) = ACUM
   20    CONTINUE
   25 CONTINUE
C                                  CENTER RANGE
   30 NF = N+1-KK1+KMID
      DO 45 I=KMID,NF
         ITE = I
         CALL DVCPU(ITE,KK1,KMID,C,A,X,.5*(X(I+1)+X(I)))
         I1 = I-1
         II = I1-KMID
         IT = I1*M
         DO 40 L=1,M
            ACUM = 0.
            DO 35 J=1,KK1
   35       ACUM = ACUM+C(J)*Y((II+J)*M+L)
            S(IT+L) = ACUM
   40    CONTINUE
   45 CONTINUE
C                                  RIGHT BOUNDARY
      KMIDP1 = KMID+1
      II = N-KK
      II1 = II-1
      DO 60 I=KMIDP1,KK
         IAD = II+I
         ITE = I
         CALL DVCPU(IAD,KK1,ITE,C,A,X,.5*(X(IAD+1)+X(IAD)))
         IT = (II1+I)*M
         DO 55 L=1,M
            ACUM = 0.
            DO 50 J=1,KK1
   50       ACUM = ACUM+C(J)*Y((II1+J)*M+L)
            S(IT+L) = ACUM
   55    CONTINUE
   60 CONTINUE
C                                  REGULAR EXIT
      RETURN
   65 IERROR = 1
      RETURN
      END

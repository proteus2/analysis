C   IMSL ROUTINE NAME   - EHESSC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCC
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHESSC (AR,AI,K,L,N,IA,ID)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,N,IA
      REAL               AR(IA,1),AI(IA,1),ID(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LA,KP1,M,I,J,MM1,MP1
      REAL               XR,XI,YR,YI,T1(2),T2(2),ZERO
      COMPLEX            X,Y
      EQUIVALENCE        (X,T1(1)),(T1(1),XR),(T1(2),XI),(Y,T2(1)),
     1                   (T2(1),YR),(T2(2),YI)
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      LA=L-1
      KP1=K+1
      IF (LA .LT. KP1) GO TO 45
      DO 40 M=KP1,LA
         I=M
         XR=ZERO
         XI=ZERO
         DO 5 J=M,L
            IF (ABS(AR(J,M-1))+ABS(AI(J,M-1)) .LE. ABS(XR)+ABS(XI))
     1      GO TO 5
            XR=AR(J,M-1)
            XI=AI(J,M-1)
            I=J
    5    CONTINUE
         ID(M)=I
         IF (I .EQ. M) GO TO 20
C                                  INTERCHANGE ROWS AND COLUMNS OF
C                                    ARRAYS AR AND AI
         MM1=M-1
         DO 10 J=MM1,N
            YR=AR(I,J)
            AR(I,J)=AR(M,J)
            AR(M,J)=YR
            YI=AI(I,J)
            AI(I,J)=AI(M,J)
            AI(M,J)=YI
   10    CONTINUE
         DO 15 J=1,L
            YR=AR(J,I)
            AR(J,I)=AR(J,M)
            AR(J,M)=YR
            YI=AI(J,I)
            AI(J,I)=AI(J,M)
            AI(J,M)=YI
   15    CONTINUE
C                                  END INTERCHANGE
   20    IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 40
         MP1=M+1
         DO 35 I=MP1,L
            YR=AR(I,M-1)
            YI=AI(I,M-1)
            IF (YR .EQ. ZERO .AND. YI .EQ. ZERO) GO TO 35
            Y=Y/X
            AR(I,M-1)=YR
            AI(I,M-1)=YI
            DO 25 J=M,N
               AR(I,J)=AR(I,J)-YR*AR(M,J)+YI*AI(M,J)
               AI(I,J)=AI(I,J)-YR*AI(M,J)-YI*AR(M,J)
   25       CONTINUE
            DO 30 J=1,L
               AR(J,M)=AR(J,M)+YR*AR(J,I)-YI*AI(J,I)
               AI(J,M)=AI(J,M)+YR*AI(J,I)+YI*AR(J,I)
   30       CONTINUE
   35    CONTINUE
   40 CONTINUE
   45 RETURN
      END

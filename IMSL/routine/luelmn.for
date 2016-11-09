C   IMSL ROUTINE NAME   - LUELMN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE LEQT2F
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUELMN (A,IA,N,B,APVT,X)
C
      DIMENSION          A(IA,1),B(1),APVT(1),X(1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  SOLVE LY = B FOR Y
      DO 5 I=1,N
    5 X(I) = B(I)
      IW = 0
      DO 20 I=1,N
         IP = APVT(I)
         SUM = X(IP)
         X(IP) = X(I)
         IF (IW .EQ. 0) GO TO 15
         IM1 = I-1
         DO 10 J=IW,IM1
            SUM = SUM-A(I,J)*X(J)
   10    CONTINUE
         GO TO 20
   15    IF (SUM .NE. 0.) IW = I
   20 X(I) = SUM
C                                  SOLVE UX = Y FOR X
      DO 30 IB=1,N
         I = N+1-IB
         IP1 = I+1
         SUM = X(I)
         IF (IP1 .GT. N) GO TO 30
         DO 25 J=IP1,N
            SUM = SUM-A(I,J)*X(J)
   25   CONTINUE
   30 X(I) = SUM/A(I,I)
      RETURN
      END

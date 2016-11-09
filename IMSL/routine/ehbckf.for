C   IMSL ROUTINE NAME   - EHBCKF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRF
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHBCKF (Z,H,D,N,MM,IZH,K,L)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MM,IZH,K,L
      REAL               Z(IZH,1),H(IZH,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LM2,KI,LTEMP,M,MA,MP2,I,J
      REAL               G,ZERO
      DATA               ZERO/0.0/
C                                  ADAPTED FROM EISPACK ROUTINE ORTBAK
C                                  FIRST EXECUTABLE STATEMENT
      LM2=L-2
      IF(LM2.LT.K) GO TO 9005
      LTEMP=LM2+K
      DO 30 KI=K,LM2
         M=LTEMP-KI
         MA=M+1
         IF(H(MA,M).EQ.ZERO) GO TO 30
         MP2=M+2
         IF(MP2.GT.L) GO TO 10
         DO 5 I=MP2,L
            D(I)=H(I,M)
    5    CONTINUE
   10    IF(MA.GT.L) GO TO 30
         DO 25 J=1,MM
            G=ZERO
            DO 15 I=MA,L
               G=G+D(I)*Z(I,J)
   15       CONTINUE
C                                  DOUBLE DIVISION AVOIDS POSSIBLE
C                                  UNDERFLOW
            G = (G/D(MA))/H(MA,M)
            DO 20 I=MA,L
               Z(I,J)=Z(I,J)+G*D(I)
   20       CONTINUE
   25    CONTINUE
   30 CONTINUE
 9005 RETURN
      END

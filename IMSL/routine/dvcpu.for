C   IMSL ROUTINE NAME   - DVCPU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DVCPR
C
C   REQD. IMSL ROUTINES - NONE
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPU (I0,N,NP,C,BB,X,XBAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            I0,N,NP
      REAL               C(1),BB(1),X(1),XBAR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JM1,J,KM1,K,LL,N1,NN
      REAL               ALF(50),HI
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I=1,N
    5 C(I) = BB(I)
      DO 10 I=1,N
         HI = X(I0+1)-X(I0)
         ALF(I) = (X(I0-NP+I)-XBAR)/HI
   10 CONTINUE
      NN = N-1
      N1 = N+1
      DO 20 I=1,NN
         LL = N-I
         DO 15 J=1,LL
            K = N1-J
            C(K) = C(K)-ALF(I)*C(K-1)
   15    CONTINUE
   20 CONTINUE
      DO 30 I=1,NN
         K = N-I
         KM1 = K+1
         DO 25 J=KM1,N
            C(J) = C(J)/(ALF(J)-ALF(J-K))
            JM1 = J-1
            C(JM1) = C(JM1)-C(J)
   25    CONTINUE
   30 CONTINUE
C
      RETURN
      END

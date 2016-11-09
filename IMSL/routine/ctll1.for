C   IMSL ROUTINE NAME   - CTLL1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE CTLLF
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CTLL1  (NVAR,X,Y,LOCX,LOCY,NVAL,IWK,IWK1)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NVAR,LOCX,LOCY,NVAL(1),IWK(1),IWK1(1)
      REAL               X(1),Y(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IND,I,J,K,LOCU,L,NVARP1,N
C                                  FIRST EXECUTABLE STATEMENT
      IWK1(1) = 1
      NVARP1 = NVAR+1
      DO 5 K=1,NVAR
         L = IWK(K)
         IF (L.EQ.0) GO TO 10
         IWK1(K+1) = IWK1(K)*NVAL(L)
    5 CONTINUE
C                                  FIND NUMBER OF VARIABLES IN
C                                    CONFIGURATION
      K = NVAR+1
   10 N = K-1
      LOCU = LOCY+IWK1(K)-1
      DO 15 J=LOCY,LOCU
         Y(J) = 0.0
   15 CONTINUE
C                                  INITIALIZE COORDINATES
      DO 20 K=1,NVAR
         IND = NVARP1+K
         IWK1(IND) = 0
   20 CONTINUE
C                                  FIND LOCATIONS IN TABLES
      I = LOCX
   25 J = LOCY
      DO 30 K=1,N
         L = IWK(K)
         J = J+IWK1(NVARP1+L)*IWK1(K)
   30 CONTINUE
      Y(J) = Y(J)+X(I)
C                                  UPDATE COORDINATES
      I = I+1
      DO 35 K=1,NVAR
         IND = NVARP1+K
         IWK1(IND) = IWK1(IND)+1
         IF (IWK1(IND).LT.NVAL(K)) GO TO 25
         IWK1(IND) = 0.0
   35 CONTINUE
      RETURN
      END

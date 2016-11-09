C   IMSL ROUTINE NAME   - CTLL2
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
      SUBROUTINE CTLL2  (NVAR,X,Y,Z,LOCX,LOCY,LOCZ,NVAL,IWK,D,IWK1)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NVAR,LOCX,LOCY,LOCZ,NVAL(1),IWK(1),IWK1(1)
      REAL               X(1),Y(1),Z(1),D
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IND,I,J,K,L,NVARP1,N
      REAL               E
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
C                                  TEST SIZE OF DEVIATION
      L = IWK1(K)
      J = LOCY
      K = LOCZ
      DO 15 I=1,L
         E = ABS(Z(K)-Y(J))
         IF (E.GT.D) D = E
         J = J+1
         K = K+1
   15 CONTINUE
C                                  INITIALIZE COORDINATES
      DO 20 K=1,NVAR
         IND = NVARP1+K
         IWK1(IND) = 0
   20 CONTINUE
      I = LOCX
C                                  PERFORM ADJUSTMENT
   25 J = 0
      DO 30 K=1,N
         L = IWK(K)
         J = J+IWK1(NVARP1+L)*IWK1(K)
   30 CONTINUE
      K = J+LOCZ
      J = J+LOCY
C                                  NOTE THAT Y(J) SHOULD BE
C                                    NON-NEGATIVE
      IF (Y(J).LE.0.0) X(I) = 0.0
      IF (Y(J).GT.0.0) X(I) = X(I)*Z(K)/Y(J)
C                                  UPDATE COORDINATES
      I = I+1
      DO 35 K=1,NVAR
         IND = NVARP1+K
         IWK1(IND) = IWK1(IND)+1
         IF (IWK1(IND).LT.NVAL(K)) GO TO 25
         IWK1(IND) = 0
   35 CONTINUE
      RETURN
      END

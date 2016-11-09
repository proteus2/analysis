C   IMSL ROUTINE NAME   - DTPTE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DTPTB
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DTPTE  (FCNI,FCNJ,STEMP,NN,X,Y,YPRIME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NN
      REAL               X,Y(NN),YPRIME(NN),STEMP(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N
      INTEGER            IJ,IK,I,J,KJ,K,N2
      REAL               ASGN
      COMMON /DTPTZ /    ASGN,N
C                                  FIRST EXECUTABLE STATEMENT
      IF (NN.GT.N) GO TO 5
      CALL FCNI(NN,X,Y,YPRIME)
      RETURN
    5 N2=N*N+1
      CALL FCNI(N,X,Y(N2),YPRIME(N2))
      CALL FCNJ(N,X,Y(N2),YPRIME)
      DO 20 I=1,N
         DO 10 J=1,N
            KJ=(J-1)*N
            IK=I-N
            STEMP(J)=0.
            DO 10 K=1,N
               KJ=KJ+1
               IK=IK+N
               STEMP(J)=STEMP(J)+YPRIME(IK)*Y(KJ)
   10    CONTINUE
         DO 15 J=1,N
            IJ=(J-1)*N+I
            YPRIME(IJ)=ASGN*STEMP(J)
   15    CONTINUE
   20 CONTINUE
      RETURN
      END

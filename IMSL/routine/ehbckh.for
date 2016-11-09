C   IMSL ROUTINE NAME   - EHBCKH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCH
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHBCKH (AR,AI,TAU,N,ZR,ZI,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IZ
      REAL               AR(1),AI(1),TAU(2,1),ZR(IZ,1),ZI(IZ,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,K,NR,L,NRM1,INX1,INX2,K1
      REAL               DELTA,ZERO,ALPHA1,ALPHA2
      DATA               ZERO/0.0/
C                                  TRANSFORM THE EIGENVECTORS OF THE
C                                    REAL SYMMETRIC TRIDIAGONAL MATRIX
C                                    TO THOSE OF THE HERMITIAN TRIDIA-
C                                    GONAL MATRIX
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 J=1,N
         DO 5 K=1,N
            ZI(J,K)=-ZR(J,K)*TAU(2,J)
            ZR(J,K)=ZR(J,K)*TAU(1,J)
    5 CONTINUE
      IF (N .LE. 2) GO TO 30
C                                  RECOVER THE HOUSEHOLDER MATRICES IN
C                                    REVERSE ORDER
      DO 25 L=3,N
         NR=N-L+2
         NRM1=NR-1
         INX1=(NR*(NRM1))/2+NR
         INX2=INX1-1
         IF (AI(INX1) .EQ. ZERO) GO TO 25
         DELTA=AI(INX1)* SQRT(AR(INX2)**2+AI(INX2)**2)
         DO 20 J=1,N
            ALPHA1=ZERO
            ALPHA2=ZERO
            DO 10 K=NR,N
               K1=(K*(K-1))/2+NRM1
               ALPHA1=ALPHA1+AR(K1)*ZR(K,J)+AI(K1)*ZI(K,J)
               ALPHA2=ALPHA2-AI(K1)*ZR(K,J)+AR(K1)*ZI(K,J)
   10       CONTINUE
            ALPHA1=ALPHA1/DELTA
            ALPHA2=ALPHA2/DELTA
            DO 15 K=NR,N
               K1=(K*(K-1))/2+NRM1
               ZR(K,J)=ZR(K,J)-AR(K1)*ALPHA1+AI(K1)*ALPHA2
               ZI(K,J)=ZI(K,J)-AR(K1)*ALPHA2-AI(K1)*ALPHA1
   15       CONTINUE
   20    CONTINUE
   25 CONTINUE
   30 RETURN
      END

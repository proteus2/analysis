C   IMSL ROUTINE NAME   - EBBCKC
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
      SUBROUTINE EBBCKC (ZR,ZI,N,IZ,K,L,M,D)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IZ,K,L,M
      REAL               ZR(IZ,1),ZI(IZ,1),D(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,II,KK
      REAL               S
C                                  FIRST EXECUTABLE STATEMENT
      IF (L .EQ. K) GO TO 15
      DO 10 I = K,L
         S = D(I)
C                                  LEFT HAND EIGENVECTORS ARE BACK
C                                    TRANSFORMED IF THE ABOVE
C                                    STATEMENT IS REPLACED BY S=1.0/D(I)
         DO 5 J = 1,M
            ZR(I,J) = ZR(I,J)*S
            ZI(I,J) = ZI(I,J)*S
    5    CONTINUE
   10 CONTINUE
C                                  DO 25 I=K-1,1,-1 AND
C                                    DO I=L+1,N,1
   15 DO 25 II = 1,N
         I = II
         IF (I .GE. K .AND. I .LE. L) GO TO 25
         IF (I .LT. K) I = K-II
         KK = D(I)
         IF (KK .EQ. I) GO TO 25
         DO 20 J = 1,M
            S = ZR(I,J)
            ZR(I,J) = ZR(KK,J)
            ZR(KK,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(KK,J)
            ZI(KK,J) = S
   20    CONTINUE
   25 CONTINUE
      RETURN
      END

C   IMSL ROUTINE NAME   - EHOBKS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHOBKS (A,N,M1,M2,Z,IZ)
C
      DIMENSION          A(1),Z(IZ,1)
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .EQ. 1) GO TO 30
      DO 25 I=2,N
         L = I-1
         IA = (I*L)/2
         H = A(IA+I)
         IF (H.EQ.0.) GO TO 25
C                                  DERIVES EIGENVECTORS M1 TO M2 OF
C                                  THE ORIGINAL MATRIX FROM EIGENVECTORS
C                                  M1 TO M2 OF THE SYMMETRIC
C                                  TRIDIAGONAL MATRIX
         DO 20 J = M1,M2
            S = 0.0
            DO 10 K = 1,L
               S = S+A(IA+K)*Z(K,J)
   10       CONTINUE
            S = S/H
            DO 15 K=1,L
               Z(K,J) = Z(K,J)-S*A(IA+K)
   15       CONTINUE
   20    CONTINUE
   25 CONTINUE
   30 RETURN
      END


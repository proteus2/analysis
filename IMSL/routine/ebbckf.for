C   IMSL ROUTINE NAME   - EBBCKF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
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
      SUBROUTINE EBBCKF (D,Z,K,L,MM,N,IZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,MM,N,IZ
      REAL               Z(IZ,1),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,KM1,II,JJ,LP1
      REAL               S
C                                  COLUMN SCALE Z BY APPROPRIATE D VALUE
C                                  FIRST EXECUTABLE STATEMENT
      IF (L.EQ.0) GO TO 15
      DO 10 I=K,L
         S = D(I)
         DO 5 J=1,MM
            Z(I,J) = Z(I,J)*S
    5    CONTINUE
   10 CONTINUE
C                                  INTERCHANGE ROWS IF PERMUTATIONS
C                                    OCCURRED IN EBALAF
   15 IF (K.EQ.1) GO TO 30
      KM1 = K-1
      DO 25 I=1,KM1
         II = K-I
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 25
         DO 20 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   20    CONTINUE
   25 CONTINUE
   30 IF (L.EQ.N) GO TO 45
      LP1 = L+1
      DO 40 II=LP1,N
         JJ = D(II)
         IF (II.EQ.JJ) GO TO 40
         DO 35 J=1,MM
            S = Z(II,J)
            Z(II,J) = Z(JJ,J)
            Z(JJ,J) = S
   35    CONTINUE
   40 CONTINUE
   45 RETURN
      END

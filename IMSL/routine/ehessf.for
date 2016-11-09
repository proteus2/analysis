C   IMSL ROUTINE NAME   - EHESSF
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
      SUBROUTINE EHESSF (A,K,L,N,IA,D)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,N,IA
      REAL               A(IA,N),D(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            LA,KP1,M,I,MP,II,J,JJ
      REAL               F,G,H,SCALE,ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      LA = L - 1
      KP1 = K + 1
      IF (LA .LT. KP1) GO TO 50
      DO 45 M = KP1, LA
         H = ZERO
         D(M) = ZERO
         SCALE = ZERO
C                                  SCALE COLUMN
         DO 5 I = M, L
            SCALE = SCALE + ABS(A(I,M-1))
    5    CONTINUE
         IF (SCALE .EQ. ZERO ) GO TO 45
         MP = M + L
C                                  DO 10 I=L,M,-1
         DO 10 II = M, L
            I = MP - II
            D(I) = A(I,M-1) / SCALE
            H = H + D(I) * D(I)
   10    CONTINUE
         G = -SIGN(SQRT(H),D(M))
         H = H - D(M) * G
         D(M) = D(M) - G
C                                  FORM (I-(U*UT)/H) * A
         DO 25 J = M,N
            F = ZERO
C                                  DO 15 I=L,M,-1
            DO 15 II = M, L
               I = MP - II
               F = F + D(I) * A(I,J)
   15       CONTINUE
            F = F / H
            DO 20 I = M, L
               A(I,J) = A(I,J) - F * D(I)
   20       CONTINUE
   25    CONTINUE
C                                  FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
         DO 40 I = 1,L
            F = ZERO
C                                  DO 30 J=L,M,-1
            DO 30 JJ = M, L
               J = MP - JJ
               F = F + D(J) * A(I,J)
   30       CONTINUE
            F = F / H
            DO 35 J = M, L
               A(I,J) = A(I,J) - F * D(J)
   35       CONTINUE
   40    CONTINUE
         D(M) = SCALE * D(M)
         A(M,M-1) = SCALE * G
   45 CONTINUE
   50 RETURN
      END

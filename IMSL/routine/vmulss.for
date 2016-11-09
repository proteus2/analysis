C   IMSL ROUTINE NAME   - VMULSS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (SYMMETRIC STORAGE
C                           MODE)
C
C   USAGE               - CALL VMULSS (A,B,N,C,IC)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING AN
C                           N BY N SYMMETRIC MATRIX STORED IN SYMMETRIC
C                           STORAGE MODE.
C                B      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING AN
C                           N BY N SYMMETRIC MATRIX STORED IN SYMMETRIC
C                           STORAGE MODE.
C                N      - ORDER OF MATRICES A, B, AND C. (INPUT)
C                C      - N BY N MATRIX IN FULL STORAGE MODE CONTAINING
C                           THE PRODUCT C = A*B. (OUTPUT)
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VMULSS (A,B,N,C,IC)
C
      REAL               A(1),B(1),C(IC,N)
      DOUBLE PRECISION   TEMP
C                                  FIRST EXECUTABLE STATEMENT
      LI=1
         DO 30 IR=1,N
         LJ=1
            DO 25 JC = 1,N
            I=LI
            J=LJ
            TEMP=0.0D0
               DO 20 K=1,N
               TEMP=DBLE(A(I))*DBLE(B(J))+TEMP
               IF (K .LT. IR) GO TO 5
               I=I+K
               GO TO 10
    5          I=I+1
   10          IF (K .LT. JC) GO TO 15
               J=J+K
               GO TO 20
   15          J=J+1
   20          CONTINUE
            C(IR,JC) = TEMP
            LJ = LJ+JC
   25       CONTINUE
         LI=LI+IR
   30    CONTINUE
      RETURN
      END

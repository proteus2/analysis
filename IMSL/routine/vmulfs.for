C   IMSL ROUTINE NAME   - VMULFS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (FULL BY SYMMETRIC
C                           MATRICES)
C
C   USAGE               - CALL VMULFS (A,B,L,M,IA,C,IC)
C
C   ARGUMENTS    A      - L BY M MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                B      - M BY M SYMMETRIC MATRIX STORED IN SYMMETRIC
C                           MODE. (INPUT)
C                L      - NUMBER OF ROWS IN A. (INPUT)
C                M      - NUMBER OF COLUMNS IN A (SAME AS NUMBER OF
C                           ROWS IN B). (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - L BY M MATRIX CONTAINING THE PRODUCT C = A*B.
C                           (OUTPUT)
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
      SUBROUTINE VMULFS (A,B,L,M,IA,C,IC)
C
      REAL               A(IA,M),B(1),C(IC,M)
      DOUBLE PRECISION   TEMP
C                                  FIRST EXECUTABLE STATEMENT
         DO 15 I=1,L
         LI=1
            DO 15 J=1,M
            LS=LI
            TEMP=0.0D0
               DO 10 K=1,M
               TEMP=TEMP+DBLE(A(I,K))*DBLE(B(LS))
               IF (K .GE. J) GO TO 5
               LS=LS+1
               GO TO 10
    5          LS=LS+K
   10          CONTINUE
            C(I,J)=TEMP
   15       LI=LI+J
      RETURN
      END

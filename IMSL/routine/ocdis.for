C   IMSL ROUTINE NAME   - OCDIS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - PAIRWISE EUCLIDEAN DISTANCES BETWEEN THE
C                           COLUMNS OF A MATRIX
C
C   USAGE               - CALL OCDIS (X,IX,N,M,WK,XDIS)
C
C   ARGUMENTS    X      - INPUT N BY M MATRIX WHOSE COLUMNS ARE USED
C                           TO COMPUTE THE PAIRWISE EUCLIDEAN
C                           DISTANCES
C                IX     - INPUT ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM
C                N      - INPUT NUMBER OF ROWS OF MATRIX X
C                M      - INPUT NUMBER OF COLUMNS OF MATRIX X
C                WK     - WORK VECTOR OF LENGTH N
C                XDIS   - OUTPUT MATRIX OF ORDER M CONTAINING THE
C                           EUCLIDEAN N-SPACE DISTANCES BETWEEN
C                           PAIRS OF COLUMNS OF MATRIX X. XDIS IS
C                           IN SYMMETRIC STORAGE MODE (A VECTOR OF
C                           LENGTH (M*(M+1))/2). THE DISTANCE BETWEEN
C                           COLUMN J AND COLUMN K IS XDIS(L) WHERE
C                           L = J*(J-1)/2 + K FOR K LESS THAN J.
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2
C                       - DOUBLE/VBLA=DNRM2
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OCDIS  (X,IX,N,M,WK,XDIS)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,N,M
      REAL               WK(N),X(IX,M),XDIS(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,INCX,J,K,L,IEND
      REAL               SNRM2
C                                  FIRST EXECUTABLE STATEMENT
      XDIS(1) = 0.0
      IF(M.EQ.1) GO TO 20
      INCX = 1
      L = 2
      IEND = 1
      DO 15 K=2,M
         DO 10 I=1,IEND
            DO 5 J=1,N
               WK(J) = X(J,K)-X(J,I)
    5       CONTINUE
            XDIS(L) = SNRM2(N,WK,INCX)
            L = L+1
   10    CONTINUE
         XDIS(L) = 0.0
         L = L+1
         IEND = K
   15 CONTINUE
   20 RETURN
      END

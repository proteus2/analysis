C   IMSL ROUTINE NAME   - LUELMF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ELIMINATION PART OF SOLUTION OF AX=B
C                           (FULL STORAGE MODE)
C
C   USAGE               - CALL LUELMF (A,B,IPVT,N,IA,X)
C
C   ARGUMENTS    A      - A = LU (THE RESULT COMPUTED IN THE IMSL
C                           ROUTINE LUDATF) WHERE L IS A LOWER
C                           TRIANGULAR MATRIX WITH ONES ON THE MAIN
C                           DIAGONAL. U IS UPPER TRIANGULAR. L AND U
C                           ARE STORED AS A SINGLE MATRIX A AND THE
C                           UNIT DIAGONAL OF L IS NOT STORED. (INPUT)
C                B      - B IS A VECTOR OF LENGTH N ON THE RIGHT HAND
C                           SIDE OF THE EQUATION AX=B. (INPUT)
C                IPVT   - THE PERMUTATION MATRIX RETURNED FROM THE
C                           IMSL ROUTINE LUDATF, STORED AS AN N LENGTH
C                           VECTOR. (INPUT)
C                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
C                IA     - ROW DIMENSION OF A EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                X      - THE RESULT X. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUELMF (A,B,IPVT,N,IA,X)
C
      DIMENSION          A(IA,1),B(1),IPVT(1),X(1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  SOLVE LY = B FOR Y
      DO 5 I=1,N
    5 X(I) = B(I)
      IW = 0
      DO 20 I=1,N
         IP = IPVT(I)
         SUM = X(IP)
         X(IP) = X(I)
         IF (IW .EQ. 0) GO TO 15
         IM1 = I-1
         DO 10 J=IW,IM1
            SUM = SUM-A(I,J)*X(J)
   10    CONTINUE
         GO TO 20
   15    IF (SUM .NE. 0.) IW = I
   20 X(I) = SUM
C                                  SOLVE UX = Y FOR X
      DO 30 IB=1,N
         I = N+1-IB
         IP1 = I+1
         SUM = X(I)
         IF (IP1 .GT. N) GO TO 30
         DO 25 J=IP1,N
            SUM = SUM-A(I,J)*X(J)
   25   CONTINUE
   30 X(I) = SUM/A(I,I)
      RETURN
      END

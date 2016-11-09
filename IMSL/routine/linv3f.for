C   IMSL ROUTINE NAME   - LINV3F
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - IN PLACE INVERSE, EQUATION SOLUTION, AND/OR
C                           DETERMINANT EVALUATION - FULL STORAGE MODE
C
C   USAGE               - CALL LINV3F (A,B,IJOB,N,IA,D1,D2,WKAREA,IER)
C
C   ARGUMENTS    A      - INPUT/OUTPUT MATRIX OF DIMENSION N BY N. SEE
C                           PARAMETER IJOB.
C                B      - INPUT/OUTPUT VECTOR OF LENGTH N WHEN IJOB =
C                           2 OR 3. OTHERWISE, B IS NOT USED.
C                         ON INPUT, B CONTAINS THE RIGHT HAND SIDE OF
C                           OF THE EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION X REPLACES B.
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES:
C                           I = 1, INVERT MATRIX A. A IS REPLACED BY
C                             ITS INVERSE.
C                           I = 2, SOLVE THE EQUATION AX = B. A IS
C                             REPLACED BY THE LU DECOMPOSITION OF A
C                             ROWWISE PERMUTATION OF A, WHERE U IS
C                             UPPER TRIANGULAR AND L IS LOWER
C                             TRIANGULAR WITH UNIT DIAGONAL.
C                             THE UNIT DIAGONAL OF L IS NOT STORED.
C                           I = 3, SOLVE AX = B AND INVERT MATRIX A.
C                             A IS REPLACED BY ITS INVERSE.
C                           I = 4, COMPUTE THE DETERMINANT OF A.
C                             A IS REPLACED BY THE LU DECOMPOSITION
C                             OF A ROWWISE PERMUTATION OF A.
C                N      - ORDER OF A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                D1     - INPUT/OUTPUT. IF THE D1 AND D2 COMPONENTS
C                D2         OF DETERMINANT(A) = D1*2**D2 ARE DESIRED,
C                           INPUT D1.GE.0. OTHERWISE, INPUT D1.LT.0.
C                           D2 IS NEVER INPUT.
C                WKAREA - WORK AREA OF LENGTH AT LEAST 2*N FOR IJOB=1
C                           OR IJOB=3.
C                         WORK AREA OF LENGTH AT LEAST N FOR IJOB=2
C                           OR IJOB=4.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING WITH FIX
C                           IER = 65 INDICATES THAT IJOB WAS LESS THAN
C                             1 OR GREATER THAN 4. IJOB IS ASSUMED TO
C                             BE 4.
C                         TERMINAL ERROR
C                           IER = 130 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - LUDATN,LUELMN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LINV3F (A,B,IJOB,N,IA,D1,D2,WKAREA,IER)
C
      REAL               A(IA,1),B(1),WKAREA(1),C1,C2,D1,D2,WA,ZERO,
     *                   ONE,SUM,C
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  LU DECOMPOSITION OF A
      CALL LUDATN (A,IA,N,A,IA,0,C1,C2,WKAREA,WKAREA,WA,IER)
      IF (D1 .LT. ZERO .AND. IJOB .GE. 1 .AND. IJOB .LT. 4) GO TO 5
      D1 = C1
      D2 = C2
    5 IF (IER .GE. 128) GO TO 60
      IF (IJOB .LE. 0 .OR. IJOB .GT. 4) GO TO 55
C                                  SOLVE AX = B
      IF (IJOB .EQ. 2 .OR . IJOB .EQ. 3) CALL LUELMN (A,IA,N,B,WKAREA,B)
      IF (IJOB .NE. 1 .AND. IJOB .NE. 3) GO TO 9005
C                                  MATRIX INVERSION
      A(N,N) = ONE/A(N,N)
      NM1 = N-1
      IF (NM1 .LT. 1) GO TO 9005
      DO 40 II=1,NM1
         L = N-II
         M = L+1
         DO 15 I=M,N
            SUM = ZERO
            DO 10 K=M,N
               SUM = SUM-A(I,K)*A(K,L)
   10       CONTINUE
            WKAREA(N+I) = SUM
   15    CONTINUE
         DO 20 I=M,N
            A(I,L) = WKAREA(N+I)
   20    CONTINUE
         DO 30 J=L,N
            SUM = ZERO
            IF (J .EQ. L) SUM = ONE
            DO 25 K=M,N
               SUM = SUM-A(L,K)*A(K,J)
   25       CONTINUE
            WKAREA(N+J) = SUM/A(L,L)
   30    CONTINUE
         DO 35 J=L,N
            A(L,J) = WKAREA(N+J)
   35    CONTINUE
   40 CONTINUE
C                                  PERMUTE COLUMNS OF A INVERSE
      DO 50 I=1,N
         J = N-I+1
         JP = WKAREA(J)
         IF (J .EQ. JP) GO TO 50
         DO 45 K=1,N
            C = A(K,JP)
            A(K,JP) = A(K,J)
            A(K,J) = C
   45    CONTINUE
   50 CONTINUE
      GO TO 9005
   55 CONTINUE
C                                  WARNING WITH FIX - IJOB WAS SET
C                                  INCORRECTLY
      IER = 65
      GO TO 9000
C                                  TERMINAL ERROR - MATRIX A IS
C                                  ALGORITHMICALLY SINGULAR
   60 IER = 130
 9000 CONTINUE
      CALL UERTST(IER,6HLINV3F)
 9005 RETURN
      END

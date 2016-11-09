C   IMSL ROUTINE NAME   - LEQT2P
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATIONS SOLUTION - POSITIVE DEFINITE
C                           MATRIX - SYMMETRIC STORAGE MODE - HIGH
C                           ACCURACY SOLUTION
C
C   USAGE               - CALL LEQT2P (A,M,N,IB,B,IDGT,D1,D2,WKAREA,
C                           IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE POSITIVE DEFINITE SYMMETRIC COEFFICIENT
C                           MATRIX OF THE EQUATION AX = B. A IS STORED
C                           IN SYMMETRIC STORAGE MODE.
C                M      - NUMBER OF RIGHT HAND SIDES. (INPUT)
C                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           THE M RIGHT HAND SIDES OF THE EQUATION
C                           AX = B.
C                         ON OUTPUT, THE N BY M MATRIX OF SOLUTIONS
C                           REPLACES B.
C                IDGT   - THE NUMBER OF DIGITS IN THE SOLUTION WHICH
C                           WERE UNCHANGED AFTER THE FIRST ITERATIVE
C                           IMPROVEMENT. (OUTPUT)
C                D1     - COMPONENTS OF THE DETERMINANT OF A.
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                WKAREA - WORK AREA OF DIMENSION GREATER THAN OR
C                           EQUAL TO N*(N+1)/2+2N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ORIGINAL
C                             MATRIX IS ALGORITHMICALLY NOT POSITIVE
C                             DEFINITE.
C                           IER = 130 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE. MATRIX
C                             A IS TOO ILL-CONDITIONED.
C
C   REQD. IMSL ROUTINES - SINGLE/LUDECP,LUELMP,LUREFP,UERTST,UGETIO
C                       - DOUBLE/LUDECP,LUELMP,LUREFP,UERTST,UGETIO,
C                           VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQT2P (A,M,N,IB,B,IDGT,D1,D2,WKAREA,IER)
C
      DIMENSION          A(1),B(IB,1),WKAREA(1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
      J = N+1
      K = N+J
C                                  DECOMPOSE A
      CALL LUDECP(A,WKAREA(K),N,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 10 I = 1,M
C                                  PERFORMS THE ELIMINATION PART OF
C                                  AX = B
         CALL LUELMP (WKAREA(K),B(1,I),N,WKAREA)
C                                  REFINEMENT OF SOLUTION TO AX = B
         CALL LUREFP (A,B(1,I),WKAREA(K),N,WKAREA,IDGT,WKAREA(J),IER)
         DO 5 II = 1,N
            B(II,I)=WKAREA(II)
    5    CONTINUE
         IF (IER .NE. 0) GO TO 15
   10 CONTINUE
      GO TO 9005
   15 IER = 130
 9000 CONTINUE
      CALL UERTST (IER,6HLEQT2P)
 9005 RETURN
      END

C   IMSL ROUTINE NAME   - LEQT1P
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE DEFINITE
C                           MATRIX - SYMMETRIC STORAGE MODE - SPACE
C                           ECONOMIZER SOLUTION
C
C   USAGE               - CALL LEQT1P (A,M,N,B,IB,IDGT,D1,D2,IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING THE
C                           N BY N COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B. A IS A POSITIVE DEFINITE SYMMETRIC
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
C                         ON OUTPUT, A IS REPLACED BY THE LOWER
C                           TRIANGULAR MATRIX L WHERE A = L*L-TRANSPOSE.
C                           L IS STORED IN SYMMETRIC STORAGE MODE WITH
C                           THE DIAGONAL ELEMENTS OF L IN RECIPROCAL
C                           FORM.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
C                B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           THE RIGHT-HAND SIDES OF THE EQUATION
C                           AX = B.
C                         ON OUTPUT, THE N BY M SOLUTION MATRIX X
C                           REPLACES B.
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT
C                           TO IDGT DECIMAL DIGITS. (CURRENTLY NOT USED)
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE INPUT MATRIX
C                             A IS ALGORITHMICALLY NOT POSITIVE
C                             DEFINITE. (SEE THE CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - LUDECP,LUELMP,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQT1P (A,M,N,B,IB,IDGT,D1,D2,IER)
C
      DIMENSION          A(1),B(IB,1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
C                                  DECOMPOSE A
      CALL LUDECP (A,A,N,D1,D2,IER)
      IF (IER.NE.0) GO TO 9000
C                                  PERFORM ELIMINATION
      DO 5 I = 1,M
         CALL LUELMP (A,B(1,I),N,B(1,I))
    5 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HLEQT1P)
 9005 RETURN
      END

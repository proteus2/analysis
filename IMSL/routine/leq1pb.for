C   IMSL ROUTINE NAME   - LEQ1PB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE
C                           DEFINITE SYMMETRIC BAND MATRIX - BAND
C                           SYMMETRIC STORAGE MODE - SPACE ECONOMIZER
C                           SOLUTION
C
C   USAGE               - CALL LEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
C
C   ARGUMENTS    A      - THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B, WHERE A IS ASSUMED TO BE AN N X N
C                           POSITIVE DEFINITE BAND SYMMETRIC MATRIX. A
C                           IS STORED IN BAND SYMMETRIC STORAGE MODE
C                           AND THEREFORE HAS DIMENSION N BY (NC+1).
C                           (INPUT)
C                         ON OUTPUT, A IS REPLACED BY L WHERE
C                           A = L*L-TRANSPOSE. L IS A LOWER BAND
C                           MATRIX STORED IN BAND FORM AND THEREFORE
C                           HAS DIMENSION N BY (NC+1). NOTE THAT THE
C                           DIAGONAL ELEMENTS OF L ARE STORED IN
C                           RECIPROCAL FORM. (OUTPUT)
C                N      - ORDER OF MATRIX A AND NUMBER OF ROWS IN B.
C                           (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS OF A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT MATRIX OF DIMENSION N X M CONTAINING
C                           THE M RIGHT-HAND SIDES OF THE EQUATION
C                           AX = B.
C                         ON OUTPUT, THE N X M SOLUTION MATRIX X
C                           REPLACES B.
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT
C                           TO IDGT DECIMAL DIGITS. (INPUT - CURRENTLY
C                           NOT USED)
C                D1     - COMPONENTS OF THE DETERMINANT OF A.
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C
C   REQD. IMSL ROUTINES - LUDAPB,LUELPB,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
C
      REAL               A(IA,1),B(IB,1),D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  DECOMPOSITION OF MATRIX A INTO
C                                  L*L-TRANSPOSE
      CALL LUDAPB(A,N,NC,IA,A,IA,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 5 I = 1,M
C                                  SOLUTION OF AX = B
         CALL LUELPB(A,B(1,I),N,NC,IA,B(1,I))
    5 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLEQ1PB)
 9005 RETURN
      END

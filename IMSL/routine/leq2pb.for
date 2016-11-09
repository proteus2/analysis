C   IMSL ROUTINE NAME   - LEQ2PB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE DEFINITE
C                           BAND SYMMETRIC MATRIX - BAND SYMMETRIC
C                           STORAGE MODE - HIGH ACCURACY SOLUTION
C
C   USAGE               - CALL LEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,
C                           WK,IER)
C
C   ARGUMENTS    A      - THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B, WHERE A IS ASSUMED TO BE AN N X N
C                           POSITIVE DEFINITE BAND SYMMETRIC MATRIX. A
C                           IS STORED IN BAND SYMMETRIC STORAGE MODE
C                           AND THEREFORE HAS DIMENSION N BY (NC+1).
C                           (INPUT)
C                N      - ORDER OF A AND THE NUMBER OF ROWS IN B.
C                           (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - INPUT MATRIX OF DIMENSION N X M CONTAINING
C                           THE M RIGHT HAND SIDES OF THE EQUATION
C                           AX = B.
C                         ON OUTPUT, THE N X M MATRIX OF SOLUTIONS
C                           REPLACES B.
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN
C                           B). (INPUT)
C                IDGT   - THE NUMBER OF DIGITS IN THE ANSWER WHICH
C                           WERE UNCHANGED AFTER THE FIRST ITERATIVE
C                           IMPROVEMENT. (OUTPUT)
C                D1     - COMPONENTS OF THE DETERMINANT OF A.
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                WK     - VECTOR OR MATRIX HAVING N*(NC+3) LOCATIONS
C                           USED AS WORK STORAGE.
C                         ON RETURN, THE FIRST N*(NC+1) LOCATIONS OF WK
C                           WILL CONTAIN THE MATRIX L STORED IN BAND
C                           STORAGE MODE WHERE A = L*L-TRANSPOSE.
C                           NOTE THAT THE DIAGONAL ELEMENTS OF L ARE
C                           STORED IN RECIPROCAL FORM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C                           IER = 130 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE. THE
C                             MATRIX A IS TOO ILL-CONDITIONED.
C
C   REQD. IMSL ROUTINES - SINGLE/LUDAPB,LUELPB,LUREPB,UERTST,UGETIO
C                       - DOUBLE/LUDAPB,LUELPB,LUREPB,UERTST,UGETIO,
C                           VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
C
      REAL               A(IA,1),B(IB,1),WK(N,1),D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      K = NC+2
      K1 = K+1
      CALL LUDAPB(A,N,NC,IA,WK,N,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
         DO 10 I = 1,M
         CALL LUELPB(WK,B(1,I),N,NC,N,WK(1,K))
         CALL LUREPB(A,N,NC,IA,WK,N,B(1,I),WK(1,K),IDGT,WK(1,K1),IER)
            DO 5 J = 1,N
            B(J,I) = WK(J,K)
    5       CONTINUE
         IF (IER .NE. 0) GO TO 15
   10    CONTINUE
      GO TO 9005
   15 IER = 130
 9000 CONTINUE
      CALL UERTST(IER,6HLEQ2PB)
 9005 RETURN
      END

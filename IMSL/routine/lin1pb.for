C   IMSL ROUTINE NAME   - LIN1PB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INVERSION OF A MATRIX - POSITIVE DEFINITE BAND
C                           SYMMETRIC MATRIX - BAND SYMMETRIC STORAGE
C                           MODE - SPACE ECONOMIZER SOLUTION
C
C   USAGE               - CALL LIN1PB (A,N,NC,IA,AINV,IDGT,D1,D2,IER)
C
C   ARGUMENTS    A      - N BY N POSITIVE DEFINITE BAND SYMMETRIC
C                           MATRIX TO BE INVERTED. A IS STORED IN
C                           BAND SYMMETRIC STORAGE MODE AND THEREFORE
C                           HAS DIMENSION N BY (NC+1). (INPUT)
C                         ON OUTPUT, A IS REPLACED BY THE LOWER BAND
C                           MATRIX L WHERE A = L*L-TRANSPOSE. L IS
C                           STORED IN BAND STORAGE MODE WITH THE
C                           DIAGONAL ELEMENTS OF L STORED IN RECIPROCAL
C                           FORM. (OUTPUT)
C                N      - ORDER OF A AND AINV. (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                AINV   - THE N BY N INVERSE OF A. AINV IS STORED IN
C                           SYMMETRIC STORAGE MODE. AINV IS A VECTOR OF
C                           LENGTH N*(N+1)/2. (OUTPUT)
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT
C                           TO IDGT DECIMAL DIGITS (CURRENTLY NOT USED).
C                D1     - COMPONENTS OF THE DETERMINANT OF A.
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT THE MATRIX A IS
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
      SUBROUTINE LIN1PB (A,N,NC,IA,AINV,IDGT,D1,D2,IER)
C
      REAL               A(IA,1),AINV(1),ZERO,ONE,D1,D2
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      L = 1
      N1 = N-1
      LN = N
      CALL LUDAPB(A,N,NC,IA,A,IA,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 10 I = 1,N
         DO 5 J = L,LN
            AINV(J) = ZERO
    5    CONTINUE
         AINV(L+I-1) = ONE
         CALL LUELPB(A,AINV(L),N,NC,IA,AINV(L))
         L = L+I
         LN = L+N1
   10 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HLIN1PB)
 9005 RETURN
      END

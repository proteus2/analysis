C   IMSL ROUTINE NAME   - LINV1P
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INVERSION OF MATRIX - POSITIVE DEFINITE -
C                           SYMMETRIC STORAGE MODE - SPACE ECONOMIZER
C                           SOLUTION
C
C   USAGE               - CALL LINV1P (A,N,AINV,IDGT,D1,D2,IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE N BY N POSITIVE DEFINITE SYMMETRIC
C                           MATRIX TO BE INVERTED. A IS STORED IN
C                           SYMMETRIC STORAGE MODE.
C                         ON OUTPUT, A IS REPLACED BY THE LOWER
C                           TRIANGULAR MATRIX L WHERE A = L*L-TRANSPOSE.
C                           L IS STORED IN SYMMETRIC STORAGE MODE WITH
C                           THE DIAGONAL ELEMENTS OF L STORED IN
C                           RECIPROCAL FORM.
C                N      - ORDER OF A. (INPUT)
C                AINV   - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE N BY N INVERSE OF A. AINV IS STORED IN
C                           SYMMETRIC STORAGE MODE.
C                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT
C                           TO IDGT DECIMAL DIGITS. (CURRENTLY NOT USED)
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ORIGINAL
C                             MATRIX A IS ALGORITHMICALLY NOT POSITIVE
C                             DEFINITE.  (SEE THE CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - LUDECP,LUELMP,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LINV1P (A,N,AINV,IDGT,D1,D2,IER)
C
      DIMENSION          A(1),AINV(1)
      DATA               ZERO,ONE/0.0,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      L = 1
      N1 = N-1
      LN = N
C                                  DECOMPOSE A
      CALL LUDECP(A,A,N,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 10 I = 1,N
         DO 5 J = L,LN
            AINV(J) = ZERO
    5    CONTINUE
         AINV(L+I-1) = ONE
C                                  OBTAIN THE SOLUTION
         CALL LUELMP(A,AINV(L),N,AINV(L))
         L = L+I
         LN = L+N1
   10 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HLINV1P)
 9005 RETURN
      END

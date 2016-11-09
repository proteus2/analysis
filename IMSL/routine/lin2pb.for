C   IMSL ROUTINE NAME   - LIN2PB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INVERSION OF MATRIX - POSITIVE DEFINITE BAND
C                           SYMMETRIC MATRIX - BAND SYMMETRIC STORAGE
C                           MODE - HIGH ACCURACY SOLUTION
C
C   USAGE               - CALL LIN2PB (A,N,NC,IA,AINV,IDGT,D1,D2,WK,
C                           IER)
C
C   ARGUMENTS    A      - N BY N POSITIVE DEFINITE BAND SYMMETRIC
C                           MATRIX TO BE INVERTED. A IS STORED IN BAND
C                           SYMMETRIC STORAGE MODE AND THEREFORE HAS
C                           DIMENSION N BY (NC+1). (INPUT)
C                N      - ORDER OF A AND AINV. (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                AINV   - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE N BY N INVERSE OF MATRIX A. AINV IS
C                           STORED IN SYMMETRIC STORAGE MODE.
C                IDGT   - THE APPROXIMATE NUMBER OF DIGITS IN THE
C                           ANSWER WHICH WERE UNCHANGED AFTER
C                           IMPROVEMENT. (OUTPUT)
C                D1     - COMPONENTS OF THE DETERMINANT OF A.
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                WK     - VECTOR OF LENGTH N*(NC+3) USED AS WORK
C                           STORAGE. ON OUTPUT THE FIRST N*(NC+1)
C                           LOCATIONS OF WK WILL CONTAIN THE
C                           MATRIX L IN BAND STORAGE MODE WHERE A =
C                           L*L-TRANSPOSE (THE DIAGONAL ELEMENTS OF L
C                           ARE STORED IN RECIPROCAL FORM).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           N = 129 INDICATES THAT THE MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C                             (SEE THE CHAPTER L PRELUDE).
C                           N = 130 INDICATES THAT ITERATIVE IMPROVE-
C                             MENT FAILED TO CONVERGE. THE MATRIX A
C                             IS TOO ILL-CONDITIONED.
C
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
      SUBROUTINE LIN2PB (A,N,NC,IA,AINV,IDGT,D1,D2,WK,IER)
C
      REAL               A(IA,1),AINV(1),WK(N,1),ONE,ZERO,D1,D2
      DATA               ZERO/0./,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      L = 1
      K1 = NC+2
      K2 = NC+3
C                                  DECOMPOSE A INTO L*L-TRANSPOSE
      CALL LUDAPB(A,N,NC,IA,WK,N,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 10 I = 1,N
         DO 5 J = 1,N
            WK(J,K1) = ZERO
    5    CONTINUE
         WK(I,K1) = ONE
C                                  OBTAIN THE INVERSE
         CALL LUELPB(WK,WK(1,K1),N,NC,N,AINV(L))
C                                  REFINE THE INVERSE
         CALL LUREPB(A,N,NC,IA,WK,N,WK(1,K1),AINV(L),IDGT,WK(1,K2),IER)
         IF (IER .NE. 0) GO TO 15
         L = L+I
   10 CONTINUE
      GO TO 9005
   15 IER = 130
 9000 CONTINUE
      CALL UERTST(IER,6HLIN2PB)
 9005 RETURN
      END

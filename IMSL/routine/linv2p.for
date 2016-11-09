C   IMSL ROUTINE NAME   - LINV2P
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INVERSION OF A MATRIX - POSITIVE DEFINITE -
C                           SYMMETRIC STORAGE MODE - HIGH ACCURACY
C                           SOLUTION
C
C   USAGE               - CALL LINV2P (A,N,AINV,IDGT,D1,D2,WKAREA,IER)
C
C   ARGUMENTS    A      - N BY N POSITIVE DEFINITE SYMMETRIC MATRIX TO
C                           BE INVERTED. A IS STORED IN SYMMETRIC
C                           STORAGE MODE. (INPUT)
C                N      - ORDER OF A. (INPUT)
C                AINV   - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE INVERSE OF MATRIX A. AINV IS STORED IN
C                           SYMMETRIC STORAGE MODE.
C                IDGT   - THE APPROXIMATE NUMBER OF DIGITS IN THE
C                           ANSWER WHICH WERE UNCHANGED AFTER
C                           IMPROVEMENT. (OUTPUT)
C                D1     - COMPONENTS OF THE DETERMINANT OF A.
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                WKAREA - WORK AREA OF DIMENSION GREATER THAN OR EQUAL
C                           TO (N*(N+1))/2+2*N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C                             (SEE THE CHAPTER L PRELUDE).
C                           IER = 130 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE.
C                             MATRIX A IS TOO ILL-CONDITIONED.
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
      SUBROUTINE LINV2P (A,N,AINV,IDGT,D1,D2,WKAREA,IER)
C
      REAL               A(1),AINV(1),WKAREA(1),D1,D2,ZERO,ONE
      DATA               ZERO,ONE/0.0,1.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER = 0
      L = 1
      K1 = N+1
      K2 = K1+N
C                                  DECOMPOSE A
      CALL LUDECP(A,WKAREA(K2),N,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 10 I = 1,N
         DO 5 J = 1,N
            WKAREA(J) = ZERO
    5    CONTINUE
         WKAREA(I) = ONE
         CALL LUELMP(WKAREA(K2),WKAREA,N,AINV(L))
         CALL LUREFP(A,WKAREA,WKAREA(K2),N,AINV(L),IDGT,WKAREA(K1),IER)
         IF (IER .NE. 0) GO TO 15
         L = L+I
   10 CONTINUE
      GO TO 9005
   15 IER = 130
 9000 CONTINUE
      CALL UERTST (IER,6HLINV2P)
 9005 RETURN
      END
 
R; T=0.02/0.20 22:20:08

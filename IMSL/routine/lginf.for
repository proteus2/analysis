C   IMSL ROUTINE NAME   - LGINF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERALIZED INVERSE OF A REAL MATRIX
C
C   USAGE               - CALL LGINF (A,IA,M,N,TOL,AINV,IAINV,S,WK,IER)
C
C   ARGUMENTS    A      - M BY N MATRIX. (INPUT)  A IS DESTROYED.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                           A IS USED BY LGINF AS WORK STORAGE FOR AN
C                           N BY N MATRIX. THEREFORE, IA MUST BE
C                           GREATER THAN OR EQUAL TO MAX(M,N).
C                M      - NUMBER OF ROWS IN A. (INPUT)
C                N      - NUMBER OF COLUMNS IN A. (INPUT)
C                TOL    - TOLERANCE PARAMETER. (INPUT)
C                           IF TOL IS LESS THAN OR EQUAL TO ZERO ON
C                           INPUT, LGINF COMPUTES THE GENERALIZED
C                           INVERSE OF A.
C                           IF TOL IS GREATER THAN ZERO ON INPUT, LGINF
C                           COMPUTES THE GENERALIZED INVERSE OF A MATRIX
C                           CLOSE TO A, BUT HAVING CONDITION NUMBER LESS
C                           THAN 1.0/TOL.
C                AINV   - N BY M MATRIX. (OUTPUT)
C                           AINV CONTAINS THE GENERALIZED INVERSE OF A
C                           OR THE GENERALIZED INVERSE OF A MATRIX
C                           CLOSE TO A. SEE TOL ARGUMENT DESCRIPTION.
C                IAINV  - ROW DIMENSION OF MATRIX AINV EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                S      - VECTOR OF LENGTH N. S CONTAINS THE ORDERED
C                           SINGULAR VALUES OF A.  S(1) .GE. S(2),...,
C                           .GE. S(N) .GE. 0. (OUTPUT)
C                WK     - WORK VECTOR OF LENGTH 2N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT MATRIX A IS NOT
C                             FULL RANK OR VERY ILL-CONDITIONED. SMALL
C                             SINGULAR VALUES MAY NOT BE VERY ACCURATE.
C                           IER=34 INDICATES THAT EITHER N.LE.0 OR
C                             M.LE.0.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT CONVERGENCE WAS
C                             NOT OBTAINED BY LSVDB AND COMPUTATION
C                             WAS DISCONTINUED.
C
C   REQD. IMSL ROUTINES - LGING,LSVDB,LSVG1,LSVG2,VHS12,UERSET,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LGINF  (A,IA,M,N,TOL,AINV,IAINV,S,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IAINV,IER,M,N
      REAL               A(IA,N),S(N),TOL,WK(N,2),AINV(IAINV,M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L
      REAL               T,TOLL
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      DO 10 I=1,N
         DO 5 J=1,M
            AINV(I,J)=0.0
            IF (I.EQ.J) AINV(I,J)=1.0
    5    CONTINUE
   10 CONTINUE
      CALL LGING(A,IA,M,N,AINV,IAINV,M,S,WK,IER)
      IF (IER.GT.128) GO TO 9000
      L = MIN0(M,N)
C                                  U**T HAS BEEN PLACED IN AINV
      TOLL = S(1)*AMAX1(TOL,0.0)
C                                  COMPUTE Q**(+) U**T
      DO 30 I=1,L
         IF (S(I).LE.TOLL) GO TO 20
         DO 15 J=1,M
   15    AINV(I,J)=AINV(I,J)/S(I)
         GO TO 30
   20    DO 25 J=1,M
   25    AINV(I,J)=0.0
   30 CONTINUE
C                                  COMPUTE V Q**(+) U**T
      DO 50 J=1,M
         DO 40 I=1,N
            T = 0.0
            DO 35 K=1,L
   35       T = T+A(I,K)*AINV(K,J)
            WK(I,1)=T
   40    CONTINUE
         DO 45 I=1,N
   45    AINV(I,J)=WK(I,1)
   50 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLGINF )
 9005 RETURN
      END

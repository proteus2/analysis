MODULE BVPShared
!
   USE Base, ONLY: i4, r8
!
   IMPLICIT NONE
!
   PRIVATE
!
!  COMMON /mchprs/
!
   REAL(r8),    PUBLIC :: flmin, flmax, epsmch
!
!  COMMON /algprs/
!
   LOGICAL,     PUBLIC :: pdebug, use_c, comp_c
   REAL(r8),    PUBLIC :: uval0
   INTEGER(i4), PUBLIC :: nminit, iprint, idum
!
!  COMMON /consts/
!   
   REAL(r8),    PUBLIC :: alp1, alp2, alp3, bet0, bet2, bet3, bet4
   REAL(r8),    PUBLIC :: a1, b1, c1, d1, e2, f2, c2, d2, b2, a2
   REAL(r8),    PUBLIC :: p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4
   REAL(r8),    PUBLIC :: a5, b5, c5, d5, e5, f5, a6, b6, c6
!
!  COMMON /cons1/
!
   REAL(r8),    PUBLIC :: a21a, a22a, a23a, a24a, a31a, a32a, a33a, a34a
   REAL(r8),    PUBLIC :: c1a, c2a, c16a, c26a, c123a, c223a, c14a, c24a
!
!  COMMON /cons2/
!
   REAL(r8),    PUBLIC :: a21b, a22b, a23b, a24b, a25b, a31b, a32b, a34b
   REAL(r8),    PUBLIC :: a35b, a41b, a42b, a43b, a44b, a45b
   REAL(r8),    PUBLIC :: b1b, b2b, b3b, c1b, c2b, c3b, c16b, c26b, c36b
   REAL(r8),    PUBLIC :: c123b, c223b, c323b, c14b, c24b, c34b
!
!  COMMON /flags/
!
   INTEGER(i4), PUBLIC :: ifinal, iback, iprec
!
   PUBLIC :: ABDNRM, DASUM, DONEST
   PUBLIC :: SPRT, MPRT, COLROW, COLROW1, INVERSE, CRDCMP, CRDCMP1, CRSLVE
   PUBLIC :: LUFAC, LUSOL, DCOPY, DAXPY, DDOT, DSCAL, DSWAP, IDAMAX, DLOAD
   PUBLIC :: MAXPY, MATCOP, MTLOAD, MTLOAD1, MSSQ, DSSQ
!
CONTAINS
!
!-------------------------------------------------------------------------------
!  FUNCTION ABDNRM
!-------------------------------------------------------------------------------
!
   FUNCTION ABDNRM(nbloks, ntop, nbot, novrlp, nrwblk, nclblk, top, a, bot) &
      RESULT(retval)
!
!     ABDNRM IS USED IN CONJUNCTION WITH DONEST TO COMPUTE THE
!     CONDITION NUMBER OF AN ALMOST BLOCK DIAGONAL MATRIX LIKE
!     THE ONES HANDLED BY COLROW. [SEE COMMENTS IN COLROW, CRDCMP]
!
!     THE BLAS ARE REQUIRED BY ABDNRM. SPECIFICALLY, DASUM IS USED.
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: nbloks
      INTEGER(i4), INTENT(IN)    :: ntop
      INTEGER(i4), INTENT(IN)    :: nbot
      INTEGER(i4), INTENT(IN)    :: novrlp
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      REAL(r8),    INTENT(IN)    :: top(ntop,*)
      REAL(r8),    INTENT(IN)    :: a(nrwblk,nclblk,*)
      REAL(r8),    INTENT(INOUT) :: bot(nbot,*)
!
      INTEGER(i4) :: j,k
      REAL(r8) :: temp, retval
!!!
!
      temp = 0.0D0
!
!     FIRST, GO OVER THE COLUMNS OF TOP AND THE FIRST BLOCK:
!
      DO j = 1, novrlp
         temp = MAX(temp, DASUM(ntop,top(1,j),1) + DASUM(nrwblk,a(1,j,1),1))
      END DO
      DO k = 1, nbloks-1
!
!        IN EACH BLOCK: FIRST, THE COLUMNS FROM THE KTH BLOCK ALONE:
!
         DO j = novrlp+1, nrwblk
            temp = MAX(temp, DASUM(nrwblk,a(1,j,k),1))
         END DO
!
!        NOW, TH COLUMNS WHICH INTERSECT BOTH THE KTH AND (K+1)ST BLOCKS.
!
         DO j = nrwblk+1, nclblk
            temp = MAX(temp, DASUM(nrwblk,a(1,j,k),1) +  &
                             DASUM(nrwblk,a(1,j-nrwblk,k+1),1))
         END DO
      END DO
!
!     FINALLY, THE COLUMNS OF THE LAST BLOCK WHICH DO NOT OVERLAP WITH ANYTHING.
!
      DO j = novrlp+1, nrwblk
         temp = MAX(temp, DASUM(nrwblk,a(1,j,nbloks),1))
      END DO
!
!     THOSE COLUMNS OVERLAPPING WITH BOTH THE LAST BLOCK AND THE BOTTOM BLOCK.
!
      DO j = nrwblk+1, nclblk
         temp = MAX(temp, DASUM(nrwblk,a(1,j,nbloks),1) +  &
                          DASUM(nbot,bot(1,j-nrwblk),1))
      END DO
!
      retval = temp   ! abdnrm = temp
!
      RETURN
   END FUNCTION ABDNRM
!
!-------------------------------------------------------------------------------
!  FUNCTION DASUM
!-------------------------------------------------------------------------------
!
   FUNCTION DASUM(n, dx, incx) RESULT(retval)
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN) :: n, incx
      REAL(r8),    INTENT(IN) :: dx(*)  ! dx(*)
!
      REAL(r8) :: dtemp, retval
      INTEGER(i4) :: i, m, mp1, nincx
!!!
!
      retval = 0.0D0
      dtemp = 0.0D0
      IF (n <= 0 .OR. incx <= 0) RETURN
      IF (incx == 1) GO TO 20
!
!     code for increment not equal to 1
!
      nincx = n*incx
      DO i = 1, nincx, incx
        dtemp = dtemp + DABS(dx(i))
      END DO
      retval = dtemp
      RETURN
!
!     code for increment equal to 1
!
!     clean-up loop
!
20    m = MOD(n,6)
      IF (m == 0) GO TO 40
      DO i = 1, m
         dtemp = dtemp + DABS(dx(i))
      END DO
      IF (n < 6) GO TO 60
40    mp1 = m + 1
      DO i = mp1, n, 6
         dtemp = dtemp + DABS(dx(i)) + DABS(dx(i + 1)) + DABS(dx(i + 2))  &
                       + DABS(dx(i + 3)) + DABS(dx(i + 4)) + DABS(dx(i + 5))
      END DO
60    retval = dtemp
!
      RETURN
   END FUNCTION DASUM
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DONEST
!-------------------------------------------------------------------------------
!
   SUBROUTINE DONEST(n, v, x, isgn, est, kase)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: v(n), x(n)
      INTEGER(i4), INTENT(INOUT) :: isgn(n)
      REAL(r8),    INTENT(INOUT) :: est
      INTEGER(i4), INTENT(INOUT) :: kase
!
!     DONEST ESTIMATES THE 1-NORM OF A SQUARE, REAL(r8) MATRIX  A.
!     REVERSE COMMUNICATION IS USED FOR EVALUATING MATRIX-VECTOR PRODUCTS.
!
!     ON ENTRY
!
!        N       INTEGER(i4)
!                THE ORDER OF THE MATRIX.  N .GE. 1.
!
!        ISGN    INTEGER(i4)(N)
!                USED AS WORKSPACE.
!
!        KASE    INTEGER(i4)
!                = 0.
!
!     ON INTERMEDIATE RETURNS
!
!        KASE    = 1 OR 2.
!
!        X       REAL(r8)(N)
!                MUST BE OVERWRITTEN BY
!
!                     A*X,             IF KASE=1,
!                     TRANSPOSE(A)*X,  IF KASE=2,
!
!                AND DONEST MUST BE RE-CALLED, WITH ALL THE OTHER
!                PARAMETERS UNCHANGED.
!
!     ON FINAL RETURN
!
!        KASE    = 0.
!
!        EST     REAL(r8)
!                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
!
!        V       REAL(r8)(N)
!                = A*W,   WHERE  EST = NORM(V)/NORM(W)
!                         (W  IS NOT RETURNED).
!
!     THIS VERSION DATED MARCH 16, 1988.
!     NICK HIGHAM, UNIVERSITY OF MANCHESTER.
!
!     MODIFIED FOR REAL(r8) ON JUNE 11, 1996.
!
!     REFERENCE
!     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
!     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS
!     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
!     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
!
!     SUBROUTINES AND FUNCTIONS
!     BLAS     IDAMAX, DASUM, DCOPY
!     GENERIC  ABS, NINT, FLOAT, SIGN
!
      INTEGER(i4), PARAMETER :: itmax = 5
!
      REAL(r8), PARAMETER :: zero = 0.0D0
      REAL(r8), PARAMETER :: one = 1.0D0
      REAL(r8), PARAMETER :: two = 2.0D0
!
!     INTERNAL VARIABLES
!
      INTEGER(i4) :: i, iter, j, jlast, jump
      REAL(r8) :: altsgn, estold, temp
!
      SAVE
!!!
!
      IF (kase == 0) THEN
         DO i = 1, n
            x(i) = one / FLOAT(n)
         END DO
         kase = 1
         jump = 1
         RETURN
      END IF
!
      SELECT CASE ( jump )
         CASE (    1)
            GO TO 100
         CASE (    2)
            GO TO  200
         CASE (    3)
            GO TO  300
         CASE (    4)
            GO TO  400
         CASE (    5)
            GO TO  500
      END SELECT
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
100   CONTINUE
!
      IF (n == 1) THEN
         v(1) = x(1)
         est = ABS(v(1))
!        ... QUIT
         GO TO 510
      END IF
      est = DASUM(n,x,1)
!
      DO i = 1, n
         x(i) = SIGN(one,x(i))
         isgn(i) = NINT(x(i))
      END DO
      kase = 2
      jump = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
200   CONTINUE
!
      j = IDAMAX(n,x,1)
      iter = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
220   CONTINUE
!
      DO i = 1, n
         x(i) = zero
      END DO
      x(j) = one
      kase = 1
      jump = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
300   CONTINUE
!
      CALL DCOPY(n,x,1,v,1)
      estold = est
      est = DASUM(n,v,1)
      DO i = 1,n
         IF ( NINT( SIGN(one,x(i)) ) /= isgn(i) ) GO TO 320
      END DO
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 410
!
320   CONTINUE
!     TEST FOR CYCLING.
      IF (est <= estold) GO TO 410
!
      DO i = 1, n
         x(i) = SIGN(one,x(i))
         isgn(i) = NINT(x(i))
      END DO
      kase = 2
      jump = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
400   CONTINUE
!
      jlast = j
      j = IDAMAX(n,x,1)
      IF ((x(jlast) /= ABS(x(j))) .AND. (iter < itmax)) THEN
         iter = iter + 1
         GO TO 220
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
410   CONTINUE
!
      altsgn = one
      DO i = 1, n
         x(i) = altsgn * (one + FLOAT(i-1)/FLOAT(n-1))
         altsgn = -altsgn
      END DO
      kase = 1
      jump = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
500   CONTINUE
!
      temp = two*dasum(n,x,1)/FLOAT(3*n)
      IF (temp > est) THEN
         CALL DCOPY(n,x,1,v,1)
         est = temp
      END IF
!
510   kase = 0
!
      RETURN
   END SUBROUTINE DONEST
!

!-------------------------------------------------------------------------------
!  SUBROUTINE SPRT
!-------------------------------------------------------------------------------
!
   SUBROUTINE SPRT(n, array)
!
!  sprt prints an array.
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN) :: n
      REAL(r8),    INTENT(IN) :: array(n)
!
      INTEGER(i4) :: i
!!!
!
      WRITE(6,900) (array(i), i=1,n)
900   FORMAT(1H ,(7(1PE11.3)))
!
      RETURN
   END SUBROUTINE SPRT
!
!-------------------------------------------------------------------------------
!  SUBROUTINE MPRT
!-------------------------------------------------------------------------------
!
   SUBROUTINE MPRT(nrowd, nrow, ncol, array)
!
!  mprt prints a matrix.
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: nrowd, nrow, ncol
      REAL(r8),    INTENT(INOUT) :: array(nrowd, ncol)
!
      INTEGER(i4) :: i, j
!!!
      DO i = 1, nrow
         WRITE(6,900) i,(array(i,j), j=1,ncol)
      END DO
900   FORMAT(1H ,i5,(6(1PE11.3)))
!
      RETURN
   END SUBROUTINE MPRT
!
!-------------------------------------------------------------------------------
!
!  THE AUGUST 27 1992 VERSION OF COLROW IN WHICH X IS NO LONGER REQUIRED,
!  WITH THE SOLUTION BEING RETURNED IN B, THE RIGHT HAND SIDE. IN ADDITION,
!  ALL VARIABLES ARE EXPLICITLY DECLARED. A PARAMETER "JOB" IS INCLUDED,
!  TO SPECIFY WHICH OF A.X = B OR TRANSPOSE(A).X = B IS TO BE SOLVED.
!
!  THIS PROGRAM SOLVES ONE OF THE LINEAR SYSTEMS A*X = B OR TRANSPOSE(A)*X = B,
!  WHERE A IS AN ALMOST BLOCK DIAGONAL MATRIX OF THE FORM
!
!               TOPBLK
!               ARRAY(1)
!                     ARRAY(2)
!                          .
!                             .
!                                .
!                                   .
!                                    ARRAY(NBLOKS)
!                                           BOTBLK
!
!  WHERE
!           TOPBLK IS  NRWTOP  BY NOVRLP
!           ARRAY(K), K=1,NBLOKS, ARE NRWBLK BY NRWBLK+NOVRLP
!           BOTBLK IS NRWBOT BY NOVRLP,
!  AND
!           NOVRLP = NRWTOP + NRWBOT
!  WITH
!           NOVRLP.LE.NRWBLK .
!
!  THE LINEAR SYSTEM IS OF ORDER  N = NBLOKS*NRWBLK + NOVRLP.
!
!  THE METHOD IMPLEMENTED IS BASED ON GAUSS ELIMINATION WITH
!  ALTERNATE ROW AND COLUMN ELIMINATION WITH PARTIAL PIVOTING,
!  WHICH PRODUCES A STABLE DECOMPOSITION OF THE MATRIX  A
!  WITHOUT INTRODUCING FILL-IN.
!
!-------------------------------------------------------------------------------
!
!  TO OBTAIN A SINGLE PRECISION VERSION OF THIS PACKAGE, REMOVE
!  ALL REAL(r8) STATEMENTS.  THERE IS ONE SUCH STATEMENT
!  IN C O L R O W, THREE IN C R D C M P, AND TWO IN C R S O L V.
!  IN ADDITION, REFERENCES TO BUILT-IN FUNCTIONS DABS AND DMAX1
!  MUST BE REPLACED BY ABS AND AMAX1, RESPECTIVELY.  DABS OCCURS
!  NINE TIMES, IN C R D C M P.  DMAX1 OCCURS FOUR TIMES, IN
!  C R D C M P.  FINALLY, ZERO IS INITIALISED TO 0.D0 IN A
!  DATA STATEMENT IN C R D C M P.  THIS MUST BE REPLACED BY:
!               DATA ZERO/0.0/
!
!-------------------------------------------------------------------------------
!
!               *****  PARAMETERS  *****
!
!       *** ON ENTRY ...
!
!               N      - INTEGER(i4)
!                         THE ORDER OF THE LINEAR SYSTEM,
!                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
!
!               TOPBLK - REAL(r8)(NRWTOP,NOVRLP)
!                         THE FIRST BLOCK OF THE ALMOST BLOCK
!                         DIAGONAL MATRIX A
!
!               NRWTOP - INTEGER(i4)
!                         NUMBER OF ROWS IN THE BLOCK TOPBLK
!
!               NOVRLP - INTEGER(i4)
!                         THE NUMBER OF COLUMNS IN WHICH SUCC-
!                         ESSIVE BLOCKS OVERLAP, WHERE
!                                NOVRLP = NRWTOP + NRWBOT
!
!               ARRAY  - REAL(r8)(NRWBLK,NCLBLK,NBLOKS)
!                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
!                         BY NCLBLK BLOCK OF THE MATRIX A
!
!               NRWBLK - INTEGER(i4)
!                         NUMBER OF ROWS IN K-TH BLOCK
!
!               NCLBLK - INTEGER(i4)
!                         NUMBER OF COLUMNS IN K-TH BLOCK
!
!               NBLOKS - INTEGER(i4)
!                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
!                         THE MATRIX A
!
!               BOTBLK - REAL(r8)(NRWBOT,NOVRLP)
!                         THE LAST BLOCK OF THE MATRIX A
!
!               NRWBOT - INTEGER(i4)
!                         NUMBER OF ROWS IN THE BLOCK BOTBLK
!
!                PIVOT - INTEGER(i4)(N)
!                         WORK SPACE
!
!                    B - REAL(r8)(N)
!                         THE RIGHT HAND SIDE VECTOR
!
!               JOB    - INTEGER(i4), INDICATING:
!                      = 0: SOLVE A*X = B;
!                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
!
!       *** ON RETURN  ...
!
!               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
!                        DESIRED DECOMPOSITION OF THE MATRIX A
!                        (IF IFLAG = 0)
!
!                PIVOT - INTEGER(i4)(N)
!                         RECORDS THE PIVOTING INDICES DETER-
!                         MINED IN THE DECOMPOSITION
!
!                    B - REAL(r8)(N)
!                         THE SOLUTION VECTOR (IF IFLAG = 0)
!
!               IFLAG  - INTEGER(i4)
!                         =  1, IF INPUT PARAMETERS ARE INVALID
!                         = -1, IF MATRIX IS SINGULAR
!                         =  0, OTHERWISE
!
!-------------------------------------------------------------------------------
!
!               *****  AUXILIARY PROGRAMS  *****
!
!       CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
!    *     BOTBLK,NRWBOT,PIVOT,IFLAG)
!            - DECOMPOSES THE MATRIX  A  USING MODIFIED
!              ALTERNATE ROW AND COLUMN ELIMINATON WITH
!              PARTIAL PIVOTING, AND IS USED FOR THIS
!              PURPOSE IN C O L R O W.
!              THE ARGUMENTS ARE AS IN C O L R O W.
!
!       CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
!    *     BOTBLK,NRWBOT,PIVOT,B,JOB)
!            - SOLVES THE SYSTEM A*X = B ONCE A IS DECOMPOSED.
!              THE ARGUMENTS ARE ALL AS IN C O L R O W.
!
!-------------------------------------------------------------------------------
!
!       THE SUBROUTINE  C O L R O W  AUTOMATICALLY SOLVES THE
!  INPUT SYSTEM WHEN IFLAG=0.  C O L R O W  IS CALLED ONLY ONCE
!  FOR A GIVEN SYSTEM. THE SOLUTION FOR A SEQUENCE OF P RIGHT
!  HAND SIDES CAN BE OBTAINED BY ONE CALL TO  C O L R O W  AND
!  P-1 CALLS TO CRSLVE ONLY. SINCE THE ARRAYS TOPBLK,ARRAY,
!  BOTBLK AND PIVOT CONTAIN THE DECOMPOSITION OF THE GIVEN
!  COEFFICIENT MATRIX AND PIVOTING INFORMATION ON RETURN FROM
!  C O L R O W , THEY MUST NOT BE ALTERED BETWEEN SUCCESSIVE
!  CALLS TO CRSLVE WITH THE SAME LEFT HAND SIDES. FOR THE
!  SAME REASON, IF THE USER WISHES TO SAVE THE COEFFICIENT
!  MATRIX, THE ARRAYS TOPBLK,ARRAY,BOTBLK MUST BE COPIED
!  BEFORE A CALL TO  C O L R O W .
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE COLROW(n, topblk, nrwtop, novrlp, array, nrwblk,  &
                     nclblk, nbloks, botblk, nrwbot, pivot, b, iflag, job)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(INOUT) :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(INOUT) :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(INOUT) :: pivot(*)
      REAL(r8),    INTENT(INOUT) :: b(*)
      INTEGER(i4), INTENT(INOUT) :: iflag
      INTEGER(i4), INTENT(INOUT) :: job
!
!     INTEGER(i4) :: idamax, i, ifail
!!!
!
!     NOUT = 6, DO THE FACTORIZATION USING CRDCMP:
!
      CALL CRDCMP(n, topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                  botblk, nrwbot, pivot, iflag)
!
      IF (iflag /= 0) RETURN
!
!     *****************solving the linear system********************
!
      job = 0
!
      CALL CRSLVE(topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                  botblk, nrwbot, pivot, b, job)
!
      RETURN
   END SUBROUTINE COLROW
!!!
!!!
   SUBROUTINE COLROW1(n, topblk, nrwtop, novrlp, array, nrwblk,  &
                      nclblk, nbloks, botblk, nrwbot, pivot, b, iflag, job)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(INOUT) :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(INOUT) :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(INOUT) :: pivot(*)
      REAL(r8),    INTENT(INOUT) :: b(*)
      INTEGER(i4), INTENT(INOUT) :: iflag
      INTEGER(i4), INTENT(INOUT) :: job
!
!     INTEGER(i4) :: idamax, i, ifail
!!!
!
!     NOUT = 6, DO THE FACTORIZATION USING CRDCMP:
!
      CALL CRDCMP1(n, topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                  botblk, nrwbot, pivot, iflag)
!
      IF (iflag /= 0) RETURN
!
!     *****************solving the linear system********************
!
      job = 0
!
      CALL CRSLVE(topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                  botblk, nrwbot, pivot, b, job)
!
      RETURN
   END SUBROUTINE COLROW1

!-------------------------------------------------------------------------------
!  SUBROUTINE INVERSE
!-------------------------------------------------------------------------------
!
   SUBROUTINE INVERSE(n, topblk, nrwtop, novrlp, array, nrwblk, nclblk,  &
                      nbloks, botblk, nrwbot, pivot, inmat)
!
!     INVERSE COMPUTES THE INVERSE OF A MATRIX
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(INOUT) :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(INOUT) :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(INOUT) :: pivot(*)
      REAL(r8),    INTENT(INOUT) :: inmat(n,n)
!
      REAL(r8) :: work(n)
      INTEGER(i4) :: k, l, j
!!!
!
      DO k = 1, n
         DO l = 1, n
            work(l) = 0.0D0
            IF (l == k) work(l) = 1.0D0
         END DO
         CALL CRSLVE(topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                     botblk, nrwbot, pivot, work, 0)
         DO l = 1, n
            inmat(l,k) = work(l)
         END DO
      END DO
!
      RETURN
   END SUBROUTINE inverse
!
!-------------------------------------------------------------------------------
!
!  C R D C M P DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
!  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
!  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
!  TOPBLK, ARRAY, AND BOTBLK.
!
!-------------------------------------------------------------------------------
!
!               *****  PARAMETERS  *****
!
!       *** ON ENTRY ...
!
!               N      - INTEGER(i4)
!                         THE ORDER OF THE LINEAR SYSTEM,
!                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
!
!               TOPBLK - REAL(r8)(NRWTOP,NOVRLP)
!                         THE FIRST BLOCK OF THE ALMOST BLOCK
!                         DIAGONAL MATRIX A TO BE DECOMPOSED
!
!               NRWTOP - INTEGER(i4)
!                         NUMBER OF ROWS IN THE BLOCK TOPBLK
!
!               NOVRLP - INTEGER(i4)
!                         THE NUMBER OF COLUMNS IN WHICH SUCC-
!                         ESSIVE BLOCKS OVERLAP, WHERE
!                                NOVRLP = NRWTOP + NRWBOT
!
!               ARRAY  - REAL(r8)(NRWBLK,NCLBLK,NBLOKS)
!                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
!                         BY NCLBLK BLOCK OF THE MATRIX A
!
!               NRWBLK - INTEGER(i4)
!                         NUMBER OF ROWS IN K-TH BLOCK
!
!               NCLBLK - INTEGER(i4)
!                         NUMBER OF COLUMNS IN K-TH BLOCK
!
!               NBLOKS - INTEGER(i4)
!                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
!                         THE MATRIX A
!
!               BOTBLK - REAL(r8)(NRWBOT,NOVRLP)
!                         THE LAST BLOCK OF THE MATRIX A
!
!               NRWBOT - INTEGER(i4)
!                         NUMBER OF ROWS IN THE BLOCK BOTBLK
!
!                PIVOT - INTEGER(i4)(N)
!                         WORK SPACE
!       *** ON RETURN  ...
!
!               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
!                        DESIRED DECOMPOSITION OF THE MATRIX A
!                        (IF IFLAG = 0)
!
!                PIVOT - INTEGER(i4)(N)
!                         RECORDS THE PIVOTING INDICES DETER-
!                         MINED IN THE DECOMPOSITION
!
!               IFLAG  - INTEGER(i4)
!                         =  1, IF INPUT PARAMETERS ARE INVALID
!                         = -1, IF MATRIX IS SINGULAR
!                         =  0, OTHERWISE
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE CRDCMP(n, topblk, nrwtop, novrlp, array, nrwblk, nclblk,  &
                     nbloks, botblk, nrwbot, pivot, iflag)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(INOUT) :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(INOUT) :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(INOUT) :: pivot(*)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
      REAL(r8) :: zero = 0.0D+0
      REAL(r8) :: colmax, pivmax, tempiv, swap, colpiv, colmlt
      REAL(r8) :: rowmax, rowpiv, rowmlt
      INTEGER(i4) :: nrwtp1, nrowel, nrwel1, nvrlp0
      INTEGER(i4) :: iplus1, iplusn, ipvt, jplus1, jminn, loop, incrj
      INTEGER(i4) :: incr, kplus1, incrn, irwblk, ipvblk, jrwblk
      INTEGER(i4) :: i, j, k, l
!
!     DEFINE THE CONSDTANTS USED THROUGHOUT
!
      iflag = 0
      pivmax = zero
      nrwtp1 = nrwtop+1
      nrowel = nrwblk-nrwtop
      nrwel1 = nrowel+1
      nvrlp0 = novrlp-1
!
!     CHECK VALIDITY OF THE INPUT PARAMETERS....
!     IF PARAMETERS ARE INVALID THEN TERMINATE AT 10; ELSE CONTINUE AT 100.
!
      IF (n /= nbloks*nrwblk+novrlp) GO TO 10
      IF (novrlp /= nrwtop+nrwbot) GO TO 10
      IF (nclblk /= novrlp+nrwblk) GO TO 10
      IF (novrlp > nrwblk) GO TO 10
!
!     PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
!
      GO TO 20
10    CONTINUE
!
!     PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
!
      iflag = 1
      RETURN
20    CONTINUE
!
!     FIRST, IN TOPBLK....
!     APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN PIVOTING ....
!
      DO i = 1, nrwtop
!
         iplus1 = i+1
!
!        DETERMINE COLUMN PIVOT AND PIVOT INDEX
!
         ipvt = i
         colmax = DABS(topblk(i,i))
         DO j = iplus1, novrlp
            tempiv = DABS(topblk(i,j))
            IF (tempiv <= colmax) CYCLE
            ipvt = j
            colmax = tempiv
         END DO
!
!        TEST FOR SINGULARITY:
!        IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
         IF (pivmax+colmax == pivmax) GO TO 410
         pivmax = DMAX1(colmax,pivmax)
!
!        IF NECESSARY INTERCHANGE COLUMNS
!
         pivot(i) = ipvt
         IF (ipvt == i) GO TO 60
         DO l = i, nrwtop
            swap = topblk(l,ipvt)
            topblk(l,ipvt) = topblk(l,i)
            topblk(l,i) = swap
         END DO
         DO l = 1, nrwblk
            swap = array(l,ipvt,1)
            array(l,ipvt,1) = array(l,i,1)
            array(l,i,1) = swap
         END DO
60       CONTINUE
!
!        COMPUTE MULTIPLIERS AND PERFORM COLUMN ELIMINATION
!
         colpiv = topblk(i,i)
         DO j = iplus1, novrlp
            colmlt = topblk(i,j)/colpiv
            topblk(i,j) = colmlt
            IF (iplus1 > nrwtop) GO TO 80
            DO l = iplus1, nrwtop
               topblk(l,j) = topblk(l,j)-colmlt*topblk(l,i)
            END DO
80          CONTINUE
            DO l = 1, nrwblk
               array(l,j,1) = array(l,j,1)-colmlt*array(l,i,1)
            END DO
         END DO
!
      END DO
!
!     IN EACH BLOCK ARRAY(,,K)....
!
      incr = 0
!
      DO k = 1, nbloks
!
         kplus1 = k+1
!
!        FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH ROW PIVOTING....
!
         DO j = nrwtp1, nrwblk
!
            jplus1 = j+1
            jminn = j-nrwtop
!
!           DETERMINE ROW PIVOT AND PIVOT INDEX
!
            ipvt = jminn
            rowmax = DABS(array(jminn,j,k))
            loop = jminn+1
            DO i = loop, nrwblk
               tempiv = DABS(array(i,j,k))
               IF (tempiv <= rowmax) CYCLE
               ipvt = i
               rowmax = tempiv
            END DO
!
!           TEST FOR SINGULARITY:
!           IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
            IF (pivmax+rowmax == pivmax) GO TO 410
            pivmax = DMAX1(rowmax,pivmax)
!
!           IF NECESSARY INTERCHANGE ROWS
!
            incrj = incr+j
            pivot(incrj) = incr+ipvt+nrwtop
            IF (ipvt == jminn) GO TO 140
            DO l = j, nclblk
               swap = array(ipvt,l,k)
               array(ipvt,l,k) = array(jminn,l,k)
               array(jminn,l,k) = swap
            END DO
140         CONTINUE
!
!           COMPUTE MULTIPLERS
!
            rowpiv = array(jminn,j,k)
            DO i = loop, nrwblk
               array(i,j,k) = array(i,j,k)/rowpiv
            END DO
!
!           PERFORM ROW ELIMINATION WITH COLUMN INDEXING
!
            DO l = jplus1, nclblk
               rowmlt = array(jminn,l,k)
               DO i = loop, nrwblk
                  array(i,l,k) = array(i,l,k)-rowmlt*array(i,j,k)
               END DO
            END DO
!
         END DO
!
!        NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN PIVOTING....
!
         DO i = nrwel1, nrwblk
!
            iplusn = i+nrwtop
            iplus1 = i+1
!
!           DETERMINE COLUMN PIVOT AND PIVOT INDEX
!
            ipvt = iplusn
            colmax = DABS(array(i,ipvt,k))
            loop = iplusn+1
            DO j = loop, nclblk
               tempiv = DABS(array(i,j,k))
               IF (tempiv <= colmax) CYCLE
               ipvt = j
               colmax = tempiv
            END DO
!
!           TEST FOR SINGULARITY:
!           IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
            IF (pivmax+colmax == pivmax) GO TO 410
            pivmax = DMAX1(colmax,pivmax)
!
!           IF NECESSARY INTERCHANGE COLUMNS
!
            incrn = incr+iplusn
            pivot(incrn) = incr+ipvt
            irwblk = iplusn-nrwblk
            IF (ipvt == iplusn) GO TO 240
            DO l = i, nrwblk
               swap = array(l,ipvt,k)
               array(l,ipvt,k) = array(l,iplusn,k)
               array(l,iplusn,k) = swap
            END DO
            ipvblk = ipvt-nrwblk
            IF (k == nbloks) GO TO 220
            DO l = 1, nrwblk
               swap = array(l,ipvblk,kplus1)
               array(l,ipvblk,kplus1) = array(l,irwblk,kplus1)
               array(l,irwblk,kplus1) = swap
            END DO
            GO TO 240
220         CONTINUE
            DO l = 1, nrwbot
               swap = botblk(l,ipvblk)
               botblk(l,ipvblk) = botblk(l,irwblk)
               botblk(l,irwblk) = swap
            END DO
240         CONTINUE
!
!           COMPUTE MULTIPLIERS AND PERFORM COLUMN ELIMINATION
!
            colpiv = array(i,iplusn,k)
            DO j = loop, nclblk
               colmlt = array(i,j,k)/colpiv
               array(i,j,k) = colmlt
               IF (i == nrwblk) GO TO 260
               DO l = iplus1, nrwblk
                  array(l,j,k) = array(l,j,k)-colmlt*array(l,iplusn,k)
               END DO
260            CONTINUE
               jrwblk = j-nrwblk
               IF (k == nbloks) GO TO 280
               DO l = 1, nrwblk
                  array(l,jrwblk,kplus1) = array(l,jrwblk,kplus1)-colmlt  &
                                          *array(l,irwblk,kplus1)
               END DO
               CYCLE
280            CONTINUE
               DO l = 1, nrwbot
                  botblk(l,jrwblk) = botblk(l,jrwblk)-colmlt*botblk(l, irwblk)
               END DO
            END DO
!
         END DO
!
         incr = incr+nrwblk
!
      END DO
!
!     FINALLY, IN BOTBLK....
!     APPLY NRWBOT ROW ELIMINATIONS WITH ROW PIVOTING....
!     IF BOT HAS JUST ONE ROW GO TO 500
!
      IF (nrwbot == 1) GO TO 400
!
      DO j = nrwtp1, nvrlp0
!
         jplus1 = j+1
         jminn = j-nrwtop
!
!        DETERMINE ROW PIVOT AND PIVOT INDEX
!
         ipvt = jminn
         rowmax = DABS(botblk(jminn,j))
         loop = jminn+1
         DO i = loop, nrwbot
            tempiv = DABS(botblk(i,j))
            IF (tempiv <= rowmax) CYCLE
            ipvt = i
            rowmax = tempiv
         END DO
!
!        TEST FOR SINGULARITY:
!        IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
         IF (pivmax+rowmax == pivmax) GO TO 410
         pivmax = DMAX1(rowmax,pivmax)
!
!        IF NECESSARY INTERCHANGE ROWS
!
         incrj = incr+j
         pivot(incrj) = incr+ipvt+nrwtop
         IF (ipvt == jminn) GO TO 350
         DO l = j, novrlp
            swap = botblk(ipvt,l)
            botblk(ipvt,l) = botblk(jminn,l)
            botblk(jminn,l) = swap
         END DO
350      CONTINUE
!
!        COMPUTE MULTIPLIERS
!
         rowpiv = botblk(jminn,j)
         DO i = loop, nrwbot
            botblk(i,j) = botblk(i,j)/rowpiv
         END DO
!
!        PERFORM ROW ELIMINATION WITH COLUMN INDEXING
!
         DO l = jplus1, novrlp
            rowmlt = botblk(jminn,l)
            DO i = loop, nrwbot
               botblk(i,l) = botblk(i,l)-rowmlt*botblk(i,j)
            END DO
         END DO
!
      END DO
!
400   CONTINUE
!
!     DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
!
      IF (pivmax+DABS(botblk(nrwbot,novrlp)) /= pivmax) RETURN
!
!     MATRIX IS SINGULAR - SET IFLAG = - 1.
!     TERMINATE AT 1000.
!
410   CONTINUE
!
      iflag = -1
!
      RETURN
   END SUBROUTINE CRDCMP
!!!
!!!
   SUBROUTINE CRDCMP1(n, topblk, nrwtop, novrlp, array, nrwblk, nclblk,  &
                      nbloks, botblk, nrwbot, pivot, iflag)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(INOUT) :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(INOUT) :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(INOUT) :: pivot(*)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
      REAL(r8) :: zero = 0.0D+0
      REAL(r8) :: colmax, pivmax, tempiv, swap, colpiv, colmlt
      REAL(r8) :: rowmax, rowpiv, rowmlt, pivtol
      INTEGER(i4) :: nrwtp1, nrowel, nrwel1, nvrlp0
      INTEGER(i4) :: iplus1, iplusn, ipvt, jplus1, jminn, loop, incrj
      INTEGER(i4) :: incr, kplus1, incrn, irwblk, ipvblk, jrwblk
      INTEGER(i4) :: i, j, k, l
!
!     DEFINE THE CONSDTANTS USED THROUGHOUT
!
      iflag = 0
      pivmax = zero
      nrwtp1 = nrwtop+1
      nrowel = nrwblk-nrwtop
      nrwel1 = nrowel+1
      nvrlp0 = novrlp-1
      pivtol = 10.0d+0*epsmch
!
!     CHECK VALIDITY OF THE INPUT PARAMETERS....
!     IF PARAMETERS ARE INVALID THEN TERMINATE AT 10; ELSE CONTINUE AT 100.
!
      IF (n /= nbloks*nrwblk+novrlp) GO TO 10
      IF (novrlp /= nrwtop+nrwbot) GO TO 10
      IF (nclblk /= novrlp+nrwblk) GO TO 10
      IF (novrlp > nrwblk) GO TO 10
!
!     PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
!
      GO TO 20
10    CONTINUE
!
!     PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
!
      iflag = 1
      RETURN
20    CONTINUE
!
!     FIRST, IN TOPBLK....
!     APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN PIVOTING ....
!
      DO i = 1, nrwtop
!
         iplus1 = i+1
!
!        DETERMINE COLUMN PIVOT AND PIVOT INDEX
!
         ipvt = i
         colmax = DABS(topblk(i,i))
         DO j = iplus1, novrlp
            tempiv = DABS(topblk(i,j))
            IF (tempiv <= colmax) CYCLE
            ipvt = j
            colmax = tempiv
         END DO
!
!        TEST FOR SINGULARITY:
!        IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
         IF (colmax <= pivtol) THEN
            iflag = -1
            RETURN
         END IF
!
         IF (pivmax+colmax == pivmax) GO TO 410
         pivmax = DMAX1(colmax,pivmax)
!
!        IF NECESSARY INTERCHANGE COLUMNS
!
         pivot(i) = ipvt
         IF (ipvt == i) GO TO 60
         DO l = i, nrwtop
            swap = topblk(l,ipvt)
            topblk(l,ipvt) = topblk(l,i)
            topblk(l,i) = swap
         END DO
         DO l = 1, nrwblk
            swap = array(l,ipvt,1)
            array(l,ipvt,1) = array(l,i,1)
            array(l,i,1) = swap
         END DO
60       CONTINUE
!
!        COMPUTE MULTIPLIERS AND PERFORM COLUMN ELIMINATION
!
         colpiv = topblk(i,i)
         DO j = iplus1, novrlp
            colmlt = topblk(i,j)/colpiv
            topblk(i,j) = colmlt
            IF (iplus1 > nrwtop) GO TO 80
            DO l = iplus1, nrwtop
               topblk(l,j) = topblk(l,j)-colmlt*topblk(l,i)
            END DO
80          CONTINUE
            DO l = 1, nrwblk
               array(l,j,1) = array(l,j,1)-colmlt*array(l,i,1)
            END DO
         END DO
!
      END DO
!
!     IN EACH BLOCK ARRAY(,,K)....
!
      incr = 0
!
      DO k = 1, nbloks
!
         kplus1 = k+1
!
!        FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH ROW PIVOTING....
!
         DO j = nrwtp1, nrwblk
!
            jplus1 = j+1
            jminn = j-nrwtop
!
!           DETERMINE ROW PIVOT AND PIVOT INDEX
!
            ipvt = jminn
            rowmax = DABS(array(jminn,j,k))
            loop = jminn+1
            DO i = loop, nrwblk
               tempiv = DABS(array(i,j,k))
               IF (tempiv <= rowmax) CYCLE
               ipvt = i
               rowmax = tempiv
            END DO
!
!           TEST FOR SINGULARITY:
!           IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
            IF (rowmax <= pivtol) THEN
               iflag = -1
               RETURN
            END IF
!
            IF (pivmax+rowmax == pivmax) GO TO 410
            pivmax = DMAX1(rowmax,pivmax)
!
!           IF NECESSARY INTERCHANGE ROWS
!
            incrj = incr+j
            pivot(incrj) = incr+ipvt+nrwtop
            IF (ipvt == jminn) GO TO 140
            DO l = j, nclblk
               swap = array(ipvt,l,k)
               array(ipvt,l,k) = array(jminn,l,k)
               array(jminn,l,k) = swap
            END DO
140         CONTINUE
!
!           COMPUTE MULTIPLERS
!
            rowpiv = array(jminn,j,k)
            DO i = loop, nrwblk
               array(i,j,k) = array(i,j,k)/rowpiv
            END DO
!
!           PERFORM ROW ELIMINATION WITH COLUMN INDEXING
!
            DO l = jplus1, nclblk
               rowmlt = array(jminn,l,k)
               DO i = loop, nrwblk
                  array(i,l,k) = array(i,l,k)-rowmlt*array(i,j,k)
               END DO
            END DO
!
         END DO
!
!        NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN PIVOTING....
!
         DO i = nrwel1, nrwblk
!
            iplusn = i+nrwtop
            iplus1 = i+1
!
!           DETERMINE COLUMN PIVOT AND PIVOT INDEX
!
            ipvt = iplusn
            colmax = DABS(array(i,ipvt,k))
            loop = iplusn+1
            DO j = loop, nclblk
               tempiv = DABS(array(i,j,k))
               IF (tempiv <= colmax) CYCLE
               ipvt = j
               colmax = tempiv
            END DO
!
!           TEST FOR SINGULARITY:
!           IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
            IF (colmax <= pivtol) THEN
               iflag = -1
               RETURN
            END IF
!
            IF (pivmax+colmax == pivmax) GO TO 410
            pivmax = DMAX1(colmax,pivmax)
!
!           IF NECESSARY INTERCHANGE COLUMNS
!
            incrn = incr+iplusn
            pivot(incrn) = incr+ipvt
            irwblk = iplusn-nrwblk
            IF (ipvt == iplusn) GO TO 240
            DO l = i, nrwblk
               swap = array(l,ipvt,k)
               array(l,ipvt,k) = array(l,iplusn,k)
               array(l,iplusn,k) = swap
            END DO
            ipvblk = ipvt-nrwblk
            IF (k == nbloks) GO TO 220
            DO l = 1, nrwblk
               swap = array(l,ipvblk,kplus1)
               array(l,ipvblk,kplus1) = array(l,irwblk,kplus1)
               array(l,irwblk,kplus1) = swap
            END DO
            GO TO 240
220         CONTINUE
            DO l = 1, nrwbot
               swap = botblk(l,ipvblk)
               botblk(l,ipvblk) = botblk(l,irwblk)
               botblk(l,irwblk) = swap
            END DO
240         CONTINUE
!
!           COMPUTE MULTIPLIERS AND PERFORM COLUMN ELIMINATION
!
            colpiv = array(i,iplusn,k)
            DO j = loop, nclblk
               colmlt = array(i,j,k)/colpiv
               array(i,j,k) = colmlt
               IF (i == nrwblk) GO TO 260
               DO l = iplus1, nrwblk
                  array(l,j,k) = array(l,j,k)-colmlt*array(l,iplusn,k)
               END DO
260            CONTINUE
               jrwblk = j-nrwblk
               IF (k == nbloks) GO TO 280
               DO l = 1, nrwblk
                  array(l,jrwblk,kplus1) = array(l,jrwblk,kplus1)-colmlt  &
                                          *array(l,irwblk,kplus1)
               END DO
               CYCLE
280            CONTINUE
               DO l = 1, nrwbot
                  botblk(l,jrwblk) = botblk(l,jrwblk)-colmlt*botblk(l, irwblk)
               END DO
            END DO
!
         END DO
!
         incr = incr+nrwblk
!
      END DO
!
!     FINALLY, IN BOTBLK....
!     APPLY NRWBOT ROW ELIMINATIONS WITH ROW PIVOTING....
!     IF BOT HAS JUST ONE ROW GO TO 500
!
      IF (nrwbot == 1) GO TO 400
!
      DO j = nrwtp1, nvrlp0
!
         jplus1 = j+1
         jminn = j-nrwtop
!
!        DETERMINE ROW PIVOT AND PIVOT INDEX
!
         ipvt = jminn
         rowmax = DABS(botblk(jminn,j))
         loop = jminn+1
         DO i = loop, nrwbot
            tempiv = DABS(botblk(i,j))
            IF (tempiv <= rowmax) CYCLE
            ipvt = i
            rowmax = tempiv
         END DO
!
!        TEST FOR SINGULARITY:
!        IF SINGULAR THEN TERMINATE AT 1000; ELSE CONTINUE.
!
         IF (rowmax <= pivtol) THEN
            iflag = -1
            RETURN
         END IF
!
         IF (pivmax+rowmax == pivmax) GO TO 410
         pivmax = DMAX1(rowmax,pivmax)
!
!        IF NECESSARY INTERCHANGE ROWS
!
         incrj = incr+j
         pivot(incrj) = incr+ipvt+nrwtop
         IF (ipvt == jminn) GO TO 350
         DO l = j, novrlp
            swap = botblk(ipvt,l)
            botblk(ipvt,l) = botblk(jminn,l)
            botblk(jminn,l) = swap
         END DO
350      CONTINUE
!
!        COMPUTE MULTIPLIERS
!
         rowpiv = botblk(jminn,j)
         DO i = loop, nrwbot
            botblk(i,j) = botblk(i,j)/rowpiv
         END DO
!
!        PERFORM ROW ELIMINATION WITH COLUMN INDEXING
!
         DO l = jplus1, novrlp
            rowmlt = botblk(jminn,l)
            DO i = loop, nrwbot
               botblk(i,l) = botblk(i,l)-rowmlt*botblk(i,j)
            END DO
         END DO
!
      END DO
!
400   CONTINUE
!
!     DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
!
      IF (DABS(botblk(nrwbot,novrlp)) > pivtol) RETURN
      IF (pivmax+DABS(botblk(nrwbot,novrlp)) /= pivmax) RETURN
!
!     MATRIX IS SINGULAR - SET IFLAG = - 1.
!     TERMINATE AT 1000.
!
410   CONTINUE
!
      iflag = -1
!
      RETURN
   END SUBROUTINE CRDCMP1
!
!-------------------------------------------------------------------------------
!
!  C R S L V E  SOLVES THE LINEAR SYSTEM  A*X = B
!  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
!
!-------------------------------------------------------------------------------
!
!               *****  PARAMETERS  *****
!
!       *** ON ENTRY  ...
!
!               TOPBLK - REAL(r8)(NRWTOP,NOVRLP)
!                         OUTPUT FROM  C R D C M P
!
!               NOVRLP - INTEGER(i4)
!                         THE NUMBER OF COLUMNS IN WHICH SUCC-
!                         ESSIVE BLOCKS OVERLAP, WHERE
!                                NOVRLP = NRWTOP + NRWBOT
!
!               NRWTOP - INTEGER(i4)
!                         NUMBER OF ROWS IN THE BLOCK TOPBLK
!
!               ARRAY  - REAL(r8)(NRWBLK,NCLBLK,NBLOKS)
!                         OUTPUT FROM  C R D C M P
!
!               NRWBLK - INTEGER(i4)
!                         NUMBER OF ROWS IN K-TH BLOCK
!
!               NCLBLK - INTEGER(i4)
!                         NUMBER OF COLUMNS IN K-TH BLOCK
!
!               NBLOKS - INTEGER(i4)
!                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
!                         THE MATRIX A
!
!               BOTBLK - REAL(r8)(NRWBOT,NOVRLP)
!                         OUTPUT FROM  C R D C M P
!
!               NRWBOT - INTEGER(i4)
!                         NUMBER OF ROWS IN THE BLOCK BOTBLK
!
!                PIVOT - INTEGER(i4)(N)
!                         THE PIVOT VECTOR FROM  C R D C M P
!
!                    B - REAL(r8)(N)
!                         THE RIGHT HAND SIDE VECTOR
!
!               JOB    - INTEGER(i4), INDICATING:
!                      = 0: SOLVE A*X = B;
!                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
!
!       *** ON RETURN  ...
!
!                    B - REAL(r8)(N)
!                         THE SOLUTION VECTOR
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE CRSLVE(topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                     botblk, nrwbot, pivot, b, job)
!
      IMPLICIT NONE
!
      REAL(r8),    INTENT(IN)    :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(IN)    :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(IN)    :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(IN)    :: pivot(*)
      REAL(r8),    INTENT(INOUT) :: b(*)
      INTEGER(i4), INTENT(IN)    :: job
!
      REAL(r8) :: dotprd,bj,xincrj,bincrj,swap,bi
      INTEGER(i4) :: nrwtp1,nrwbk1,nvrlp1,nrwbt1,nrowel,nvrlp0,nblks1
      INTEGER(i4) :: nbktop,j,i,loop,incr,incrj,incri,jpivot,jrwtop
      INTEGER(i4) :: ll,l1,iplusn,incrn,nrwtp0,nrwel1,k,incrtp,nrwbtl
      INTEGER(i4) :: ipvtn,nrwell,ipvti,l
!!!
!
!     DEFINE THE CONSTANTS USED THROUGHOUT
!
      nrwtp1 = nrwtop+1
      nrwbk1 = nrwblk+1
      nvrlp1 = novrlp+1
      nrwtp0 = nrwtop-1
      nrwbt1 = nrwbot+1
      nrowel = nrwblk-nrwtop
      nrwel1 = nrowel+1
      nvrlp0 = novrlp-1
      nblks1 = nbloks+1
      nbktop = nrwblk+nrwtop
!
!     IF JOB IS NON-ZERO, TRANSFER TO THE SECTION DEALING WITH
!     TRANSPOSE(A)*X = B.
!
      IF ( job /= 0 ) GO TO 530
!
!     FORWARD RECURSION
!     FIRST, IN TOPBLK....  FORWARD SOLUTION
!
      DO j = 1,nrwtop
         b(j) = b(j) / topblk(j,j)
         IF (j == nrwtop) GO TO 120
         bj = -b(j)
         loop = j+1
         DO i = loop,nrwtop
            b(i) = b(i) + topblk(i,j) * bj
         END DO
120      CONTINUE
      END DO
!
!     IN EACH BLOCK ARRAY(,,K)....
!
      incr = 0
!
      DO k = 1,nbloks
!
         incrtp = incr+nrwtop
!
!        FORWARD MODIFICATION
!
         DO j = 1,nrwtop
            incrj = incr+j
            xincrj = -b(incrj)
            DO i = 1,nrwblk
               incri = incrtp+i
               b(incri) = b(incri) + array(i,j,k) * xincrj
            END DO
         END DO
!
!        FORWARD ELIMINATION
!
         DO j = nrwtp1, nrwblk
            incrj = incr+j
            jpivot = pivot(incrj)
            IF (jpivot == incrj) GO TO 225
            swap = b(incrj)
            b(incrj) = b(jpivot)
            b(jpivot) = swap
225         CONTINUE
            bincrj = -b(incrj)
            loop = j-nrwtp0
            DO i = loop,nrwblk
               incri = incrtp+i
               b(incri) = b(incri) + array(i,j,k) * bincrj
            END DO
         END DO
!
!        FORWARD SOLUTION
!
         DO j = nrwbk1, nbktop
            incrj = incr+j
            jrwtop = j -nrwtop
            b(incrj) = b(incrj)/array(jrwtop,j,k)
            IF (j == nbktop) GO TO 260
            xincrj = -b(incrj)
            loop = j-nrwtp0
            DO i = loop,nrwblk
               incri = incrtp+i
               b(incri) = b(incri) + array(i,j,k) * xincrj
            END DO
260         CONTINUE
         END DO
!
         incr = incr+nrwblk
!
      END DO
!
!     FINALLY, IN BOTBLK....     FORWARD MODIFICATION
!
      incrtp = incr+nrwtop
      DO j = 1, nrwtop
         incrj = incr+j
         xincrj = -b(incrj)
         DO i = 1, nrwbot
            incri = incrtp+i
            b(incri) = b(incri)+botblk(i,j)*xincrj
         END DO
      END DO
!
!     FORWARD ELIMINATION
!
      IF (nrwbot == 1) GO TO 350
!
      DO j = nrwtp1, nvrlp0
         incrj = incr+j
         jpivot = pivot(incrj)
         IF (jpivot == incrj) GO TO 325
         swap = b(incrj)
         b(incrj) = b(jpivot)
         b(jpivot) = swap
325      CONTINUE
         bincrj = -b(incrj)
         loop = j-nrwtp0
         DO i = loop, nrwbot
            incri = incrtp+i
            b(incri) = b(incri)+botblk(i,j)*bincrj
         END DO
      END DO
350   CONTINUE
!
!     BACKWARD RECURSION
!     FIRST IN BOTBLK....  BACKWARD SOLUTION
!
      DO ll = 1,nrwbot
         j = nvrlp1-ll
         incrj = incr+j
         nrwbtl = nrwbt1-ll
         b(incrj) = b(incrj) / botblk(nrwbtl,j)
         IF (ll == nrwbot) GO TO 420
         xincrj = -b(incrj)
         loop = nrwbot-ll
         DO i = 1,loop
            incri = incrtp+i
            b(incri) = b(incri) + botblk(i,j) * xincrj
         END DO
420      CONTINUE
      END DO
!
!     THEN IN EACH BLOCK ARRAY(,,K)....
!
      DO l = 1,nbloks
!
!        BACKWARD ELIMINATION
!
         k = nblks1-l
         incr = incr-nrwblk
         DO l1 = nrwel1, nrwblk
            i = nrwblk+nrwel1-l1
            iplusn = i+nrwtop
            loop = iplusn+1
            incrn = incr+iplusn
            dotprd = b(incrn)
            DO j = loop, nclblk
               incrj = incr+j
               dotprd = dotprd-array(i,j,k)*b(incrj)
            END DO
            b(incrn) = dotprd
            ipvtn = pivot(incrn)
            IF (incrn == ipvtn) GO TO 445
            swap = b(incrn)
            b(incrn) = b(ipvtn)
            b(ipvtn) = swap
445         CONTINUE
         END DO
!
!        BACKWARD MODIFICATION
!
         incrtp = incr+nrwtop
         DO j = nrwbk1,nclblk
            incrj = incr+j
            xincrj = -b(incrj)
            DO i = 1,nrowel
               incri = incrtp+i
               b(incri) = b(incri)+array(i,j,k)*xincrj
            END DO
         END DO
!
!        BACKWARD SOLUTION
!
         DO ll = 1,nrowel
            j = nrwbk1-ll
            incrj = incr+j
            nrwell = nrwel1-ll
            b(incrj) = b(incrj)/array(nrwell,j,k)
            IF(ll == nrowel)GO TO 470
            xincrj = -b(incrj)
            loop = nrowel-ll
            DO i = 1,loop
               incri = incrtp+i
               b(incri) = b(incri)+array(i,j,k)*xincrj
            END DO
470         CONTINUE
         END DO
!
      END DO
!
!     IN TOPBLK FINISH WITH.... BACKWARD ELIMINATION
!
      DO l = 1,nrwtop
         i = nrwtp1-l
         loop = i+1
         dotprd = b(i)
         DO j = loop,novrlp
            dotprd = dotprd-topblk(i,j)*b(j)
         END DO
         b(i) = dotprd
         ipvti = pivot(i)
         IF(i == ipvti)GO TO 515
         swap = b(i)
         b(i) = b(ipvti)
         b(ipvti) = swap
 515     CONTINUE
      END DO
!
!     RETURN FROM THE SOLUTION OF A.X = B.
!
      RETURN
!
!     IF JOB IS NON-ZERO, SOLVE TRANSPOSE(A)*X = B:
!
530   CONTINUE
!
!     FIRST, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U).
!
      DO i = 1,nrwtop
         ipvti = pivot(i)
         IF ( i /= ipvti ) THEN
            swap = b(i)
            b(i) = b(ipvti)
            b(ipvti) = swap
         END IF
         bi = -b(i)
         loop = i+1
         DO j = loop,novrlp
            b(j) = b(j) + bi*topblk(i,j)
         END DO
      END DO
!
!     IN EACH BLOCK, K = 1,..,NBLOKS:
!
      incr = nrwtop
!
      DO k = 1,nbloks
!
!        FIRST, THE FORWARD SOLUTION.
!
         DO j = 1,nrowel
            incrj = incr + j
            DO i = 1,j-1
               b(incrj) = b(incrj) - array(i,nrwtop+j,k)*b(incr+i)
            END DO
            b(incrj) = b(incrj)/array(j,nrwtop+j,k)
         END DO
!
!        FORWARD MODIFICATION.
!
         DO i = 1,novrlp
            incri = incr + nrowel + i
            loop = nrwblk + i
            DO j = 1,nrowel
               incrj = incr + j
               b(incri) = b(incri) - array(j,loop,k)*b(incrj)
            END DO
         END DO
!
!        NOW, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U). THIS
!        CORRESPONDS TO THE LOOP 540 ABOVE.
!
         incr = incr + nrowel
!
         DO i = 1,nrwtop
            incri = incr + i
            ipvti = pivot(incri)
            IF ( incri /= ipvti ) THEN
               swap = b(incri)
               b(incri) = b(ipvti)
               b(ipvti) = swap
            END IF
            loop = nrowel + i
            bi = -b(incri)
            DO j = i+1,novrlp
               incrj = incr+j
               l = nrwblk + j
               b(incrj) = b(incrj) + bi*array(loop,l,k)
            END DO
         END DO
!
         incr = incr + nrwtop
!
      END DO
!
!     FINALLY, FINISH WITH NRWBOT SOLUTIONS:
!
      DO j = 1,nrwbot
         incrj = incr + j
         DO i = 1,j-1
            b(incrj) = b(incrj) - botblk(i,j+nrwtop)*b(incr+i)
         END DO
         b(incrj) = b(incrj)/botblk(j,j+nrwtop)
      END DO
!
!     NOW, THE BACKWARD PASS:
!     FIRST, BACKWARD SOLUTION IN BOTBLK:
!
      incrj = incr + nrwbot
      DO j = 1,nrwbot-1
         incrj = incrj - 1
         DO i = nrwbot-j+1,nrwbot
            incri = incr + i
            b(incrj) = b(incrj) - botblk(i,novrlp-j)*b(incri)
         END DO
         IF ( incrj /= pivot(incrj) ) THEN
            swap = b(incrj)
            b(incrj) = b(pivot(incrj))
            b(pivot(incrj)) = swap
         END IF
      END DO
!
!     NOW DO THE DEFERRED OPERATIONS IN BOTBLOK:
!
      DO j = 1,nrwtop
         incrj = incr - j + 1
         DO i = 1,nrwbot
            incri = incr + i
            b(incrj) = b(incrj) - botblk(i,nrwtp1-j)*b(incri)
         END DO
      END DO
!
!     NOW, IN EACH BLOCK, K = NBLOKS,..,1:
!
      DO k = nbloks,1,-1
!
!        FIRST, THE BACKSUBSTITUIONS:
!
         DO j = 1,nrwtop
            incrj = incr - j + 1
            loop = nbktop - j + 1
            DO i = 1,j-1
               incri = incr - i + 1
               b(incrj) = b(incrj) - array(nrwblk-i+1,loop,k)*b(incri)
            END DO
            b(incrj) = b(incrj)/array(nrwblk-j+1,loop,k)
         END DO
!
!        THEN THE BACKWARD SOLUTION IN THE KTH BLOCK:
!
         DO j = 1,nrowel
            incrj = incr - nrwtop -j + 1
            DO i = 1,j+nrwtop-1
               incri = incrj + i
               b(incrj) = b(incrj) -  &
                          array(nrwblk-nrwtop-j+1+i,nrwblk-j+1,k)*b(incri)
            END DO
            IF ( incrj /= pivot(incrj) ) THEN
               swap = b(incrj)
               b(incrj) = b(pivot(incrj))
               b(pivot(incrj)) = swap
            END IF
         END DO
!
!        NOW, THE DEFERRED OPERATIONS ON B:
!
         incr = incr - nrwblk
         DO j = 1,nrwtop
            incrj = incr + j - nrwtop
            DO i = 1,nrwblk
               incri = incr + i
               b(incrj) = b(incrj) - array(i,j,k)*b(incri)
            END DO
         END DO
!
      END DO
!
!     FINALLY, THE LAST SET OF BACK-SUBSTITUTIONS IN TOPBLK:
!
      DO j = 1,nrwtop
         incrj = nrwtop -j + 1
         DO i = incrj+1,nrwtop
            b(incrj) = b(incrj) - topblk(i,incrj)*b(i)
         END DO
         b(incrj) = b(incrj)/topblk(incrj,incrj)
      END DO
!
!     RETURN FROM THE SOLUTION OF A-TRANSPOSE.X = B
!
      RETURN
   END SUBROUTINE CRSLVE
!
!-------------------------------------------------------------------------------
!  SUBROUTINE LUFAC
!-------------------------------------------------------------------------------
!
   SUBROUTINE LUFAC(n, ndim, a, ip, ier)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n, ndim
      REAL(r8),    INTENT(INOUT) :: a(ndim,n)
      INTEGER(i4), INTENT(INOUT) :: ip(n)
      INTEGER(i4), INTENT(INOUT) :: ier
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
!
      REAL(r8) :: piv
      INTEGER(i4) :: j, k, ipiv
!
!     The subroutine lufac is a very simple code to compute the
!     LU decomposition (with partial pivoting) of the n by n matrix a.
!
!     The LU factors are overwritten on a.  The integer array ip
!     reflects the pairwise interchanges performed.  note that ip(k)
!     therefore does not give the index in the original array of
!     the k-th pivot.
!
!     On exit, the error flag ier is zero when no zero pivots are
!     encountered.  Otherwise, ier is equal to the index of the
!     step at which a zero pivot occurred.
!!!
!
      ier = 0
      ip(n) = 0
!
!     Begin loop over columns 1 through n-1.  k is the current column index.
!
      DO k = 1, n-1
!
!        Find the row index ipiv of the element of largest magnitude in column k
!
         ipiv = k-1 + IDAMAX(n-k+1, a(k,k), 1)
         piv = a(ipiv,k)
         IF (piv == zero) THEN
            ier = k
            RETURN
         END IF
         ip(k) = ipiv
!
!        Perform interchanges if necessary.
!
         IF (ipiv /= k) THEN
            CALL DSWAP(n-k+1, a(ipiv,k), ndim, a(k,k), ndim)
         END IF
!
!        Save the (negative) multipliers in the subdiagonal elements of column k
!
         CALL DSCAL(n-k, (-one/piv), a(k+1,k), 1)
!
!        Update the remaining matrix.  Note that a(i,k) now contains
!        the negative multipliers.
!
         DO j = k+1, n
            CALL DAXPY(n-k, a(k,j), a(k+1,k), 1, a(k+1,j), 1)
         END DO
!
!        End of loop over columns.
!
      END DO
!
      IF (a(n,n) == zero) ier = n
!
      RETURN
   END SUBROUTINE LUFAC
!
!-------------------------------------------------------------------------------
!  SUBROUTINE LUSOL
!-------------------------------------------------------------------------------
!
   SUBROUTINE LUSOL(n, ndim, a, ip, b, x)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n, ndim
      REAL(r8),    INTENT(IN)    :: a(ndim,n)
      INTEGER(i4), INTENT(IN)    :: ip(n)
      REAL(r8),    INTENT(INOUT) :: b(n)
      REAL(r8),    INTENT(INOUT) :: x(n)
!
      REAL(r8) :: tem
      INTEGER(i4) :: k, kb, ipiv
!
!     blas:  daxpy, dcopy
!
!     The subroutine lusol is a simple-minded routine to solve a
!     linear system whose LU factors have been computed by lufac.
!     On entry, the matrix a should contain the LU factors, and
!     ip should contain the interchange array constructed by lufac.
!
!     Copy the right-hand side b into x.
!
!!!
!
      CALL DCOPY(n, b, 1, x, 1)
!
!     Forward solution with l (unit lower-triangular factor), which
!     is stored in the strict lower triangle of a.
!
      DO k = 1, n-1
         ipiv = ip(k)
         IF (ipiv /= k) THEN
            tem = x(ipiv)
            x(ipiv) = x(k)
            x(k) = tem
         END IF
         CALL DAXPY(n-k, x(k), a(k+1,k), 1, x(k+1), 1)
      END DO
!
!     Backward solution with u (upper-triangular factor), which is stored
!     in the upper triangle of a.
!
      DO kb = n, 1, -1
         x(kb) = x(kb)/a(kb,kb)
         CALL DAXPY(kb-1, (-x(kb)), a(1,kb), 1, x(1), 1)
      END DO
!
      RETURN
   END SUBROUTINE LUSOL
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DCOPY
!-------------------------------------------------------------------------------
!
   SUBROUTINE DCOPY(n, x, incx, y, incy)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n, incx, incy
      REAL(r8),    INTENT(IN)    :: x(*)      ! x( * )
      REAL(r8),    INTENT(INOUT) :: y(*)      ! y( * )
!
!     dcopy  performs the operation
!
!     y := x
!
!  nag fortran 77 version of the blas routine dcopy .
!  nag fortran 77 o( n ) basic linear algebra routine.
!
!  -- written on 26-november-1982.
!     sven hammarling, nag central office.
!
      INTEGER(i4) :: i, ix, iy
!!!
!
      IF( n < 1 ) RETURN
!
      IF ((incx == incy) .AND. (incy > 0)) THEN
         DO iy = 1, 1 + (n-1) * incy, incy
            y(iy) = x(iy)
         END DO
      ELSE
         IF (incx >= 0) THEN
            ix = 1
         ELSE
            ix = 1 - (n-1) * incx
         END IF
         IF (incy > 0) THEN
            DO iy = 1, 1 + (n-1) * incy, incy
               y(iy) = x(ix)
               ix = ix + incx
            END DO
         ELSE
            iy = 1 - (n-1) * incy
            DO i = 1, n
               y(iy) = x(ix)
               iy = iy + incy
               ix = ix + incx
            END DO
         END IF
      END IF
!
      RETURN
   END SUBROUTINE DCOPY
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DAXPY
!-------------------------------------------------------------------------------
!
   SUBROUTINE DAXPY( n, alpha, x, incx, y, incy )
!
      INTEGER(i4), INTENT(IN)    :: n, incx, incy
      REAL(r8),    INTENT(IN)    :: alpha
      REAL(r8),    INTENT(IN)    :: x(*)         ! x( * )
      REAL(r8),    INTENT(INOUT) :: y(*)         ! y( * )
!
!     daxpy  performs the operation
!
!        y := alpha*x + y
!
!     modified nag fortran 77 version of the blas routine daxpy .
!
!  -- written on 3-september-1982.
!     sven hammarling, nag central office.
!
      REAL(r8), PARAMETER :: zero  = 0.0+0
      INTEGER(i4) :: i, ix, iy
!!!
!
      IF (n < 1) RETURN
      IF (alpha == zero) RETURN
!
      IF ((incx == incy) .AND. (incx > 0)) THEN
         DO ix = 1, 1+(n-1)*incx, incx
            y(ix) = alpha*x(ix) + y(ix)
         END DO
      ELSE
         IF (incy >= 0) THEN
            iy = 1
         ELSE
            iy = 1-(n-1)*incy
         END IF
         IF (incx > 0) THEN
            DO ix = 1, 1+(n-1)*incx, incx
               y(iy) = alpha*x(ix) + y(iy)
               iy = iy + incy
            END DO
         ELSE
            ix = 1-(n-1)*incx
            DO i = 1, n
               y(iy) = alpha*x(ix) + y(iy)
               ix = ix + incx
               iy = iy + incy
            END DO
         END IF
      END IF
!
      RETURN
   END SUBROUTINE DAXPY
!
!------------------------------------------------------------------------------- 
!  FUNCTION DDOT
!------------------------------------------------------------------------------- 
!
   FUNCTION DDOT(n, x, incx, y, incy) RESULT(retval)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN) :: n, incx, incy
      REAL(r8),    INTENT(IN) :: x(*)      ! x( * )
      REAL(r8),    INTENT(IN) :: y(*)      ! y( * )
!
!     ddot   returns the value
!     ddot   = x'y
!     modified nag fortran 77 version of the blas routine ddot  .
!  -- written on 21-september-1982.
!     sven hammarling, nag central office.
!
      REAL(r8), PARAMETER :: zero  = 0.0D+0
      REAL(r8) :: retval, sum
      INTEGER(i4) :: i, ix, iy
!!!
!
      sum = zero
      IF (n >= 1) THEN
         IF ((incx == incy) .AND. (incx > 0)) THEN
            DO ix = 1, 1+(n-1)*incx, incx
               sum = sum + x(ix)*y(ix)
            END DO
         ELSE
            IF (incy >= 0) THEN
               iy = 1
            ELSE
               iy = 1-(n-1)*incy
            END IF
            IF (incx > 0) THEN
               DO ix = 1, 1+(n-1)*incx, incx
                  sum = sum + x(ix)*y(iy)
                  iy = iy + incy
               END DO
            ELSE
               ix = 1-(n-1)*incx
               DO i = 1, n
                  sum = sum + x(ix)*y(iy)
                  ix = ix + incx
                  iy = iy + incy
               END DO
            END IF
         END IF
      END IF
      retval = sum
!
      RETURN
   END FUNCTION DDOT
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DSCAL
!-------------------------------------------------------------------------------
!
   SUBROUTINE DSCAL(n, alpha, x, incx)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n, incx
      REAL(r8),    INTENT(IN)    :: alpha
      REAL(r8),    INTENT(INOUT) :: x(*)         ! x( * )
!
!     dscal performs the operation
!     x := alpha*x
!     modified nag fortran 77 version of the blas routine dscal .
!  -- written on 26-november-1982.
!     sven hammarling, nag central office.
!
      REAL(r8), PARAMETER :: one   = 1.0D+0
      REAL(r8), PARAMETER :: zero  = 0.0D+0
      INTEGER(i4) :: ix
!!!
!
      IF (n >= 1) THEN
         IF (alpha == zero) THEN
            DO ix = 1, 1+(n-1)*incx, incx
               x(ix) = zero
            END DO
         ELSE IF (alpha == (-one)) THEN
            DO ix = 1, 1+(n-1)*incx, incx
               x(ix) = -x(ix)
            END DO
         ELSE IF (alpha /= one) THEN
            DO ix = 1, 1+(n-1)*incx, incx
               x(ix) = alpha*x(ix)
            END DO
         END IF
      END IF
!
      RETURN
   END SUBROUTINE DSCAL
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DSWAP
!-------------------------------------------------------------------------------
!
   SUBROUTINE DSWAP(n, x, incx, y, incy)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n, incx, incy
      REAL(r8),    INTENT(INOUT) :: x(*)      ! x( * )
      REAL(r8),    INTENT(INOUT) :: y(*)      ! y( * )
!
!     dswap  performs the operations
!     temp := x,   x := y,   y := temp.
!     modified nag fortran 77 version of the blas routine dswap .
!  -- written on 26-november-1982.
!     sven hammarling, nag central office.
!
      REAL(r8) :: temp
      INTEGER(i4) :: i, ix, iy
!!!
!
      IF (n < 1) RETURN
      IF ((incx == incy) .AND. (incy > 0)) THEN
         DO iy = 1, 1+(n-1)*incy, incy
            temp = x(iy)
            x(iy) = y(iy)
            y(iy) = temp
         END DO
      ELSE
         IF (incx >= 0) THEN
            ix = 1
         ELSE
            ix = 1-(n-1)*incx
         END IF
         IF (incy > 0) THEN
            DO iy = 1, 1+(n-1)*incy, incy
               temp  = x(ix)
               x(ix) = y(iy)
               y(iy) = temp
               ix = ix + incx
            END DO
         ELSE
            iy = 1-(n-1)*incy
            DO i = 1, n
               temp = x(ix)
               x(ix) = y(iy)
               y(iy) = temp
               iy = iy + incy
               ix = ix + incx
            END DO
         END IF
      END IF
!
      RETURN
   END SUBROUTINE DSWAP
!
!-------------------------------------------------------------------------------
!  FUNCTION IDAMAX
!-------------------------------------------------------------------------------
!
   FUNCTION IDAMAX(n, x, incx) RESULT(retval)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN) :: n, incx
      REAL(r8),    INTENT(IN) :: x(*)    ! x( * )
!
!     idamax returns the smallest value of i such that
!
!               abs(x(i)) = max(abs(x(j)))
!                            j
!
!     nag fortran 77 version of the blas routine idamax.
!     nag fortran 77 o( n ) basic linear algebra routine.
!
!  -- written on 31-may-1983.
!     sven hammarling, nag central office.
!
      REAL(r8) :: xmax
      INTEGER(i4) :: i, imax, ix, retval
!!!
!
      IF( n < 1 )THEN
         retval = 0    ! idamax = 0
         RETURN
      END IF
!
      imax = 1
      IF (n > 1)THEN
         xmax = ABS(x(1))
         ix  = 1
         DO i = 2, n
            ix = ix + incx
            IF (xmax < ABS(x(ix))) THEN
               xmax = ABS(x(ix))
               imax = i
            END IF
         END DO
      END IF
!
      retval = imax   ! idamax = imax
!
      RETURN
   END FUNCTION IDAMAX
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DLOAD
!-------------------------------------------------------------------------------
!
   SUBROUTINE DLOAD(n, const, x, incx)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n, incx
      REAL(r8),    INTENT(IN)    :: const
      REAL(r8),    INTENT(INOUT) :: x(*)
!
!     dload  performs the operation
!     x = const*e,   e' = ( 1  1 ... 1 ).
!     nag fortran 77 o( n ) basic linear algebra routine.
!  -- written on 22-september-1983.
!     sven hammarling, nag central office.
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      INTEGER(i4) :: ix
!!!
!
      IF (n < 1) RETURN
      IF (const /= zero) THEN
         DO ix = 1, 1+(n-1)*incx, incx
            x(ix) = const
         END DO
      ELSE
         DO ix = 1, 1+(n-1)*incx, incx
            x(ix) = zero
         END DO
      END IF
!
      RETURN
   END SUBROUTINE DLOAD
!
!-------------------------------------------------------------------------------
!  SUBROUTINE MAXPY
!-------------------------------------------------------------------------------
!
   SUBROUTINE MAXPY(nrow, ncol, alpha, xmat, nrowy, ymat)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: nrow, ncol
      REAL(r8),    INTENT(IN)    :: alpha
      REAL(r8),    INTENT(IN)    :: xmat(nrow, ncol)
      INTEGER(i4), INTENT(IN)    :: nrowy
      REAL(r8),    INTENT(INOUT) :: ymat(nrowy, ncol)
!
!     Subroutine maxpy takes as input the scalar alpha and two matrices,
!     xmat and ymat.  xmat has declared row dimension nrow, and
!     ymat has declared row dimension nrowy, but both are
!     conceptually nrow by ncol.
!     On output, (new ymat) is alpha*xmat+ (old ymat), by analogy
!     with the vector blas routine saxpy.
!
      INTEGER(i4) :: i, j
!
      DO j = 1, ncol
         DO i = 1, nrow
            ymat(i,j) = ymat(i,j) + alpha*xmat(i,j)
         END DO
      END DO
!
      RETURN
   END SUBROUTINE MAXPY
!
!-------------------------------------------------------------------------------
!  SUBROUTINE MATCOP
!-------------------------------------------------------------------------------
!
   SUBROUTINE MATCOP(nrow1, nrow2, nrow, ncol, xmat1, xmat2)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: nrow1, nrow2
      INTEGER(i4), INTENT(IN)    :: nrow, ncol
      REAL(r8),    INTENT(IN)    :: xmat1(nrow1, ncol)
      REAL(r8),    INTENT(INOUT) :: xmat2(nrow2, ncol)
!
!     Given 2 matrices xmat1 and xmat2, where xmat1 has declared
!     row dimension nrow1, xmat2 has declared row dimension nrow2,
!     and both have column dimension ncol, the routine matcop copies
!     rows 1 through nrow, and columns 1 through ncol from xmat1 into
!     xmat2.
!
      INTEGER(i4) :: i, j
!!!
!
      IF (nrow <= 0 .OR. ncol <= 0) RETURN
      DO j = 1, ncol
         DO i = 1, nrow
            xmat2(i,j) = xmat1(i,j)
         END DO
      END DO
!
      RETURN
   END SUBROUTINE MATCOP
!
!-------------------------------------------------------------------------------
!  SUBROUTINE MTLOAD
!-------------------------------------------------------------------------------
!
   SUBROUTINE MTLOAD(nrow, ncol, const, nrowx, xmat)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: nrow, ncol
      REAL(r8),    INTENT(IN)    :: const
      INTEGER(i4), INTENT(IN)    :: nrowx
      REAL(r8),    INTENT(INOUT) :: xmat(nrowx, ncol)
!
!     mtload sets elements 1 through nrow, 1 through ncol, of the
!     matrix xmat (whose declared row dimension is nrowx) to the
!     scalar value const.
!
      INTEGER(i4) :: i, j
!!!
!
      IF (nrow <= 0 .OR. ncol <= 0) RETURN
!
      DO j = 1, ncol
         DO i = 1, nrow
            xmat(i,j) = const
         END DO
      END DO
!
      RETURN
   END SUBROUTINE MTLOAD
!
!-------------------------------------------------------------------------------
!  SUBROUTINE MTLOAD1
!-------------------------------------------------------------------------------
!
   SUBROUTINE MTLOAD1(nrow, ncol, xx, nrowx, xmat)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: nrow, ncol
      INTEGER(i4), INTENT(IN)    :: nrowx
      REAL(r8),    INTENT(INOUT) :: xx(*)
      REAL(r8),    INTENT(INOUT) :: xmat(nrowx, ncol)
!
!     mtload sets elements 1 through nrow, 1 through ncol, of the
!     matrix xmat (whose declared row dimension is nrowx) to the
!     scalar value const.
!
      INTEGER(i4) :: i, j
!!!
!
      IF (nrow <= 0 .OR. ncol <= 0) RETURN
!
      DO j = 1, ncol
         DO i = 1, nrow
            IF (i == 1) THEN
               xmat(i,j) = 2.D0 * xx(i-1) - 1.D0
            ELSE IF (i == 2) THEN
               xmat(i,j) = 2.D0
            ELSE
               xmat(i,j) = 0.D0
            END IF
         END DO
      END DO
!
      RETURN
   END SUBROUTINE MTLOAD1
!
!-------------------------------------------------------------------------------
!  SUBROUTINE MSSQ
!-------------------------------------------------------------------------------
!
   SUBROUTINE MSSQ(nrow, ncol, xmat, scale, sumsq)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: nrow, ncol
      REAL(r8),    INTENT(IN)    :: xmat(nrow,*)
      REAL(r8),    INTENT(INOUT) :: scale, sumsq
!
!     Given the nrow by ncol matrix xmat, mssq returns values
!     scale and sumsq such that (scale**2) * sumsq = sum of squares of
!     xmat(i,j), where  scale = max  abs(xmat(i,j)).
!     mssq is a stripped-down matrix version of the blas routine sssq.
!
      REAL(r8), PARAMETER :: one   = 1.0D+0
      REAL(r8), PARAMETER :: zero  = 0.0D+0
      REAL(r8) :: absxij
      INTEGER(i4) :: i, j
!!!
!
      scale = zero
      sumsq = one
      IF (nrow >= 1 .AND. ncol >= 1) THEN
         DO i = 1, nrow
            DO j = 1, ncol
               IF (xmat(i,j) /= zero) THEN
                  absxij = ABS(xmat(i,j))
                  IF (scale < absxij) THEN
                     sumsq = one + sumsq* (scale/absxij)**2
                     scale = absxij
                  ELSE
                     sumsq = sumsq + (absxij/scale)**2
                  END IF
               END IF
            END DO
         END DO
      END IF
!
      RETURN
   END SUBROUTINE mssq
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DSSQ
!-------------------------------------------------------------------------------
!
   SUBROUTINE DSSQ(n, x, incx, scale, sumsq)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(IN)    :: x(*)
      INTEGER(i4), INTENT(IN)    :: incx
      REAL(r8),    INTENT(INOUT) :: scale, sumsq
!
!     Given the n-vector x, dssq returns values scale and sumsq such that
!     (scale**2) * sumsq = sum of squares of x(i),
!     where  scale = max  abs(x(i)).
!
      REAL(r8), PARAMETER :: one   = 1.0D+0
      REAL(r8), PARAMETER :: zero  = 0.0D+0
      REAL(r8) :: absxi
      INTEGER(i4) :: ix
!!!
!
      scale = zero
      sumsq = one
      IF (n >= 1) THEN
         DO ix = 1, 1+(n-1)*incx, incx
            IF (x(ix) /= zero) THEN
               absxi = ABS(x(ix))
               IF (scale < absxi) THEN
                  sumsq = one + sumsq*(scale/absxi)**2
                  scale = absxi
               ELSE
                  sumsq = sumsq + (absxi/scale)**2
               END IF
            END IF
         END DO
      END IF
!
      RETURN
   END SUBROUTINE DSSQ
!
END

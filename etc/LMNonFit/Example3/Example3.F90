PROGRAM Example3
!
!-------------------------------------------------------------------------------
! 
!  Code converted using TO_F90 by Alan Miller
!  Date: 1999-12-09  Time: 23:25:02
!
!  This program tests codes for the least-squares solution of m nonlinear
!  equations in n variables. It consists of a driver and an interface 
!  subroutine FCN. The driver reads in data, calls the nonlinear least-squares
!  solver, and finally prints out information on the performance of the solver.
!  This is only a sample driver, many other drivers are possible. The interface
!  subroutine FCN is necessary to take into account the forms of calling
!  sequences used by the function and jacobian subroutines in the various
!  nonlinear least-squares solvers.
!
!  SUBPROGRAMS CALLED
!
!  USER-SUPPLIED ...... FCN
!
!  MINPACK-SUPPLIED ... DPMPAR,ENORM,INITPT,LMDER1,SSQFCN

!  FORTRAN-SUPPLIED ... SQRT
!
!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!-------------------------------------------------------------------------------
!
   USE Base,         ONLY: i4, r4, r8
   USE CommonRefnum
   USE LMNonFit,     ONLY: LMDIF1, ENORM
   USE LMExtern,     ONLY: INITPT, FCN, SSQFCN
!
   IMPLICIT NONE
!
   REAL(r8) :: fnm(60), fvec(65), x(40)
   REAL(r8) :: factor, fnorm1, fnorm2, tol
   INTEGER(i4) :: iwa(40), ma(60), na(60), nf(60), nj(60), np(60), nx(60)
   INTEGER(i4) :: i, ic, info, k, m, n, ntries

!  LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.
!  LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.

   INTEGER(i4), PARAMETER :: NREAD = 5, NWRITE = 6
   REAL(r8), PARAMETER :: one = 1.0_r8, ten = 10.0_r8
!!!
!
   tol = SQRT(EPSILON(one))
   ic = 0
10 READ(NREAD,50) nprob, n, m, ntries
   IF (nprob <= 0) GO TO 30
!
   factor = one
!
   DO k = 1, ntries
      ic = ic + 1
      CALL INITPT(n, x, nprob, factor)
      CALL SSQFCN(m, n, x, fvec, nprob)
      fnorm1 = ENORM(m, fvec)
      WRITE(NWRITE,60) nprob, n, m
      nfev = 0
      njev = 0
      CALL LMDIF1(m, n, x, fvec, tol, info, iwa)
      CALL SSQFCN(m, n, x, fvec, nprob)
      fnorm2 = ENORM(m, fvec)
      np(ic) = nprob
      na(ic) = n
      ma(ic) = m
      nf(ic) = nfev
      njev = njev/n
      nj(ic) = njev
      nx(ic) = info
      fnm(ic) = fnorm2
      WRITE(NWRITE,70) fnorm1, fnorm2, nfev, njev, info, x(1:n)
      factor = ten * factor
   END DO
   GO TO 10

30 WRITE(NWRITE,80) ic
   WRITE(NWRITE,90)
   DO i=1,ic
      WRITE(nwrite,100) np(i), na(i), ma(i), nf(i), nj(i), nx(i), fnm(i)
   END DO

50 FORMAT (4I5)
60 FORMAT (////'      PROBLEM', i5, '      DIMENSIONS', 2I5//)
70 FORMAT ('      INITIAL L2 NORM OF THE RESIDUALS', g15.7//  &
           '      FINAL L2 NORM OF THE RESIDUALS  ', g15.7//  &
           '      NUMBER OF FUNCTION EVALUATIONS  ', i10//  &
           '      NUMBER OF JACOBIAN EVALUATIONS  ', i10//  &
           '      EXIT PARAMETER', t39, i10//  &
           '      FINAL APPROXIMATE SOLUTION'// (t6, 5g15.7))
80 FORMAT (' SUMMARY OF ', i3, ' CALLS TO LMDER1'/)
90 FORMAT (' NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'/)
100 FORMAT (3I5, 3I6, ' ', g15.7)
!
END PROGRAM Example3

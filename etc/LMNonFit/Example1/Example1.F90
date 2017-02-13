PROGRAM Example1
!
!  Fit the 4-parameter logistic curve:
!     Y = A + C / [1 + exp{-B(X - D)}]
!  by unweighted least squares.
!  Latest revision - 25 October 2001
!
   USE LMNonFit,   ONLY: LMDER1
   USE LMExtern,   ONLY: FCN
   USE DataStorage
!
   IMPLICIT NONE
!
   REAL(r8), ALLOCATABLE :: fvec(:), fjac(:,:)
   INTEGER(i4), PARAMETER :: n = 4                 ! The number of parameters
   INTEGER(i4) :: info, iostatus, ipvt(n), m
   REAL(r8) :: p(n), tol = 1.0E-04_r8
!
!  Read in the data, counting the number of cases.
!
   OPEN(UNIT=8, FILE='lgstic4.dat', STATUS='OLD')
   m = 1
   DO
      READ(8, *, IOSTAT=iostatus) x(m), y(m)
      IF (iostatus < 0) EXIT
      IF (iostatus > 0) CYCLE    ! Skips cases with errors - e.g. a heading
      m = m + 1
   END DO
   m = m - 1
   ALLOCATE(fvec(m), fjac(m,n))
!
!  Set starting values for parameters.
!  A = p(1), B = p(2), C = p(3), D = p(4)
!
   p(1) = MINVAL( y(1:m) )
   p(3) = MAXVAL( y(1:m) ) - p(1)
   p(2) = 2.0_r8 / ( MAXVAL( x(1:m) ) - MINVAL( x(1:m) ) )
   p(4) = 0.5_r8 * ( MAXVAL( x(1:m) ) + MINVAL( x(1:m) ) )
!
   CALL LMDER1(m, n, p, fvec, fjac, tol, info, ipvt)
!
   SELECT CASE (info)
      CASE (:-1)
         WRITE(*, *) 'Users FCN returned INFO = ', -info
      CASE (0)
         WRITE(*, *) 'Improper values for input parameters'
      CASE (1:3)
         WRITE(*, *) 'Convergence'
         WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
      CASE (4)
         WRITE(*, *) 'Residuals orthogonal to the Jacobian'
         WRITE(*, *) 'There may be an error in FCN'
         WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
      CASE (5)
         WRITE(*, *) 'Too many calls to FCN'
         WRITE(*, *) 'Either slow convergence, or an error in FCN'
         WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
      CASE (6:7)
         WRITE(*, *) 'TOL was set too small'
         WRITE(*, '(a, 4g13.5)') ' Final values of A, B, C, D: ', p
      CASE DEFAULT
         WRITE(*, *) 'INFO =', info, ' ???'
   END SELECT
!
END PROGRAM Example1

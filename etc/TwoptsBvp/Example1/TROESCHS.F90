PROGRAM TROESCHS
!
   USE Base,      ONLY: i4, r8
   USE BVPShared, ONLY: pdebug, use_c, comp_c,  &
                        uval0, nminit, iprint, idum,  &
                        MTLOAD
   USE BVPExtern, ONLY: INITU, FSUB, DFSUB, GSUB, DGSUB
   USE BVPC,      ONLY: TWPBVPC_INIT, TWPBVPC
!
   IMPLICIT NONE
!
!     INTEGER LWRKFL, LWRKIN, NUDIM, NXXDIM, LISERIES, ind,
!    + IWRK, NTOL, NLBC, NMSH, NFXPNT, NCOMP, LTOL, NMAX, IPAR,
!    + IFLBVP, NMINIT, IDUM, IPRINT
!     DOUBLE PRECISION PI, TOL, ETOL, WRK, UVAL0, ALEFT, ARIGHT,
!    + FIXPNT, XX, U, CKAPPA1, GAMMA1, CKAPPA, RPAR

   INTEGER(i4), PARAMETER :: NUDIM = 2, NXXDIM = 100000
   INTEGER(i4), PARAMETER :: LWRKFL = 250000
   INTEGER(i4), PARAMETER :: LWRKIN = 50000
   INTEGER(i4), PARAMETER :: LISERIES = 10
!
   REAL(r8) :: wrk(LWRKFL)
   REAL(r8) :: u(NUDIM,NXXDIM), xx(NXXDIM)
   REAL(r8) :: tol(NUDIM)
   REAL(r8) :: rpar(1), fixpnt(1)
   INTEGER(i4) :: ltol(NUDIM), ipar(1)
   INTEGER(i4) :: iwrk(LWRKIN)
   INTEGER(i4) :: ind, nfxpnt
!
   REAL(r8)    :: etol, aleft, aright, ckappa1, gamma1, ckappa
   INTEGER(i4) :: ncomp, ntol, nlbc, nmsh, nmax, iflbvp
   LOGICAL     :: linear, giveu, givmsh
!!!
!      
   WRITE(6,'(A)') 'PROGRAM TROESCHS'
!
   CALL TWPBVPC_INIT
!
!  If you do not like to use the conditioning in the mesh selection
!  use_c  = .false.
!
!  If you do not like to compute the conditioning parameters
!  comp_c = .false.     
!
   use_c  =  .true.
   comp_c =  .true.
!     
!  RPAR (real) and IPAR (integer) are arrays that are passed to
!  the F function, G function, and their Jacobians.
!  RPAR is typically used as a parameter in the problem formulation, and
!  IPAR is normally used to select a problem from a set.
!
!  For our problem:
!  RPAR(1) = mu in Troesch's, increasing mu increases the problem stiffness      
!
   rpar(1) = 5._r8
      
!  WRK is the floating point workspace and IWRK
!  is the integer work space. We set these to zero here,
!  to ensure consistent behaviour with different compilers.
!
   DO ind = 1, LWRKFL
      WRK(ind) = 0.D0
   END DO
   DO ind = 1, LWRKIN
      IWRK(ind) = 0
   END DO
!
!  ETOL is the required error tolerance of the solution.
!  Decreasing ETOL will give a more accurate solution, but
!  more mesh points are likely to be required and the code
!  will take longer to finish.
!
   etol = 1.e-4_r8
! 
!  IPRINT controls whether or not the solver prints out
!  verbose output. A value of -1 disables these diagnostics
!
   iprint = 1
      
!  NTOL is the number of components that have to satisfy
!  the prescribed tolerance.
!  LTOL is the component that needs to be tested.
!  Most of the time one will set NTOL to the number of system components,
!  LTOL(i) to component i, and set the tolerance for each component TOL to be equal.
!
   ntol = 2     
   tol(1) = etol
   DO ind = 1, ntol
      ltol(ind) = ind
      TOL(ind) = tol(1)
   END DO      
! 
!  ALEFT and ARIGHT are the values of x at the left
!  and right boundaries.
!
   ALEFT  = 0._r8
   ARIGHT = 1._r8
!      
!  The number of components in the system
!
   ncomp = 2      
!
!  The number of boundary conditions at the left, the
!  number of the right being equal to NCOMP-NLBC
!
   nlbc = 1
! 
!  nmsh is the number of initial points we supply to
!  the solver in the form of a guess, we set this
!  to zero to indicate that we have no initial guess
!
   nmsh  = 0
!
!  fixed points are x values which are required by the
!  user to be part of the returned mesh.
!  nfxpnt is the number of fixed points, which we set to be zero
!
   nfxpnt= 0
!
!  The initial approximation is equal to UVAL0      
!
   uval0 = 0.D0  
!
!  the problem is nonlinear so we specify .false. here
!
   linear = .FALSE.
! 
!  we do not supply any initial guesses, so again we
!  choose .false.
!
   giveu  = .FALSE.
!
!  No initial mesh (a set of x values) are given, so
!  this is .false. too
!
   givmsh = .FALSE.
!
!  PDEBUG controls the low level solver diagnostics,
!  most of the time this should be set to .false.
!
   pdebug = .TRUE.
! 
   CALL TWPBVPC(ncomp, nlbc, aleft, aright, nfxpnt, fixpnt, ntol,  &
                ltol, tol, linear, givmsh, giveu, nmsh, nxxdim,  &
                xx, nudim, u, nmax, lwrkfl, wrk, lwrkin, iwrk,  &
                ckappa1, gamma1, ckappa, rpar, ipar, iflbvp)
!
!  When returning from TWPBVPC, one should immediately
!  check IFLBVP, to see whether or not the solver was
!  succesful
!     
   IF (iflbvp < 0) THEN
      WRITE(6,*) 'The code failed to converge!'
      STOP
   END IF
!     
   WRITE(6,*) 'Number of mesh points = ', nmsh
   WRITE(6,*) 'Dumping (x,y) mesh now...'
!      
!  the solution x values are stored in XX, the Y
!  values are stored in u.
!  u(i,j) refers to component i of point j in the mesh.      
!
   DO ind = 1 ,nmsh
      WRITE(6,900) xx(ind), u(1,ind)
   END DO

900 FORMAT(G14.7,G14.7)
!
END

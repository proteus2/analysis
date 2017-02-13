PROGRAM ELASTICA
!
   USE Base,      ONLY: i4, r8
   USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                        nminit, iprint, idum
   USE BVPExtern, ONLY: FSUB, DFSUB, GSUB, DGSUB
   USE BVPC,      ONLY: TWPBVPC_INIT, TWPBVPC
!
   IMPLICIT NONE
!
!
!     INTEGER LWRKFL, LWRKIN, NUDIM, NXXDIM, LISERIES, IPRINT, ind,
!    + IWRK, NTOL, NLBC, NMSH, NFXPNT, NCOMP, LTOL, NMAX, IPAR,
!    + IFLBVP, NMINIT, IDUM
!     DOUBLE PRECISION TOL, ETOL, WRK, UVAL0, ALEFT, ARIGHT,
!    + FIXPNT, XX, U, CKAPPA1, GAMMA1, CKAPPA, RPAR

   INTEGER(i4), PARAMETER :: LWRKFL=250000
   INTEGER(i4), PARAMETER :: LWRKIN=50000
   INTEGER(i4), PARAMETER :: NUDIM=5
   INTEGER(i4), PARAMETER :: NXXDIM=100000
   INTEGER(i4), PARAMETER :: LISERIES=10
!
   REAL(r8) :: u(NUDIM,NXXDIM), xx(NXXDIM)
   REAL(r8) :: tol(NUDIM), fixpnt(1), rpar(1)
   REAL(r8) :: wrk(LWRKFL), aleft, aright, etol, ckappa1, ckappa, gamma1
   INTEGER(i4) :: ltol(NUDIM), ipar(1)
   INTEGER(i4) :: iwrk(LWRKIN), ncomp, ntol, nmax, iflbvp, ind, nlbc, nfxpnt, nmsh
   LOGICAL :: linear, giveu, givmsh
!!!
!
   WRITE(*,'(/,A,/)') 'PROGRAM ELASTICA'
!
   CALL TWPBVPC_INIT
!
!  RPAR (real) and IPAR (integer) are arrays that are passed to
!  the F function, G function, and their Jacobians.
!  RPAR is typically used as a parameter in the problem formulation, and
!  IPAR is normally used to select a problem from a set.
!
!  We have no parameters in our problem, so IPAR and RPAR are left alone.
!    
!  If you do not like to use the conditioning in the mesh selection
!  USE_C  = .false.
!  If you do not like to compute the conditioning parameters
!  COMP_C = .false.     
!
   use_c  =  .true.
   comp_c =  .true.
!                
!  WRK is the floating point workspace and IWRK
!  is the integer work space. We set these to zero here,
!  to ensure consistent behaviour with different compilers.      
!
   DO ind = 1, LWRKFL
      WRK(ind) = 0d0
   END DO
   DO ind = 1, LWRKIN
      IWRK(ind) = 0
   END DO
!
!  ALEFT and ARIGHT are the values of x at the left
!  and right boundaries.      
!
   aleft  = 0.0_r8
   aright = 5.e-1_r8       
!
!  ETOL is the required error tolerance of the solution.
!  Decreasing ETOL will give a more accurate solution, but
!  more mesh points are likely to be required and the code
!  will take longer to finish.
!
   etol = 1.e-12_r8
!
!  NTOL is the number of components that have to satisfy
!  the prescribed tolerance.
!  LTOL is the component that needs to be tested.
!  Most of the time one will set NTOL to the number of system components,
!  LTOL(i) to component i, and set the tolerance for each component TOL to be equal.
!
   ntol = 5
   tol(1) = etol
   DO ind = 1, ntol
      ltol(ind) = ind
      tol(ind) = tol(1)
   END DO      
!
!  IPRINT controls whether or not the solver prints out
!  verbose output. A value of -1 disables these diagnostics      
!
   iprint = 1      
!
!  The number of components in the system      
!
   ncomp = 5
!
!  The number of boundary conditions at the left, the
!  number of the right being equal to NCOMP-NLBC
!
   nlbc = 3
!
!  nmsh is the number of initial points we supply to
!  the solver in the form of a guess, we set this
!  to zero to indicate that we have no initial guess      
!
   nmsh = 16
!
!  fixed points are x values which are required by the
!  user to be part of the returned mesh.
!  NFXPNT is the number of fixed points, which we set to be zero            
!
   nfxpnt = 0
!
!  The initial approximation is equal to UVAL0      
!
   uval0 = 0.D0  
!
!  the problem is nonlinear so we specify .false. here
   linear = .FALSE.
!      
!  we do not supply any initial guesses, so again we
!  choose .false.      
   giveu  = .FALSE.
      
!  No initial mesh (a set of x values) are given, so
!  this is .false. too      
   givmsh = .FALSE.
      
!  PDEBUG controls the low level solver diagnostics,
!  most of the time this should be set to .false.      
!
   pdebug = .TRUE.
!           
   CALL TWPBVPC(ncomp, nlbc, aleft, aright, nfxpnt, fixpnt, ntol,  &
                ltol, tol, linear, givmsh, giveu, nmsh, nxxdim,  &
                xx, NUDIM, u, nmax, LWRKFL, wrk, LWRKIN, iwrk,  &
                ckappa1, gamma1, ckappa, rpar, ipar, iflbvp) 
     
!     When returning from TWPBVPC, one should immediately
!     check IFLBVP, to see whether or not the solver was
!     succesful    
      IF (iflbvp .LT. 0) THEN
         WRITE(6,*) 'The code failed to converge!'
         STOP
      END IF
     
!     the solution x values are stored in XX, the Y
!     values are stored in U.
!     U(i,j) refers to component i of point j in the mesh.           
     
      WRITE(6,*) 'Number of mesh points = ', NMSH
      WRITE(6,*) 'P = ', U(5,1)
      WRITE(6,*) 'Theta0 = ', U(3,1)
      
      WRITE(6,*) 'Dumping (x,y) mesh now...'
      
      DO IND=1,NMSH
         WRITE(6,900) U(1,IND), U(2,IND)
900    FORMAT(F10.6, G14.6)
      END DO
      
END

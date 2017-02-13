      PROGRAM  TROESCHS
      IMPLICIT NONE

      WRITE(*,'(/,A,/)') 'PROGRAM TROESCHS'
      CALL RUNPROB()
   
      STOP
      END


C  This routine must be provided to reset u after re-meshing 
C  for linear problems or for nonlinear problems
C  when interpolation of the old solution is not used.

      SUBROUTINE initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
      IMPLICIT NONE
      INTEGER ncomp, nmsh, nudim, rpar, ipar, nminit, iprint, idum
      DOUBLE PRECISION xx, u, uval0
      DIMENSION rpar(*),ipar(*)
      DIMENSION xx(*), u(nudim, *)
      
      LOGICAL pdebug, use_c, comp_c
      COMMON/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

C  This version sets all elements of u to the constant uval0,
C  using a presupplied routine mtload.
      CALL mtload(ncomp, nmsh, uval0, nudim, u)

      RETURN
      END     
      
      
      SUBROUTINE FSUB(NCOMP,X,Z,F,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER NCOMP, IPAR
      DOUBLE PRECISION F, Z, RPAR, X, MU
      DIMENSION Z(*),F(*)
      DIMENSION RPAR(*), IPAR(*)
      
      MU = RPAR(1)
      
      F(1)=Z(2)
      F(2)=MU*sinh(mU*Z(1))

      RETURN
      END      
      
      SUBROUTINE DFSUB(NCOMP,X,Z,DF,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER NCOMP, IPAR, I, J
      DOUBLE PRECISION X, Z, DF, RPAR, MU
      DIMENSION Z(*),DF(NCOMP,*)
      DIMENSION RPAR(*), IPAR(*)
      
      MU = RPAR(1)
      
C     dF1/dZ1
      DF(1,1)=0.D0
C     dF1/dZ2
      DF(1,2)=1.D0
      
C     dF2/dZ1
      DF(2,1)=MU*MU*cosh(MU*Z(1))
C     dF2/dZ2
      DF(2,2)=0.D0
      
      RETURN
      END
      
      SUBROUTINE GSUB(I,NCOMP,Z,G,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER I, NCOMP, IPAR
      DOUBLE PRECISION Z, RPAR, G      
      DIMENSION Z(*)
      DIMENSION RPAR(*), IPAR(*)
      
C     I == the boundary condition "number".
C     The conditions at the left are enumerated first,
C     then the ones at the right.
C     The order of left or right conditions does not matter,
C     but we must be consistent when defining the jacobian
C     of the boundary conditions!

C     We have specified 1 left bc, and 2 bcs total.
C     This means that:
C     BC(1) = the condition at the "left"
C     BC(2) = the condition at the "right"

      IF (I.EQ.1) G=Z(1)
      IF (I.EQ.2) G=Z(1)-1.D0
      
      RETURN
      END
    
      SUBROUTINE DGSUB(I,NCOMP,Z,DG,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER I, NCOMP, IPAR
      DOUBLE PRECISION Z, DG, RPAR
      DIMENSION Z(*),DG(*)
      DIMENSION RPAR(*), IPAR(*)
      
C     I == the boundary condition "number".
C     The conditions at the left are enumerated first,
C     then the ones at the right.
C     The order of left or right conditions does not matter,
C     but we must be consistent when defining the jacobian
C     of the boundary conditions!
                 
      IF (I.EQ.1) THEN
C     dG1/dZ1      
         DG(1) = 1.D0
C     dG1/dZ2
         DG(2) = 0.D0
      END IF
            
      IF (I.EQ.2) THEN
C     dG2/dZ1
         DG(1)=1.D0      
C     dG2/dZ2         
         DG(2)=0.D0
      END IF
     
      RETURN
      END 

      SUBROUTINE RUNPROB()
      IMPLICIT NONE
      
      INTEGER LWRKFL, LWRKIN, NUDIM, NXXDIM, LISERIES, ind,
     + IWRK, NTOL, NLBC, NMSH, NFXPNT, NCOMP, LTOL, NMAX, IPAR,
     + IFLBVP, NMINIT, IDUM, IPRINT
      DOUBLE PRECISION PI, TOL, ETOL, WRK, UVAL0, ALEFT, ARIGHT,
     + FIXPNT, XX, U, CKAPPA1, GAMMA1, CKAPPA, RPAR

      PARAMETER(LWRKFL=250000)
      PARAMETER(LWRKIN=50000)
      DIMENSION WRK(LWRKFL),IWRK(LWRKIN)
      PARAMETER(NUDIM=2,NXXDIM=100000)
      DIMENSION U(NUDIM,NXXDIM), XX(NXXDIM)
      DIMENSION TOL(NUDIM), LTOL(NUDIM), FIXPNT(1)
      DIMENSION RPAR(1), IPAR(1)
      PARAMETER (LISERIES=10)
      LOGICAL   LINEAR, GIVEU, GIVMSH
      EXTERNAL  TWPBVPC, FSUB, DFSUB, GSUB, DGSUB
      
      LOGICAL PDEBUG, use_c, comp_c
c this common need to be defined in order to run TWPBVPC successfully
c give information about some parameters      
      COMMON /ALGPRS/ NMINIT,PDEBUG,IPRINT,IDUM,UVAL0,use_c,comp_c
      
C     If you do not like to use the conditioning in the mesh selection
C     USE_C  = .false.
C     If you do not like to compute the conditioning parameters
C     COMP_C = .false.     

      USE_C  =  .true.
      COMP_C =  .true.
      
C     RPAR (real) and IPAR (integer) are arrays that are passed to
C     the F function, G function, and their Jacobians.
C     RPAR is typically used as a parameter in the problem formulation, and
C     IPAR is normally used to select a problem from a set.
C
C     For our problem:
C     RPAR(1) = mu in Troesch's, increasing mu increases the problem stiffness      
      RPAR(1) = 5.D0
      
C     WRK is the floating point workspace and IWRK
C     is the integer work space. We set these to zero here,
C     to ensure consistent behaviour with different compilers.
      DO ind=1,LWRKFL
         WRK(ind) = 0.D0
      ENDDO
      DO ind=1,LWRKIN
         IWRK(ind) = 0
      ENDDO
      
      
C     ETOL is the required error tolerance of the solution.
C     Decreasing ETOL will give a more accurate solution, but
C     more mesh points are likely to be required and the code
C     will take longer to finish.
      ETOL = 1.D-4
      
      
C     IPRINT controls whether or not the solver prints out
C     verbose output. A value of -1 disables these diagnostics
      IPRINT = 1
      
C     NTOL is the number of components that have to satisfy
C     the prescribed tolerance.
C     LTOL is the component that needs to be tested.
C     Most of the time one will set NTOL to the number of system components,
C     LTOL(i) to component i, and set the tolerance for each component TOL to be equal.
      NTOL = 2     
      TOL(1) = ETOL
      DO ind=1,NTOL
        LTOL(ind)=ind
        TOL(ind) =TOL(1)
      ENDDO      
      
C     ALEFT and ARIGHT are the values of x at the left
C     and right boundaries.
      ALEFT  = 0.D0
      ARIGHT = 1.D0
      
C     The number of components in the system
      NCOMP = 2      
C     The number of boundary conditions at the left, the
C     number of the right being equal to NCOMP-NLBC
      NLBC=1
          
C     NMSH is the number of initial points we supply to
C     the solver in the form of a guess, we set this
C     to zero to indicate that we have no initial guess
      NMSH  = 0

C     fixed points are x values which are required by the
C     user to be part of the returned mesh.
C     NFXPNT is the number of fixed points, which we set to be zero
      NFXPNT= 0
C     The initial approximation is equal to UVAL0      
      UVAL0 = 0.D0  

C     the problem is nonlinear so we specify .false. here
      LINEAR = .FALSE.
      
C     we do not supply any initial guesses, so again we
C     choose .false.
      GIVEU  = .FALSE.
      
C     No initial mesh (a set of x values) are given, so
C     this is .false. too
      GIVMSH = .FALSE.
      
C     PDEBUG controls the low level solver diagnostics,
C     most of the time this should be set to .false.
      PDEBUG = .TRUE.
         
      CALL TWPBVPC(NCOMP,NLBC,ALEFT,ARIGHT,NFXPNT,FIXPNT,NTOL,
     +            LTOL,TOL,LINEAR,GIVMSH,GIVEU,NMSH,NXXDIM,
     +            XX,NUDIM,U,NMAX,LWRKFL,WRK,LWRKIN,IWRK,
     +            FSUB,DFSUB,GSUB,DGSUB,
     +            ckappa1,gamma1,ckappa,rpar,ipar,IFLBVP)

C     When returning from TWPBVPC, one should immediately
C     check IFLBVP, to see whether or not the solver was
C     succesful
     
      IF (IFLBVP .LT. 0) THEN
         WRITE(6,*) 'The code failed to converge!'
         RETURN
      END IF
     
      WRITE(6,*) 'Number of mesh points = ', NMSH
      WRITE(6,*) 'Dumping (x,y) mesh now...'
      
C     the solution x values are stored in XX, the Y
C     values are stored in U.
C     U(i,j) refers to component i of point j in the mesh.      
      
      DO IND=1,NMSH
         WRITE(6,900) XX(IND), U(1,IND)
  900    FORMAT (G14.8,G14.8)
      END DO
      
      RETURN
      END


MODULE BVPExtern
!
   USE Base,      ONLY: i4, r8
   USE BVPShared, ONLY: pdebug, use_c, comp_c,  &
                        uval0, nminit, iprint, idum,  &
                        MTLOAD
!
   IMPLICIT NONE
!
   PRIVATE
   PUBLIC  :: INITU, FSUB, DFSUB, GSUB, DGSUB
!
CONTAINS
!
   SUBROUTINE INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
!
      USE Base,      ONLY: i4, r8
      USE BVPShared, ONLY: pdebug, use_c, comp_c,  &
                           uval0, nminit, iprint, idum,  &
                           MTLOAD
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN) :: ncomp, nudim
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(*), u(ncomp,*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!!!
!
!     This version sets all elements of u to the constant uval0,
!     using a presupplied routine mtload.
!
      CALL MTLOAD(ncomp, nmsh, uval0, nudim, u)
!
      RETURN
   END SUBROUTINE INITU
!!!
!!!
   SUBROUTINE FSUB(ncomp, x, z, f, rpar, ipar)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      REAL(r8),    INTENT(IN)    :: x, z(*)
      REAL(r8),    INTENT(INOUT) :: f(*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8) :: mu
!!!
!
      mu = rpar(1)
! 
      f(1) = z(2)
      f(2) = mu*SINH(mu*z(1))
!
      RETURN
   END SUBROUTINE FSUB
!!!
!!!
   SUBROUTINE DFSUB(ncomp, x, z, df, rpar, ipar)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      REAL(r8),    INTENT(IN)    :: x, z(*)
      REAL(r8),    INTENT(INOUT) :: df(ncomp,*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8) :: mu
!!!
!
      mu = rpar(1)
!
!     dF1/dZ1
      df(1,1)=0.D0
!     dF1/dZ2
      df(1,2)=1.D0
!     dF2/dZ1
      df(2,1)=mu*mu*COSH(mu*z(1))
!     dF2/dZ2
      df(2,2)=0.D0
!
      RETURN
   END SUBROUTINE DFSUB
!!!
!!!
!!!
!!!
   SUBROUTINE GSUB(i, ncomp, z, g, rpar, ipar)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: i, ncomp
      REAL(r8),    INTENT(IN)    :: z(*)
      REAL(r8),    INTENT(INOUT) :: g, rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
!     I == the boundary condition "number".
!     The conditions at the left are enumerated first,
!     then the ones at the right.
!     The order of left or right conditions does not matter,
!     but we must be consistent when defining the jacobian
!     of the boundary conditions!
!
!     We have specified 1 left bc, and 2 bcs total.
!     This means that:
!     BC(1) = the condition at the "left"
!     BC(2) = the condition at the "right"
!
!!!
!
      IF (I == 1) g=z(1)
      IF (I == 2) g=z(1)-1.D0
!
      RETURN
   END SUBROUTINE GSUB
!!!
!!!
   SUBROUTINE DGSUB(I,NCOMP,Z,DG,RPAR,IPAR)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: i, ncomp
      REAL(r8),    INTENT(IN)    :: z(*)
      REAL(r8),    INTENT(INOUT) :: dg(*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
!     I == the boundary condition "number".
!     The conditions at the left are enumerated first,
!     then the ones at the right.
!     The order of left or right conditions does not matter,
!     but we must be consistent when defining the jacobian
!     of the boundary conditions!
!
!!!
!
      IF (I == 1) THEN
!     dG1/dZ1      
         dg(1) = 1.D0
!     dG1/dZ2
         dg(2) = 0.D0
      END IF
!            
      IF (I == 2) THEN
!     dG2/dZ1
         dg(1)=1.D0
!     dG2/dZ2         
         dg(2)=0.D0
      END IF
!
      RETURN
   END SUBROUTINE DGSUB
!
END

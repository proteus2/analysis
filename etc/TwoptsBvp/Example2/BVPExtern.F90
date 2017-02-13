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
!!!
!
      f(1)=COS(z(3))
      f(2)=SIN(z(3))
      f(3)=z(4)
      f(4)=z(5)*COS(z(3))
      f(5)=0._r8
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
!!!
!
!     dF1/dZ3
      df(1,3)=-SIN(z(3))
!     dF2/dZ3      
      df(2,3)=COS(z(3))
!     dF3/dZ4
      df(3,4)=1.0D0
!     dF4/dZ3
      df(4,3)=-z(5)*SIN(z(3))
!     dF4/dZ4
      df(4,4)=1.0D0
!     dF4/dZ5
      df(4,5)=COS(z(3))
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
!     We have specified 3 left bcs, and 5 bcs total.
!     This means that:
!     BC(1) = x(0) = 0
!     BC(2) = y(0) = 0
!     BC(3) = kappa(0) = 0
!     BC(4) = y(0.5) = 0
!     BC(5) = phi(0.5) = -pi/2
!
!!!
!
      IF (i == 1) g = z(1)
      IF (i == 2) g = z(2)
      IF (i == 3) g = z(4)
      IF (i == 4) g = z(2)
      IF (i == 5) g = z(3)+1.5707963267948966192313216916397514D0
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
!     We have specified 3 left bcs, and 5 bcs total.
!     This means that:
!     BC(1) = x(0) = 0
!     BC(2) = y(0) = 0
!     BC(3) = kappa(0) = 0
!     BC(4) = y(0.5) = 0
!     BC(5) = phi(0.5) = -pi/2      
!
      dg(1) = 0.0_r8
      dg(2) = 0.0_r8
      dg(3) = 0.0_r8
      dg(4) = 0.0_r8
      dg(5) = 0.0_r8
!
!     dG1/dZ1
      IF (i == 1) dg(1)=1.0_r8
!     dG2/dZ2
      IF (i == 2) dg(2)=1.0_r8
!     dG3/dZ4
      IF (i == 3) dg(4)=1.0_r8
!     dG4/dZ2
      IF (i == 4) dg(2)=1.0_r8
!     dG5/dZ3
      IF (i == 5) dg(3)=1.0_r8
!
      RETURN
   END SUBROUTINE DGSUB
!
END

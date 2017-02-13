MODULE LMExtern
!
   USE Base, ONLY: i4, r4, r8
!
   IMPLICIT NONE
!
   PRIVATE
!
   INTERFACE FCN
      MODULE PROCEDURE FCNJAC
      MODULE PROCEDURE FCNNOJAC
   END INTERFACE
!
   PUBLIC :: FCN
!
CONTAINS
!
   SUBROUTINE FCNJAC(m, n, x, fvec, fjac, iflag)
!
      USE Base, ONLY: i4, r4, r8
!
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN)    :: m, n
      REAL(r8),    INTENT(IN)    :: x(:)
      REAL(r8),    INTENT(INOUT) :: fvec(:)
      REAL(r8),    INTENT(OUT)   :: fjac(:,:)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
!     Local variables
!
      INTEGER(i4) :: i
!!!
!
      RETURN
   END SUBROUTINE FCNJAC
!
   SUBROUTINE FCNNOJAC(m, n, x, fvec, iflag)
!
      USE Base, ONLY: i4, r4, r8
!
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN)    :: m, n
      REAL(r8),    INTENT(IN)    :: x(:)
      REAL(r8),    INTENT(INOUT) :: fvec(:)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
!     Local variables
!
      INTEGER(i4) :: i
!!!
!
      RETURN
   END SUBROUTINE FCNNOJAC
!
END

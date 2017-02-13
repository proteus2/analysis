MODULE LMExtern
!
   USE Base,     ONLY: i4, r4, r8
   USE DataStorage
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
   SUBROUTINE FCNJAC(m, n, p, fvec, fjac, iflag)
!
!     Calculate either residuals or the Jacobian matrix.
!     A = p(1), B = p(2), C = p(3), D = p(4)
!     m = no. of cases, n = no. of parameters (4)
!
      USE Base,        ONLY: i4, r4, r8
      USE DataStorage
!
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN)    :: m, n
      REAL(r8),    INTENT(IN)    :: p(:)
      REAL(r8),    INTENT(INOUT) :: fvec(:)
      REAL(r8),    INTENT(OUT)  :: fjac(:,:)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
!     Local variables
!
      REAL(r8), PARAMETER :: one = 1.0_r8
      REAL(r8)            :: expntl, temp
      INTEGER(i4)         :: i
!!!
!
      IF (iflag == 1) THEN
         fvec = y(1:m) - p(1) - p(3)/(one + EXP(-p(2)*(x(1:m) - p(4))))
      ELSE IF (iflag == 2) THEN
         fjac(1:m,1) = -one
         DO i = 1, m
            expntl = EXP(-p(2)*(x(i) - p(4)))
            temp = one / (one + expntl)
            fjac(i,2) = - p(3) * temp**2 * (x(i) - p(4)) * expntl
            fjac(i,3) = - temp
            fjac(i,4) =   p(3) * temp**2 * expntl * p(2)
         END DO
      END IF
!
      RETURN
   END SUBROUTINE FCNJAC
!
   SUBROUTINE FCNNOJAC(m, n, p, fvec, iflag)
!
      USE Base,        ONLY: i4, r4, r8
!
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN)    :: m, n
      REAL(r8),    INTENT(IN)    :: p(:)
      REAL(r8),    INTENT(INOUT) :: fvec(:)
      INTEGER(i4), INTENT(INOUT) :: iflag
!
!     Local variables
!
      INTEGER(i4)         :: i
!!!
!
      RETURN
   END SUBROUTINE FCNNOJAC
!
END

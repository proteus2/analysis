!-------------------------------------------------------------------------------
!
!  MODULE Base
!
!> @brief
!> - Declare type parameters
!
!  REVISION HISTORY
!
!  16AUG2015 : In-Sun Song : Reviewed
!
!-------------------------------------------------------------------------------
!
MODULE Base
!
   IMPLICIT NONE
!
   PRIVATE
!
   INTEGER, PARAMETER, PUBLIC :: r8 = SELECTED_REAL_KIND(12,60)
   INTEGER, PARAMETER, PUBLIC :: r4 = SELECTED_REAL_KIND( 6)
   INTEGER, PARAMETER, PUBLIC :: rn = KIND(1.0)   ! native real
   INTEGER, PARAMETER, PUBLIC :: i8 = SELECTED_INT_KIND(13)
   INTEGER, PARAMETER, PUBLIC :: i4 = SELECTED_INT_KIND( 6)
   INTEGER, PARAMETER, PUBLIC :: in = KIND(1)     ! native integer
!
!  INTEGER, PARAMETER, PUBLIC :: r16 = 8
   INTEGER, PARAMETER, PUBLIC :: r16 = SELECTED_REAL_KIND(12)
!
END MODULE Base

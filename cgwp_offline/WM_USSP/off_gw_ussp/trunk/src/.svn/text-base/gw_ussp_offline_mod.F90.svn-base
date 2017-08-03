! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENSE.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Global standard conversions

MODULE conversions_mod

! Description:
!       Model and code section invariant physical constants

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

IMPLICIT NONE

! 
REAL, PARAMETER :: pi                  = 3.14159265358979323846
!
! Conversion factor degrees to radians
REAL, PARAMETER :: pi_over_180         = pi/180.0
!
! zerodegc is a conversion between degrees centigrade and kelvin
REAL, PARAMETER :: zerodegc            = 273.15

END MODULE conversions_mod
! ------------------------
! --------------
MODULE  parkind1

!
! Description:
!   Dummy module to replace the DrHook library. Defines data types
!   which would otherwise be declared by the DrHook library.
!
! Code Owner: Please refer to the utils/off_gw_ussp file README.txt
! This file belongs in section: Dummy libraries
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.


IMPLICIT NONE

INTEGER, PARAMETER :: jpim = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: jprb = SELECTED_REAL_KIND(13,300)

END MODULE parkind1
! -----------------
! ------------
MODULE yomhook

!
! Description:
!   Dummy module to replace the DrHook library.
!
! Code Owner: Please refer to the utils/off_gw_ussp file README.txt
! This file belongs in section: Dummy libraries
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.

USE parkind1, ONLY: jpim, jprb
IMPLICIT NONE

LOGICAL, PARAMETER :: lhook = .FALSE.

INTERFACE dr_hook
MODULE PROCEDURE dr_hook
END INTERFACE dr_hook

CONTAINS

SUBROUTINE dr_hook(NAME,code,handle)
IMPLICIT NONE

!Arguments
CHARACTER(LEN=*),   INTENT(IN)    :: NAME
INTEGER(KIND=jpim), INTENT(IN)    :: code
REAL(KIND=jprb),    INTENT(INOUT) :: handle

!Nothing to do

RETURN

END SUBROUTINE dr_hook
END MODULE yomhook
! ----------------

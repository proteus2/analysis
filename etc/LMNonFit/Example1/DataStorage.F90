MODULE DataStorage

! The data is stored in this module.
! This is equivalent to a COMMON area in old Fortran.

   USE Base, ONLY: i4, r4, r8
!
   IMPLICIT NONE
!
   REAL(r8), SAVE :: x(100), y(100)
!
END MODULE DataStorage
!

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENSE.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! 
! Subroutine Interface:

MODULE runvnamelist_off_mod

IMPLICIT NONE
!
! purpose:    This module provides runtime variables required for offline
!             running of the core subroutine that calculates the 
!             vertically propagated flux of pseudomomentum due to gravity
!             waves as parametrized by the Ultra Simple Spectral gravity 
!             wave parametrization
!
! Code Owner: Please refer to the utils/off_gw_ussp file README.txt
! This file belongs in section: Gravity Wave Drag
!
! code description:
!   language: fortran 90
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

! Specify field dimension variable inputs
!
! Hardwire maximum value of grid points in row
INTEGER, PARAMETER  :: max_row_length = 192
!
! Hardwire maximum value of rows
INTEGER, PARAMETER  :: max_rows = 145
!
! Hardwire maximum value of levels to permit namelist definition
INTEGER, PARAMETER  :: max_model_lev = 200
!
! Hardwire maximum number of directions
INTEGER, PARAMETER  :: max_idir = 16

! Hardwire Fortran unit for input namelist file
INTEGER, PARAMETER :: In_UNIT = 16
!
! Number of grid points in row
INTEGER  :: row_length
!
! Number of rows
INTEGER  :: rows
!
! Number of model levels
INTEGER  :: model_levels
!
! Max number of directions, typically four.
INTEGER  :: idir
!
! Level for gw launch - in principle this could be launchlev(i,j)
INTEGER  :: launchlev

! Desired Fortran unit for output file
INTEGER :: IUNIT
!
! Cos(phi_j) - azimuthal direction specify consistent with idir
REAL            :: cosazim(max_idir)
!
! Sin(phi_j) - azimuthal direction specify consistent with idir
REAL            :: sinazim(max_idir)
!
! Specify latitude values for calculating Sin_theta_lat
! Interval e.g. 180. / (rows -1)
REAL    ::  delta_phi
!
! Zero point for j = 1, N angle is base_phi + (j*delta_phi)
REAL    ::  base_phi
!
! Factor enhancement for invariant global wave launch amplitude
REAL    :: ussp_launch_factor
!
! Characteristic (spectrum peak) wavelength (m)
REAL    :: wavelstar
!
! Factor enhancement for conversion from convective rain to GW source flux
REAL    :: cgw_scale_factor
!
! Flag to indicate CGW flux calculation should be carried out
LOGICAL :: L_add_cgw


NAMELIST /RUNVARS/                                                      &
      row_length, rows, model_levels, idir, launchlev,                  &
      ussp_launch_factor, wavelstar, L_add_cgw, cgw_scale_factor,       &
      base_phi, delta_phi, cosazim, sinazim, IUNIT


CONTAINS

SUBROUTINE set_runvnamelist_off()

IMPLICIT NONE

! ---------------------------------------------------------------------+--------
! Evaluate default value for parameters (SETCONA)
! ---------------------------------------------------------------------+--------
delta_phi = 180. / (rows - 1.)
base_phi  = (-90. - delta_phi)


END SUBROUTINE set_runvnamelist_off
! ---------------------------------
! ---------------------------------------
SUBROUTINE read_runvnamelist_off(unit_in)

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: ErrorStatus

! Set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 6
INTEGER, PARAMETER :: n_real = 5
INTEGER, PARAMETER :: n_log = 1

READ(UNIT=unit_in, NML=RUNVARS, IOSTAT=ErrorStatus)

IF (ErrorStatus /= 0) THEN
  WRITE (0, *) 'Error ', ErrorStatus, ' occurred on opening input file'&
                       , unit_in
END IF

END SUBROUTINE read_runvnamelist_off
! ----------------------------------

END MODULE runvnamelist_off_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENSE.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Wrapper program for off-line running of the Non-orographic Gravity
!  Wave (Ultra-Simple Spectral Parametrization) Scheme Core.
PROGRAM GW_ussp_core_offline
!
! purpose: This program enables offline running of the core subroutine 
!          that calculates the vertically propagated flux of
!          pseudomomentum due to gravity waves as parametrized by the
!          Ultra Simple Spectral gravity wave Parametrization
!         (originally Warner and McIntyre, since reengineered for use
!          in the UM and extended).
!
! Code Owner: Please refer to the utils/off_gw_ussp file README.txt
! This file belongs in section: Gravity Wave Drag
!
! code description:
!   language: fortran 90
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

USE conversions_mod, ONLY: pi, pi_over_180
USE runvnamelist_off_mod, ONLY: max_row_length, max_rows,              &
          max_model_lev, max_idir, In_UNIT, row_length, rows,          &
          model_levels, idir, launchlev, ussp_launch_factor, wavelstar,&
          L_add_cgw, cgw_scale_factor, base_phi, delta_phi,            &
          cosazim, sinazim, IUNIT, RUNVARS,                            &
          read_runvnamelist_off, set_runvnamelist_off
! ---------------------------------------------------------------------+--------
! Subroutines defined from modules
! ----------------------------------------------------------------------+-------
USE gw_ussp_core_mod, ONLY: gw_ussp_core
!     
! ----------------------------------------------------------------------+-------

IMPLICIT NONE
! ---------------------------------------------------------------------+--------
! Assorted parameters.
! ---------------------------------------------------------------------+--------
! Unit number on which to read in data arrays
INTEGER, PARAMETER :: ifunit = 17
!
! omega is Angular speed of Earth's rotation (set in SETCONA)
! For rotating earth value set to 2*pi/siderial day (23h56m04s)
REAL, PARAMETER :: omega     = 7.292116E-5 
REAL, PARAMETER :: two_omega = 2 * omega
!
! Strength coefficient constant for Launch spectrum (CCL / A0)
REAL, PARAMETER :: ccl0              = 3.41910625e-9
!
! Scale conversion factor from convective rain to GW flux (mPa)
REAL, PARAMETER :: cor2flux          = 7.20507e-04
!
! Minimum allowed non-zero value of precipitation (0.1 mm / day)
REAL, PARAMETER ::  rppnmin           = 1.0 / 1.1574e-06
!
! Translate between latitude (-90. to +90.) and (0. to 180.)
REAL, PARAMETER ::  offset_phi        = 90.0
!
! Switch to set fluxes to zero at Upper Boundary
LOGICAL, PARAMETER :: L_ussp_opaque  = .TRUE.
! ----------------------------------------------------------------------+-------
!     Local variables (dynamic arrays) used in GW_USSP
! ----------------------------------------------------------------------+-------
!
! Grid point latitudes
REAL    :: sin_theta_latitude(max_row_length, max_rows)
!
! Rho (kg m^-3)
REAL    :: rho(max_row_length, max_rows, max_model_lev)
!
! U wind (m s^-1)
REAL    :: u(max_row_length, max_rows, max_model_lev)
!
! V wind (m s^-1)
REAL    :: v(max_row_length, max_rows, max_model_lev)
!
! Buoyancy [Brunt Vaisala] Frequency (rad s^-1)
REAL    :: nbv(max_row_length, max_rows, max_model_lev)
!
! Total Precipitation for CGW source
REAL    :: total_ppn(max_row_length, max_rows)
!
! Heights of input levels (m)
REAL    :: zlev(max_row_length, max_rows, max_model_lev)
!
! Component of wind in phi_jdir direction (m s^-1)
REAL    :: udotk(max_row_length, max_rows, max_model_lev, max_idir)
!
! Pseudomomentum flux integrated over azimuthal sector (kg m^-1 s^-2)
REAL    :: fptot(max_row_length, max_rows, max_model_lev, max_idir)

! ---------------------------------------------------------------------+--------
! Assorted scalars.
! ---------------------------------------------------------------------+--------
! Return code from I/O operations
INTEGER :: IOCode
!
! Longitude index
INTEGER         :: i
!
! Latitude index
INTEGER         :: j
!
! Level index
INTEGER         :: k
!
! Azimuthal direction index
INTEGER :: jdir
!
! Loop terminators (for compatibility with online code)
INTEGER         :: i_start, i_end, j_start, j_end, k_start, k_end
!
! Array dimensions from input file
INTEGER :: idim(3)
!
! ==Main Block==-------------------------------------------------------+--------
!
! Set default values for namelist RUNVARS
! ---------------------------------------
! row_length = 192
! rows       = 144
! model_levels   = 85
! idir       = 4
! launchlev  = 19
! ussp_launch_factor = 1.0
! wavelstar  = 4300.0
! L_add_cgw  = .FALSE.
! cgw_scale_factor = 1.0
!
CALL read_runvnamelist_off(In_UNIT)

WRITE(6,RUNVARS)
!
! --------------------------------------------------------
! Calculation from atm_fields_bounds_mod (ENDGAME version)
! --------------------------------------------------------
!     dimensions of theta-field:
i_start           = 1
i_end             = row_length
j_start           = 1
j_end             = rows
k_start           = 1
k_end             = model_levels

! ---------------------------------------------------------------------+--------
! Evaluated field (SETCONA) Sin_theta_latitude
! ---------------------------------------------------------------------+--------
DO j = j_start, j_end
  DO i = i_start, i_end
    sin_theta_latitude(i, j) = SIN(pi_over_180 *                       &
                                   (offset_phi+base_phi+(j*delta_phi)))
    total_ppn(i, j) = 0.0
  END DO
END DO
!
! ---------------------------------------------------------------------+--------
! Read input fields required to run gw_ussp_core offline
! L_add_cgw = T : calculate variable CGW launch flux
! ---------------------------------------------------------------------+--------
OPEN(ifunit, STATUS='old', IOSTAT=IOCode)
IF (IOCode /= 0) THEN
  WRITE (0, *) 'Error ', IOCode, ' occurred on opening input file',    &
                         ifunit
  GOTO 999
END IF

READ(ifunit, *) idim
READ(ifunit, *) u(i_start:i_end,j_start:j_end,k_start:k_end)
READ(ifunit, *) v(i_start:i_end,j_start:j_end,k_start:k_end)
READ(ifunit, *) rho(i_start:i_end,j_start:j_end,k_start:k_end)
READ(ifunit, *) nbv(i_start:i_end,j_start:j_end,k_start:k_end)
READ(ifunit, *) zlev(i_start:i_end,j_start:j_end,k_start:k_end)
IF (L_add_cgw) THEN
  READ(ifunit, *) total_ppn(i_start:i_end,j_start:j_end)
END IF
CLOSE(ifunit)
!
DO k = k_start, k_end
  DO j = j_start, j_end
    DO i = i_start, i_end
!     rho(i, j, k) = 0.0
!     nbv(i, j, k) = 0.0
      DO jdir = 1, idir
        udotk(i,j,k,jdir) = (u(i, j, k) * cosazim(jdir)) +             &
                            (v(i, j, k) * sinazim(jdir))
      END DO
    END DO
  END DO
END DO
!
! ---------------------------------------------------------------------+--------
! Run core calculation of the Ultra-Simple Spectral Parametrization for 
! non-orographic gravity waves
! L_add_cgw = F : calculate standard USSP isotropic GW launch flux
! L_add_cgw = T : calculate variable CGW launch flux
! ---------------------------------------------------------------------+--------
CALL gw_ussp_core(rows, row_length, idir, launchlev, i_start, i_end,   &
     j_start, j_end, k_start, k_end, ussp_launch_factor, wavelstar,    &
     cgw_scale_factor, two_omega,                                      &
     sin_theta_latitude(i_start:i_end,j_start:j_end),                  &
     rho(i_start:i_end,j_start:j_end,k_start:k_end),                   &
     nbv(i_start:i_end,j_start:j_end,k_start:k_end),                   &
     udotk(i_start:i_end,j_start:j_end,k_start:k_end,1:idir),          &
     total_ppn(i_start:i_end,j_start:j_end),                           &
     fptot(i_start:i_end,j_start:j_end,k_start:k_end,1:idir),          &
     L_ussp_opaque, L_add_cgw)
!
! ----------------------------------------------------------------------+-------
! Sending output to Unit=IUNIT
! ----------------------------------------------------------------------+-------
WRITE (6, *) 'Now sending output.'
WRITE (6, *) fptot(i_start:i_end,j_start:j_end,k_start:k_end,1:idir)
!
OPEN(IUNIT, FORM='unformatted', STATUS='new')
WRITE (IUNIT) fptot(i_start:i_end,j_start:j_end,k_start:k_end,1:idir)
! ----------------------------------------------------------------------+-------
!

 999  CONTINUE  ! error exit GO TOs come here

CLOSE (IUNIT)

STOP('Offline code completed.')

END PROGRAM GW_ussp_core_offline

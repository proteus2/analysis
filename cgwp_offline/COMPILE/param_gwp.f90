MODULE param_gwp

!     nphi    - Number of wave-propagation directions considered
!
!     nc      - Number of positive phase speeds in the discrete spectrum
!               Total number is nc*2+1 (including zero and negatives).
!
!     c_max   - Maximum phase speed in the spectrum
!
!     cfactor - Conversion factor for magnitude of the cloud-top
!               momentum flux
!               (See the appendix in Song et al. (2007, JAS).)
!
!     phi_deg - Wave-propagation directions considered
!               They must be in [0,180).
!               45, 135 (deg) are chosen by Choi and Chun (2011, JAS).

  implicit none

  integer, public ::  nc, nphi
  real   , public ::  dc

  real, dimension(:), allocatable, public ::  phi_deg, c_phase, phi_deg2
  real, dimension(:), allocatable, public ::  cosphi, sinphi
  real, dimension(:), allocatable, public ::  c_m05dc

  real, private, parameter :: pi = 3.14159265358979323846

!:::  CGW source scheme (Song and Chun, 2005, JAS) :::::::::::::::::::::

  ! See Song et al., 2007, JAS for the values.
  real, parameter, public ::  hscale = 5.e3
  real, parameter, public ::  tscale = 1200.
  real, parameter, public ::  lt     = tscale
  real, parameter, public ::  ah     = pi*1.e4*tscale*tscale

  ! Criterion of the maximum subgrid-scale heating rate [K/s]
  real, parameter, public ::  schm_c = 0.0001/86400.

  real, public ::  cfactor

!:::  Warner and McIntyre scheme  ::::::::::::::::::::::::::::::::::::::

  ! See UMDP34.
  ! The code is written assuming (s,t) = (1,3).
  real, public ::  mstar_wm, p_wm, beta_wm
  real, public ::  s_wm, t_wm   ! not used

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CONTAINS


SUBROUTINE set_spec_param

  integer ::  ic

  allocate( c_phase(-nc:nc) )
  do ic=-nc, nc
    c_phase(ic) = float(ic)*dc
  enddo

  allocate( c_m05dc(-nc:nc+1) )
  c_m05dc(-nc:nc) = c_phase(:) - 0.5*dc
  c_m05dc(nc+1) = c_phase(nc) + 0.5*dc

  allocate( cosphi(nphi), sinphi(nphi) )
  cosphi(:) = cos(phi_deg(:)*(pi/180.))
  sinphi(:) = sin(phi_deg(:)*(pi/180.))

  allocate( phi_deg2(nphi*2) )
  phi_deg2(:nphi  ) = phi_deg(:)
  phi_deg2(nphi+1:) = phi_deg(:) + 180.

END subroutine set_spec_param

SUBROUTINE get_wm_default

  mstar_wm = 2.0*pi/4.3e3
  p_wm = 5.0/3.0
  s_wm = 1.0
  t_wm = 3.0
  beta_wm = 2.0   ! should be checked !

END subroutine get_wm_default

SUBROUTINE get_wm_hg2cgwp

  call get_wm_default

  beta_wm = 0.8

END subroutine get_wm_hg2cgwp

END module param_gwp


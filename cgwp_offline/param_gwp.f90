MODULE param_gwp

  implicit none

  ! parameters for Warner and McIntyre scheme (see UMDP34.)
  ! The code is written assuming (s,t) = (1,3).
  real, public ::  mstar_wm, p_wm, beta_wm
  real, public ::  s_wm, t_wm   ! not used

  real, private, parameter :: pi = 3.14159265358979323846

CONTAINS


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


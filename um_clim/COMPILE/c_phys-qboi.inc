  real, parameter :: g = 9.806
  real, parameter :: rd = 287.0
  real, parameter :: cp = 1004.0
  real, parameter :: kappa = rd/cp
  real, parameter :: p0 = 1000.e2
!  real, parameter :: h_scale = 6.6e3   ! tropical strato.
  real, parameter :: h_scale = 7.0e3   ! extratropics; QBOi
  real, parameter :: rho_s = p0/g/h_scale
  real, parameter :: r_earth = 6371.e3
  real, parameter :: lapse_rate_sfc = 6.5e-3
  real, parameter :: omega_earth = 7.292116e-5
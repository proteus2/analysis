MODULE subr_common

  implicit none

  real, private, parameter :: pi = 3.14159265358979323846
  real, private, parameter :: deg2rad = pi/180.0

CONTAINS


SUBROUTINE basic_u_phi(u,v,phi_deg, u_phi)

  real, dimension(:,:), intent(in) ::  u, v
  real, dimension(:)  , intent(in) ::  phi_deg

  real, dimension(:,:,:), intent(out) ::  u_phi

  real, dimension(size(phi_deg)) ::  phi_rad
  integer ::  nphi, iphi

  nphi = size(phi_deg)

  phi_rad(:) = phi_deg(:)*deg2rad

  do iphi=1, nphi
    u_phi(:,:,iphi) = u(:,:)*cos(phi_rad(iphi)) +                        &
                      v(:,:)*sin(phi_rad(iphi))
  enddo

END subroutine basic_u_phi

END module subr_common


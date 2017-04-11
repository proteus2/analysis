MODULE subr_common

  implicit none

  real, private, parameter :: pi = 3.14159265358979323846
  real, private, parameter :: deg2rad = pi/180.0

  contains


SUBROUTINE basic_u_phi(u,v,phi_deg,u_phi)

  real, dimension(:,:), intent(in) ::  u, v
  real, dimension(:)  , intent(in) ::  phi_deg

  real, dimension(:,:,:), intent(out) ::  u_phi

  real, dimension(size(phi_deg)) ::  cosphi, sinphi
  integer ::  nphi, iphi

  nphi = size(phi_deg)

  cosphi(:) = cos(phi_deg(:)*deg2rad)
  sinphi(:) = sin(phi_deg(:)*deg2rad)

  do iphi=1, nphi
    u_phi(:,:,iphi) = u(:,:)*cosphi(iphi) +                &
                      v(:,:)*sinphi(iphi)
  enddo

END subroutine basic_u_phi

SUBROUTINE ind_c_critlayer(c_phase,ic1,u_phi, ic_cl)

  real, dimension(:),   intent(in) ::  c_phase
  integer,              intent(in) ::  ic1
  real, dimension(:,:), intent(in) ::  u_phi

!  real, dimension(size(u_phi)), intent(out) ::  ic_cl
  real, dimension(size(u_phi,dim=1),size(u_phi,dim=2)), intent(out) ::  ic_cl

  integer ::  ncol, nphi, l, iphi
  real, dimension(size(c_phase)) ::  c_int

  ncol = size(u_phi, dim=1)
  nphi = size(u_phi, dim=2)

  do iphi=1, nphi
  do l=1, ncol
    c_int(:) = c_phase(:) - u_phi(l,iphi)
    ic_cl(l,iphi) = minloc(abs(c_int),1)+(ic1-1)
  enddo
  enddo

END subroutine ind_c_critlayer

END module subr_common


MODULE T_G_eqn
! w_zz + m2(z)*w = 0, where W = Re{ w*exp[i(kx-ft)] }/sqrt[rho(z)]

  implicit none

  real,    parameter ::  vr_null = -999.
  complex, parameter ::  vc_null = (-999.,-999.)

  integer ::  nz_tg, m2n_bvp
  real    ::  eps_tg, m_scl

  real, dimension(:), allocatable ::  z_nd, m2_nd, m2_zz_nd

  contains

SUBROUTINE t_g_eqn_ivp(nz,m2,zt,z1,w1,wp1, w,wz,kf)

  use nr, only: spline, odeint, stifbs

  implicit none

  integer,             intent(in) ::  nz
  real,                intent(in) ::  zt, z1
  complex,             intent(in) ::  w1, wp1
  real, dimension(nz), intent(in) ::  m2

  complex, dimension(nz), intent(out) ::  w, wz
  real,    dimension(nz), intent(out) ::  kf

  integer ::  k, k1
  real    ::  dz_nd, tmp, yr(2), yi(2), tmp1d(nz)
  complex ::  tmpc1d(nz*2)

  integer ::  ijob, ier
  complex ::  a(0:nz+1,0:nz+1), b(0:nz+1), warr((nz+2)*(nz+4))
  real    ::  wkarr(0:nz+1)

  interface
    SUBROUTINE derivs(x,y,dydx)
    USE nrtype
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: y
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
    END SUBROUTINE derivs
  end interface
!
! setup
!
  if ( eps_tg <= 0. .or. eps_tg > 0.1 ) then
    print*, 'eps_tg should be defined in the main program.'  ;  STOP
  end if

  if ( .not. allocated(z_nd) ) then
    nz_tg = nz
    m_scl = sqrt(sum(abs(m2(:)))/float(nz))
    dz_nd = zt/float(nz-1)*m_scl
    allocate( z_nd(nz), m2_nd(nz), m2_zz_nd(nz) )
    do k=1, nz
      z_nd(k) = float(k-1)*dz_nd
    enddo
    m2_nd(:) = m2(:)/(m_scl*m_scl)
    call spline(z_nd,m2_nd,0.,0.,m2_zz_nd)
  end if
!
! initial condition
!
  k1 = min(nz, nint(z1/(zt/float(nz-1)))+1)
  w(k1) = w1
  wz(k1) = wp1/m_scl
!
! solver
!
  if ( minval(m2) < 0. .and. m2n_bvp == 1 ) then
    call t_g_eqn_bvp(3,nz,m2,zt,z1,w1,z1,wp1, w,wz,kf)
    deallocate( z_nd, m2_nd, m2_zz_nd )
    RETURN
  end if

  do k=k1+1, nz
    yr = real ( (/ w(k-1), wz(k-1) /) )
    yi = aimag( (/ w(k-1), wz(k-1) /) )
    call odeint(yr,z_nd(k-1),z_nd(k),eps_tg,dz_nd,0.,derivs,stifbs)
    call odeint(yi,z_nd(k-1),z_nd(k),eps_tg,dz_nd,0.,derivs,stifbs)
    w(k) = cmplx(yr(1),yi(1))  ;  wz(k) = cmplx(yr(2),yi(2))
  enddo
  do k=k1-1, 1, -1
    yr = real ( (/ w(k+1), wz(k+1) /) )
    yi = aimag( (/ w(k+1), wz(k+1) /) )
    call odeint(yr,z_nd(k+1),z_nd(k),eps_tg,dz_nd,0.,derivs,stifbs)
    call odeint(yi,z_nd(k+1),z_nd(k),eps_tg,dz_nd,0.,derivs,stifbs)
    w(k) = cmplx(yr(1),yi(1))  ;  wz(k) = cmplx(yr(2),yi(2))
  enddo
  wz(:) = wz(:)*m_scl
!
! momentum flux
!
  kf(:) = -0.5*aimag(wz(:)*conjg(w(:)))

  deallocate( z_nd, m2_nd, m2_zz_nd )

  RETURN

END subroutine t_g_eqn_ivp


SUBROUTINE t_g_eqn_bvp(bdy,nz,m2,zt,z1,w1,z2,w2, w,wz,kf)

  implicit none

  integer,             intent(in) ::  bdy
  ! 1, Dirichlet BDY at two specified pts        :  w(z1), w(z2)
  ! 2, Dirichlet BDY at k1 and Neumann BDY at k2 :  w(z1), w'(z2)
  ! 3, Cauchy BDY at k1 (initial condition)      :  w(z1), w'(z1)
  integer,             intent(in) ::  nz
  real,                intent(in) ::  zt, z1, z2
  complex,             intent(in) ::  w1, w2
  real, dimension(nz), intent(in) ::  m2

  complex, dimension(nz), intent(out) ::  w, wz
  real,    dimension(nz), intent(out) ::  kf

  integer ::  k, k1, k2
  real    ::  dz, tmp

  integer ::  ijob, ier
  complex ::  a(0:nz+1,0:nz+1), b(0:nz+1), warr((nz+2)*(nz+4))
  real    ::  wkarr(0:nz+1)
!
! constants
!
  dz = zt/float(nz-1)
!
! matrix for linear eqn set (2nd-order differencing)
!
  a(:,:) = 0.
  tmp = dz*dz
  do k=1, nz
    a(k,k) = m2(k)*tmp - 2.
    a(k,k-1) = 1.
    a(k,k+1) = 1.
  enddo
  b(:) = 0.
!
! boundary condition
!
  k1 = min(nz, nint(z1/dz)+1)
  k2 = min(nz, nint(z2/dz)+1)
  select case ( bdy )
  case ( 1 )
    a(0  ,k1) = 1.  ;  b(0  ) = w1
    a(nz+1,k2) = 1.  ;  b(nz+1) = w2
  case ( 2 )
    a(0,k1) = 1.  ;  b(0) = w1
    a(nz+1,k2-1) = -1.  ;  a(nz+1,k2+1) = 1.  ;  b(nz+1) = w2*(2.*dz)
  case ( 3 )
    a(0,k1) = 1.  ;  b(0) = w1
    a(nz+1,k1-1) = -1.  ;  a(nz+1,k1+1) = 1.  ;  b(nz+1) = w2*(2.*dz)
  end select
!
! solver
!
  ijob = 0
  call leq2c(a,nz+2,nz+2,b,1,nz+2,ijob,warr,wkarr,ier)
!
! w, momentum flux
!
  w(:) = b(1:nz)
  do k=1, nz
    wz(k) = (b(k+1)-b(k-1))/(2.*dz)
  enddo
  kf(:) = -0.5*aimag(wz(:)*conjg(w(:)))

  RETURN

END subroutine t_g_eqn_bvp


SUBROUTINE wave_decomp(nz,dz,m2,w,wz, w_pm,w_nm,kf_pm,kf_nm,r)
! Assumption: 1) The WKB assumption is used (i.e., m' = 0)
!             2) m2 > 0, everywhere
! This code is made with the assumption that m2 > 0 everywhere.

  implicit none

  integer,                intent(in) ::  nz
  real,                   intent(in) ::  dz
  real,    dimension(nz), intent(in) ::  m2
  complex, dimension(nz), intent(in) ::  w
  complex, dimension(:),  intent(in) ::  wz

  complex, dimension(nz), intent(out) ::  w_pm, w_nm
  real,    dimension(nz), intent(out) ::  kf_pm, kf_nm, r

  integer ::  k
  complex ::  der_z(nz)

  complex, parameter ::  ai = (0.,1.)

  if (size(wz) == nz) then
    der_z(:) = wz(:)
  else if (size(wz) == 2) then
    do k=2, nz-1
      der_z(k) = (w(k+1)-w(k-1))/(2.*dz)
    enddo
    der_z(1 ) = wz(1)
    der_z(nz) = wz(2)
  else
    print*, 'Error: size of wz in wave_decomp'  ;  STOP
  end if

  r(:) = 0.
  w_pm(:) = 0.  ;  w_nm(:) = 0.
  where (m2(:) > 0.)
    w_pm(:) = 0.5*( w(:) - ai/sqrt(m2(:))*der_z(:) )
    w_nm(:) = w(:) - w_pm(:)
    r(:) = abs(w_pm(:)/w_nm(:))
  end where

  kf_pm(:) = 0.  ;  kf_nm(:) = 0.
  do k=2, nz-1
    der_z(k) = (w_pm(k+1)-w_pm(k-1))/(2.*dz)
    kf_pm(k) = -0.5*aimag(der_z(k)*conjg(w_pm(k)))
    der_z(k) = (w_nm(k+1)-w_nm(k-1))/(2.*dz)
    kf_nm(k) = -0.5*aimag(der_z(k)*conjg(w_nm(k)))
  enddo
  where (m2(2:nz-1) <= 0.)
    kf_pm(1:nz-2) = 0.
    kf_pm(3:nz  ) = 0.
    kf_nm(1:nz-2) = 0.
    kf_nm(3:nz  ) = 0.
  end where
  kf_pm(1 ) = kf_pm(2   )
  kf_pm(nz) = kf_pm(nz-1)
  kf_nm(1 ) = kf_nm(2   )
  kf_nm(nz) = kf_nm(nz-1)

  RETURN

END subroutine wave_decomp


END module T_G_eqn


SUBROUTINE jacobn(x,y,dfdx,dfdy)

  use nrtype
  use T_G_eqn, only: nz_tg, z_nd, m2_nd, m2_zz_nd
  use nr, only: splint, locate

  implicit none

  real(sp),                 intent(in)  :: x
  real(sp), dimension(:),   intent(in)  :: y
  real(sp), dimension(:),   intent(out) :: dfdx
  real(sp), dimension(:,:), intent(out) :: dfdy

  integer ::  klo, khi
  real ::  mm, mm_z, h, a, b

  klo = max(min(locate(z_nd,x),nz_tg-1),1)
  khi = klo+1
  h   = z_nd(khi)-z_nd(klo)
  a   = (z_nd(khi)-x)/h
  b   = (x-z_nd(klo))/h

  mm   = splint(z_nd,m2_nd,m2_zz_nd,x)
  mm_z = (m2_nd(khi)-m2_nd(klo))/h - h/6.* &
         ( (3.*(a*a)-1.)*m2_zz_nd(klo) -   &
           (3.*(b*b)-1.)*m2_zz_nd(khi) )

  dfdx(:) = 0.  ;  dfdy(:,:) = 0.
  dfdx(2)   = -mm_z*y(1)
  dfdy(1,2) = 1.
  dfdy(2,1) = -mm

END subroutine jacobn

SUBROUTINE derivs(x,y,dydx)

  use nrtype
  use T_G_eqn, only: z_nd, m2_nd, m2_zz_nd
  use nr, only: splint

  implicit none

  real(sp),               intent(in)  :: x
  real(sp), dimension(:), intent(in)  :: y
  real(sp), dimension(:), intent(out) :: dydx

  real ::  mm

  mm = splint(z_nd,m2_nd,m2_zz_nd,x)

  dydx(:) = 0.
  dydx(1) = y(2)
  dydx(2) = -mm*y(1)

END SUBROUTINE derivs

SUBROUTINE derivs_2nd(x,y,dydx)  ! for bsstep_2nd.f90

  use nrtype
  use T_G_eqn, only: z_nd, m2_nd, m2_zz_nd
  use nr, only: splint

  implicit none

  real(sp),               intent(in)  :: x
  real(sp), dimension(:), intent(in)  :: y
  real(sp), dimension(:), intent(out) :: dydx

  real ::  mm

  mm = splint(z_nd,m2_nd,m2_zz_nd,x)

  dydx(:) = 0.
  dydx(1) = -mm*y(1)
  dydx(2) = dydx(1)

END SUBROUTINE derivs_2nd


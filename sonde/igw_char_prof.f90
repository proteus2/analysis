MODULE igw_char_prof

  implicit none

  real, parameter, private ::  pi = 3.14159265358979323846
  real, parameter, private ::  rad2deg = 180./pi
  real, parameter, private ::  deg2rad = pi/180.


CONTAINS


SUBROUTINE waveno_m_avg(u_prt,v_prt,dz,wf,m0)
!-----------------------------------------------------------------------
! u_prt : u'
! v_prt : v'
! dz    : vertical spacing of the profile
! wf    : window function
!-----------------------------------------------------------------------
! m0    : mean vertical wavenumber
!-----------------------------------------------------------------------

  use fft

  real, dimension(:), intent(in)  ::  u_prt, v_prt, wf
  real,               intent(in)  ::  dz
  real,               intent(out) ::  m0

  integer                           ::  n
  real,    dimension(size(u_prt))   ::  m, windowed
  complex, dimension(size(u_prt))   ::  fc
  real,    dimension(size(u_prt)/2) ::  wgt

  n = size(u_prt)

  m(:) = waveno_fft(n,dz)

  wgt(:) = 0.

  windowed(:) = (u_prt(:) - sum(u_prt)/float(n))*wf(:)
  call fft1d_f(windowed,fc)
  wgt(2:n/2) = wgt(2:n/2) + real(fc(2:n/2)*conjg(fc(2:n/2)))

  windowed(:) = (v_prt(:) - sum(v_prt)/float(n))*wf(:)
  call fft1d_f(windowed,fc)
  wgt(2:n/2) = wgt(2:n/2) + real(fc(2:n/2)*conjg(fc(2:n/2)))

  ! 2*delta wavelength not considered: fc(n/2+1)

  m0 = sum(m(2:n/2)*wgt(2:n/2))/sum(wgt(2:n/2))

END subroutine waveno_m_avg

SUBROUTINE phase_dir(u_prt,v_prt,t_n,ph_dir)
!-----------------------------------------------------------------------
! u_prt  : u'
! v_prt  : v'
! t_n    : T'/T0 (T0 : basic-state T)
!-----------------------------------------------------------------------
! ph_dir : phase propagation direction [deg]
!-----------------------------------------------------------------------

  real, dimension(:), intent(in)  ::  u_prt, v_prt, t_n
  real,               intent(out) ::  ph_dir

  real                       ::  x, y
  real, dimension(size(t_n)) ::  t_n_hil

  call hilbert_transform(t_n,t_n_hil)
   
  x = sum(u_prt(:)*t_n_hil(:))
  y = sum(v_prt(:)*t_n_hil(:))

  if ( x == 0. .and. y == 0. ) then
    ph_dir = -999.  ;  RETURN
  end if

  if (abs(x) > abs(y)) then
    ph_dir = atan(y/x)*rad2deg
  else
    ph_dir = sign(90.-atan(abs(x/y))*rad2deg, x*y)
    if (x == 0.)  ph_dir = sign(90.,y)
  end if
  ! -180 <= ph_dir <= 180 [deg]

  if (x < 0.)  ph_dir = ph_dir + 180.
  ! -180 <= ph_dir < 360 [deg]

  if (ph_dir < 0.)  ph_dir = ph_dir + 360.
  ! 0 <= ph_dir < 360 [deg]

END subroutine phase_dir

SUBROUTINE stokes_param(u_prt,v_prt,i,d,p,q,i_m,d_m,p_m,q_m)
!-----------------------------------------------------------------------
! u_prt : u'
! v_prt : v'
!-----------------------------------------------------------------------
! i     : Stokes parameter I
! d     : Stokes parameter D
! p     : Stokes parameter P
! q     : Stokes parameter Q
! i_m   : Stokes parameter I(m), optional
! d_m   : Stokes parameter D(m), optional
! p_m   : Stokes parameter P(m), optional
! q_m   : Stokes parameter Q(m), optional
!-----------------------------------------------------------------------

  use fft

  real, dimension(:), intent(in)  ::  u_prt, v_prt
  real,               intent(out) ::  i, d, p, q
  real, dimension(size(u_prt)/2), optional, intent(out) ::  i_m, d_m, p_m, q_m

  integer                         ::  n
  real, dimension(size(u_prt)/2)  ::  im, dm, pm, qm
  complex, dimension(size(u_prt)) ::  fc_u, fc_v

  n = size(u_prt)/2

  call fft1d_f(u_prt,fc_u)
  call fft1d_f(v_prt,fc_v)

  fc_u(n+1) = 0.  ;  fc_v(n+1) = 0.
  ! 2*delta wavelength not considered

  im(:) = 2.*(real(fc_u(2:n+1)*conjg(fc_u(2:n+1)) +  &
                   fc_v(2:n+1)*conjg(fc_v(2:n+1))))
  dm(:) = 2.*(real(fc_u(2:n+1)*conjg(fc_u(2:n+1)) -  &
                   fc_v(2:n+1)*conjg(fc_v(2:n+1))))
  pm(:) = 4.*real(conjg(fc_u(2:n+1))*fc_v(2:n+1))
  qm(:) = 4.*aimag(conjg(fc_u(2:n+1))*fc_v(2:n+1))

  if ( present(i_m) )  i_m(:) = im(:)
  if ( present(d_m) )  d_m(:) = dm(:)
  if ( present(p_m) )  p_m(:) = pm(:)
  if ( present(q_m) )  q_m(:) = qm(:)

  i = sum(im)
  d = sum(dm)
  p = sum(pm)
  q = sum(qm)

END subroutine stokes_param

SUBROUTINE degree_polar(i,d,p,q,d_pol)
!-----------------------------------------------------------------------
! i     : Stokes parameter I
! d     : Stokes parameter D
! p     : Stokes parameter P
! q     : Stokes parameter Q
!-----------------------------------------------------------------------
! d_pol : degree of polarization
!-----------------------------------------------------------------------

  real, intent(in)  ::  i, d, p, q
  real, intent(out) ::  d_pol

  d_pol = sqrt(d*d + p*p + q*q)/i

END subroutine degree_polar

SUBROUTINE axial_ratio(d,p,q, r_ax)
!-----------------------------------------------------------------------
! d    : Stokes parameter D
! p    : Stokes parameter P
! q    : Stokes parameter Q
!-----------------------------------------------------------------------
! r_ax : ellipse axial ratio of the wind perturbation hodograph
!-----------------------------------------------------------------------

  real, intent(in)  ::  d, p, q
  real, intent(out) ::  r_ax

  r_ax = tan(0.5*asin(abs(q)/sqrt(d*d + p*p + q*q)))

END subroutine axial_ratio

SUBROUTINE intr_freq(d,p,q,ph_dir,f,u0s,v0s,n0,ome_i,r_ax)
!-----------------------------------------------------------------------
! d      : Stokes parameter D
! p      : Stokes parameter P
! q      : Stokes parameter Q
! ph_dir : phase propagation direction [deg]
! f      : Coriolis parameter
! u0s    : vertical-mean shear of the background u
! v0s    : vertical-mean shear of the background v
! n0     : vertical-mean background Brunt-Vaisala frequency
!-----------------------------------------------------------------------
! ome_i  : intrinsic frequency
! r_ax   : ellipse axial ratio of the wind perturbation hodograph,
!          optional
!-----------------------------------------------------------------------

  real,           intent(in)  ::  d, p, q, ph_dir, f, u0s, v0s, n0
  real,           intent(out) ::  ome_i
  real, optional, intent(out) ::  r_ax

  real ::  rax, shear_tr

  call axial_ratio(d,p,q,rax)
  if ( present(r_ax) )  r_ax = rax

  if ( u0s == 0. .and. v0s == 0. ) then
    ome_i = f/rax
  else
    ! shear projected to the direction transverse to ph_dir
    shear_tr = -u0s*sin(ph_dir*deg2rad) + v0s*cos(ph_dir*deg2rad)
    ome_i = f/(rax - shear_tr/n0)
  end if

END subroutine intr_freq

SUBROUTINE hilbert_transform(x,y)

  use fft

  real, dimension(:), intent(in)  ::  x
  real, dimension(:), intent(out) ::  y

  integer                     ::  n
  complex, dimension(size(x)) ::  fc

  n = size(x)

  call fft1d_f(x,fc)

  fc(2:n/2) = fc(2:n/2)*(0.,-1.)
  fc(n/2+1:n-n/2+1) = 0.
  fc(n-n/2+2:) = fc(n-n/2+2:)*(0.,1.)

  call fft1d_b(fc,y)

END subroutine hilbert_transform

END module igw_char_prof


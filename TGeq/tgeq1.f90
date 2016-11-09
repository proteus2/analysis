PROGRAM Taylor_Goldstein_eqn

  use T_G_eqn
  use nr, only: spline, splint
  use netio

  implicit none

  real,    parameter ::  zt = 25.e3
  integer, parameter ::  nz = 250+1
  real,    parameter ::  lh = 500.e3
  integer, parameter ::  bdy = 3
    ! 1, Dirichlet BDY at two specified pts        :  w(z1), w(z2)
    ! 2, Dirichlet BDY at k1 and Neumann BDY at k2 :  w(z1), w'(z2)
    ! 3, Cauchy BDY at k1 (initial condition)      :  w(z1), w'(z1)
  real,    parameter ::  z1 = zt, z2 = zt
  complex            ::  w1, w2, wp, r_coef
  data  w2 /vc_null/, wp /vc_null/, r_coef /vc_null/
  parameter  ( w1 = (1.,0.) )
!  parameter  ( w2 = (0.,0.) )
!  parameter  ( wp = (0,-6.3e-4) )
!  parameter  ( r_coef = (0.,0.) )
  real,    parameter ::  eps_err = 1.e-6

  integer, parameter ::  nphi = 72
  integer, parameter ::  nabs = 21

  real,    dimension(nz) ::  m2, z, wr, wi, wr_pm, wi_pm, wr_nm, wi_nm
  real,    dimension(nz) ::  r, f, f_pm, f_nm
  complex, dimension(nz) ::  w, wz, w_pm, w_nm

  real, dimension(nabs)      ::  r_coef_abs
  real, dimension(nphi)      ::  r_coef_phi
  real, dimension(nabs,nphi) ::  r2d, f2d_pm, f2d_nm

  integer ::  i,j,k, k1
  real    ::  dz, kwn
  complex ::  w1a, wpa, w2a, w_norm

  integer ::  ncid, nzi
  real    ::  phi, cp, hrho
  logical ::  lerr
  real, dimension(:), allocatable ::  zi, m2i, ubi, vbi, nb2i, tmpcs

  complex, parameter ::  ai = (0.,1.)

  eps_tg  = eps_err
  m2n_bvp = 0

  kwn = 2.*acos(-1.)/lh
!
! axis
!
  dz = zt/float(nz-1)
  do k=1, nz
    z(k) = float(k-1)*dz
  enddo

  do i=1, nabs
    r_coef_abs(i) = float(i-1)/float(nabs-1)
  enddo
  do j=1, nphi
    r_coef_phi(j) = float(j-1)/float(nphi)*360.
  enddo
!
! m^2 profile
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  cp  = 10.    ! horizontal phase speed
  phi = 0.     ! wave-propagation direction [rad]
  hrho = 7.e3  ! density scale height
  call opennc('background/uvn_background2_p1_region.nc',ncid)
  call dilen(ncid,'height_m',nzi,lerr)
  allocate( zi(nzi), m2i(nzi), tmpcs(nzi) )
  allocate( ubi(nzi), vbi(nzi), nb2i(nzi) )
  call get1d(ncid,'height_m',nzi,zi)
  call get1d(ncid,'ubavg',nzi,ubi)
  call get1d(ncid,'vbavg',nzi,vbi)
  call get1d(ncid,'nbavg',nzi,nb2i)
  nb2i(:) = nb2i(:)*nb2i(:)
  call closenc(ncid)
  zi(:) = zi(:) - zi(1)
  !! set z_bottom = 0.
  m2i(:) = nb2i(:)/(cp-(ubi(:)*cos(phi)+vbi(:)*sin(phi)))**2 - kwn*kwn &
           - 0.25/(hrho*hrho)
  deallocate( ubi, vbi, nb2i )
  call spline(zi,m2i,0.,0.,tmpcs)
  do k=1, nz
    m2(k) = splint(zi,m2i,tmpcs,z(k))
  enddo
  deallocate( zi, m2i, tmpcs )
!  m2(:) = ( 0.5*(9.+6.3)*1.e-4 - 0.5*(9.-6.3)*1.e-4*tanh(z(:)/1.e3-14.) )**2
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  m2(1 ) = m2(2   )
  m2(nz) = m2(nz-1)
!
! == loop ==============================================================

  R_PHI:  DO j=1, nphi
  R_ABS:  DO i=1, nabs
print*,i,j


  r_coef = r_coef_abs(i)*exp(ai*acos(-1.)/180.*r_coef_phi(j))
!
! boundary condition
!
  w1a = w1
  if (bdy == 1) then
    w2a = w2
  else
    wpa = wp
    if (wp == vc_null) then  ! wp is not defined
      if ( bdy == 3 .or. ( z1 == z2 ) ) then
        k1 = min(nz, nint(z1/dz)+1)
        if (m2(k1) > 0.) then
          if (r_coef /= vc_null) then
            wpa = w1*ai*(-1.)*sqrt(m2(k1))*(1.-r_coef)/(1.+r_coef)
          else
            print*, 'Either wp or r_coef should be defined for the B.C.'
            STOP
          end if
        else
          wpa = w1*(-1.)*sqrt(-m2(k1))
        end if
      else
        print*, 'wp(z2) should be defined for the B.C.'  ;  STOP
      end if
    end if
  end if
!
! w, momentum flux
!
  if (bdy == 3) then
    call t_g_eqn_ivp(nz,m2,zt,z1,w1a,wpa, w,wz,f)
  else
    if (bdy == 2)  w2a = wpa
    call t_g_eqn_bvp(bdy,nz,m2,zt,z1,w1a,z2,w2a, w,wz,f)
  end if
  f(:) = f(:)/kwn
!
! w, momentum flux for positive/negative m, reflection coef.
!
  call wave_decomp(nz,dz,m2,w,wz, w_pm,w_nm,f_pm,f_nm,r)
  f_pm(:) = f_pm(:)/kwn
  f_nm(:) = f_nm(:)/kwn
!
! normalizing by w_nm(1)
!
  w_norm = w_nm(1)
  w   (:) = w   (:)/w_norm
  w_pm(:) = w_pm(:)/w_norm
  w_nm(:) = w_nm(:)/w_norm
  f   (:) = f   (:)/abs(w_norm)**2
  f_pm(:) = f_pm(:)/abs(w_norm)**2
  f_nm(:) = f_nm(:)/abs(w_norm)**2
!
! output
!
  wr   (:) = real (w   (:))
  wi   (:) = aimag(w   (:))
  wr_pm(:) = real (w_pm(:))
  wi_pm(:) = aimag(w_pm(:))
  wr_nm(:) = real (w_nm(:))
  wi_nm(:) = aimag(w_nm(:))

  r2d(i,j) = r(1)
  f2d_pm(i,j) = f_pm(1)
  f2d_nm(i,j) = f_nm(1)


  ENDDO  R_ABS
  ENDDO  R_PHI

! ======================================================================
!
! dump
!
  call out1d('res.nc',11,(/'m2   ','wr','wi','r', &
             'wr_pm','wi_pm','wr_nm','wi_nm','f','f_pm','f_nm'/), &
             (/m2,wr,wi,r,wr_pm,wi_pm,wr_nm,wi_nm,f,f_pm,f_nm/), &
             'z',nz,z/1.e3,'T-G eqn solution')

  call out2d('res2.nc',3,(/'r   ','f_pm','f_nm'/),(/r2d,f2d_pm,f2d_nm/), &
             'R_abs',nabs,r_coef_abs,'R_phi',nphi,r_coef_phi, &
             'T-G eqn solution')


  STOP

END program Taylor_Goldstein_eqn


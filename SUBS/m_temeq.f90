MODULE temeq

!  Modules used here
!    utiletc
!    deriv                 
!    avg

  implicit none

  logical               ::  uniform_y = .FALSE., uniform_z = .FALSE.
  logical, dimension(2) ::  l_ybdy = (/.FALSE.,.FALSE./), &
                            l_zbdy = (/.FALSE.,.FALSE./)


  CONTAINS

! density approx.
! hydrostaticity
! z/zp coord.

!=======================================================================

SUBROUTINE utend_tem_hydp(nx,ny,nz,nt,lat,p,u,v,ome,t,outopt,         &
                          vmr,wmr,vadv_u,wadv_u,cor,epdrag,epfy,epfz, &
                          vadv_t,wadv_t,eddy_t,                       &
                          epfy1,epfy2,epfz1,epfz2,epdy1,epdy2,epdz1,epdz2,&
                          epdrags,epfys,epfzs, &
                          epfy1s,epfy2s,epfz1s,epfz2s,epdy1s,epdy2s,epdz1s,epdz2s)

! h_scale affects only wmr and epf (height coord.-dependent var.s).
! rho_s affects only epfy and epfz

! 1: v,w
! 2: epf, epf_stationary
! 3: v,w,adv_u,cor,epf
! 4: v,w,adv_t,t_chi
! 9: all (except epf_stationary)

  use avg
  use deriv
  use utiletc

  implicit none

  integer,                      intent(in) ::  nx, ny, nz, nt
  integer,                      intent(in) ::  outopt
  real, dimension(ny),          intent(in) ::  lat
  real, dimension(nz),          intent(in) ::  p
  real, dimension(nx,ny,nz,nt), intent(in) ::  u, v, ome, t

  real, dimension(ny,nz,nt), intent(inout) ::  vmr, wmr
  real, dimension(ny,nz,nt), intent(inout) ::  vadv_u, wadv_u, cor
  real, dimension(ny,nz,nt), intent(inout) ::  epdrag, epfy, epfz
  real, dimension(ny,nz,nt), intent(inout) ::  epdrags, epfys, epfzs
  real, dimension(ny,nz,nt), intent(inout) ::  vadv_t, wadv_t, eddy_t
  real, dimension(ny,nz,nt), intent(inout) ::  epfy1,epfy2,epfz1,epfz2
  real, dimension(ny,nz,nt), intent(inout) ::  epdy1,epdy2,epdz1,epdz2
  real, dimension(ny,nz,nt), intent(inout) ::  epfy1s,epfy2s,epfz1s,epfz2s
  real, dimension(ny,nz,nt), intent(inout) ::  epdy1s,epdy2s,epdz1s,epdz2s

  integer                      ::  j,k,n
  real, dimension(ny)          ::  y
  real, dimension(nz)          ::  zp, tempz
  real, dimension(ny,nz,nt)    ::  rho0, exner
  real, dimension(ny,nz,nt)    ::  um, vm, wm, ptm, temp3d
  real, dimension(ny,nz,nt)    ::  divy_um, dumdz, dptmdy, dptmdz
  real, dimension(ny,nz,nt)    ::  rhovum, rhowum, vptm, wptm, psi, chi
  real, dimension(ny,nz,nt)    ::  st1,st2,st3,st4,st5,st6,st7,st8,st9,st10,st11
  real, dimension(nx,ny,nz,nt) ::  uprt, vprt, wprt, ptprt, temp4d
  real, dimension(nx,ny,nz,1 ) ::  u_s, v_s, w_s, pt_s

  include 'c_math.inc'
  include 'c_phys.inc'


print*, 'perturbation'
! calculate zonal mean and perturbation

  y(:) = r_earth*lat(:)*deg2rad

  do k=1, nz
    zp(k) = -h_scale*log(p(k)/p0)
  enddo

  call zm_prt(nx,ny,nz,nt,u,1., um,uprt)
  call zm_prt(nx,ny,nz,nt,v,1., vm,vprt)

  do n=1, nt
  do k=1, nz
    temp4d(:,:,k,n) = -ome(:,:,k,n)*h_scale/p(k)
  enddo
  enddo
  call zm_prt(nx,ny,nz,nt,temp4d,1., wm,wprt)

  tempz(:) = (p(:)/p0)**kappa
  do n=1, nt
  do k=1, nz
    temp4d(:,:,k,n) = t(:,:,k,n)/tempz(k)
  enddo
  enddo
  call zm_prt(nx,ny,nz,nt,temp4d,1., ptm,ptprt)

  do n=1, nt
  do k=1, nz
    exner(:,k,n) = tempz(k)
  enddo
  enddo

  if (outopt == 2) then
    call avg_d4((/nx,ny,nz,nt/),uprt ,1.,u_s )
    call avg_d4((/nx,ny,nz,nt/),vprt ,1.,v_s )
    call avg_d4((/nx,ny,nz,nt/),wprt ,1.,w_s )
    call avg_d4((/nx,ny,nz,nt/),ptprt,1.,pt_s)
  end if

print*, 'grad'
! calculate the gradient of um, ptm

  call div_lat(1,ny,nz,nt,uniform_y,lat,um,l_ybdy,0., divy_um)
  call deriv1d((/ny,nz,nt,1/),2,uniform_z,zp,um,l_zbdy,0., dumdz)
  call deriv1d_exp((/ny,nz,nt,1/),2,uniform_z,zp,ptm,h_scale/kappa, &
                   l_zbdy,0., dptmdz)

  if ( outopt == 4 .or. outopt == 9 )                &
     call deriv1d((/ny,nz,nt,1/),1,uniform_y,y,ptm,l_ybdy,0., dptmdy)

print*, 'flx'
! calculate zonal-mean wave flux

  tempz(:) = rho_s*exp(-zp(:)/h_scale)
  do k=1, nz
    rho0(:,k,:) = tempz(k)
  enddo

  call zmflx(nx,ny,nz,nt,vprt,uprt,1., rhovum)
  call zmflx(nx,ny,nz,nt,wprt,uprt,1., rhowum)
  rhovum(:,:,:) = rho0(:,:,:)*rhovum(:,:,:)
  rhowum(:,:,:) = rho0(:,:,:)*rhowum(:,:,:)

  call zmflx(nx,ny,nz,nt,vprt,ptprt,1., vptm)
  if ( outopt == 4 .or. outopt == 9 )  &
     call zmflx(nx,ny,nz,nt,wprt,ptprt,1., wptm)

print*, 'psi,chi'
! calculate psi and chi in hydro. eqn.

  call psi_tem_hyd(ny,nz,nt,rho0,dptmdz,vptm, psi)

  if ( outopt == 4 .or. outopt == 9 )  &
     call chi_tem_hyd(ny,nz,nt,rho0,dptmdy,dptmdz,vptm,wptm, chi)

! calculate residual meridional velocities and adv. / Cor.

  if (outopt /= 2) then

    call res_mer_vel(ny,nz,nt,lat,zp,rho0,vm,wm,psi, vmr,wmr)
    ! save wm at eddy_t
    if (outopt == 1)  eddy_t(:,:,:) = wm(:,:,:)

    if ( outopt == 3 .or. outopt == 9 )                    &
       call utend_mer(ny,nz,nt,lat,vmr,wmr,divy_um,dumdz, &
                      vadv_u,wadv_u,cor)

    if ( outopt == 4 .or. outopt == 9 )                                 &
       call ttend_hyd(ny,nz,nt,zp,rho0,exner,vmr,wmr,dptmdy,dptmdz,chi, &
                      vadv_t,wadv_t,eddy_t)

  end if

print*, 'epf'
! calculate drag by E-P flux

  if ( outopt == 2 .or. outopt == 3 .or. outopt == 9 )                &
     call epflx(ny,nz,nt,lat,zp,rho0,divy_um,dumdz,rhovum,rhowum,psi, &
                epfy,epfz,epdrag,                                     &
                epfy1,epfy2,epfz1,epfz2,epdy1,epdy2,epdz1,epdz2)

  if ( outopt == 2 ) then
    call zmflx(nx,ny,nz,1,v_s,u_s,1., rhovum(:,:,1))
    call zmflx(nx,ny,nz,1,w_s,u_s,1., rhowum(:,:,1))
    rhovum(:,:,1) = rho0(:,:,1)*rhovum(:,:,1)
    rhowum(:,:,1) = rho0(:,:,1)*rhowum(:,:,1)
    call zmflx(nx,ny,nz,1,v_s,pt_s,1., vptm(:,:,1))

    do n=2, nt
      rhovum(:,:,n) = rhovum(:,:,1)
      rhowum(:,:,n) = rhowum(:,:,1)
      vptm  (:,:,n) = vptm  (:,:,1)
    enddo

    call psi_tem_hyd(ny,nz,nt,rho0,dptmdz,vptm, psi)

    call epflx(ny,nz,nt,lat,zp,rho0,divy_um,dumdz,rhovum,rhowum,psi, &
               epfys,epfzs,epdrags,                                  &
               epfy1s,epfy2s,epfz1s,epfz2s,epdy1s,epdy2s,epdz1s,epdz2s)
  end if

  RETURN

END subroutine utend_tem_hydp

!=======================================================================

SUBROUTINE psi_tem_hyd(ny,nz,nt,rho0,dptmdz,vptm, psi)

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny,nz,nt), intent(in) ::  rho0
  real, dimension(ny,nz,nt), intent(in) ::  dptmdz, vptm

  real, dimension(ny,nz,nt), intent(out) ::  psi

  integer ::  n


  psi(:,:,:) = rho0(:,:,:)*vptm(:,:,:)/dptmdz(:,:,:)

  RETURN

END subroutine psi_tem_hyd

!=======================================================================

SUBROUTINE psi_tem_nonhyd(ny,nz,nt,rho0,dptmdy,dptmdz,vptm,wptm, psi)

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny,nz,nt), intent(in) ::  rho0
  real, dimension(ny,nz,nt), intent(in) ::  dptmdy, dptmdz, vptm, wptm

  real, dimension(ny,nz,nt), intent(out) ::  psi

  integer ::  n


  psi(:,:,:) = rho0(:,:,:)*( vptm(:,:,:)*dptmdz(:,:,:) - &
                             wptm(:,:,:)*dptmdy(:,:,:) ) &
                       / ( dptmdy(:,:,:)*dptmdy(:,:,:) + &
                           dptmdz(:,:,:)*dptmdz(:,:,:) )

  RETURN

END subroutine psi_tem_nonhyd

!=======================================================================

SUBROUTINE chi_tem_hyd(ny,nz,nt,rho0,dptmdy,dptmdz,vptm,wptm, chi)

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny,nz,nt), intent(in) ::  rho0
  real, dimension(ny,nz,nt), intent(in) ::  dptmdy, dptmdz, vptm, wptm

  real, dimension(ny,nz,nt), intent(out) ::  chi

  integer ::  n


  chi(:,:,:) = rho0(:,:,:)*( vptm(:,:,:)*dptmdy(:,:,:) + &
                             wptm(:,:,:)*dptmdz(:,:,:) ) &
                       / ( dptmdz(:,:,:)*dptmdz(:,:,:) )

  RETURN

END subroutine chi_tem_hyd

!=======================================================================

SUBROUTINE chi_tem_nonhyd(ny,nz,nt,rho0,dptmdy,dptmdz,vptm,wptm, chi)

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny,nz,nt), intent(in) ::  rho0
  real, dimension(ny,nz,nt), intent(in) ::  dptmdy, dptmdz, vptm, wptm

  real, dimension(ny,nz,nt), intent(out) ::  chi

  integer ::  n


  chi(:,:,:) = rho0(:,:,:)*( vptm(:,:,:)*dptmdy(:,:,:) + &
                             wptm(:,:,:)*dptmdz(:,:,:) ) &
                       / ( dptmdy(:,:,:)*dptmdy(:,:,:) + &
                           dptmdz(:,:,:)*dptmdz(:,:,:) )

  RETURN

END subroutine chi_tem_nonhyd

!=======================================================================

SUBROUTINE res_mer_vel(ny,nz,nt,lat,z,rho0,vm,wm,psi, vmr,wmr)

  use deriv

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  z
  real, dimension(ny,nz,nt), intent(in) ::  rho0
  real, dimension(ny,nz,nt), intent(in) ::  vm, wm, psi

  real, dimension(ny,nz,nt), intent(out) ::  vmr, wmr

  integer                   ::  n
  real, dimension(ny,nz,nt) ::  temp


  call deriv1d((/ny,nz,nt,1/),2,uniform_z,z,psi,l_zbdy,0., temp)
  vmr(:,:,:) = vm(:,:,:) - temp(:,:,:)/rho0(:,:,:)

  call div_lat(1,ny,nz,nt,uniform_y,lat,psi,l_ybdy,0., temp)
  wmr(:,:,:) = wm(:,:,:) + temp(:,:,:)/rho0(:,:,:)

  RETURN

END subroutine res_mer_vel

!=======================================================================

SUBROUTINE utend_mer(ny,nz,nt,lat,vmr,wmr,divy_um,dumdz, vadv,wadv,cor)

  use utiletc

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(ny,nz,nt), intent(in) ::  vmr, wmr, divy_um, dumdz

  real, dimension(ny,nz,nt), intent(out) ::  vadv, wadv, cor

  integer ::  k,n
  real    ::  tempy(ny), f(ny,nz)


  vadv(:,:,:) = -vmr(:,:,:)*divy_um(:,:,:)
  wadv(:,:,:) = -wmr(:,:,:)*dumdz  (:,:,:)

  call f_coriolis(ny,lat, tempy)
  do k=1, nz
    f(:,k) = tempy(:)
  enddo
  do n=1, nt
    cor(:,:,n) = f(:,:)*vmr(:,:,n)
  enddo

  RETURN

END subroutine utend_mer

!=======================================================================

SUBROUTINE ttend_hyd(ny,nz,nt,zp,rho0,exner,vmr,wmr,dptmdy,dptmdz,chi, &
                     vadv,wadv,eddy)

  use deriv

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(nz),       intent(in) ::  zp
  real, dimension(ny,nz,nt), intent(in) ::  rho0, exner
  real, dimension(ny,nz,nt), intent(in) ::  vmr, wmr, dptmdy, dptmdz
  real, dimension(ny,nz,nt), intent(in) ::  chi

  real, dimension(ny,nz,nt), intent(out) ::  vadv, wadv, eddy

  integer ::  k,n
  real    ::  temp3d(ny,nz,nt)


  vadv(:,:,:) = -vmr(:,:,:)*dptmdy(:,:,:)*exner(:,:,:)
  wadv(:,:,:) = -wmr(:,:,:)*dptmdz(:,:,:)*exner(:,:,:)

  temp3d(:,:,:) = chi(:,:,:)*dptmdz(:,:,:)
  call deriv1d((/ny,nz,nt,1/),2,uniform_z,zp,temp3d,l_zbdy,0., eddy)
  eddy(:,:,:) = -eddy(:,:,:)/rho0(:,:,:)*exner(:,:,:)

  RETURN

END subroutine ttend_hyd

!=======================================================================

SUBROUTINE ttend_nonhyd(ny,nz,nt,lat,zp,rho0,exner,vmr,wmr,dptmdy, &
                        dptmdz,chi, vadv,wadv,eddy)

  use deriv

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  zp
  real, dimension(ny,nz,nt), intent(in) ::  rho0, exner
  real, dimension(ny,nz,nt), intent(in) ::  vmr, wmr, dptmdy, dptmdz
  real, dimension(ny,nz,nt), intent(in) ::  chi

  real, dimension(ny,nz,nt), intent(out) ::  vadv, wadv, eddy

  integer ::  k,n
  real    ::  temp2d(ny,nz), temp3d1(ny,nz,nt), temp3d2(ny,nz,nt)


  vadv(:,:,:) = -vmr(:,:,:)*dptmdy(:,:,:)*exner(:,:,:)
  wadv(:,:,:) = -wmr(:,:,:)*dptmdz(:,:,:)*exner(:,:,:)

  temp3d1(:,:,:) = chi(:,:,:)*dptmdy(:,:,:)
  temp3d2(:,:,:) = chi(:,:,:)*dptmdz(:,:,:)
  call div_lat(1,ny,nz,nt,uniform_y,lat,temp3d1,l_ybdy,0., eddy)
  call deriv1d((/ny,nz,nt,1/),2,uniform_z,zp,temp3d2,l_zbdy,0., temp3d1)
  eddy(:,:,:) = -(eddy(:,:,:)+temp3d1(:,:,:))/rho0(:,:,:)*exner(:,:,:)

  RETURN

END subroutine ttend_nonhyd

!=======================================================================

SUBROUTINE epflx(ny,nz,nt,lat,z,rho0,divy_um,dumdz, &
                 rhovum,rhowum,psi,epfy,epfz,epdrag,&
                 epfy1,epfy2,epfz1,epfz2,epdy1,epdy2,epdz1,epdz2)
  use utiletc
  use deriv

  implicit none

  integer,                   intent(in) ::  ny, nz, nt
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  z
  real, dimension(ny,nz,nt), intent(in) ::  rho0
  real, dimension(ny,nz,nt), intent(in) ::  divy_um, dumdz
  real, dimension(ny,nz,nt), intent(in) ::  rhovum, rhowum, psi

  real, dimension(ny,nz,nt), intent(out) ::  epfy, epfz, epdrag
  real, dimension(ny,nz,nt), intent(out) ::  epfy1,epfy2,epfz1,epfz2
  real, dimension(ny,nz,nt), intent(out) ::  epdy1,epdy2,epdz1,epdz2

  integer                   ::  k,n, iy1, iy2
  real, dimension(ny)       ::  cosphi
  real, dimension(ny,nz,nt) ::  temp3d1, temp3d2
  real                      ::  tempy(ny), tempyz(ny,nz)

  include 'c_phys.inc'


  call f_coriolis(ny,lat, tempy)
  do n=1, nt
  do k=1, nz
    temp3d1(:,k,n) = tempy(:)
  enddo
  enddo
  call coslat(ny,lat, cosphi)

  ! r_earth*cos(phi)
  tempy(:) = r_earth*cosphi(:)
  do n=1, nt
  do k=1, nz
    temp3d2(:,k,n) = tempy(:)
  enddo
  enddo

  ! Fy
  epfy1(:,:,:) = -temp3d2(:,:,:)*rhovum(:,:,:)
  epfy2(:,:,:) = -temp3d2(:,:,:)*(-psi(:,:,:)*dumdz(:,:,:))

  epfy = epfy1 + epfy2

  ! Fz
  epfz1(:,:,:) = -temp3d2(:,:,:)*rhowum(:,:,:)
  epfz2(:,:,:) = -temp3d2(:,:,:)*(psi(:,:,:)*(divy_um(:,:,:)-temp3d1(:,:,:)))

  epfz = epfz1 + epfz2

  ! rho0*r_earth*cos(phi)
  temp3d2(:,:,:) = rho0(:,:,:)*temp3d2(:,:,:)

  ! Drag by EP-div. [m/s/day]
  !  = Div[F] / [rho0*r_earth*cos(phi)]
  call div_lat(1,ny,nz,nt,uniform_y,lat,epfy1,l_ybdy,0., epdy1)
  call div_lat(1,ny,nz,nt,uniform_y,lat,epfy2,l_ybdy,0., epdy2)
  call deriv1d((/ny,nz,nt,1/),2,uniform_z,z,epfz1,l_zbdy,0., epdz1)
  call deriv1d((/ny,nz,nt,1/),2,uniform_z,z,epfz2,l_zbdy,0., epdz2)

  if ( l_ybdy(1) .and. l_ybdy(2) ) then
    epdy1(:,:,:) = epdy1(:,:,:)/temp3d2(:,:,:) * 86400.
    epdy2(:,:,:) = epdy2(:,:,:)/temp3d2(:,:,:) * 86400.
    epdz1(:,:,:) = epdz1(:,:,:)/temp3d2(:,:,:) * 86400.
    epdz2(:,:,:) = epdz2(:,:,:)/temp3d2(:,:,:) * 86400.
  else
    epdy1(1,:,:) = 0.  ;  epdy1(ny,:,:) = 0.
    epdy2(1,:,:) = 0.  ;  epdy2(ny,:,:) = 0.
    epdz1(1,:,:) = 0.  ;  epdz1(ny,:,:) = 0.
    epdz2(1,:,:) = 0.  ;  epdz2(ny,:,:) = 0.
    iy1 = 2  ;  iy2 = ny-1
    if ( l_ybdy(1) )  iy1 = 1
    if ( l_ybdy(2) )  iy2 = ny
    epdy1(iy1:iy2,:,:) = epdy1(iy1:iy2,:,:)/temp3d2(iy1:iy2,:,:) * 86400.
    epdy2(iy1:iy2,:,:) = epdy2(iy1:iy2,:,:)/temp3d2(iy1:iy2,:,:) * 86400.
    epdz1(iy1:iy2,:,:) = epdz1(iy1:iy2,:,:)/temp3d2(iy1:iy2,:,:) * 86400.
    epdz2(iy1:iy2,:,:) = epdz2(iy1:iy2,:,:)/temp3d2(iy1:iy2,:,:) * 86400.
  end if
  epdrag = epdy1 + epdy2 + epdz1 + epdz2

  RETURN

END subroutine epflx

!=======================================================================

SUBROUTINE wmr_dynbal_hydp(ny,nz,nt,dt,lat,p,um,ndr,drag,ilat_avg, &
                           lato,nyo,wmr_d)

! following Randel et al. (2002)

  use avg
  use deriv
  use integ
  use utiletc

  implicit none

  integer,                       intent(in) ::  ny, nz, nt, ndr
  integer,                       intent(in) ::  ilat_avg(2), nyo
  real,                          intent(in) ::  dt
  real, dimension(ny),           intent(in) ::  lat
  real, dimension(nz),           intent(in) ::  p
  real, dimension(ny,nz,nt),     intent(in) ::  um
  real, dimension(ny,nz,nt,ndr), intent(in) ::  drag

  real, dimension(nyo),          intent(out) ::  lato
  real, dimension(nyo,nz,ndr+2), intent(out) ::  wmr_d

  integer                   ::  j,k,n,m, k1, kn, ki, j1, j2
  real                      ::  ztop, ypole1, ypole2
  real, dimension(ny)       ::  y, f, cosphi, rsinphi
  real, dimension(nz)       ::  zp, rho0s, rrdsinphi
  real, dimension(ny)       ::  trid_a, trid_b, trid_c, trid_d
  real, dimension(nz)       ::  tempz1, tempz2
  real, dimension(0:ny+1)   ::  fn_str2
  real, dimension(ny,nz)    ::  rho0, r_fac
  real, dimension(nyo,nz)   ::  temp2d
  real, dimension(ny,nz,nt) ::  f_hat, dumdzdfh, dumdt, temp3d
  real, dimension(ny,nz,nt) ::  fn_str
  real, dimension(3,ny)     ::  dcoef
  real*8                    ::  dcoef8(3)
  logical                   ::  l_zrev

  include 'c_math.inc'
  include 'c_phys.inc'


  j1 = ilat_avg(1)  ;  j2 = ilat_avg(2)

  if (p(1) < p(2))  l_zrev = .TRUE.

! latr, zp, cosphir, rsinphir, rho0r, rho0o, f_hatr, dumdtr

  do k=1, nz
    zp(k) = -h_scale*log(p(k)/p0)
  enddo
  ztop = maxval(zp)

  y(:) = r_earth*lat(:)*deg2rad
  ypole1 = r_earth*90.*deg2rad
  ypole1 = sign(ypole1,y(1))
  ypole2 = -ypole1

  call coslat(ny,lat, cosphi)
  rsinphi(:) = r_earth*sin(lat(:)*deg2rad)

  rho0s(:) = rho_s*exp(-zp(:)/h_scale)

  ! for hydrostatic eqn.
  do k=1, nz
    rho0(:,k) = rho0s(k)
  enddo
  r_fac(:,:) = 1.

!  ! for rho0(y,z) in the nonhydrostatic eqn.
!  do k=1, nz
!    r_fac(:,k) = rho0(:,k) / rho0s(k)
!  enddo

  rrdsinphi(:) = 0.
  do j=j1, j2-1
    rrdsinphi(:) = rrdsinphi(:) + 0.5*(r_fac(j,:)+r_fac(j+1,:))* &
                   (rsinphi(j+1)-rsinphi(j))
  enddo


  call f_coriolis(ny,lat, f)
  call div_lat(1,ny,nz,nt,uniform_y,lat,um,l_ybdy,0., temp3d)
  do n=1, nt
  do k=1, nz
    f_hat(:,k,n) = f(:) - temp3d(:,k,n)
  enddo
  enddo

  call deriv1d((/ny,nz,nt,1/),2,uniform_z,zp,um,l_zbdy,0., dumdzdfh)
  dumdzdfh(:,:,:) = dumdzdfh(:,:,:) / f_hat(:,:,:)

  call deriv1d_uni((/ny,nz,nt,1/),3,dt,um,(/.FALSE.,.FALSE./),0., dumdt)

  dcoef(:,:) = 0.
  do j=2, ny-1
    call fdcoef(2,3,dble(y(j)),dble(y(j-1:j+1)), dcoef8)
    dcoef(:,j) = dcoef8(:)
  enddo
  if (abs(lat(1)) /= 90.) then
    call fdcoef(2,3,dble(y(1)),dble((/ypole1,y(1),y(2)/)), dcoef8)
    dcoef(:,1) = dcoef8(:)
  end if
  if (abs(lat(ny)) /= 90.) then
    call fdcoef(2,3,dble(y(ny)),dble((/y(ny-1),y(ny),ypole2/)), dcoef8)
    dcoef(:,ny) = dcoef8(:)
  end if

  dcoef(1,j1) = -1./(y(j1)-y(j1-1))
  dcoef(2,j1) =  1./(y(j1)-y(j1-1))
  dcoef(3,j1) =  0.

  dcoef(1,j2) =  0.
  dcoef(2,j2) = -1./(y(j2+1)-y(j2))
  dcoef(3,j2) =  1./(y(j2+1)-y(j2))

dcoef = 0.

  if ( l_zrev ) then
    k1 = nz
    kn = 2
    ki = 1
  else
    k1 = 1
    kn = nz-1
    ki = -1
  end if

  do k=kn, k1, ki
    tempz1(k) = 2./(zp(k)-zp(k-ki)) - 1./h_scale
    tempz2(k) = 2./(zp(k)-zp(k-ki)) + 1./h_scale
  enddo

  wmr_d(:,:,1) = 0.


! for each forcing

  N_DR:  DO m=2, ndr+2

  if (m == 2) then
    do n=1, nt
      temp3d(:,:,n) = -dumdt(:,:,n)*r_fac(:,:)/f_hat(:,:,n)
    enddo
  else
    do n=1, nt
      temp3d(:,:,n) = drag(:,:,n,m-2)*r_fac(:,:)/f_hat(:,:,n)
    enddo
  end if

  fn_str(:,:,:) = 0.
  fn_str2(0) = 0.  ;  fn_str2(ny+1) = 0.

  do n=2, nt-1
  do k=kn, k1, ki

    trid_a(:) = dcoef(1,:)*dumdzdfh(:,k,n)
    trid_b(:) = dcoef(2,:)*dumdzdfh(:,k,n) + tempz1(k)
    trid_c(:) = dcoef(3,:)*dumdzdfh(:,k,n)

    fn_str2(1:ny) = fn_str(:,k-ki,n)
    do j=1, ny
      trid_d(j) = -sum( dcoef(:,j)*fn_str2(j-1:j+1) )* &
                  dumdzdfh(j,k-k1,n)
    enddo
    trid_d(:) = trid_d(:) + tempz2(k)*fn_str(:,k-ki,n) + &
                cosphi(:)*(temp3d(:,k-ki,n)+temp3d(:,k,n))

    ! boundary condition
    trid_a(1 ) = 0.
    trid_c(ny) = 0.

    ! solve
    call tridag(j1,trid_a(1:j1),trid_b(1:j1),trid_c(1:j1), &
                trid_d(1:j1))
    fn_str(1:j1,k,n) = trid_d(1:j1)

    call tridag(ny-j2+1,trid_a(j2:ny),trid_b(j2:ny),trid_c(j2:ny), &
                trid_d(j2:ny))
    fn_str(j2:ny,k,n) = trid_d(j2:ny)

  enddo
  enddo

  call avg_d3((/j1,nz,nt-2,1/),fn_str(1:j1,:,2:nt-1),1., temp2d(1:j1,:))

  call avg_d3((/nyo-j1,nz,nt-2,1/),fn_str(j2:ny,:,2:nt-1),1., &
              temp2d(j1+1:nyo,:))

  call deriv1d((/j1,nz,1,1/),1,.False.,rsinphi(1:j1), &
               temp2d(1:j1,:),l_ybdy,0., wmr_d(1:j1,:,m))

  call deriv1d((/nyo-j1,nz,1,1/),1,.False.,rsinphi(j2:ny), &
               temp2d(j1+1:nyo,:),l_ybdy,0., wmr_d(j1+1:nyo,:,m))

  wmr_d(1:j1-1  ,:,m) = wmr_d(1:j1-1  ,:,m) / r_fac(1:j1-1 ,:)
  wmr_d(j1+2:nyo,:,m) = wmr_d(j1+2:nyo,:,m) / r_fac(j2+1:ny,:)

  wmr_d(j1  ,:,m) = (temp2d(j1+1,:)-temp2d(j1,:))/rrdsinphi(:)
  wmr_d(j1+1,:,m) = wmr_d(j1,:,m)

  wmr_d(:,:,1) = wmr_d(:,:,1) + wmr_d(:,:,m)


  ENDDO  N_DR


  lato(1   :j1 ) = lat(1 :j1)
  lato(j1+1:nyo) = lat(j2:ny)

  RETURN

END subroutine wmr_dynbal_hydp

!=======================================================================

END module temeq


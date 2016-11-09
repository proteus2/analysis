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
                          vadv_t,wadv_t,eddy_t)

! h_scale affects only wmr and epf (height coord.-dependent var.s).
! rho_s affects only epfy and epfz

! 1: v,w
! 2: epf
! 3: v,w,adv_u,cor,epf
! 4: v,w,adv_t,t_chi
! 9: all

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
  real, dimension(ny,nz,nt), intent(inout) ::  vadv_t, wadv_t, eddy_t

  integer                      ::  j,k,n
  integer                      ::  revy, revz
  real, dimension(ny)          ::  lata, ya
  real, dimension(nz)          ::  zp, tempz
  real, dimension(ny,nz,nt)    ::  rho0, exner
  real, dimension(ny,nz,nt)    ::  um, vm, wm, ptm, temp3d
  real, dimension(ny,nz,nt)    ::  divy_um, dumdz, dptmdy, dptmdz
  real, dimension(ny,nz,nt)    ::  rhovum, rhowum, vptm, wptm, psi, chi
  real, dimension(nx,ny,nz,nt) ::  uprt, vprt, wprt, ptprt, temp4d

  include 'c_math.inc'
  include 'c_phys.inc'


  revy = 0  ;  revz = 0
  if (lat(1) > lat(2))  revy = 1
  if (p  (1) < p  (2))  revz = 1
print*, 'reorder'

! reorder variables / calculate zonal mean and perturbation
! --> lata, zp, um,vm,wm,ptm, uprt,vprt,wprt,ptprt

  call reorder((/ny,1,1,1/),(/revy,0,0,0/),lat, lata)
  ya(:) = r_earth*lata(:)*deg2rad

  do k=1, nz
    zp(k) = -h_scale*log(p(k)/p0)
  enddo
  call reorder((/nz,1,1,1/),(/revz,0,0,0/),zp, zp)

  call reorder((/nx,ny,nz,nt/),(/0,revy,revz,0/),u, temp4d)
  call zm_prt(nx,ny,nz,nt,temp4d,1., um,uprt)

  call reorder((/nx,ny,nz,nt/),(/0,revy,revz,0/),v, temp4d)
  call zm_prt(nx,ny,nz,nt,temp4d,1., vm,vprt)

  do n=1, nt
  do k=1, nz
    temp4d(:,:,k,n) = -ome(:,:,k,n)*h_scale/p(k)
  enddo
  enddo
  call reorder((/nx,ny,nz,nt/),(/0,revy,revz,0/),temp4d, temp4d)
  call zm_prt(nx,ny,nz,nt,temp4d,1., wm,wprt)

  tempz(:) = (p(:)/p0)**kappa
  do n=1, nt
  do k=1, nz
    temp4d(:,:,k,n) = t(:,:,k,n)/tempz(k)
  enddo
  enddo
  call reorder((/nx,ny,nz,nt/),(/0,revy,revz,0/),temp4d, temp4d)
  call zm_prt(nx,ny,nz,nt,temp4d,1., ptm,ptprt)

  call reorder((/nz,1,1,1/),(/revz,0,0,0/),tempz, tempz)
  do n=1, nt
  do k=1, nz
    exner(:,k,n) = tempz(k)
  enddo
  enddo
print*, 'grad'

! calculate the gradient of um, ptm

  call div_lat(1,ny,nz,nt,uniform_y,lata,um,l_ybdy,0., divy_um)
  call deriv1d((/ny,nz,nt,1/),2,uniform_z,zp,um,l_zbdy,0., dumdz)
  call deriv1d_exp((/ny,nz,nt,1/),2,uniform_z,zp,ptm,h_scale/kappa, &
                   l_zbdy,0., dptmdz)

  if ( outopt == 4 .or. outopt == 9 )                &
     call deriv1d((/ny,nz,nt,1/),1,uniform_y,ya,ptm,l_ybdy,0., dptmdy)
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

    call res_mer_vel(ny,nz,nt,lata,zp,rho0,vm,wm,psi, vmr,wmr)
    ! save wm at eddy_t
    if (outopt == 1)  eddy_t(:,:,:) = wm(:,:,:)

    if ( outopt == 3 .or. outopt == 9 )                    &
       call utend_mer(ny,nz,nt,lata,vmr,wmr,divy_um,dumdz, &
                      vadv_u,wadv_u,cor)

    if ( outopt == 4 .or. outopt == 9 )                                 &
       call ttend_hyd(ny,nz,nt,zp,rho0,exner,vmr,wmr,dptmdy,dptmdz,chi, &
                      vadv_t,wadv_t,eddy_t)

  end if
print*, 'epf'

! calculate drag by E-P flux

  if ( outopt == 2 .or. outopt == 3 .or. outopt == 9 )                 &
     call epflx(ny,nz,nt,lata,zp,rho0,divy_um,dumdz,rhovum,rhowum,psi, &
                epfy,epfz,epdrag)
print*,'last reorder'

! reorder output variables

  if (outopt /= 2) then
    call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),vmr, vmr)
    call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),wmr, wmr)
    if (outopt == 1) then
      call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),eddy_t, eddy_t)
    end if
    if ( outopt == 3 .or. outopt == 9 ) then 
      call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),vadv_u, vadv_u)
      call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),wadv_u, wadv_u)
      call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),cor   , cor   )
    end if
    if ( outopt == 4 .or. outopt == 9 ) then
      call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),vadv_t, vadv_t)
      call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),wadv_t, wadv_t)
      call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),eddy_t, eddy_t)
    end if
  end if
  if ( outopt == 2 .or. outopt == 3 .or. outopt == 9 ) then
    call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),epfy  , epfy  )
    call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),epfz  , epfz  )
    call reorder((/ny,nz,nt,1/),(/revy,revz,0,0/),epdrag, epdrag)
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
                 rhovum,rhowum,psi,epfy,epfz,epdrag)

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
  epfy(:,:,:) = -temp3d2(:,:,:)*( rhovum(:,:,:) - psi(:,:,:)* &
                dumdz(:,:,:) )

  ! Fz
  epfz(:,:,:) = -temp3d2(:,:,:)*( rhowum(:,:,:) + psi(:,:,:)* &
                (divy_um(:,:,:)-temp3d1(:,:,:)) )

  ! rho0*r_earth*cos(phi)
  temp3d2(:,:,:) = rho0(:,:,:)*temp3d2(:,:,:)

  ! Drag by EP-div. [m/s/day]
  !  = Div[F] / [rho0*r_earth*cos(phi)]
  call div_lat(1,ny,nz,nt,uniform_y,lat,epfy,l_ybdy,0., epdrag)
  call deriv1d((/ny,nz,nt,1/),2,uniform_z,z,epfz,l_zbdy,0., temp3d1)
  if ( l_ybdy(1) .and. l_ybdy(2) ) then
    epdrag(:,:,:) = (epdrag(:,:,:)+temp3d1(:,:,:))/temp3d2(:,:,:)
  else
    epdrag(1 ,:,:) = 0.
    epdrag(ny,:,:) = 0.
    iy1 = 2  ;  iy2 = ny-1
    if ( l_ybdy(1) )  iy1 = 1
    if ( l_ybdy(2) )  iy2 = ny
    epdrag(iy1:iy2,:,:) = (epdrag(iy1:iy2,:,:)+temp3d1(iy1:iy2,:,:)) &
                         / temp3d2(iy1:iy2,:,:)
  end if
  epdrag(:,:,:) = epdrag(:,:,:) * 86400.

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
  real                      ::  ztop
  real, dimension(ny)       ::  y, f, cosphi, rsinphi
  real, dimension(nz)       ::  zp, rho0s, rrdsinphi
  real, dimension(ny)       ::  trid_a, trid_b, trid_c, trid_d
  real, dimension(nz)       ::  tempz1, tempz2
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

  do j=2, ny-1
    call fdcoef(2,3,dble(y(j)),dble(y(j-1:j+1)), dcoef8)
    dcoef(:,j) = dcoef8(:)
  enddo

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

  do n=2, nt-1
  do k=kn, k1, ki

    trid_a(:) = dcoef(1,:)*dumdzdfh(:,k,n)
    trid_b(:) = dcoef(2,:)*dumdzdfh(:,k,n) + tempz1(k)
    trid_c(:) = dcoef(3,:)*dumdzdfh(:,k,n)

    do j=2, ny-1
      trid_d(j) = -sum( dcoef(:,j)*fn_str(j-1:j+1,k-ki,n) )* &
                  dumdzdfh(j,k-k1,n)
    enddo
    trid_d(:) = trid_d(:) + tempz2(k)*fn_str(:,k-ki,n) + &
                cosphi(:)*(temp3d(:,k-ki,n)+temp3d(:,k,n))

    ! boundary condition
    trid_a(2   ) = 0.
    trid_c(ny-1) = 0.

    ! solve
    call tridag(j1-1,trid_a(2:j1),trid_b(2:j1),trid_c(2:j1), &
                trid_d(2:j1))
    fn_str(2:j1,k,n) = trid_d(2:j1)

    call tridag(ny-j2,trid_a(j2:ny-1),trid_b(j2:ny-1),trid_c(j2:ny-1), &
                trid_d(j2:ny-1))
    fn_str(j2:ny-1,k,n) = trid_d(j2:ny-1)

  enddo
  enddo

  call avg_d3((/j1,nz,nt-2,1/),fn_str(1:j1,:,2:nt-1),1., temp2d(1:j1,:))

  call avg_d3((/nyo-j1,nz,nt-2,1/),fn_str(j2:ny,:,2:nt-1),1., &
              temp2d(j1+1:nyo,:))

  call deriv1d((/j1,nz,1,1/),1,uniform_y,rsinphi(1:j1), &
               temp2d(1:j1,:),l_ybdy,0., wmr_d(1:j1,:,m))

  call deriv1d((/nyo-j1,nz,1,1/),1,uniform_y,rsinphi(j2:ny), &
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


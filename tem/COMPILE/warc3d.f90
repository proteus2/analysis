MODULE warc3d

  implicit none

  private ::  set_gridvar3d_p
  private ::  gradx_2nd, grady_2nd, gradz_2nd_irr
  private ::  missing_bdy

  character(len=64), dimension(11) ::  varname_warc_p =                  &
      (/'v_res  ','w_res  ','cor    ','uadv_y ','uadv_z ',               &
        'epd    ','epd_y  ','epd_z  ','f_y    ','f_z    ',               &
        'u_force'/)
  character(len=64), dimension(10) ::  varname_warc_qg =                 &
      (/'v_res  ','w_res  ','cor    ','epd    ','epd_y  ',               &
        'epd_z  ','f_y    ','f_z    ','u_force','tadv_z '/)

  integer,                           private ::  ny_pre, nz_pre
  real, dimension(:),   allocatable, private ::  lat_pre, ht_pre, zp
  real, dimension(:,:), allocatable, private ::  coslat, f, rho0

  real, parameter, private ::  g = 9.80665, rd = 287.05, cp = 1005.
  real, parameter, private ::  kappa = rd/cp
  real, parameter, private ::  a_earth = 6371229.
  real, parameter, private ::  pi = 3.14159265358979323846
  real, parameter, private ::  ome_earth = 7.292116e-5
  real, parameter, private ::  deg2rad = pi/180.


  CONTAINS


SUBROUTINE warc_hydro_p(                                                 &
     nx,ny,nz,lat,p,u,v,w,pt,h_scale,missv,                              &
     vres,wres,cor,uadvy,uadvz,epd,epd_y,epd_z,fy,fz,u_force )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  u, v, w, pt

  real, dimension(ny,nz), intent(out) ::  vres, wres, fy, fz
  real, dimension(ny,nz), intent(out) ::  cor, uadvy, uadvz,             &
                                          epd, epd_y, epd_z
  real, dimension(ny,nz), intent(out) ::  u_force

  real, dimension(ny,nz) ::  um, ptm
  real, dimension(ny,nz) ::  phi, rvpt, rwpt, rvu, rwu
  real, dimension(ny,nz) ::  dptmdy, dptmdz, divy_um, dumdz
  real, dimension(ny,nz) ::  temp

  real, dimension(:,:,:), allocatable ::  prt, v_prt, w_prt

  integer ::  j,k

  call set_gridvar3d_p(ny,nz,lat,p,h_scale)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), v_prt(nx,ny,nz), w_prt(nx,ny,nz) )

  um  (:,:) = sum(u , dim=1)/float(nx)
  vres(:,:) = sum(v , dim=1)/float(nx)
  wres(:,:) = sum(w , dim=1)/float(nx)
  ptm (:,:) = sum(pt, dim=1)/float(nx)

  v_prt(:,:,:) = v(:,:,:) - spread(vres(:,:),1,nx)
  w_prt(:,:,:) = w(:,:,:) - spread(wres(:,:),1,nx)

  prt(:,:,:) = pt(:,:,:) - spread(ptm(:,:),1,nx)

  rvpt(:,:) = sum(v_prt*prt, dim=1)/float(nx)
  rwpt(:,:) = sum(w_prt*prt, dim=1)/float(nx)

  prt(:,:,:) = u(:,:,:) - spread(um(:,:),1,nx)

  rvu(:,:) = sum(v_prt*prt, dim=1)/float(nx)
  rwu(:,:) = sum(w_prt*prt, dim=1)/float(nx)

  deallocate( prt, v_prt, w_prt )

  rvpt(:,:) = rho0(:,:)*rvpt(:,:)
  rwpt(:,:) = rho0(:,:)*rwpt(:,:)
  rvu (:,:) = rho0(:,:)*rvu (:,:)
  rwu (:,:) = rho0(:,:)*rwu (:,:)

  ! grad, pt
!  call grady_2nd(ny,nz,ptm,lat, dptmdy)

  temp(:,:) = log(ptm(:,:))
  call gradz_2nd_irr(1,ny,nz,temp,zp, dptmdz)
  dptmdz(:,:) = ptm(:,:)*dptmdz(:,:)

  ! grad, u
  call grady_2nd(1,ny,nz,um*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(1,ny,nz,um,zp, dumdz)

  ! phi
  phi(:,:) = rvpt(:,:)/dptmdz(:,:)
!!  call missing_bdy0(ny,nz,phi,missv,0,0,1,0)

  ! residual mean meridional circulation
  call gradz_2nd_irr(1,ny,nz,phi,zp, temp)
  vres(:,:) = vres(:,:) - temp(:,:)/rho0(:,:)
  call grady_2nd(1,ny,nz,phi*coslat,lat, temp)
  do k=1, nz
  do j=2, ny-1
    wres(j,k) = wres(j,k) + temp(j,k)/rho0(j,k)/coslat(j,k)
  enddo
  enddo
  call missing_bdy0(ny,nz,vres,missv,0,0,1,0)
  call missing_bdy0(ny,nz,wres,missv,1,1,1,0)

  ! adv
  uadvy(:,:) = -vres(:,:)*divy_um(:,:) * 86400.
  uadvz(:,:) = -wres(:,:)*dumdz  (:,:) * 86400.
  call missing_bdy0(ny,nz,uadvy,missv,1,1,1,0)
  call missing_bdy0(ny,nz,uadvz,missv,1,1,1,0)

  ! coriolis term
  cor(:,:) = vres(:,:)*f(:,:) * 86400.
  call missing_bdy0(ny,nz,cor,missv,0,0,1,0)

  ! epf, epd
  fy(:,:) = -a_earth*coslat(:,:)*(rvu(:,:) - phi(:,:)*dumdz(:,:))
  fz(:,:) = -a_earth*coslat(:,:)*(rwu(:,:) + phi(:,:)*(divy_um(:,:)-f(:,:)))
  call missing_bdy0(ny,nz,fy,missv,0,0,1,0)
!mv  call missing_bdy0(ny,nz,fz,missv,1,1,1,0)

  call grady_2nd(1,ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) * 86400.
  enddo
  enddo
  call missing_bdy0(ny,nz,epd_y,missv,1,1,1,0)

  call gradz_2nd_irr(1,ny,nz,fz,zp, epd_z)
  epd_z(:,:) = epd_z(:,:)/(rho0(:,:)*a_earth*coslat(:,:)) * 86400.
  call missing_bdy0(ny,nz,epd_z,missv,1,1,1,0)

  call missing_bdy0(ny,nz,fz,missv,1,1,1,0)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy0(ny,nz,epd,missv,1,1,1,0)

  ! u_force
  u_force(:,:) = cor(:,:) + uadvy(:,:) + uadvz(:,:) + epd(:,:)
  call missing_bdy0(ny,nz,u_force,missv,1,1,1,0)

END subroutine warc_hydro_p

!::::  Decomposition of fields  ::::::::::::::::::::::::::::::::::::::::
!  U = U_stationary + U_transient  : time mean + departure
!  U_stationary = U0 (y,z) + [U] (x,y,z)
!  U_transient = U' = U_bar' (y,z,t) + U'' (x,y,z,t)
!  ( zonal-mean U = U0 + U_bar' )
!  Here, the stationary wave perturbation is defined as [U],
!  and the transient wave perturbation is as U'.
!  Note that the former does not have a zonally symmetric component,
!  while the latter has (i.e., U_bar').
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

SUBROUTINE warc_s_qg(                                                    &
     nx,ny,nz,lat,p,u,v,pt,gp,dlon,h_scale,missv,                        &
     u0,epd,epd_x,epd_y,epd_z,fx,fy,fz )
! Plumb (1985, JAS, Eq. 7.1)

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  dlon, h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  u, v, pt, gp

  real, dimension(ny,nz),    intent(out) ::  u0
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  epd, epd_x, epd_y, epd_z

  real, dimension(ny,nz) ::  pt0, gp0
  real, dimension(ny,nz) ::  dpt0dz
  real, dimension(ny,nz) ::  temp

  real, dimension(:,:,:), allocatable ::  prt, u_prt, v_prt

  integer ::  j,k

  call set_gridvar3d_p(ny,nz,lat,p,h_scale)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), u_prt(nx,ny,nz), v_prt(nx,ny,nz) )

  u0  (:,:) = sum(u , dim=1)/float(nx)
  temp(:,:) = sum(v , dim=1)/float(nx)
  pt0 (:,:) = sum(pt, dim=1)/float(nx)
  gp0 (:,:) = sum(gp, dim=1)/float(nx)

  u_prt(:,:,:) = u(:,:,:) - spread(u0  (:,:),1,nx)
  v_prt(:,:,:) = v(:,:,:) - spread(temp(:,:),1,nx)

  temp(:,:) = log(pt0(:,:))
  call gradz_2nd_irr(1,ny,nz,temp,zp, dpt0dz)
  dpt0dz(:,:) = pt0(:,:)*dpt0dz(:,:)

  prt(:,:,:) = (gp(:,:,:) - spread(gp0(:,:),1,nx))/spread(f(:,:),1,nx)
  
  call gradx_2nd(nx,ny,nz,v_prt*prt,lat,dlon, fx)
  fx(:,:,:) = v_prt(:,:,:)*v_prt(:,:,:) - 0.5*fx(:,:,:)

  call gradx_2nd(nx,ny,nz,u_prt*prt,lat,dlon, fy)
  fy(:,:,:) = -(v_prt(:,:,:)*u_prt(:,:,:)) + 0.5*fy(:,:,:)

  u_prt(:,:,:) = prt(:,:,:)  ! save [gp]/f in u_prt, temporally

  prt(:,:,:) = pt(:,:,:) - spread(pt0(:,:),1,nx)

  call gradx_2nd(nx,ny,nz,prt*u_prt,lat,dlon, fz)
  fz(:,:,:) = (v_prt(:,:,:)*prt(:,:,:)) - 0.5*fz(:,:,:)
  fz(:,:,:) = fz(:,:,:)*spread(f(:,:)/dpt0dz(:,:),1,nx)

  deallocate( u_prt, v_prt )

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, epd_x)
  epd_x(:,:,:) = epd_x(:,:,:) * 86400.
!  call missing_bdy(nx,ny,nz,epd_x,missv,0,0,0,0)

  fx(:,:,:) = fx(:,:,:)*spread(rho0(:,:)*coslat(:,:),1,nx)
!  call missing_bdy(nx,ny,nz,fx,missv,0,0,0,0)

  prt(:,:,:) = spread(coslat(:,:),1,nx)

  fy(:,:,:) = fy(:,:,:)*prt(:,:,:)

  call grady_2nd(nx,ny,nz,fy*prt,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(:,j,k) = epd_y(:,j,k)/(prt(:,j,k)*prt(:,j,k)) * 86400.
  enddo
  enddo
  call missing_bdy(nx,ny,nz,epd_y,missv,1,1,0,0)

  prt(:,:,:) = spread(rho0(:,:),1,nx)

  fy(:,:,:) = fy(:,:,:)*prt(:,:,:)
!  call missing_bdy(nx,ny,nz,fy,missv,0,0,0,0)

  fz(:,:,:) = fz(:,:,:)*prt(:,:,:)

  call gradz_2nd_irr(nx,ny,nz,fz,zp, epd_z)
  epd_z(:,:,:) = epd_z(:,:,:)/prt(:,:,:) * 86400.
  call missing_bdy(nx,ny,nz,epd_z,missv,0,0,1,1)

  fz(:,:,:) = fz(:,:,:)*spread(coslat(:,:),1,nx)
  call missing_bdy(nx,ny,nz,fz,missv,0,0,0,0)

  deallocate( prt )

  epd(:,:,:) = epd_x(:,:,:) + epd_y(:,:,:) + epd_z(:,:,:)
  call missing_bdy(nx,ny,nz,epd,missv,1,1,1,1)

END subroutine warc_s_qg

SUBROUTINE warc_s_qg_gp(                                                 &
     nx,ny,nz,lat,p,gp,pt0,dlon,h_scale,missv,                           &
     u0,epd,epd_x,epd_y,epd_z,fx,fy,fz )
! Plumb (1985, JAS, Eq. 5.7)

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  dlon, h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(ny,nz),    intent(in) ::  pt0
  real, dimension(nx,ny,nz), intent(in) ::  gp

  real, dimension(ny,nz),    intent(out) ::  u0
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  epd, epd_x, epd_y, epd_z

  real, dimension(ny,nz) ::  gp0
  real, dimension(ny,nz) ::  dpt0dz
  real, dimension(ny,nz) ::  temp

  real, dimension(:,:,:), allocatable ::  prt, prt_x, prt_a

  integer ::  j,k

  call set_gridvar3d_p(ny,nz,lat,p,h_scale)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), prt_x(nx,ny,nz), prt_a(nx,ny,nz) )

  gp0(:,:) = sum(gp, dim=1)/float(nx)

  call grady_2nd(1,ny,nz,-gp0/f,lat, u0)

  temp(:,:) = log(pt0(:,:))
  call gradz_2nd_irr(1,ny,nz,temp,zp, dpt0dz)
  dpt0dz(:,:) = pt0(:,:)*dpt0dz(:,:)

  prt(:,:,:) = (gp(:,:,:) - spread(gp0(:,:),1,nx))/spread(f(:,:),1,nx)

  call gradx_2nd(nx,ny,nz,prt,lat,dlon, prt_x)
  call gradx_2nd(nx,ny,nz,prt_x,lat,dlon, fx)
  fx(:,:,:) = 0.5*(prt_x(:,:,:)*prt_x(:,:,:) - prt(:,:,:)*fx(:,:,:))

  call grady_2nd(nx,ny,nz,prt,lat, prt_a)
  call gradx_2nd(nx,ny,nz,prt_a,lat,dlon, fy)
!  call grady_2nd(nx,ny,nz,prt_x,lat, fy)
  fy(:,:,:) = 0.5*(prt_x(:,:,:)*prt_a(:,:,:) - prt(:,:,:)*fy(:,:,:))

  call gradz_2nd_irr(nx,ny,nz,prt,zp, prt_a)
  call gradx_2nd(nx,ny,nz,prt_a,lat,dlon, fz)
!  call gradz_2nd_irr(nx,ny,nz,prt_x,zp, fz)
  fz(:,:,:) = 0.5*(prt_x(:,:,:)*prt_a(:,:,:) - prt(:,:,:)*fz(:,:,:))

  temp(:,:) = f(:,:)*f(:,:)/dpt0dz(:,:)*                                 &
              spread((h_scale/rd)*(1.e5/p(:))**kappa,1,ny)
  fz(:,:,:) = fz(:,:,:)*spread(temp(:,:),1,nx)

  deallocate( prt_x, prt_a )

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, epd_x)
  epd_x(:,:,:) = epd_x(:,:,:) * 86400.
!  call missing_bdy(nx,ny,nz,epd_x,missv,0,0,0,0)

  fx(:,:,:) = fx(:,:,:)*spread(rho0(:,:)*coslat(:,:),1,nx)
!  call missing_bdy(nx,ny,nz,fx,missv,0,0,0,0)

  prt(:,:,:) = spread(coslat(:,:),1,nx)

  fy(:,:,:) = fy(:,:,:)*prt(:,:,:)

  call grady_2nd(nx,ny,nz,fy*prt,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(:,j,k) = epd_y(:,j,k)/(prt(:,j,k)*prt(:,j,k)) * 86400.
  enddo
  enddo
  call missing_bdy(nx,ny,nz,epd_y,missv,2,2,0,0)

  prt(:,:,:) = spread(rho0(:,:),1,nx)

  fy(:,:,:) = fy(:,:,:)*prt(:,:,:)
  call missing_bdy(nx,ny,nz,fy,missv,1,1,0,0)

  fz(:,:,:) = fz(:,:,:)*prt(:,:,:)

  call gradz_2nd_irr(nx,ny,nz,fz,zp, epd_z)
  epd_z(:,:,:) = epd_z(:,:,:)/prt(:,:,:) * 86400.
  call missing_bdy(nx,ny,nz,epd_z,missv,0,0,1,1)

  fz(:,:,:) = fz(:,:,:)*spread(coslat(:,:),1,nx)
  call missing_bdy(nx,ny,nz,fz,missv,0,0,0,0)

  deallocate( prt )

  epd(:,:,:) = epd_x(:,:,:) + epd_y(:,:,:) + epd_z(:,:,:)
  call missing_bdy(nx,ny,nz,epd,missv,1,1,1,1)

END subroutine warc_s_qg_gp

SUBROUTINE set_gridvar3d_p(ny,nz,lat,p,h_scale)

  integer,             intent(in) ::  ny, nz
  real,                intent(in) ::  h_scale
  real, dimension(ny), intent(in) ::  lat
  real, dimension(nz), intent(in) ::  p

  if ( allocated(lat_pre) ) then
    if ( ny == ny_pre .and. nz == nz_pre ) then
      if ( all(lat == lat_pre) .and. all(p == ht_pre) )  RETURN
    end if
    deallocate( lat_pre, ht_pre, zp, coslat, f, rho0 )
  end if

  allocate( zp(nz), coslat(ny,nz), f(ny,nz), rho0(ny,nz) )
  allocate( lat_pre(ny), ht_pre(nz) )

  zp(:) = -h_scale*(log(p(:)/1.e5))

  coslat(:,:) = spread(cos(lat(:)*deg2rad),2,nz)
  if (abs(lat(1 )) == 90.)  coslat(1 ,:) = 0.
  if (abs(lat(ny)) == 90.)  coslat(ny,:) = 0.
  f   (:,:) = spread(2.*ome_earth*sin(lat(:)*deg2rad),2,nz)
  rho0(:,:) = spread(p(:)/g/h_scale                  ,1,ny)

  ny_pre = ny  ;  nz_pre = nz
  lat_pre(:) = lat(:)
  ht_pre (:) = p  (:)

END subroutine set_gridvar3d_p

SUBROUTINE gradx_2nd(nx,ny,nz,var,lat,dlon, gradx)

  integer,                   intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real,                      intent(in)  ::  dlon
  real, dimension(nx,ny,nz), intent(out) ::  gradx

  integer ::  j,k
  real    ::  inv_2dx(ny)

  inv_2dx(:) = 0.5/(dlon*deg2rad*a_earth)/cos(lat(:)*deg2rad)
  if (abs(lat(1 )) == 90.)  inv_2dx(1 ) = 0.
  if (abs(lat(ny)) == 90.)  inv_2dx(ny) = 0.

  gradx(2:nx-1,:,:) = var(3:nx,:,:) - var(1:nx-2,:,:)
  gradx(1 ,:,:) = var(2,:,:) - var(nx  ,:,:)
  gradx(nx,:,:) = var(1,:,:) - var(nx-1,:,:)

  do k=1, nz
  do j=1, ny
    gradx(:,j,k) = gradx(:,j,k)*inv_2dx(j)
  enddo
  enddo

END subroutine gradx_2nd

SUBROUTINE grady_2nd(nx,ny,nz,var,lat, grady)

  integer,                   intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real, dimension(nx,ny,nz), intent(out) ::  grady

  integer ::  j,k
  real    ::  inv_2dy(ny)

  inv_2dy(2:ny-1) = 1./((lat(3:ny)-lat(1:ny-2))*deg2rad*a_earth)

  do k=1, nz
  do j=2, ny-1
    grady(:,j,k) = (var(:,j+1,k)-var(:,j-1,k))*inv_2dy(j)
  enddo
  enddo
  grady(:,1 ,:) = 0.
  grady(:,ny,:) = 0.

END subroutine grady_2nd

SUBROUTINE gradz_2nd_irr(nx,ny,nz,var,z, gradz)

  integer,                   intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(nz),       intent(in)  ::  z
  real, dimension(nx,ny,nz), intent(out) ::  gradz

  integer ::  j,k
  real    ::  coef(3,nz), inv_dz_1, inv_dz_n

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coef(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coef(:,k))
  enddo
  inv_dz_1 = 1./(z(2 ) - z(1   ))
  inv_dz_n = 1./(z(nz) - z(nz-1))

  do k=2, nz-1
    gradz(:,:,k) = coef(1,k)*var(:,:,k-1) +  &
                   coef(2,k)*var(:,:,k  ) +  &
                   coef(3,k)*var(:,:,k+1)
  enddo

  gradz(:,:,1 ) = (var(:,:,2 ) - var(:,:,1   ))*inv_dz_1
  gradz(:,:,nz) = (var(:,:,nz) - var(:,:,nz-1))*inv_dz_n

END subroutine gradz_2nd_irr

SUBROUTINE missing_bdy(nx,ny,nz,var,missv,nm_y1,nm_y2,nm_z1,nm_z2)

  integer, intent(in) ::  nx, ny, nz, nm_y1, nm_y2, nm_z1, nm_z2
  real,    intent(in) ::  missv

  real, dimension(nx,ny,nz), intent(inout) ::  var

  if (nm_y1 > 0)  var(:,1         :nm_y1,          :     ) = missv
  if (nm_y2 > 0)  var(:,ny+1-nm_y2:ny   ,          :     ) = missv
  if (nm_z1 > 0)  var(:,          :     ,1         :nm_z1) = missv
  if (nm_z2 > 0)  var(:,          :     ,nz+1-nm_z2:nz   ) = missv

END subroutine missing_bdy

SUBROUTINE missing_bdy0(ny,nz,var,missv,nm_y1,nm_y2,nm_z1,nm_z2)

  integer, intent(in) ::  ny, nz, nm_y1, nm_y2, nm_z1, nm_z2
  real,    intent(in) ::  missv

  real, dimension(ny,nz), intent(inout) ::  var

  if (nm_y1 > 0)  var(1         :nm_y1,          :     ) = missv
  if (nm_y2 > 0)  var(ny+1-nm_y2:ny   ,          :     ) = missv
  if (nm_z1 > 0)  var(          :     ,1         :nm_z1) = missv
  if (nm_z2 > 0)  var(          :     ,nz+1-nm_z2:nz   ) = missv

END subroutine missing_bdy0

END module warc3d


MODULE tem3d

  implicit none

  private ::  set_gridvar3d_p, div_flx_unwgt, flx_wgt
  private ::  gradx_2nd, grady_2nd, gradz_2nd_irr, gradxx_2nd
  private ::  missing_bdy

  character(len=64), dimension(11) ::  varname_tem3d_p =                 &
      (/'v_res  ','w_res  ','cor    ','uadv_y ','uadv_z ',               &
        'epd    ','epd_y  ','epd_z  ','f_y    ','f_z    ',               &
        'u_force'/)
!  character(len=64), dimension(10) ::  varname_tem3d_qg =                &
!      (/'v_res  ','w_res  ','cor    ','epd    ','epd_y  ',               &
!        'epd_z  ','f_y    ','f_z    ','u_force','tadv_z '/)
  character(len=64), dimension(8) ::  varname_tem3d_qg =                 &
      (/'epd    ','epd_x  ','epd_y  ','epd_z  ','f_x    ',               &
        'f_y    ','f_z    ','U0     '/)

  integer,                           private ::  ny_pre, nz_pre
  integer,                           private ::  j5s, j5e
  real, dimension(:),   allocatable, private ::  lat_pre, ht_pre, zp
  real, dimension(:,:), allocatable, private ::  coslat, f, rho0

  real, parameter, private ::  g = 9.80665, rd = 287.05, cp = 1005.
  real, parameter, private ::  kappa = rd/cp
  real, parameter, private ::  a_earth = 6371229.
  real, parameter, private ::  pi = 3.14159265358979323846
  real, parameter, private ::  ome_earth = 7.292116e-5
  real, parameter, private ::  deg2rad = pi/180.

  interface grady_2nd
    module procedure grady_2nd_3d, grady_2nd_2d
  end interface
  interface gradz_2nd_irr
    module procedure gradz_2nd_irr_3d, gradz_2nd_irr_2d
  end interface
  interface missing_bdy
    module procedure missing_bdy_3d, missing_bdy_2d
  end interface


  CONTAINS


!SUBROUTINE tem3d_hydro_p(                                                &
!END subroutine tem3d_hydro_p

!::::  Decomposition of fields  ::::::::::::::::::::::::::::::::::::::::
!  U = U_stationary + U_transient  : time mean + departure
!  U_stationary = U0 (y,z) + [U] (x,y,z)
!  U_transient = U' = U_bar' (y,z,t) + U'' (x,y,z,t)
!  ( zonal-mean U = U0 + U_bar' )
!  Here, the stationary wave perturbation is defined as [U],
!!!  and the transient wave perturbation is as U'.
!!  Note that the former does not have a zonally symmetric component,
!!  while the latter has (i.e., U_bar').
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

SUBROUTINE tem3d_s_qg(                                                   &
     nx,ny,nz,lat,p,u,v,t,gp,dlon,h_scale,missv,                         &
     divf,divf_x,divf_y,divf_z,fx,fy,fz,u0 )
! Plumb (1985, JAS, Eq. 7.1)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon, h_scale, missv
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nz)      , intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  u, v, t, gp

  real, dimension(ny,nz)   , intent(out) ::  u0
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  divf, divf_x, divf_y, divf_z

  real, dimension(ny,nz) ::  t0, gp0
  real, dimension(ny,nz) ::  s_t0

  real, dimension(:,:,:), allocatable ::  prt, u_prt, v_prt

  integer ::  j,k

  call set_gridvar3d_p(ny,nz,lat,p,h_scale)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), u_prt(nx,ny,nz), v_prt(nx,ny,nz) )

  u0  (:,:) = sum(u , dim=1)/float(nx)
  t0  (:,:) = sum(t , dim=1)/float(nx)
  s_t0(:,:) = sum(v , dim=1)/float(nx)   ! rent s_t0
  gp0 (:,:) = sum(gp, dim=1)/float(nx)

  u_prt(:,:,:) = u(:,:,:) - spread(u0  (:,:),1,nx)
  v_prt(:,:,:) = v(:,:,:) - spread(s_t0(:,:),1,nx)

  call gradz_2nd_irr(ny,nz,t0,zp, s_t0)
  s_t0(:,:) = s_t0(:,:) + t0(:,:)*(kappa/h_scale)

  prt(:,:,:) = (gp(:,:,:) - spread(gp0(:,:),1,nx))/spread(f(:,:),1,nx)
  if ( minval(abs(lat)) == 0. ) then
    j = sum(minloc(abs(lat)))
    prt(:,j,:) = 0.5*(prt(:,j-1,:) + prt(:,j+1,:))
  end if
  
  call gradx_2nd(nx,ny,nz,v_prt*prt,lat,dlon, fx)
  fx(:,:,:) = v_prt(:,:,:)*v_prt(:,:,:) - 0.5*fx(:,:,:)

  call gradx_2nd(nx,ny,nz,u_prt*prt,lat,dlon, fy)
  fy(:,:,:) = -(v_prt(:,:,:)*u_prt(:,:,:)) + 0.5*fy(:,:,:)

  u_prt(:,:,:) = t(:,:,:) - spread(t0(:,:),1,nx)   ! rent u_prt

  call gradx_2nd(nx,ny,nz,u_prt*prt,lat,dlon, fz)
  fz(:,:,:) = v_prt(:,:,:)*u_prt(:,:,:) - 0.5*fz(:,:,:)
  fz(:,:,:) = fz(:,:,:)*spread(f(:,:)/s_t0(:,:),1,nx)

  deallocate( prt, u_prt, v_prt )

  fy(:,:,:) = fy(:,:,:)*spread(coslat(:,:),1,nx)
  fz(:,:,:) = fz(:,:,:)*spread(rho0  (:,:),1,nx)

  ! F not fully weighted :  fx = fx0  /  fy = cos*fy0  /  fz = rho*fz0

  call div_flx_unwgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,                   &
                     divf_x,divf_y,divf_z,divf)
  ! div = div[F]/(rho*cos)  [m/s/day]

  call flx_wgt(ny,nz,fx,fy,fz)
  ! F = (rho*cos)*F0

  call missing_bdy(nx,ny,nz,fz    ,missv,0,0,1,1)  ! due to s_t0
  call missing_bdy(nx,ny,nz,divf_z,missv,0,0,2,2)
  call missing_bdy(nx,ny,nz,divf  ,missv,0,0,2,2)

  if (j5s /= 0) then
    divf  (:,j5s:j5e,:) = missv
    divf_x(:,j5s:j5e,:) = missv
    divf_y(:,j5s:j5e,:) = missv
    divf_z(:,j5s:j5e,:) = missv
    fx    (:,j5s:j5e,:) = missv
    fy    (:,j5s:j5e,:) = missv
    fz    (:,j5s:j5e,:) = missv
  end if

END subroutine tem3d_s_qg

SUBROUTINE tem3d_s_qg_gp(                                                &
     nx,ny,nz,lat,p,gp,t,dlon,h_scale,missv,                             &
     divf,divf_x,divf_y,divf_z,fx,fy,fz,u0 )
! Plumb (1985, JAS, Eq. 5.7)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon, h_scale, missv
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nz)      , intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  gp, t

  real, dimension(ny,nz)   , intent(out) ::  u0
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  divf, divf_x, divf_y, divf_z

  real, dimension(ny,nz) ::  gp0, t0
  real, dimension(ny,nz) ::  s_t0

  real, dimension(:,:,:), allocatable ::  prt, prt_x, prt_a

  integer ::  j,k

  call set_gridvar3d_p(ny,nz,lat,p,h_scale)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), prt_x(nx,ny,nz), prt_a(nx,ny,nz) )

  gp0(:,:) = sum(gp, dim=1)/float(nx)
  t0 (:,:) = sum(t , dim=1)/float(nx)

  call grady_2nd(ny,nz,-gp0/f,lat, u0)
  if ( minval(abs(lat)) == 0. ) then
    j = sum(minloc(abs(lat)))
    u0(j,:) = 0.5*(u0(j-1,:) + u0(j+1,:))
  end if

  call gradz_2nd_irr(ny,nz,t0,zp, s_t0)
  s_t0(:,:) = s_t0(:,:) + t0(:,:)*(kappa/h_scale)

  prt(:,:,:) = (gp(:,:,:) - spread(gp0(:,:),1,nx))/spread(f(:,:),1,nx)
  if ( minval(abs(lat)) == 0. ) then
    j = sum(minloc(abs(lat)))
    prt(:,j,:) = 0.5*(prt(:,j-1,:) + prt(:,j+1,:))
  end if

  call gradx_2nd(nx,ny,nz,prt,lat,dlon, prt_x)
!  call gradx_2nd(nx,ny,nz,prt_x,lat,dlon, fx)
  call gradxx_2nd(nx,ny,nz,prt,lat,dlon, fx)
  fx(:,:,:) = 0.5*(prt_x(:,:,:)*prt_x(:,:,:) - prt(:,:,:)*fx(:,:,:))

  call grady_2nd(nx,ny,nz,prt,lat, prt_a)
  call grady_2nd(nx,ny,nz,prt_x,lat, fy)
  fy(:,:,:) = 0.5*(prt_x(:,:,:)*prt_a(:,:,:) - prt(:,:,:)*fy(:,:,:))

  prt_a(:,:,:) = t(:,:,:) - spread(t0(:,:),1,nx)
  call gradx_2nd(nx,ny,nz,prt_a,lat,dlon, fz)
  fz(:,:,:) = 0.5*(prt_x(:,:,:)*prt_a(:,:,:) - prt(:,:,:)*fz(:,:,:))
  fz(:,:,:) = fz(:,:,:)*spread(f(:,:)/s_t0(:,:),1,nx)

  deallocate( prt, prt_x, prt_a )

  fy(:,:,:) = fy(:,:,:)*spread(coslat(:,:),1,nx)
  fz(:,:,:) = fz(:,:,:)*spread(rho0  (:,:),1,nx)

  ! F not fully weighted :  fx = fx0  /  fy = cos*fy0  /  fz = rho*fz0

  call div_flx_unwgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,                   &
                     divf_x,divf_y,divf_z,divf)
  ! div = div[F]/(rho*cos)  [m/s/day]

  call flx_wgt(ny,nz,fx,fy,fz)
  ! F = (rho*cos)*F0

  call missing_bdy(nx,ny,nz,fy    ,missv,1,1,0,0)
  call missing_bdy(nx,ny,nz,divf_y,missv,2,2,0,0)
  call missing_bdy(nx,ny,nz,fz    ,missv,0,0,1,1)  ! due to s_t0
  call missing_bdy(nx,ny,nz,divf_z,missv,0,0,2,2)
  call missing_bdy(nx,ny,nz,divf  ,missv,2,2,2,2)

  if (j5s /= 0) then
    divf  (:,j5s:j5e,:) = missv
    divf_x(:,j5s:j5e,:) = missv
    divf_y(:,j5s:j5e,:) = missv
    divf_z(:,j5s:j5e,:) = missv
    fx    (:,j5s:j5e,:) = missv
    fy    (:,j5s:j5e,:) = missv
    fz    (:,j5s:j5e,:) = missv
    u0(j5s:j5e,:) = missv
  end if

END subroutine tem3d_s_qg_gp

SUBROUTINE set_gridvar3d_p(ny,nz,lat,p,h_scale)

  integer,             intent(in) ::  ny, nz
  real,                intent(in) ::  h_scale
  real, dimension(ny), intent(in) ::  lat
  real, dimension(nz), intent(in) ::  p

  integer ::  j

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

  j5s = 0  ;  j5e = 0
  do j=1, ny
    if ( abs(lat(j)) < 5. ) then
      j5s = j  ;  EXIT
    end if
  enddo
  do j=ny, 1, -1
    if ( abs(lat(j)) < 5. ) then
      j5e = j  ;  EXIT
    end if
  enddo

  ny_pre = ny  ;  nz_pre = nz
  lat_pre(:) = lat(:)
  ht_pre (:) = p  (:)

END subroutine set_gridvar3d_p

SUBROUTINE div_flx_unwgt(nx,ny,nz,lat,dlon,missv,fx1,fyc,fzr,            &
                         div_x,div_y,div_z,div)

  integer               , intent(in) ::  nx, ny, nz
  real                  , intent(in) ::  dlon, missv
  real, dimension(:)    , intent(in) ::  lat
  real, dimension(:,:,:), intent(in) ::  fx1, fyc, fzr

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  call gradx_2nd(nx,ny,nz,fx1,lat,dlon, div_x)
  div_x(:,:,:) = div_x(:,:,:) * 86400.

  div_z(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div_z

  call grady_2nd(nx,ny,nz,fyc*div_z,lat, div_y)
  div_y(:,:,:) = div_y(:,:,:)/(div_z(:,:,:)*div_z(:,:,:)) * 86400.
  call missing_bdy(nx,ny,nz,div_y,missv,1,1,0,0)

  call gradz_2nd_irr(nx,ny,nz,fzr,zp, div_z)
  div_z(:,:,:) = div_z(:,:,:)/spread(rho0(:,:),1,nx) * 86400.
  call missing_bdy(nx,ny,nz,div_z,missv,0,0,1,1)

  div(:,:,:) = div_x(:,:,:) + div_y(:,:,:) + div_z(:,:,:)
  call missing_bdy(nx,ny,nz,div,missv,1,1,1,1)
 
END subroutine div_flx_unwgt

SUBROUTINE flx_wgt(ny,nz,fx,fy,fz)

  integer, intent(in) ::  ny, nz

  real, dimension(:,:,:), intent(inout) ::  fx, fy, fz

  integer ::  j,k

  do k=1, nz
  do j=1, ny
    fx(:,j,k) = fx(:,j,k)*(rho0(j,k)*coslat(j,k))
    fy(:,j,k) = fy(:,j,k)*rho0(j,k)
    fz(:,j,k) = fz(:,j,k)*coslat(j,k)
  enddo
  enddo

END subroutine flx_wgt

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

SUBROUTINE grady_2nd_3d(nx,ny,nz,var,lat, grady)

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

END subroutine grady_2nd_3d

SUBROUTINE gradz_2nd_irr_3d(nx,ny,nz,var,z, gradz)

  integer,                   intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(nz)      , intent(in)  ::  z
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

END subroutine gradz_2nd_irr_3d

SUBROUTINE grady_2nd_2d(ny,nz,var,lat, grady)

  integer,                intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(ny)   , intent(in)  ::  lat
  real, dimension(ny,nz), intent(out) ::  grady

  integer ::  j,k
  real    ::  inv_2dy(ny)

  inv_2dy(2:ny-1) = 1./((lat(3:ny)-lat(1:ny-2))*deg2rad*a_earth)

  do k=1, nz
  do j=2, ny-1
    grady(j,k) = (var(j+1,k)-var(j-1,k))*inv_2dy(j)
  enddo
  enddo
  grady(1 ,:) = 0.
  grady(ny,:) = 0.

END subroutine grady_2nd_2d

SUBROUTINE gradz_2nd_irr_2d(ny,nz,var,z, gradz)

  integer,                intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(nz)   , intent(in)  ::  z
  real, dimension(ny,nz), intent(out) ::  gradz

  integer ::  j,k
  real    ::  coef(3,nz), inv_dz_1, inv_dz_n

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coef(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coef(:,k))
  enddo
  inv_dz_1 = 1./(z(2 ) - z(1   ))
  inv_dz_n = 1./(z(nz) - z(nz-1))

  do k=2, nz-1
    gradz(:,k) = coef(1,k)*var(:,k-1) +  &
                 coef(2,k)*var(:,k  ) +  &
                 coef(3,k)*var(:,k+1)
  enddo

  gradz(:,1 ) = (var(:,2 ) - var(:,1   ))*inv_dz_1
  gradz(:,nz) = (var(:,nz) - var(:,nz-1))*inv_dz_n

END subroutine gradz_2nd_irr_2d

SUBROUTINE gradxx_2nd(nx,ny,nz,var,lat,dlon, gradxx)

  integer,                   intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real,                      intent(in)  ::  dlon
  real, dimension(nx,ny,nz), intent(out) ::  gradxx

  integer ::  j,k
  real    ::  inv_dx2(ny)

  inv_dx2(:) = 1.0/((dlon*deg2rad*a_earth)*cos(lat(:)*deg2rad))**2
  if (abs(lat(1 )) == 90.)  inv_dx2(1 ) = 0.
  if (abs(lat(ny)) == 90.)  inv_dx2(ny) = 0.

  gradxx(2:nx-1,:,:) = (var(1:nx-2,:,:) + var(3:nx,:,:)) -               &
                       2.*var(2:nx-1,:,:)
  gradxx(1 ,:,:) = (var(nx  ,:,:) + var(2,:,:)) - 2.*var(1 ,:,:)
  gradxx(nx,:,:) = (var(nx-1,:,:) + var(1,:,:)) - 2.*var(nx,:,:)

  do k=1, nz
  do j=1, ny
    gradxx(:,j,k) = gradxx(:,j,k)*inv_dx2(j)
  enddo
  enddo

END subroutine gradxx_2nd

SUBROUTINE missing_bdy_3d(nx,ny,nz,var,missv,nm_y1,nm_y2,nm_z1,nm_z2)

  integer, intent(in) ::  nx, ny, nz, nm_y1, nm_y2, nm_z1, nm_z2
  real,    intent(in) ::  missv

  real, dimension(nx,ny,nz), intent(inout) ::  var

  if (nm_y1 > 0)  var(:,1         :nm_y1,          :     ) = missv
  if (nm_y2 > 0)  var(:,ny+1-nm_y2:ny   ,          :     ) = missv
  if (nm_z1 > 0)  var(:,          :     ,1         :nm_z1) = missv
  if (nm_z2 > 0)  var(:,          :     ,nz+1-nm_z2:nz   ) = missv

END subroutine missing_bdy_3d

SUBROUTINE missing_bdy_2d(ny,nz,var,missv,nm_y1,nm_y2,nm_z1,nm_z2)

  integer, intent(in) ::  ny, nz, nm_y1, nm_y2, nm_z1, nm_z2
  real,    intent(in) ::  missv

  real, dimension(ny,nz), intent(inout) ::  var

  if (nm_y1 > 0)  var(1         :nm_y1,          :     ) = missv
  if (nm_y2 > 0)  var(ny+1-nm_y2:ny   ,          :     ) = missv
  if (nm_z1 > 0)  var(          :     ,1         :nm_z1) = missv
  if (nm_z2 > 0)  var(          :     ,nz+1-nm_z2:nz   ) = missv

END subroutine missing_bdy_2d

END module tem3d


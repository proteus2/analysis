MODULE tem3d

  use util,  only: gradx_2nd, grady_2nd, gradz_2nd_irr, gradxx_2nd,      &
                   get_ind_bnd, missing_bdy, gradxx_4nd
  use const_glob,  only:  g, rd, kappa, a_earth, ome_earth, deg2rad

  implicit none

  private ::  set_gridvar_p, set_gridvar_qg, div_flx_unwgt, div_flx_wgt
  private ::  g, rd, kappa, a_earth, ome_earth, deg2rad

  character(len=64), dimension(11) ::  varname_tem3d_p =                 &
      (/'v_res  ','w_res  ','cor    ','uadv_y ','uadv_z ',               &
        'epd    ','epd_y  ','epd_z  ','f_y    ','f_z    ',               &
        'u_force'/)
  character(len=64), dimension(11) ::  varname_waf3d_qg =                &
      (/'divF   ','divF_x ','divF_y ','divF_z ','f_x    ',               &
        'f_y    ','f_z    ','wad    ','q2     ','U0     ',               &
        'S      '/)
!        'beta_eff'/)

  integer,                           private ::  ny_pre, nz_pre
  integer,                           private ::  j5s, j5e
  integer,                           private ::  jns, jne, jss, jse
  logical,                           private ::  l_nh, l_sh
  real, dimension(:),   allocatable, private ::  lat_pre, ht_pre, zp
  real, dimension(:,:), allocatable, private ::  coslat, f, rho0, beta


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

SUBROUTINE waf3d_s_qg(                                                   &
     nx,ny,nz,lat,p,u,v,t,gp,dlon,h_scale,missv,                         &
     divf,divf_x,divf_y,divf_z,fx,fy,fz,wad,q2,u0,s_t0 )
! Plumb (1985, JAS, Eq. 7.1)

! E = (1/2)*[ u^2 + v^2 + (R/H/S)T^2 ]

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon, h_scale, missv
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nz)      , intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  u, v, t, gp

  real, dimension(ny,nz)   , intent(out) ::  u0, s_t0
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  divf, divf_x, divf_y, divf_z
  real, dimension(nx,ny,nz), intent(out) ::  wad, q2

  real, dimension(ny,nz)    ::  t0, gp0, v0
  real, dimension(nx,ny,nz) ::  prt, u_prt, v_prt, e_tot

  integer ::  j,k

  call set_gridvar_qg(ny,nz,lat,p,h_scale)

  u0 (:,:) = sum(u , dim=1)/float(nx)
  t0 (:,:) = sum(t , dim=1)/float(nx)
  v0 (:,:) = sum(v , dim=1)/float(nx)
  gp0(:,:) = sum(gp, dim=1)/float(nx)

  call gradz_2nd_irr(ny,nz,t0,zp, s_t0)
  s_t0(:,:) = s_t0(:,:) + t0(:,:)*(kappa/h_scale)

  u_prt(:,:,:) = u(:,:,:) - spread(u0(:,:),1,nx)
  v_prt(:,:,:) = v(:,:,:) - spread(v0(:,:),1,nx)

  prt(:,:,:) = (gp(:,:,:) - spread(gp0(:,:),1,nx))/spread(f(:,:),1,nx)
  if ( minval(abs(lat)) == 0. ) then
    j = sum(minloc(abs(lat)))
    prt(:,j,:) = 0.5*(prt(:,j-1,:) + prt(:,j+1,:))
  end if

  ! kinetic E * 2
  e_tot(:,:,:) = u_prt(:,:,:)*u_prt(:,:,:) + v_prt(:,:,:)*v_prt(:,:,:)

  call gradx_2nd(nx,ny,nz,v_prt*prt,lat,dlon, fx)
  fx(:,:,:) = v_prt(:,:,:)*v_prt(:,:,:) - 0.5*fx(:,:,:)

  call gradx_2nd(nx,ny,nz,u_prt*prt,lat,dlon, fy)
  fy(:,:,:) = -(v_prt(:,:,:)*u_prt(:,:,:)) + 0.5*fy(:,:,:)

  u_prt(:,:,:) = t(:,:,:) - spread(t0(:,:),1,nx)   ! rent u_prt

  ! total E
  e_tot(:,:,:) = 0.5*( e_tot(:,:,:) + (u_prt(:,:,:)*u_prt(:,:,:))*       &
                                   spread((rd/h_scale)/s_t0(:,:),1,nx) )

  call gradx_2nd(nx,ny,nz,u_prt*prt,lat,dlon, fz)
  fz(:,:,:) = ( v_prt(:,:,:)*u_prt(:,:,:) - 0.5*fz(:,:,:) )*             &
              spread(f(:,:)/s_t0(:,:),1,nx)

  ! F = F0 : not mass-weighted

  call div_flx_wgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,                     &
                   divf_x,divf_y,divf_z,divf)
  ! F = (rho*cos) * F0 : mass-weighted
  ! div = div[F]/(rho*cos)  [m/s/day]

  call missing_bdy(nx,ny,nz,fz    ,missv,0,0,1,1)  ! due to s_t0
  call missing_bdy(nx,ny,nz,divf_z,missv,0,0,2,2)
  call missing_bdy(nx,ny,nz,divf  ,missv,0,0,2,2)

  call missing_bdy(ny,nz,s_t0,missv,0,0,1,1)

  if (j5s /= 0) then
    divf  (:,j5s:j5e,:) = missv
    divf_x(:,j5s:j5e,:) = missv
    divf_y(:,j5s:j5e,:) = missv
    divf_z(:,j5s:j5e,:) = missv
    fx    (:,j5s:j5e,:) = missv
    fy    (:,j5s:j5e,:) = missv
    fz    (:,j5s:j5e,:) = missv
  end if

END subroutine waf3d_s_qg

SUBROUTINE waf3d_s_qg_gp(                                                &
     nx,ny,nz,lat,p,gp,t,dlon,h_scale,missv,                             &
     divf,divf_x,divf_y,divf_z,fx,fy,fz,wad,q2,u0,s_t0 )
! Plumb (1985, JAS, Eq. 5.7)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon, h_scale, missv
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nz)      , intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  gp, t

  real, dimension(ny,nz)   , intent(out) ::  u0, s_t0
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  divf, divf_x, divf_y, divf_z
  real, dimension(nx,ny,nz), intent(out) ::  wad, q2

  real, dimension(ny,nz)    ::  gp0, t0
  real, dimension(nx,ny,nz) ::  prt, u_prt, v_prt, e_tot

  integer ::  j,k

  call set_gridvar_qg(ny,nz,lat,p,h_scale)

  gp0(:,:) = sum(gp, dim=1)/float(nx)
  t0 (:,:) = sum(t , dim=1)/float(nx)

  call gradz_2nd_irr(ny,nz,t0,zp, s_t0)
  s_t0(:,:) = s_t0(:,:) + t0(:,:)*(kappa/h_scale)

  call grady_2nd(ny,nz,gp0,lat, u0)
  u0(:,:) = (-1.)*u0(:,:)/f(:,:)
  if ( minval(abs(lat)) == 0. ) then
    j = sum(minloc(abs(lat)))
    u0(j,:) = 0.5*(u0(j-1,:) + u0(j+1,:))
  end if

  prt(:,:,:) = (gp(:,:,:) - spread(gp0(:,:),1,nx))/spread(f(:,:),1,nx)
  if ( minval(abs(lat)) == 0. ) then
    j = sum(minloc(abs(lat)))
    prt(:,j,:) = 0.5*(prt(:,j-1,:) + prt(:,j+1,:))
  end if

  call gradx_2nd(nx,ny,nz,prt,lat,dlon, v_prt)
  call grady_2nd(nx,ny,nz,-prt,lat, u_prt)

  ! kinetic E * 2
  e_tot(:,:,:) = u_prt(:,:,:)*u_prt(:,:,:) + v_prt(:,:,:)*v_prt(:,:,:)

  call gradx_2nd(nx,ny,nz,v_prt,lat,dlon, fx)
!  call gradxx_2nd(nx,ny,nz,prt,lat,dlon, fx)
!  call gradxx_4nd(nx,ny,nz,prt,lat,dlon, fx)
  fx(:,:,:) = 0.5*( v_prt(:,:,:)*v_prt(:,:,:) - prt(:,:,:)*fx(:,:,:) )

  call grady_2nd(nx,ny,nz,v_prt,lat, fy(:,:,:))
  fy(:,:,:) = -0.5*( v_prt(:,:,:)*u_prt(:,:,:) + prt(:,:,:)*fy(:,:,:) )

  u_prt(:,:,:) = t(:,:,:) - spread(t0(:,:),1,nx)   ! rent u_prt

  ! total E
  e_tot(:,:,:) = 0.5*( e_tot(:,:,:) + (u_prt(:,:,:)*u_prt(:,:,:))*       &
                                   spread((rd/h_scale)/s_t0(:,:),1,nx) )

  call gradx_2nd(nx,ny,nz,u_prt,lat,dlon, fz)
  fz(:,:,:) = 0.5*( v_prt(:,:,:)*u_prt(:,:,:) - prt(:,:,:)*fz(:,:,:) )*  &
              spread(f(:,:)/s_t0(:,:),1,nx)

  ! F = F0 : not mass-weighted

  call div_flx_wgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,                     &
                   divf_x,divf_y,divf_z,divf)
  ! F = (rho*cos) * F0 : mass-weighted
  ! div = div[F]/(rho*cos)  [m/s/day]

  call missing_bdy(nx,ny,nz,fy    ,missv,1,1,0,0)
  call missing_bdy(nx,ny,nz,divf_y,missv,2,2,0,0)
  call missing_bdy(nx,ny,nz,fz    ,missv,0,0,1,1)  ! due to s_t0
  call missing_bdy(nx,ny,nz,divf_z,missv,0,0,2,2)
  call missing_bdy(nx,ny,nz,divf  ,missv,2,2,2,2)

  call missing_bdy(ny,nz,s_t0,missv,0,0,1,1)

  if (j5s /= 0) then
    divf  (:,j5s:j5e,:) = missv
    divf_x(:,j5s:j5e,:) = missv
    divf_y(:,j5s:j5e,:) = missv
    divf_z(:,j5s:j5e,:) = missv
    fx    (:,j5s:j5e,:) = missv
    fy    (:,j5s:j5e,:) = missv
    fz    (:,j5s:j5e,:) = missv
    u0   (j5s:j5e,:) = missv
  end if

END subroutine waf3d_s_qg_gp

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

SUBROUTINE set_gridvar_p(ny,nz,lat,p,h_scale)

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

  call get_ind_bnd(ny,lat,(/-5.,5./), j5s,j5e)

  ny_pre = ny  ;  nz_pre = nz
  lat_pre(:) = lat(:)
  ht_pre (:) = p  (:)

END subroutine set_gridvar_p

SUBROUTINE set_gridvar_qg(ny,nz,lat,p,h_scale)

  integer,             intent(in) ::  ny, nz
  real,                intent(in) ::  h_scale
  real, dimension(ny), intent(in) ::  lat
  real, dimension(nz), intent(in) ::  p

  if ( allocated(beta) ) then
    if ( ny == ny_pre .and. nz == nz_pre ) then
      if ( all(lat == lat_pre) .and. all(p == ht_pre) )  RETURN
    end if
    deallocate( beta )
  end if

  allocate( beta(ny,nz) )

  beta(:,:) = spread(2.*ome_earth/a_earth*cos(lat(:)*deg2rad),2,nz)

  call get_ind_bnd(ny,lat,(/20. ,90. /), jns,jne)  ! following Plumb (1985)
  call get_ind_bnd(ny,lat,(/-20.,-90./), jss,jse)

  l_nh = ( maxval(lat) > 20.  )
  l_sh = ( minval(lat) < -20. )

  call set_gridvar_p(ny,nz,lat,p,h_scale)

END subroutine set_gridvar_qg

SUBROUTINE div_flx_unwgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,               &
                         div_x,div_y,div_z,div)

  integer               , intent(in) ::  nx, ny, nz
  real                  , intent(in) ::  dlon, missv
  real, dimension(:)    , intent(in) ::  lat
  real, dimension(:,:,:), intent(in) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)
  div_x(:,:,:) = div_x(:,:,:) * 86400.

  div(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div

  call grady_2nd(nx,ny,nz,(fy*div)*div,lat, div_y)
  div_y(:,:,:) = div_y(:,:,:)/(div(:,:,:)*div(:,:,:)) * 86400.
  call missing_bdy(nx,ny,nz,div_y,missv,1,1,0,0)

  div(:,:,:) = spread(rho0(:,:),1,nx)   ! rent div

  call gradz_2nd_irr(nx,ny,nz,fz*div,zp, div_z)
  div_z(:,:,:) = div_z(:,:,:)/div(:,:,:) * 86400.
  call missing_bdy(nx,ny,nz,div_z,missv,0,0,1,1)

  div(:,:,:) = div_x(:,:,:) + div_y(:,:,:) + div_z(:,:,:)
  call missing_bdy(nx,ny,nz,div,missv,1,1,1,1)
 
END subroutine div_flx_unwgt

SUBROUTINE div_flx_wgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,                 &
                       div_x,div_y,div_z,div)

  integer           , intent(in) ::  nx, ny, nz
  real              , intent(in) ::  dlon, missv
  real, dimension(:), intent(in) ::  lat

  real, dimension(:,:,:), intent(inout) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  integer ::  j,k

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)
  div_x(:,:,:) = div_x(:,:,:) * 86400.

  div(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div

  fy(:,:,:) = fy(:,:,:)*div(:,:,:)

  call grady_2nd(nx,ny,nz,fy(:,:,:)*div,lat, div_y)
  div_y(:,:,:) = div_y(:,:,:)/(div(:,:,:)*div(:,:,:)) * 86400.
  call missing_bdy(nx,ny,nz,div_y,missv,1,1,0,0)

  div(:,:,:) = spread(rho0(:,:),1,nx)   ! rent div

  fz(:,:,:) = fz(:,:,:)*div(:,:,:)

  call gradz_2nd_irr(nx,ny,nz,fz,zp, div_z)
  div_z(:,:,:) = div_z(:,:,:)/div(:,:,:) * 86400.
  call missing_bdy(nx,ny,nz,div_z,missv,0,0,1,1)

  div(:,:,:) = div_x(:,:,:) + div_y(:,:,:) + div_z(:,:,:)
  call missing_bdy(nx,ny,nz,div,missv,1,1,1,1)
 
  do k=1, nz
  do j=1, ny
    fx(:,j,k) = fx(:,j,k)*(rho0(j,k)*coslat(j,k))
    fy(:,j,k) = fy(:,j,k)*rho0(j,k)
    fz(:,j,k) = fz(:,j,k)*coslat(j,k)
  enddo
  enddo

END subroutine div_flx_wgt

END module tem3d


MODULE tem3d

  use util,  only: gradx_2nd, grady_2nd, gradz_2nd_irr, gradxx_2nd,      &
                   get_ind_bnd, missing_bdy, gradxx_4nd
  use const_glob,  only:  g, rd, kappa, a_earth, ome_earth, deg2rad

  implicit none

  private ::  set_gridvar_p, set_gridvar_qg, div_flx, div_norm_flx0,     &
              div_norm_flx0_wgt,                                         &
              get_reference_qg, geostrophic_wind, smooth_y, smooth_xy
  private ::  g, rd, kappa, a_earth, ome_earth, deg2rad

  character(len=64), dimension(11) ::  varname_tem3d_p =                 &
      (/'v_res  ','w_res  ','cor    ','uadv_y ','uadv_z ',               &
        'epd    ','epd_y  ','epd_z  ','f_y    ','f_z    ',               &
        'u_force'/)
  character(len=64), dimension(20) ::  varname_waf3d_s_qg =              &
      (/'divF   ','divF_x ','divF_y ','divF_z ','f_x    ',               &
        'f_y    ','f_z    ',                                             &
        'divB   ','divB_x ','divB_y ','divB_z ','b_x    ',               &
        'b_y    ','b_z    ',                                             &
        'A_F    ','A_B    ','QGPV   ',                                   &
        'U0     ','S      ','beta_eff'/)
  character(len=64), dimension(22) ::  varname_waf3d_nons_qg =           &
      (/'divM   ','divM_x ','divM_y ','divM_z ','m_x    ',               &
        'm_y    ','m_z    ',                                             &
        'divF   ','divF_x ','divF_y ','divF_z ','n_x    ',               &
        'n_y    ','n_z    ',                                             &
        'A_M    ','u_res  ','v_res  ','w_res  ',                         &
        'S      ','Q      ','dQdy   ','dQdx   '/)

  integer,                           private ::  ny_pre, nz_pre
  integer,                           private ::  j5s, j5e, jeq
  integer,                           private ::  jns, jne, jss, jse
  integer,                           private ::  jtns, jtne, jtss, jtse
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
     divf,divf_x,divf_y,divf_z,fx,fy,fz,                                 &
     divb,divb_x,divb_y,divb_z,bx,by,bz,                                 &
     wadf,wadb,qgpv,                                                     &
     u0,s_t0,dq0dy, l_scaling )
! Plumb (1985, JAS, Eq. 7.1)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon, h_scale, missv
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nz)      , intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  u, v, t, gp
  logical                  , intent(in), optional ::  l_scaling

  real, dimension(ny,nz)   , intent(out) ::  u0, s_t0, dq0dy
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  divf, divf_x, divf_y, divf_z
  real, dimension(nx,ny,nz), intent(out) ::  bx, by, bz
  real, dimension(nx,ny,nz), intent(out) ::  divb, divb_x, divb_y, divb_z
  real, dimension(nx,ny,nz), intent(out) ::  wadf, wadb, qgpv

  real, dimension(ny,nz)    ::  t0, gp0, v0, dvor0dy, dt0dy, s_t00, scl
  real, dimension(nx,ny,nz) ::  prt, u_prt, v_prt, e_tot, tmp

  logical ::  l_scl
  integer ::  j,k

  call set_gridvar_qg(ny,nz,lat,p,h_scale)

  u0 (:,:) = sum(u , dim=1)/float(nx)
  t0 (:,:) = sum(t , dim=1)/float(nx)
  v0 (:,:) = sum(v , dim=1)/float(nx)
  gp0(:,:) = sum(gp, dim=1)/float(nx)

  ! S
  call gradz_2nd_irr(ny,nz,t0,zp, s_t0)
  s_t0(:,:) = s_t0(:,:) + t0(:,:)*(kappa/h_scale)

!  call get_reference_qg(s_t0, s_t00)
  call smooth_y(s_t0,7, s_t00)

  s_t0(:,:) = s_t00(:,:)

  ! effective beta
  call grady_2nd(ny,nz,-u0*coslat,lat, dt0dy)   ! rent dt0dy
  call grady_2nd(ny,nz,dt0dy/coslat,lat, dvor0dy)
  call grady_2nd(ny,nz,t0,lat, dt0dy)
  call gradz_2nd_irr(ny,nz,rho0*dt0dy/s_t0,zp, dq0dy)
  dq0dy(:,:) = beta(:,:) + dvor0dy(:,:) + f(:,:)*dq0dy(:,:)/rho0(:,:)

  call smooth_y(dq0dy(3:ny-2,:),7, s_t00(3:ny-2,:))
  dq0dy(3:ny-2,:) = s_t00(3:ny-2,:)

  l_scl = .False.
  if ( present(l_scaling) ) then
    if ( l_scaling ) then
      l_scl = .True.
      scl(:,:) = dq0dy(:,:)
    end if
  end if

  ! u', v'
  u_prt(:,:,:) = u(:,:,:) - spread(u0(:,:),1,nx)
  v_prt(:,:,:) = v(:,:,:) - spread(v0(:,:),1,nx)

  ! kinetic E * 2
  e_tot(:,:,:) = u_prt(:,:,:)*u_prt(:,:,:) + v_prt(:,:,:)*v_prt(:,:,:)

  ! relative vortivity pert.
  tmp(:,:,:) = spread(coslat(:,:),1,nx)
  call gradx_2nd(nx,ny,nz,v_prt,lat,dlon, qgpv)
  call grady_2nd(nx,ny,nz,-u_prt*tmp,lat, prt)
  qgpv(:,:,:) = qgpv(:,:,:) + prt(:,:,:)/tmp(:,:,:)

  prt(:,:,:) = gp(:,:,:) - spread(gp0(:,:),1,nx)

  ! B_x, B_y, f*G_x, f*G_y (not mass-weighted; UA and -E excluded)
  bx(:,:,:) = v_prt(:,:,:)*v_prt(:,:,:)
  by(:,:,:) = -(v_prt(:,:,:)*u_prt(:,:,:))
  call gradx_2nd(nx,ny,nz,-0.5*(v_prt*prt),lat,dlon, fx)
  call gradx_2nd(nx,ny,nz,0.5*(u_prt*prt) ,lat,dlon, fy)

  u_prt(:,:,:) = t(:,:,:) - spread(t0(:,:),1,nx)   ! rent u_prt

  ! total E
  e_tot(:,:,:) = 0.5*( e_tot(:,:,:) + (u_prt(:,:,:)*u_prt(:,:,:))*       &
                           spread((rd/h_scale)/s_t0(:,:),1,nx) )

  ! QG PV pert.
  call gradz_2nd_irr(nx,ny,nz,u_prt*spread(rho0/s_t0,1,nx),zp, tmp)
  qgpv(:,:,:) = qgpv(:,:,:) + tmp(:,:,:)*spread(f(:,:)/rho0(:,:),1,nx)

  ! wave activity density, not mass-weighted
  wadb(:,:,:) = 0.5*(qgpv(:,:,:)*qgpv(:,:,:))

  if ( l_scl ) then
    wadb(:,:,:) = wadb(:,:,:)*spread(scl(:,:)/dq0dy(:,:),1,nx)
    e_tot(:,:,:) = e_tot(:,:,:)*spread(scl(:,:),1,nx)
  else
    wadb(:,:,:) = wadb(:,:,:)/spread(dq0dy(:,:),1,nx)
  end if

  ! B_z, G_z (not mass-weighted)
  bz(:,:,:) = v_prt(:,:,:)*u_prt(:,:,:)*spread(f(:,:)/s_t0(:,:),1,nx)
  call gradx_2nd(nx,ny,nz,-0.5*(u_prt*prt),lat,dlon, fz)
  fz(:,:,:) = fz(:,:,:)/spread(s_t0(:,:),1,nx)

  ! B and G will be mass-weighted (by rho0*cos[lat]) by the call below.
  ! G_x and G_y should also be divided by f, later.

  call div_norm_flx0_wgt(nx,ny,nz,lat,dlon,missv,bx,by,bz,               &
                         divb_x,divb_y,divb_z,divb)
  call div_norm_flx0_wgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,               &
                         divf_x,divf_y,divf_z,divf)
  ! div = div[B]/(rho0*cos[lat])  [m/s/day]

  tmp(:,:,:) = spread(f(:,:),1,nx)

  ! F = B + G (B,G: UA and -E excluded)
  fx(:,:,:) = bx(:,:,:) + fx(:,:,:)/tmp(:,:,:)
  fy(:,:,:) = by(:,:,:) + fy(:,:,:)/tmp(:,:,:)
  fz(:,:,:) = bz(:,:,:) + fz(:,:,:)

  divf_x(:,:,:) = divb_x(:,:,:) + divf_x(:,:,:)/tmp(:,:,:)
  divf_y(:,:,:) = divb_y(:,:,:) + divf_y(:,:,:)/tmp(:,:,:)
  divf_z(:,:,:) = divb_z(:,:,:) + divf_z(:,:,:)

  divf(:,:,:) = divf_x(:,:,:) + divf_y(:,:,:) + divf_z(:,:,:)

  if ( l_scl ) then
    tmp(:,:,:) = spread(scl(:,:),1,nx)
    fx(:,:,:) = fx(:,:,:)*tmp(:,:,:)  ;  bx(:,:,:) = bx(:,:,:)*tmp(:,:,:)
    fy(:,:,:) = fy(:,:,:)*tmp(:,:,:)  ;  by(:,:,:) = by(:,:,:)*tmp(:,:,:)
    fz(:,:,:) = fz(:,:,:)*tmp(:,:,:)  ;  bz(:,:,:) = bz(:,:,:)*tmp(:,:,:)
    divf_x(:,:,:) = divf_x(:,:,:)*tmp(:,:,:)
    divf_y(:,:,:) = divf_y(:,:,:)*tmp(:,:,:)
    divf_z(:,:,:) = divf_z(:,:,:)*tmp(:,:,:)
    divf  (:,:,:) = divf  (:,:,:)*tmp(:,:,:)
    divb_x(:,:,:) = divb_x(:,:,:)*tmp(:,:,:)
    divb_y(:,:,:) = divb_y(:,:,:)*tmp(:,:,:)
    divb_z(:,:,:) = divb_z(:,:,:)*tmp(:,:,:)
    divb  (:,:,:) = divb  (:,:,:)*tmp(:,:,:)
  end if

  ! include UA - E into B_x
  u_prt(:,:,:) = spread(u0(:,:),1,nx)   ! rent u_prt
  prt(:,:,:) = wadb(:,:,:)*u_prt(:,:,:) - e_tot(:,:,:)   ! rent prt
  call gradx_2nd(nx,ny,nz,prt,lat,dlon, tmp)
  tmp(:,:,:) = tmp(:,:,:) * 86400.
  divb_x(:,:,:) = divb_x(:,:,:) + tmp(:,:,:)
  divb  (:,:,:) = divb  (:,:,:) + tmp(:,:,:)
  tmp(:,:,:) = spread(rho0(:,:)*coslat(:,:),1,nx)
  bx(:,:,:) = bx(:,:,:) + prt(:,:,:)*tmp(:,:,:)

  ! wave activity density
  wadb(:,:,:) = wadb(:,:,:)*tmp(:,:,:)

  wadf(:,:,:) = wadb(:,:,:) + e_tot(:,:,:)*tmp(:,:,:)/u_prt(:,:,:)

  call missing_bdy(nx,ny,nz,bx    ,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,divb_x,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,bz    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,nz,divb_z,missv,0,0,2,2)
  call missing_bdy(nx,ny,nz,divb  ,missv,2,2,2,2)

  call missing_bdy(nx,ny,nz,divf_y,missv,1,1,0,0)
  call missing_bdy(nx,ny,nz,fz    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,nz,divf_z,missv,0,0,2,2)
  call missing_bdy(nx,ny,nz,divf  ,missv,1,1,2,2)

  call missing_bdy(nx,ny,nz,wadb,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,wadf,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,qgpv,missv,2,2,2,2)

  call missing_bdy(ny,nz,s_t0 ,missv,0,0,1,1)
  call missing_bdy(ny,nz,dq0dy,missv,2,2,2,2)

  if ( l_scl ) then
    call missing_bdy(nx,ny,nz,by    ,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,divb_y,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,bz    ,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,divb_z,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,fx    ,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,divf_x,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,fy    ,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,divf_y,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,fz    ,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,divf_z,missv,2,2,2,2)
    call missing_bdy(nx,ny,nz,divf  ,missv,2,2,2,2)
  end if

  if (j5e /= 0) then
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
     divf,divf_x,divf_y,divf_z,fx,fy,fz,u0,s_t0 )
! Plumb (1985, JAS, Eq. 5.7)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon, h_scale, missv
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nz)      , intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  gp, t

  real, dimension(ny,nz)   , intent(out) ::  u0, s_t0
  real, dimension(nx,ny,nz), intent(out) ::  fx, fy, fz
  real, dimension(nx,ny,nz), intent(out) ::  divf, divf_x, divf_y, divf_z

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
  if ( jeq /= -999 ) then
    u0(jeq,:) = 0.5*(u0(jeq-1,:) + u0(jeq+1,:))
  end if

  prt(:,:,:) = (gp(:,:,:) - spread(gp0(:,:),1,nx))/spread(f(:,:),1,nx)
  if ( jeq /= -999 ) then
    prt(:,jeq,:) = 0.5*(prt(:,jeq-1,:) + prt(:,jeq+1,:))
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

  call div_norm_flx0_wgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,               &
                         divf_x,divf_y,divf_z,divf)
  ! F = (rho*cos) * F0 : mass-weighted
  ! div = div[F]/(rho*cos)  [m/s/day]

  call missing_bdy(nx,ny,nz,fy    ,missv,1,1,0,0)
  call missing_bdy(nx,ny,nz,divf_y,missv,2,2,0,0)
  call missing_bdy(nx,ny,nz,fz    ,missv,0,0,1,1)  ! due to s_t0
  call missing_bdy(nx,ny,nz,divf_z,missv,0,0,2,2)
  call missing_bdy(nx,ny,nz,divf  ,missv,2,2,2,2)

  call missing_bdy(ny,nz,s_t0,missv,0,0,1,1)

  if (j5e /= 0) then
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

SUBROUTINE waf3d_nons_qg(                                                &
     nx,ny,nz,lat,p,gp_prt,t_prt,gps,ts,us,vs,ws,dlon,h_scale,missv,     &
     divm,divm_x,divm_y,divm_z,mx,my,mz,                                 &
     divf,divf_x,divf_y,divf_z,lx,ly,lz,                                 &
     wadm,u_ar,v_ar,w_ar,                                                &
     s_ts,qs,dqy,dqx )
! Plumb (1986, JAS, Eq. 2.24)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon, h_scale, missv
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nz)      , intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  gps, ts, us, vs, ws

  real, dimension(nx,ny,nz), intent(inout) ::  gp_prt, t_prt

  real, dimension(nx,ny,nz), intent(out) ::  s_ts, qs, dqy, dqx
  real, dimension(nx,ny,nz), intent(out) ::  mx, my, mz
  real, dimension(nx,ny,nz), intent(out) ::  divm, divm_x, divm_y, divm_z
  real, dimension(nx,ny,nz), intent(out) ::  lx, ly, lz
  real, dimension(nx,ny,nz), intent(out) ::  divf, divf_x, divf_y, divf_z
  real, dimension(nx,ny,nz), intent(out) ::  wadm
  real, dimension(nx,ny,nz), intent(out) ::  u_ar, v_ar, w_ar

  real, dimension(nx,ny,nz) ::  s_ts0, dqsdh
  real, dimension(nx,ny,nz) ::  u_prt, v_prt, q_prt, us_g, vs_g
  real, dimension(nx,ny,nz) ::  u2, v2, vu, ut, vt, e_tot, tmp, tmp2
  real, dimension(ny,nz)    ::  tmp2d

  integer ::  j,k

  call set_gridvar_qg(ny,nz,lat,p,h_scale)

  ! stationary (mean) geostrophic wind
  call geostrophic_wind(nx,ny,nz,lat,dlon,gps, us_g,vs_g)

  ! stationary (mean) ageostrophic wind
  u_ar(:,:,:) = us(:,:,:) - us_g(:,:,:)
  v_ar(:,:,:) = vs(:,:,:) - vs_g(:,:,:)
  w_ar(:,:,:) = ws(:,:,:)

  ! transient (wave) geostrophic wind
  gp_prt(:,:,:) = gp_prt(:,:,:) - spread(sum(gp_prt, dim=1)/float(nx),   &
                                         1,nx)
  call geostrophic_wind(nx,ny,nz,lat,dlon,gp_prt, u_prt,v_prt)
 
  ! S
  call gradz_2nd_irr(nx,ny,nz,ts,zp, s_ts)
  s_ts(:,:,:) = s_ts(:,:,:) + ts(:,:,:)*(kappa/h_scale)

  call get_reference_qg(sum(s_ts(:,:,:), dim=1)/float(nx), tmp2d)
  s_ts0(:,:,:) = spread(tmp2d(:,:),1,nx)

  s_ts(:,:,:) = s_ts0(:,:,:)

  ! mean Q and del[Q]
  tmp(:,:,:) = spread(coslat(:,:),1,nx)
  call grady_2nd(nx,ny,nz,-us_g*tmp,lat, tmp2)
  call gradx_2nd(nx,ny,nz,vs_g,lat,dlon, qs)
  qs(:,:,:) = qs(:,:,:) + tmp2(:,:,:)/tmp(:,:,:)   ! rel. vorticity
  tmp(:,:,:) = spread(rho0(:,:),1,nx)
  u2 (:,:,:) = spread(f   (:,:),1,nx)   ! rent u2
  call grady_2nd(nx,ny,nz,qs,lat, v2)   ! rent v2 : rel. vorticity gradient
  call grady_2nd(nx,ny,nz,ts,lat, tmp2)
  call gradz_2nd_irr(nx,ny,nz,tmp*tmp2/s_ts,zp, dqy)
  dqy(:,:,:) = spread(beta(:,:),1,nx) +                                  &
               v2(:,:,:) + u2(:,:,:)*dqy(:,:,:)/tmp(:,:,:)
  call gradx_2nd(nx,ny,nz,qs,lat,dlon, v2)
  call gradx_2nd(nx,ny,nz,ts,lat,dlon, tmp2)
  call gradz_2nd_irr(nx,ny,nz,tmp*tmp2/s_ts,zp, dqx)
  dqx(:,:,:) = v2(:,:,:) + u2(:,:,:)*dqx(:,:,:)/tmp(:,:,:)
  call gradz_2nd_irr(nx,ny,nz,tmp*ts/s_ts,zp, tmp2)
  qs(:,:,:) = u2(:,:,:) + qs(:,:,:) + u2(:,:,:)*tmp2(:,:,:)/tmp(:,:,:)
!n  dqsdh(:,:,:) = sqrt( dqx(:,:,:)*dqx(:,:,:) + dqy(:,:,:)*dqy(:,:,:) )
  dqsdh(:,:,:) = 1.  !n  ! for A and M not to be normalized by |del[Q]|
  dqy(:,:,:) = dqy(:,:,:)/dqsdh(:,:,:)
  dqx(:,:,:) = dqx(:,:,:)/dqsdh(:,:,:)
 
  t_prt(:,:,:) = t_prt(:,:,:) - spread(sum(t_prt, dim=1)/float(nx),1,nx)

  u2(:,:,:) = u_prt(:,:,:)*u_prt(:,:,:)
  v2(:,:,:) = v_prt(:,:,:)*v_prt(:,:,:)
  vu(:,:,:) = v_prt(:,:,:)*u_prt(:,:,:)
  ut(:,:,:) = u_prt(:,:,:)*t_prt(:,:,:)/s_ts(:,:,:)
  vt(:,:,:) = v_prt(:,:,:)*t_prt(:,:,:)/s_ts(:,:,:)

  ! total E
  e_tot(:,:,:) = 0.5*( (u2(:,:,:) + v2(:,:,:)) +                         &
                  (t_prt(:,:,:)*t_prt(:,:,:)/s_ts(:,:,:))*(rd/h_scale) )

  ! QG PV pert.
  tmp(:,:,:) = spread(coslat(:,:),1,nx)
  call grady_2nd(nx,ny,nz,-u_prt*tmp,lat, q_prt)
  q_prt(:,:,:) = q_prt(:,:,:)/tmp(:,:,:)
  call gradx_2nd(nx,ny,nz,v_prt,lat,dlon, tmp)
  q_prt(:,:,:) = tmp(:,:,:) + q_prt(:,:,:)
  call gradz_2nd_irr(nx,ny,nz,t_prt/s_ts*spread(rho0,1,nx),zp, tmp)
  q_prt(:,:,:) = q_prt(:,:,:) + tmp(:,:,:)*spread(f(:,:)/rho0(:,:),1,nx)

  ! ageostrophic residual-mean wind: stationary wind + wave effect
  wadm(:,:,:) = (u2(:,:,:) + v2(:,:,:) - e_tot(:,:,:))/                  &
                spread(f(:,:),1,nx)   ! rent wadm
  if (jeq /= -999)  wadm(:,jeq,:) = 0.5*( wadm(:,jeq-1,:) +              &
                                          wadm(:,jeq+1,:) )
  call gradx_2nd(nx,ny,nz,ut,lat,dlon, tmp)
  w_ar(:,:,:) = w_ar(:,:,:) + tmp(:,:,:)
  call gradx_2nd(nx,ny,nz,-wadm,lat,dlon, tmp)
  v_ar(:,:,:) = v_ar(:,:,:) + tmp(:,:,:)
  tmp2(:,:,:) = spread(coslat(:,:),1,nx)
  call grady_2nd(nx,ny,nz,vt*tmp2,lat, tmp)
  w_ar(:,:,:) = w_ar(:,:,:) + tmp(:,:,:)/tmp2(:,:,:)
  call grady_2nd(nx,ny,nz,wadm*tmp2,lat, tmp)
  u_ar(:,:,:) = u_ar(:,:,:) + tmp(:,:,:)/tmp2(:,:,:)
  tmp2(:,:,:) = spread(rho0(:,:),1,nx)
  call gradz_2nd_irr(nx,ny,nz,-vt*tmp2,zp, tmp)
  v_ar(:,:,:) = v_ar(:,:,:) + tmp(:,:,:)/tmp2(:,:,:)
  call gradz_2nd_irr(nx,ny,nz,-ut*tmp2,zp, tmp)
  u_ar(:,:,:) = u_ar(:,:,:) + tmp(:,:,:)/tmp2(:,:,:)

  ! wave activity density, not mass-weighted
  wadm(:,:,:) = 0.5*(q_prt(:,:,:)*q_prt(:,:,:))/dqsdh(:,:,:)

  ! M, mass-weighted
  tmp(:,:,:) = spread(rho0(:,:)*coslat(:,:),1,nx)
  mx(:,:,:) = tmp(:,:,:)*(v2(:,:,:) - e_tot(:,:,:))
  lx(:,:,:) = tmp(:,:,:)*vu(:,:,:)
  my(:,:,:) = -lx(:,:,:)
  ly(:,:,:) = tmp(:,:,:)*(e_tot(:,:,:) - u2(:,:,:))

  tmp(:,:,:) = spread(tmp(1,:,:)*f(:,:),1,nx)
  mz(:,:,:) = vt(:,:,:)*tmp(:,:,:)
  lz(:,:,:) = ut(:,:,:)*tmp(:,:,:)

  ! div = div[M]  [kg m^-3 m/s/day]
  call div_flx(nx,ny,nz,lat,dlon,missv,mx,my,mz,                         &
               divm_x,divm_y,divm_z,divm)
  call div_flx(nx,ny,nz,lat,dlon,missv,lx,ly,lz,                         &
               divf_x,divf_y,divf_z,divf)

  ! correction
  tmp(:,:,:) = spread(rho0(:,:)*spread(sin(lat(:)*deg2rad)/a_earth,2,nz) &
                      ,1,nx)*                                            &
      ( (t_prt(:,:,:)*t_prt(:,:,:)/s_ts(:,:,:))*(rd/h_scale) ) * 86400.
  divf_y(:,:,:) = divf_y(:,:,:) + tmp(:,:,:)
  divf  (:,:,:) = divf  (:,:,:) + tmp(:,:,:)

  divm_x(:,:,:) = dqy(:,:,:)*divm_x(:,:,:) + dqx(:,:,:)*divf_x(:,:,:)
  divm_y(:,:,:) = dqy(:,:,:)*divm_y(:,:,:) + dqx(:,:,:)*divf_y(:,:,:)
  divm_z(:,:,:) = dqy(:,:,:)*divm_z(:,:,:) + dqx(:,:,:)*divf_z(:,:,:)
  divm  (:,:,:) = dqy(:,:,:)*divm  (:,:,:) + dqx(:,:,:)*divf  (:,:,:)
  mx    (:,:,:) = dqy(:,:,:)*mx    (:,:,:) + dqx(:,:,:)*lx    (:,:,:)
  my    (:,:,:) = dqy(:,:,:)*my    (:,:,:) + dqx(:,:,:)*ly    (:,:,:)
  mz    (:,:,:) = dqy(:,:,:)*mz    (:,:,:) + dqx(:,:,:)*lz    (:,:,:)

  ! include advection and mass-weight A
  call gradx_2nd(nx,ny,nz,wadm,lat,dlon, lx)
  call grady_2nd(nx,ny,nz,wadm,lat, ly)
  tmp(:,:,:) = spread(rho0(:,:)*coslat(:,:),1,nx)
  divm_x(:,:,:) = divm_x(:,:,:) + tmp(:,:,:)*(us_g(:,:,:)*lx(:,:,:))
  divm_y(:,:,:) = divm_y(:,:,:) + tmp(:,:,:)*(vs_g(:,:,:)*ly(:,:,:))
  divm  (:,:,:) = divm_x(:,:,:) + divm_y(:,:,:) + divm_z(:,:,:)
  wadm(:,:,:) = wadm(:,:,:)*tmp(:,:,:)
  mx(:,:,:) = mx(:,:,:) + us_g(:,:,:)*wadm(:,:,:)
  my(:,:,:) = my(:,:,:) + vs_g(:,:,:)*wadm(:,:,:)

!+test :  normalized as in Plumb 86
  dqsdh(:,:,:) = sqrt( dqx(:,:,:)*dqx(:,:,:) + dqy(:,:,:)*dqy(:,:,:) )
  lx = mx/dqsdh  ;  ly = my/dqsdh  ;  lz = mz/dqsdh
  call div_flx(nx,ny,nz,lat,dlon,missv,lx,ly,lz,                         &
               divf_x,divf_y,divf_z,divf)
!-test

!  mx(:,:,:) = dqy(:,:,:)*(v2(:,:,:) - e_tot(:,:,:)) + dqx(:,:,:)*vu(:,:,:)
!  my(:,:,:) = -dqx(:,:,:)*(u2(:,:,:) - e_tot(:,:,:)) - dqy(:,:,:)*vu(:,:,:)
!  lx(:,:,:) = -dqx(:,:,:)*(v2(:,:,:) - e_tot(:,:,:)) + dqy(:,:,:)*vu(:,:,:)
!  ly(:,:,:) = -dqy(:,:,:)*(u2(:,:,:) - e_tot(:,:,:)) + dqx(:,:,:)*vu(:,:,:)
!  mz(:,:,:) = (dqy(:,:,:)*vt(:,:,:) + dqx(:,:,:)*ut(:,:,:))*tmp(:,:,:)
!  lz(:,:,:) = (-dqx(:,:,:)*vt(:,:,:) + dqy(:,:,:)*ut(:,:,:))*tmp(:,:,:)
!  mx(:,:,:) = mx(:,:,:) + us_g(:,:,:)*wadm(:,:,:)
!  my(:,:,:) = my(:,:,:) + vs_g(:,:,:)*wadm(:,:,:)

!  ! div[F] = div[M + N]
!  divf_x(:,:,:) = divm_x(:,:,:) + divf_x(:,:,:)
!  divf_y(:,:,:) = divm_y(:,:,:) + divf_y(:,:,:)
!  divf_z(:,:,:) = divm_z(:,:,:) + divf_z(:,:,:)
!  divf  (:,:,:) = divm  (:,:,:) + divf  (:,:,:)

  call missing_bdy(nx,ny,nz,u_ar,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,v_ar,missv,1,1,2,2)
  call missing_bdy(nx,ny,nz,w_ar,missv,2,2,2,2)
 
  call missing_bdy(nx,ny,nz,mx    ,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,divm_x,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,my    ,missv,3,3,2,2)
  call missing_bdy(nx,ny,nz,divm_y,missv,3,3,2,2)
  call missing_bdy(nx,ny,nz,mz    ,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,divm_z,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,divm  ,missv,3,3,2,2)

  call missing_bdy(nx,ny,nz,lx    ,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,divf_x,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,ly    ,missv,3,3,2,2)
  call missing_bdy(nx,ny,nz,divf_y,missv,3,3,2,2)
  call missing_bdy(nx,ny,nz,lz    ,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,divf_z,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,divf  ,missv,3,3,2,2)

  call missing_bdy(nx,ny,nz,wadm,missv,2,2,2,2)

  call missing_bdy(nx,ny,nz,s_ts,missv,0,0,1,1)
  call missing_bdy(nx,ny,nz,qs  ,missv,2,2,2,2)
  call missing_bdy(nx,ny,nz,dqy ,missv,3,3,2,2)
  call missing_bdy(nx,ny,nz,dqx ,missv,2,2,2,2)

  if (j5e /= 0) then
    u_ar  (:,j5s:j5e,:) = missv
    v_ar  (:,j5s:j5e,:) = missv
    w_ar  (:,j5s:j5e,:) = missv
    divf  (:,j5s:j5e,:) = missv
    divf_x(:,j5s:j5e,:) = missv
    divf_y(:,j5s:j5e,:) = missv
    divf_z(:,j5s:j5e,:) = missv
    mx    (:,j5s:j5e,:) = missv
    my    (:,j5s:j5e,:) = missv
    mz    (:,j5s:j5e,:) = missv
    divm  (:,j5s:j5e,:) = missv
    divm_x(:,j5s:j5e,:) = missv
    divm_y(:,j5s:j5e,:) = missv
    divm_z(:,j5s:j5e,:) = missv
    lx    (:,j5s:j5e,:) = missv
    ly    (:,j5s:j5e,:) = missv
    lz    (:,j5s:j5e,:) = missv
    wadm  (:,j5s:j5e,:) = missv
    s_ts  (:,j5s:j5e,:) = missv
    qs    (:,j5s:j5e,:) = missv
    dqy   (:,j5s:j5e,:) = missv
    dqx   (:,j5s:j5e,:) = missv
  end if

END subroutine waf3d_nons_qg

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
  rho0(:,:) = spread(p(:)/(g*h_scale)                ,1,ny)

  call get_ind_bnd(ny,lat,(/-5.,5./), j5s,j5e)

  jeq = -999
  if ( minval(abs(lat)) == 0. )  jeq = sum(minloc(abs(lat)))
 
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
  call get_ind_bnd(ny,lat,(/1.e-6 ,20. /), jtns,jtne)
  call get_ind_bnd(ny,lat,(/-1.e-6,-20./), jtss,jtse)

  l_nh = ( maxval(lat) > 20.  )
  l_sh = ( minval(lat) < -20. )

  call set_gridvar_p(ny,nz,lat,p,h_scale)

END subroutine set_gridvar_qg

SUBROUTINE div_flx(nx,ny,nz,lat,dlon,missv,fx,fy,fz,                     &
                   div_x,div_y,div_z,div)

  integer               , intent(in) ::  nx, ny, nz
  real                  , intent(in) ::  dlon, missv
  real, dimension(:)    , intent(in) ::  lat
  real, dimension(:,:,:), intent(in) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)

  div(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div

  call grady_2nd(nx,ny,nz,fy*div,lat, div_y)
  div_y(:,:,:) = div_y(:,:,:)/div(:,:,:)
  call missing_bdy(nx,ny,nz,div_y,missv,1,1,0,0)

  call gradz_2nd_irr(nx,ny,nz,fz,zp, div_z)
  call missing_bdy(nx,ny,nz,div_z,missv,0,0,1,1)

  div_x(:,:,:) = div_x(:,:,:)*86400.
  div_y(:,:,:) = div_y(:,:,:)*86400.
  div_z(:,:,:) = div_z(:,:,:)*86400.

  div(:,:,:) = div_x(:,:,:) + div_y(:,:,:) + div_z(:,:,:)
  call missing_bdy(nx,ny,nz,div,missv,1,1,1,1)
 
END subroutine div_flx

SUBROUTINE div_norm_flx0(nx,ny,nz,lat,dlon,missv,fx,fy,fz,               &
                         div_x,div_y,div_z,div)

  integer               , intent(in) ::  nx, ny, nz
  real                  , intent(in) ::  dlon, missv
  real, dimension(:)    , intent(in) ::  lat
  real, dimension(:,:,:), intent(in) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)

  div(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div

  call grady_2nd(nx,ny,nz,(fy*div)*div,lat, div_y)
  div_y(:,:,:) = div_y(:,:,:)/(div(:,:,:)*div(:,:,:))
  call missing_bdy(nx,ny,nz,div_y,missv,1,1,0,0)

  div(:,:,:) = spread(rho0(:,:),1,nx)   ! rent div

  call gradz_2nd_irr(nx,ny,nz,fz*div,zp, div_z)
  div_z(:,:,:) = div_z(:,:,:)/div(:,:,:)
  call missing_bdy(nx,ny,nz,div_z,missv,0,0,1,1)

  div_x(:,:,:) = div_x(:,:,:)*86400.
  div_y(:,:,:) = div_y(:,:,:)*86400.
  div_z(:,:,:) = div_z(:,:,:)*86400.

  div(:,:,:) = div_x(:,:,:) + div_y(:,:,:) + div_z(:,:,:)
  call missing_bdy(nx,ny,nz,div,missv,1,1,1,1)
 
END subroutine div_norm_flx0

SUBROUTINE div_norm_flx0_wgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,           &
                             div_x,div_y,div_z,div)

  integer           , intent(in) ::  nx, ny, nz
  real              , intent(in) ::  dlon, missv
  real, dimension(:), intent(in) ::  lat

  real, dimension(:,:,:), intent(inout) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  integer ::  j,k

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)

  div(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div

  fy(:,:,:) = fy(:,:,:)*div(:,:,:)

  call grady_2nd(nx,ny,nz,fy(:,:,:)*div,lat, div_y)
  div_y(:,:,:) = div_y(:,:,:)/(div(:,:,:)*div(:,:,:))
  call missing_bdy(nx,ny,nz,div_y,missv,1,1,0,0)

  div(:,:,:) = spread(rho0(:,:),1,nx)   ! rent div

  fz(:,:,:) = fz(:,:,:)*div(:,:,:)

  call gradz_2nd_irr(nx,ny,nz,fz,zp, div_z)
  div_z(:,:,:) = div_z(:,:,:)/div(:,:,:)
  call missing_bdy(nx,ny,nz,div_z,missv,0,0,1,1)

  div_x(:,:,:) = div_x(:,:,:)*86400.
  div_y(:,:,:) = div_y(:,:,:)*86400.
  div_z(:,:,:) = div_z(:,:,:)*86400.

  div(:,:,:) = div_x(:,:,:) + div_y(:,:,:) + div_z(:,:,:)
  call missing_bdy(nx,ny,nz,div,missv,1,1,1,1)
 
  do k=1, nz
  do j=1, ny
    fx(:,j,k) = fx(:,j,k)*(rho0(j,k)*coslat(j,k))
    fy(:,j,k) = fy(:,j,k)*rho0(j,k)
    fz(:,j,k) = fz(:,j,k)*coslat(j,k)
  enddo
  enddo

END subroutine div_norm_flx0_wgt

SUBROUTINE get_reference_qg(var,ref)

  real, dimension(:,:), intent(in) ::  var

  real, dimension(size(var,1),size(var,2)), intent(out) ::  ref

  integer ::  nz
  real    ::  sumcoslat
  integer ::  k

  nz = size(var,2)

  if ( l_nh ) then
    sumcoslat = sum(coslat(jns:jne,1))
    do k=1, nz
      ref(jns:jne,k) = sum(var(jns:jne,k)*coslat(jns:jne,k), dim=1)/sumcoslat
      ref(jtns:jtne,k) = ref(jns,k)
    enddo
  end if
  if ( l_sh ) then
    sumcoslat = sum(coslat(jss:jse,1))
    do k=1, nz
      ref(jss:jse,k) = sum(var(jss:jse,k)*coslat(jss:jse,k), dim=1)/sumcoslat
      ref(jtss:jtse,k) = ref(jss,k)
    enddo
  end if
  if ( l_nh .and. .not. l_sh ) then
    do k=1, nz
      ref(:,k) = ref(jns,k)
    enddo
  else if ( l_sh .and. .not. l_nh ) then
    do k=1, nz
      ref(:,k) = ref(jss,k)
    enddo
  else
    if ( jeq /= -999 ) then
      ref(jeq,:) = 0.5*(ref(jeq-1,:) + ref(jeq+1,:))
    end if
  end if 

END subroutine get_reference_qg

SUBROUTINE geostrophic_wind(nx,ny,nz,lat,dlon,gp,ug,vg)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nx,ny,nz), intent(in) ::  gp

  real, dimension(nx,ny,nz), intent(out) ::  ug, vg

  real, dimension(nx,ny,nz) ::  tmp

  call gradx_2nd(nx,ny,nz,gp,lat,dlon, vg)
  call grady_2nd(nx,ny,nz,-gp,lat, ug)
  tmp(:,:,:) = spread(f(:,:),1,nx)
  ug(:,:,:) = ug(:,:,:)/tmp(:,:,:)
  vg(:,:,:) = vg(:,:,:)/tmp(:,:,:)
  if (jeq /= -999) then
    ug(:,jeq,:) = 0.5*(ug(:,jeq-1,:) + ug(:,jeq+1,:))
    vg(:,jeq,:) = 0.5*(vg(:,jeq-1,:) + vg(:,jeq+1,:))
  end if

END subroutine geostrophic_wind

SUBROUTINE smooth_y(var,npt,ref)

  real, dimension(:,:), intent(in) ::  var
  integer             , intent(in) ::  npt

  real, dimension(size(var,1),size(var,2)), intent(out) ::  ref

  real, dimension(size(var,1)+npt/2*2,size(var,2)) ::  varbuf
  integer ::  ny, nz, nh
  integer ::  k, jj

  ny = size(var,1)  ;  nz = size(var,2)

  nh = npt/2

  varbuf(1+nh:ny+nh,:) = var(:,:)
  do k=1, nz
    varbuf(1:nh,k) = var(1,k)
    varbuf(ny+nh+1:ny+nh+nh,k) = var(ny,k)
  enddo

  ref(:,:) = 0.
  do jj=0, nh+nh
    ref(:,:) = ref(:,:) + varbuf(1+jj:ny+jj,:)
  enddo
  ref(:,:) = ref(:,:)/float(npt)

END subroutine smooth_y

SUBROUTINE smooth_xy(var,npt,ref)

  real, dimension(:,:,:), intent(in) ::  var
  integer               , intent(in) ::  npt

  real, dimension(size(var,1),size(var,2),size(var,3)), intent(out) ::  ref

  real, dimension(size(var,1)+npt/2*2,size(var,2)+npt/2*2,size(var,3)) ::  varbuf
  integer ::  nx, ny, nz, nh
  integer ::  j,k, ii, jj

  nx = size(var,1)  ;  ny = size(var,2)  ;  nz = size(var,3)

  nh = npt/2

  varbuf(1+nh:nx+nh,1+nh:ny+nh,:) = var(:,:,:)
  do k=1, nz
    do j=1, nh
      varbuf(1+nh:nx+nh,j,k) = var(:,1,k)
    enddo
    do j=ny+nh+1, ny+nh+nh
      varbuf(1+nh:nx+nh,j,k) = var(:,ny,k)
    enddo
  enddo
  varbuf(1:nh,:,:) = varbuf(nx+1:nx+nh,:,:)
  varbuf(nx+nh+1:nx+nh+nh,:,:) = varbuf(1+nh:nh+nh,:,:)

  ref(:,:,:) = 0.
  do jj=0, nh+nh
  do ii=0, nh+nh
    ref(:,:,:) = ref(:,:,:) + varbuf(1+ii:nx+ii,1+jj:ny+jj,:)
  enddo
  enddo
  ref(:,:,:) = ref(:,:,:)/float(npt*npt)

END subroutine smooth_xy

END module tem3d


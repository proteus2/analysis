MODULE tem3d

  use util,  only: gradx_2nd, grady_2nd, gradz_2nd_irr, grady_wgt_2nd,   &
                   gradz_wgt_2nd_irr, gradxx_2nd, gradxx_4nd,            &
                   get_ind_bnd, missing_bdy, lowpass_k, filter121_y
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

  integer, parameter ::  kb_max = 5

  real, dimension(nx,ny,nz) ::  uf_prt, vf_prt, q_prt, us_g, vs_g
  real, dimension(nx,ny,nz) ::  u2f2, v2f2, vuf2, fut, fvt, f2e_k, e_p
  real, dimension(nx,ny,nz) ::  te, dqh, dqy_prt, dqx_prt, tmp
  real, dimension(ny,nz)    ::  t00

  integer ::  j,k

  call set_gridvar_qg(ny,nz,lat,p,h_scale)

  ! S
  call get_reference_qg(sum(ts(:,:,:), dim=1)/float(nx), t00)
  call gradz_2nd_irr(ny,nz,t00,zp, s_ts(1,:,:))
  s_ts(1,:,:) = s_ts(1,:,:) + t00(:,:)*(kappa/h_scale)
  s_ts(2:nx,:,:) = spread(s_ts(1,:,:),1,nx-1)

  ! T_e = T - T0(z)
  te(:,:,:) = ts(:,:,:) - spread(t00(:,:),1,nx)

  ! stationary (mean) geostrophic wind * f
  call grady_2nd(nx,ny,nz,-gps,lat, us_g)
  call gradx_2nd(nx,ny,nz,gps,lat,dlon, vs_g)

  ! mean Q and grad[Q] ; mean geostrophic wind
  call grady_wgt_2nd(nx,ny,nz,-us_g,lat,coslat(:,1), qs)
  call gradx_2nd(nx,ny,nz,vs_g,lat,dlon, tmp)
  qs(:,:,:) = tmp(:,:,:) + qs(:,:,:)   ! rel. vorticity * f
  u2f2(:,:,:) = spread(f(:,:),1,nx)   ! rent u2f2
  call grady_2nd(nx,ny,nz,qs,lat, v2f2)   ! rent v2f2 : rel. vorticity grad. * f
  call grady_2nd(nx,ny,nz,te,lat, tmp)
  call gradz_wgt_2nd_irr(nx,ny,nz,tmp/s_ts,zp,rho0(1,:), dqy)
  dqy(:,:,:) = spread(beta(:,:),1,nx) +                                  &
               v2f2(:,:,:)/u2f2(:,:,:) + u2f2(:,:,:)*dqy(:,:,:)
  call gradx_2nd(nx,ny,nz,qs,lat,dlon, v2f2)
  call gradx_2nd(nx,ny,nz,te,lat,dlon, tmp)
  call gradz_wgt_2nd_irr(nx,ny,nz,tmp/s_ts,zp,rho0(1,:), dqx)
  dqx(:,:,:) = v2f2(:,:,:)/u2f2(:,:,:) + u2f2(:,:,:)*dqx(:,:,:)
  call gradz_wgt_2nd_irr(nx,ny,nz,te/s_ts,zp,rho0(1,:), tmp)
  qs(:,:,:) = u2f2(:,:,:) + qs(:,:,:)/u2f2(:,:,:) +                      &
              u2f2(:,:,:)*tmp(:,:,:)
  us_g(:,:,:) = us_g(:,:,:)/u2f2(:,:,:)
  vs_g(:,:,:) = vs_g(:,:,:)/u2f2(:,:,:)

  ! smooth grad[Q]
  call lowpass_k(dqy,kb_max)
  call lowpass_k(dqx,kb_max)
  call filter121_y(dqy)
  call filter121_y(dqx)

  ! transient (wave) geostrophic wind * f
  gp_prt(:,:,:) = gp_prt(:,:,:) - spread(sum(gp_prt, dim=1)/float(nx),   &
                                         1,nx)
  call grady_2nd(nx,ny,nz,-gp_prt,lat, uf_prt)
  call gradx_2nd(nx,ny,nz,gp_prt,lat,dlon, vf_prt)

  t_prt(:,:,:) = t_prt(:,:,:) - spread(sum(t_prt, dim=1)/float(nx),1,nx)

  ! QG PV pert.
  call grady_wgt_2nd(nx,ny,nz,-uf_prt,lat,coslat(:,1), q_prt)
  call gradx_2nd(nx,ny,nz,vf_prt,lat,dlon, tmp)
  q_prt(:,:,:) = tmp(:,:,:) + q_prt(:,:,:)
  u2f2(:,:,:) = spread(f(:,:),1,nx)   ! rent u2f2
  call grady_2nd(nx,ny,nz,q_prt,lat, v2f2)   ! rent v2f2
  call grady_2nd(nx,ny,nz,t_prt,lat, tmp)
  call gradz_wgt_2nd_irr(nx,ny,nz,tmp/s_ts,zp,rho0(1,:), dqy_prt)
  dqy_prt(:,:,:) = v2f2(:,:,:)/u2f2(:,:,:) + u2f2(:,:,:)*dqy_prt(:,:,:)
  call gradx_2nd(nx,ny,nz,q_prt,lat,dlon, v2f2)
  call gradx_2nd(nx,ny,nz,t_prt,lat,dlon, tmp)
  call gradz_wgt_2nd_irr(nx,ny,nz,tmp/s_ts,zp,rho0(1,:), dqx_prt)
  dqx_prt(:,:,:) = v2f2(:,:,:)/u2f2(:,:,:) + u2f2(:,:,:)*dqx_prt(:,:,:)
 
  call gradz_wgt_2nd_irr(nx,ny,nz,t_prt/s_ts,zp,rho0(1,:), tmp)
  q_prt(:,:,:) = q_prt(:,:,:)/u2f2(:,:,:) + u2f2(:,:,:)*tmp(:,:,:)

  u2f2(:,:,:) = uf_prt(:,:,:)*uf_prt(:,:,:)
  v2f2(:,:,:) = vf_prt(:,:,:)*vf_prt(:,:,:)
  vuf2(:,:,:) = vf_prt(:,:,:)*uf_prt(:,:,:)

  fut(:,:,:) = uf_prt(:,:,:)*t_prt(:,:,:)/s_ts(:,:,:)
  fvt(:,:,:) = vf_prt(:,:,:)*t_prt(:,:,:)/s_ts(:,:,:)

  ! kinetic E * f^2
  f2e_k(:,:,:) = 0.5*(u2f2(:,:,:) + v2f2(:,:,:))

  ! potential E
  e_p(:,:,:) = 0.5*(rd/h_scale)*(t_prt(:,:,:)*t_prt(:,:,:)/s_ts(:,:,:))

  ! ageostrophic residual-mean wind: stationary wind + wave effect
  u_ar(:,:,:) = us(:,:,:) - us_g(:,:,:)
  v_ar(:,:,:) = vs(:,:,:) - vs_g(:,:,:)
  w_ar(:,:,:) = ws(:,:,:)
  te(:,:,:) = spread(f(:,:),1,nx)   ! rent te (discarded)
  call gradx_2nd(nx,ny,nz,fut,lat,dlon, tmp)
  w_ar(:,:,:) = w_ar(:,:,:) + tmp(:,:,:)/te(:,:,:)
  call gradx_2nd(nx,ny,nz,e_p-f2e_k/(te*te),lat,dlon, tmp)
  v_ar(:,:,:) = v_ar(:,:,:) + tmp(:,:,:)/te(:,:,:)
  call grady_wgt_2nd(nx,ny,nz,fvt,lat,coslat(:,1), tmp)
  w_ar(:,:,:) = w_ar(:,:,:) + tmp(:,:,:)/te(:,:,:)
  call grady_wgt_2nd(nx,ny,nz,-e_p,lat,coslat(:,1), tmp)
  call grady_wgt_2nd(nx,ny,nz,f2e_k,lat,coslat(:,1), wadm)   ! rent wadm
  u_ar(:,:,:) = u_ar(:,:,:) +                                            &
                ( wadm(:,:,:)/(te*te) + tmp(:,:,:) )/te(:,:,:)
  call gradz_wgt_2nd_irr(nx,ny,nz,-fvt,zp,rho0(1,:), tmp)
  v_ar(:,:,:) = v_ar(:,:,:) + tmp(:,:,:)/te(:,:,:)
  call gradz_wgt_2nd_irr(nx,ny,nz,-fut,zp,rho0(1,:), tmp)
  u_ar(:,:,:) = u_ar(:,:,:) + tmp(:,:,:)/te(:,:,:)

  ! B_2, B_1, mass-weighted
  tmp(:,:,:) = spread(rho0(:,:)*coslat(:,:),1,nx)
  te(:,:,:) = spread(f(:,:)*f(:,:),1,nx)   ! rent te (discarded)
  mx(:,:,:) = tmp(:,:,:)*( (v2f2(:,:,:) - f2e_k(:,:,:))/te(:,:,:) -      &
                           e_p(:,:,:) )
  lx(:,:,:) = tmp(:,:,:)*vuf2(:,:,:)/te(:,:,:)
  mz(:,:,:) = tmp(:,:,:)*fvt(:,:,:)
  lz(:,:,:) = tmp(:,:,:)*fut(:,:,:)
  my(:,:,:) = 0.
  ly(:,:,:) = tmp(:,:,:)*e_p(:,:,:)

  ! div[B_2], div[B_1]  [m/s/day * kg/m3]
  call div_flx(nx,ny,nz,lat,dlon,missv,mx,my,mz,                         &
               divm_x,divm_y,divm_z,divm)
  call div_flx(nx,ny,nz,lat,dlon,missv,lx,ly,lz,                         &
               divf_x,divf_y,divf_z,divf)

  ! correction of B_2_y, B_1_y and their div.
  wadm(:,:,:) = -tmp(:,:,:)*vuf2(:,:,:)   ! rent wadm
  call grady_wgt_2nd(nx,ny,nz,wadm,lat,coslat(:,1), dqh)   ! rent dqh
  dqh(:,:,:) = dqh(:,:,:)/te(:,:,:) * 86400.   ! te : f^2
  divm_y(:,:,:) = divm_y(:,:,:) + dqh(:,:,:)
  divm  (:,:,:) = divm  (:,:,:) + dqh(:,:,:)
  my(:,:,:) = my(:,:,:) + wadm(:,:,:)/te(:,:,:)
  wadm(:,:,:) = tmp(:,:,:)*(f2e_k(:,:,:) - u2f2(:,:,:))
  call grady_wgt_2nd(nx,ny,nz,wadm,lat,coslat(:,1), dqh)
  dqh(:,:,:) = dqh(:,:,:)/te(:,:,:) * 86400.
  divf_y(:,:,:) = divf_y(:,:,:) + dqh(:,:,:)
  divf  (:,:,:) = divf  (:,:,:) + dqh(:,:,:)
  ly(:,:,:) = ly(:,:,:) + wadm(:,:,:)/te(:,:,:)

  ! background wind : smooth (us_g, vs_g)
  call lowpass_k(us_g,kb_max)
  call lowpass_k(vs_g,kb_max)
  call filter121_y(us_g)
  call filter121_y(vs_g)

  ! grad[Q] (normalized or not)
!n  dqh(:,:,:) = sqrt( dqx(:,:,:)*dqx(:,:,:) + dqy(:,:,:)*dqy(:,:,:) )
  dqh(:,:,:) = 1.  !n  ! for A and M not to be normalized by |del[Q]|

  dqy(:,:,:) = dqy(:,:,:)/dqh(:,:,:)
  dqx(:,:,:) = dqx(:,:,:)/dqh(:,:,:)
!tmp=sqrt(us_g*us_g + vs_g*vs_g)
!dqy=us_g/tmp
!dqx=-vs_g/tmp

  ! M_R and div[M_R]  [m/s/day * kg/m3]
  divm_x(:,:,:) = dqy(:,:,:)*divm_x(:,:,:) + dqx(:,:,:)*divf_x(:,:,:)
  divm_y(:,:,:) = dqy(:,:,:)*divm_y(:,:,:) + dqx(:,:,:)*divf_y(:,:,:)
  divm_z(:,:,:) = dqy(:,:,:)*divm_z(:,:,:) + dqx(:,:,:)*divf_z(:,:,:)
  divm  (:,:,:) = dqy(:,:,:)*divm  (:,:,:) + dqx(:,:,:)*divf  (:,:,:)
  mx    (:,:,:) = dqy(:,:,:)*mx    (:,:,:) + dqx(:,:,:)*lx    (:,:,:)
  my    (:,:,:) = dqy(:,:,:)*my    (:,:,:) + dqx(:,:,:)*ly    (:,:,:)
  mz    (:,:,:) = dqy(:,:,:)*mz    (:,:,:) + dqx(:,:,:)*lz    (:,:,:)

  ! A : wave activity density (referred to M in P86)
  tmp(:,:,:) = spread(rho0(:,:)*coslat(:,:),1,nx)
  wadm(:,:,:) = 0.5*(q_prt(:,:,:)*q_prt(:,:,:))/dqh(:,:,:)*tmp(:,:,:)

  ! M_T ( = M_R + U_g * A )
  mx(:,:,:) = mx(:,:,:) + us_g(:,:,:)*wadm(:,:,:)
  my(:,:,:) = my(:,:,:) + vs_g(:,:,:)*wadm(:,:,:)

  ! div[M_T] ( = div[M_R] + A advection )
  tmp(:,:,:) = spread(rho0(:,:)*coslat(:,:),1,nx)*                       &
               (q_prt(:,:,:)/dqh(:,:,:))
  divm_x(:,:,:) = divm_x(:,:,:) + us_g(:,:,:)*dqx_prt(:,:,:)*tmp(:,:,:)* &
                                  86400.
  divm_y(:,:,:) = divm_y(:,:,:) + vs_g(:,:,:)*dqy_prt(:,:,:)*tmp(:,:,:)* &
                                  86400.
  divm(:,:,:) = divm_x(:,:,:) + divm_y(:,:,:) + divm_z(:,:,:)

!+test :  normalized as in Plumb 86
  dqh(:,:,:) = sqrt( dqx(:,:,:)*dqx(:,:,:) + dqy(:,:,:)*dqy(:,:,:) )
  lx = mx/dqh  ;  ly = my/dqh  ;  lz = mz/dqh
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
                   div_x,div_y,div_z,div)   ! [m/s/day * kg/m3] (rho*cos)

  integer               , intent(in) ::  nx, ny, nz
  real                  , intent(in) ::  dlon, missv
  real, dimension(:)    , intent(in) ::  lat
  real, dimension(:,:,:), intent(in) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)

  call grady_wgt_2nd(nx,ny,nz,fy,lat,coslat(:,1), div_y)
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
                         div_x,div_y,div_z,div)   ! [m/s/day]

  integer               , intent(in) ::  nx, ny, nz
  real                  , intent(in) ::  dlon, missv
  real, dimension(:)    , intent(in) ::  lat
  real, dimension(:,:,:), intent(in) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)

  div(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div

  call grady_wgt_2nd(nx,ny,nz,fy*div,lat,coslat(:,1), div_y)
  div_y(:,:,:) = div_y(:,:,:)/div(:,:,:)
  call missing_bdy(nx,ny,nz,div_y,missv,1,1,0,0)

  call gradz_wgt_2nd_irr(nx,ny,nz,fz*div,zp,rho0(1,:), div_z)
  call missing_bdy(nx,ny,nz,div_z,missv,0,0,1,1)

  div_x(:,:,:) = div_x(:,:,:)*86400.
  div_y(:,:,:) = div_y(:,:,:)*86400.
  div_z(:,:,:) = div_z(:,:,:)*86400.

  div(:,:,:) = div_x(:,:,:) + div_y(:,:,:) + div_z(:,:,:)
  call missing_bdy(nx,ny,nz,div,missv,1,1,1,1)
 
END subroutine div_norm_flx0

SUBROUTINE div_norm_flx0_wgt(nx,ny,nz,lat,dlon,missv,fx,fy,fz,           &
                             div_x,div_y,div_z,div)   ! [m/s/day]

  integer           , intent(in) ::  nx, ny, nz
  real              , intent(in) ::  dlon, missv
  real, dimension(:), intent(in) ::  lat

  real, dimension(:,:,:), intent(inout) ::  fx, fy, fz

  real, dimension(nx,ny,nz), intent(out) ::  div, div_x, div_y, div_z

  integer ::  j,k

  call gradx_2nd(nx,ny,nz,fx,lat,dlon, div_x)

  div(:,:,:) = spread(coslat(:,:),1,nx)   ! rent div

  fy(:,:,:) = fy(:,:,:)*div(:,:,:)

  call grady_wgt_2nd(nx,ny,nz,fy(:,:,:),lat,coslat(:,1), div_y)
  div_y(:,:,:) = div_y(:,:,:)/div(:,:,:)
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

SUBROUTINE geostrophic_wind_f(nx,ny,nz,lat,dlon,gp,ug,vg)

  integer                  , intent(in) ::  nx, ny, nz
  real                     , intent(in) ::  dlon
  real, dimension(ny)      , intent(in) ::  lat
  real, dimension(nx,ny,nz), intent(in) ::  gp

  real, dimension(nx,ny,nz), intent(out) ::  ug, vg

  real, dimension(nx,ny,nz) ::  tmp

  call gradx_2nd(nx,ny,nz,gp,lat,dlon, vg)
  call grady_2nd(nx,ny,nz,-gp,lat, ug)

END subroutine geostrophic_wind_f

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


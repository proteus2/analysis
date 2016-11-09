MODULE warc3d

  implicit none

  private ::  set_gridvar_p, set_gridvar_z, set_gridvar_z_cp
  private ::  grady_2nd, gradz_2nd_irr, gradz_i2m, gradz_m2i,            &
              gradzwav_m2i, gradzprt_m2i
  private ::  intp_i2m, intp_m2i, intpwav_i2m, intpwav_m2i, intpprt_i2m, &
              intpprt_m2i
  private ::  missing_bdy

  character(len=64), dimension(11) ::  varname_tem_p =                   &
      (/'v_res  ','w_res  ','cor    ','uadv_y ','uadv_z ',               &
        'epd    ','epd_y  ','epd_z  ','f_y    ','f_z    ',               &
        'u_force'/)
  character(len=64), dimension(10) ::  varname_tem_qg =                  &
      (/'v_res  ','w_res  ','cor    ','epd    ','epd_y  ',               &
        'epd_z  ','f_y    ','f_z    ','u_force','tadv_z '/)
  character(len=64), dimension(12) ::  varname_tem_z =                   &
      (/'v_res  ','w_res  ','cor    ','uadv_y ','uadv_z ',               &
        'epd    ','epd_y  ','epd_z  ','f_y    ','f_z    ',               &
        'u_force','ru_eddy'/)

  real, dimension(:,:), allocatable ::  dptmdz_from_thlev

  integer,                           private ::  ny_pre, nz_pre
  real, dimension(:),   allocatable, private ::  lat_pre, ht_pre, zp
  real, dimension(:,:), allocatable, private ::  coslat, f, rho0, rbeta, &
                                                 r_earth, rcos, rrhocos, &
                                                 r3, r3rho, r_t, r2_t,   &
                                                 r3_t

  real, parameter, private ::  g = 9.80665, rd = 287.05, cp = 1005.
  real, parameter, private ::  kappa = rd/cp
  real, parameter, private ::  a_earth = 6371229.
  real, parameter, private ::  pi = 3.14159265358979323846
  real, parameter, private ::  ome_earth = 7.292116e-5
  real, parameter, private ::  deg2rad = pi/180.

  interface gradz_i2m
    module procedure gradz_i2m_2d, gradz_i2m_3d
  end interface


  CONTAINS


SUBROUTINE tem_hydro_p(                                                  &
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

  call set_gridvar_p(ny,nz,lat,p,h_scale)

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
  call gradz_2nd_irr(ny,nz,temp,zp, dptmdz)
  dptmdz(:,:) = ptm(:,:)*dptmdz(:,:)

  ! grad, u
  call grady_2nd(ny,nz,um*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um,zp, dumdz)

  ! phi
  phi(:,:) = rvpt(:,:)/dptmdz(:,:)
!!  call missing_bdy(ny,nz,phi,missv,0,0,1,0)

  ! residual mean meridional circulation
  call gradz_2nd_irr(ny,nz,phi,zp, temp)
  vres(:,:) = vres(:,:) - temp(:,:)/rho0(:,:)
  call grady_2nd(ny,nz,phi*coslat,lat, temp)
  do k=1, nz
  do j=2, ny-1
    wres(j,k) = wres(j,k) + temp(j,k)/rho0(j,k)/coslat(j,k)
  enddo
  enddo
  call missing_bdy(ny,nz,vres,missv,0,0,1,0)
  call missing_bdy(ny,nz,wres,missv,1,1,1,0)

  ! adv
  uadvy(:,:) = -vres(:,:)*divy_um(:,:) * 86400.
  uadvz(:,:) = -wres(:,:)*dumdz  (:,:) * 86400.
  call missing_bdy(ny,nz,uadvy,missv,1,1,1,0)
  call missing_bdy(ny,nz,uadvz,missv,1,1,1,0)

  ! coriolis term
  cor(:,:) = vres(:,:)*f(:,:) * 86400.
  call missing_bdy(ny,nz,cor,missv,0,0,1,0)

  ! epf, epd
  fy(:,:) = -a_earth*coslat(:,:)*(rvu(:,:) - phi(:,:)*dumdz(:,:))
  fz(:,:) = -a_earth*coslat(:,:)*(rwu(:,:) + phi(:,:)*(divy_um(:,:)-f(:,:)))
  call missing_bdy(ny,nz,fy,missv,0,0,1,0)
!mv  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y,missv,1,1,1,0)

  call gradz_2nd_irr(ny,nz,fz,zp, epd_z)
  epd_z(:,:) = epd_z(:,:)/(rho0(:,:)*a_earth*coslat(:,:)) * 86400.
  call missing_bdy(ny,nz,epd_z,missv,1,1,1,0)

  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,1,1,1,0)

  ! u_force
  u_force(:,:) = cor(:,:) + uadvy(:,:) + uadvz(:,:) + epd(:,:)
  call missing_bdy(ny,nz,u_force,missv,1,1,1,0)

END subroutine tem_hydro_p

SUBROUTINE warc3d_qg(                                                    &
     nx,ny,nz,lat,p,u,v,pt,wm,h_scale,missv,                             &
     vres,wres,cor,epd,epd_y,epd_z,fy,fz,u_force,tadvz )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(ny,nz),    intent(in) ::  wm
  real, dimension(nx,ny,nz), intent(in) ::  u, v, pt

  real, dimension(ny,nz), intent(out) ::  vres, wres, fy, fz
  real, dimension(ny,nz), intent(out) ::  cor, epd, epd_y, epd_z
  real, dimension(ny,nz), intent(out) ::  u_force
  real, dimension(ny,nz), intent(out) ::  tadvz

  real, dimension(ny,nz) ::  um, ptm
  real, dimension(ny,nz) ::  phi, rvpt, rvu
  real, dimension(ny,nz) ::  dptmdy, dptmdz
  real, dimension(ny,nz) ::  temp

  real, dimension(:,:,:), allocatable ::  prt, v_prt

  integer ::  j,k

  call set_gridvar_p(ny,nz,lat,p,h_scale)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), v_prt(nx,ny,nz) )

  um  (:,:) = sum(u , dim=1)/float(nx)
  vres(:,:) = sum(v , dim=1)/float(nx)
  ptm (:,:) = sum(pt, dim=1)/float(nx)

  v_prt(:,:,:) = v(:,:,:) - spread(vres(:,:),1,nx)

  prt(:,:,:) = pt(:,:,:) - spread(ptm(:,:),1,nx)

  rvpt(:,:) = sum(v_prt*prt, dim=1)/float(nx)

  prt(:,:,:) = u(:,:,:) - spread(um(:,:),1,nx)

  rvu(:,:) = sum(v_prt*prt, dim=1)/float(nx)

  deallocate( prt, v_prt )

  rvpt(:,:) = rho0(:,:)*rvpt(:,:)
  rvu (:,:) = rho0(:,:)*rvu (:,:)

  ! grad, pt
!  call grady_2nd(ny,nz,ptm,lat, dptmdy)

  temp(:,:) = log(ptm(:,:))
  call gradz_2nd_irr(ny,nz,temp,zp, dptmdz)
  dptmdz(:,:) = ptm(:,:)*dptmdz(:,:)

  ! phi
  phi(:,:) = rvpt(:,:)/dptmdz(:,:)
!!  call missing_bdy(ny,nz,phi,missv,0,0,1,0)

  ! residual mean meridional circulation
  call gradz_2nd_irr(ny,nz,phi,zp, temp)
  vres(:,:) = vres(:,:) - temp(:,:)/rho0(:,:)
  call grady_2nd(ny,nz,phi*coslat,lat, temp)
  do k=1, nz
  do j=2, ny-1
    wres(j,k) = wm(j,k) + temp(j,k)/rho0(j,k)/coslat(j,k)
  enddo
  enddo
  call missing_bdy(ny,nz,vres,missv,0,0,1,0)
  call missing_bdy(ny,nz,wres,missv,1,1,1,0)

  ! coriolis term
  cor(:,:) = vres(:,:)*f(:,:) * 86400.
  call missing_bdy(ny,nz,cor,missv,0,0,1,0)

  ! epf, epd
  temp(:,:) = phi(:,:)*f(:,:)

  fy(:,:) = -a_earth*coslat(:,:)*rvu(:,:)
  fz(:,:) = a_earth*coslat(:,:)*temp(:,:)
  call missing_bdy(ny,nz,fy,missv,0,0,0,0)
  call missing_bdy(ny,nz,fz,missv,0,0,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y,missv,1,1,0,0)

  call gradz_2nd_irr(ny,nz,temp,zp, epd_z)
  epd_z(:,:) = epd_z(:,:)/rho0(:,:) * 86400.
  call missing_bdy(ny,nz,epd_z,missv,0,0,1,0)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,1,1,1,0)

  ! u_force
  u_force(:,:) = cor(:,:) + epd(:,:)
  call missing_bdy(ny,nz,u_force,missv,1,1,1,0)

  ! T_force
  tadvz(:,:) = -wres(:,:)*dptmdz(:,:) * 86400.
  call missing_bdy(ny,nz,tadvz,missv,1,1,1,0)

END subroutine warc3d_qg

SUBROUTINE tem_qg_betap_gp(                                              &
     nx,ny,nz,lat,p,gp,vm,wm,lat0,h_scale,missv,                         &
     vres,wres,cor,epd,epd_y,epd_z,fy,fz,u_force,tadvz )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  lat0, h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(ny,nz),    intent(in) ::  vm, wm
  real, dimension(nx,ny,nz), intent(in) ::  gp

  real, dimension(ny,nz), intent(out) ::  vres, wres, fy, fz
  real, dimension(ny,nz), intent(out) ::  cor, epd, epd_y, epd_z
  real, dimension(ny,nz), intent(out) ::  u_force
  real, dimension(ny,nz), intent(out) ::  tadvz

  real, dimension(ny,nz) ::  phi, rvpt, rvu
  real, dimension(ny,nz) ::  gpm, dpt0dz
  real, dimension(ny,nz) ::  temp

  real, dimension(:,:,:), allocatable ::  gpp, prt, v_prt

  integer ::  i,j,k, j0
  real    ::  f0, twodx0, twody0(ny), coef(3)
  real    ::  grad1(nz), grad2(nz), grad3(nz), coefz(nz)

  call set_gridvar_betap(ny,nz,lat,p,lat0,h_scale)
  f0 = 2.*ome_earth*sin(lat0*deg2rad)  ! beta-plane approx.

  gpm(:,:) = sum(gp, dim=1)/float(nx)

  ! obtain QG-wind and T perturbations

  allocate( gpp(nx,ny,nz), prt(nx,ny,nz), v_prt(nx,ny,nz) )

  gpp(:,:,:) = gp(:,:,:) - spread(gpm(:,:),1,nx)

  twodx0 = 2.*a_earth*cos(lat0*deg2rad)*(2.*pi/float(nx))
  twody0(2:ny-1) = a_earth*((lat(3:ny)-lat(1:ny-2))*deg2rad)

  ! v'
  v_prt(2:nx-1,:,:) = (gp(3:nx,:,:) - gp(1:nx-2,:,:))/twodx0
  v_prt(1     ,:,:) = (gp(2   ,:,:) - gp(nx    ,:,:))/twodx0
  v_prt(nx    ,:,:) = (gp(1   ,:,:) - gp(nx-1  ,:,:))/twodx0
  v_prt(:,:,:) = v_prt(:,:,:)/f0

  ! u'
  prt = 0.
  prt(:,2:ny-1,:) = (gpp(:,3:ny,:) - gpp(:,1:ny-2,:))/                   &
                    spread(spread(twody0(2:ny-1),1,nx),3,nz)
  prt(:,:,:) = (-1.)*prt(:,:,:)/f0

  rvu(:,:) = sum(v_prt*prt, dim=1)/float(nx)

  ! pt'
  prt = 0.
  do k=2, nz-1
!    call fdcoef(1,2,zp(k),zp(k-1:k+1), coef(:))
    call fdcoef_1d2o(zp(k-1:k+1), coef(:))
    prt(:,:,k) = coef(1)*gpp(:,:,k-1) + &
                 coef(2)*gpp(:,:,k  ) + &
                 coef(3)*gpp(:,:,k+1)
  enddo
  prt(:,:,1 ) = (gpp(:,:,2 ) - gpp(:,:,1   ))/(zp(2 ) - zp(1   ))
  prt(:,:,nz) = (gpp(:,:,nz) - gpp(:,:,nz-1))/(zp(nz) - zp(nz-1))
  coefz(:) = h_scale/rd*exp(kappa*zp(:)/h_scale)
  do k=1, nz
    prt(:,:,k) = prt(:,:,k)*coefz(k)
  enddo

  rvpt(:,:) = sum(v_prt*prt, dim=1)/float(nx)

  deallocate( gpp, prt, v_prt )

  rvu (:,:) = rho0(:,:)*rvu (:,:)
  rvpt(:,:) = rho0(:,:)*rvpt(:,:)
  call missing_bdy(ny,nz,rvu ,missv,1,1,0,0)
  call missing_bdy(ny,nz,rvpt,missv,0,0,1,0)

  ! grad, pt0
  j0 = minloc(abs(lat(:)-lat0),1)
  call gradz_2nd_irr(1,nz,gpm(j0,:),zp, grad1)
  call gradz_2nd_irr(1,nz,grad1,zp, grad2)
  call gradz_2nd_irr(1,nz,grad2,zp, grad3)
  dpt0dz(:,:) = spread( coefz(:)*(kappa/h_scale*grad1(:) + grad2(:)),    &
                        1,ny )
  call missing_bdy(ny,nz,dpt0dz,missv,0,0,2,1)

  ! phi
  phi(:,:) = rvpt(:,:)/dpt0dz(:,:)
  call missing_bdy(ny,nz,phi,missv,0,0,2,1)

  ! residual mean meridional circulation
  call gradz_2nd_irr(ny,nz,phi,zp, temp)
  vres(:,:) = vm(:,:) - temp(:,:)/rho0(:,:)
  call grady_2nd(ny,nz,phi,lat, temp)
  do k=1, nz
  do j=2, ny-1
    wres(j,k) = wm(j,k) + temp(j,k)/rho0(j,k)
  enddo
  enddo
  call missing_bdy(ny,nz,vres,missv,0,0,3,2)
  call missing_bdy(ny,nz,wres,missv,1,1,2,1)

  ! coriolis term
  cor(:,:) = vres(:,:)*f0 * 86400.
  call missing_bdy(ny,nz,cor,missv,0,0,3,2)

  ! epf, epd
  fy(:,:) = -rvu(:,:)
  fz(:,:) = phi(:,:)*f0
  call missing_bdy(ny,nz,fy,missv,1,1,0,0)
  call missing_bdy(ny,nz,fz,missv,1,1,2,1)

  call grady_2nd(ny,nz,fy,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/rho0(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y,missv,2,2,0,0)

  epd_z(:,:) = f0*temp(:,:)/rho0(:,:) * 86400.
  call missing_bdy(ny,nz,epd_z,missv,1,1,3,2)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,3,2)

  ! u_force
  u_force(:,:) = cor(:,:) + epd(:,:)
  call missing_bdy(ny,nz,u_force,missv,2,2,3,2)

  ! T_force
  tadvz(:,:) = -wres(:,:)*dpt0dz(:,:) * 86400.
  call missing_bdy(ny,nz,tadvz,missv,1,1,2,1)

END subroutine tem_qg_betap_gp

SUBROUTINE tem_z(                                                        &
     nx,ny,nz,lat,z,u,v,w,pt,rho,missv,                                  &
     vres,wres,cor,uadvy,uadvz,epd,epd_y,epd_z,fy,fz,u_force,            &
     rueddy )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  z
  real, dimension(nx,ny,nz), intent(in) ::  u, v, w, pt, rho
  
  real, dimension(ny,nz), intent(out) ::  vres, wres, fy, fz
  real, dimension(ny,nz), intent(out) ::  cor, uadvy, uadvz, epd,        &
                                          epd_y, epd_z, u_force
  real, dimension(ny,nz), intent(out) ::  rueddy
  
  real, dimension(ny,nz) ::  um, ptm, rhom, rvm, rwm
  real, dimension(ny,nz) ::  phi, rvpt, rwpt, rvu, rwu
  real, dimension(ny,nz) ::  dptmdy, dptmdz, divy_um, dumdz
  real, dimension(ny,nz) ::  temp, temp2
  
  real, dimension(:,:,:), allocatable ::  prt, rv_prt, rw_prt
    
  integer ::  j,k
  
  call set_gridvar_z(ny,nz,lat,z)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), rv_prt(nx,ny,nz), rw_prt(nx,ny,nz) )

  rv_prt(:,:,:) = rho(:,:,:)*v(:,:,:)
  rw_prt(:,:,:) = rho(:,:,:)*w(:,:,:)

  um  (:,:) = sum(u     , dim=1)/float(nx)
  rvm (:,:) = sum(rv_prt, dim=1)/float(nx)
  rwm (:,:) = sum(rw_prt, dim=1)/float(nx)
  ptm (:,:) = sum(pt    , dim=1)/float(nx)
  rhom(:,:) = sum(rho   , dim=1)/float(nx)

  rv_prt(:,:,:) = rv_prt(:,:,:) - spread(rvm(:,:),1,nx)
  rw_prt(:,:,:) = rw_prt(:,:,:) - spread(rwm(:,:),1,nx)

  prt(:,:,:) = pt(:,:,:) - spread(ptm(:,:),1,nx)

  rvpt(:,:) = sum(rv_prt*prt, dim=1)/float(nx)
  rwpt(:,:) = sum(rw_prt*prt, dim=1)/float(nx)

  prt(:,:,:) = u(:,:,:) - spread(um(:,:),1,nx)

  rvu(:,:) = sum(rv_prt*prt, dim=1)/float(nx)
  rwu(:,:) = sum(rw_prt*prt, dim=1)/float(nx)

  prt(:,:,:) = prt(:,:,:)*(rho(:,:,:)-spread(rhom(:,:),1,nx))

  ! -(rho'u')_bar / rho_bar
  rueddy(:,:) = (-1.)*sum(prt, dim=1)/float(nx)/rhom(:,:)

  deallocate( prt, rv_prt, rw_prt )

  ! grad, pt
  call grady_2nd(ny,nz,ptm,lat, dptmdy)

  if ( allocated(dptmdz_from_thlev) ) then
    dptmdz(:,:) = dptmdz_from_thlev(:,:)
  else
    temp(:,:) = log(ptm(:,:))
    call gradz_2nd_irr(ny,nz,temp,z, dptmdz)
    dptmdz(:,:) = ptm(:,:)*dptmdz(:,:)
  end if

  ! grad, u
  call grady_2nd(ny,nz,um*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um*r_earth,z, dumdz)
  dumdz(:,:) = dumdz(:,:)/r_earth(:,:)

  ! phi
  phi(:,:) = (rvpt(:,:)*dptmdz(:,:) - rwpt(:,:)*dptmdy(:,:))/            &
             (dptmdy(:,:)*dptmdy(:,:) + dptmdz(:,:)*dptmdz(:,:))
!  call missing_bdy(ny,nz,phi,missv,1,1,1,0)

  ! residual mean meridional circulation
  temp(:,:) = rcos(:,:)*phi(:,:)
  call gradz_2nd_irr(ny,nz,temp,z, temp2)
  do k=1, nz
  do j=2, ny-1
    vres(j,k) = (rvm(j,k) - temp2(j,k)/rcos(j,k))/rhom(j,k)
  enddo
  enddo
  call grady_2nd(ny,nz,temp,lat, temp2)
  do k=1, nz
  do j=2, ny-1
    wres(j,k) = (rwm(j,k) + temp2(j,k)/rcos(j,k))/rhom(j,k)
  enddo
  enddo
  call missing_bdy(ny,nz,vres,missv,1,1,2,1)
  call missing_bdy(ny,nz,wres,missv,2,2,1,0)

  ! adv
  uadvy(:,:) = -vres(:,:)*divy_um(:,:) * 86400.
  uadvz(:,:) = -wres(:,:)*dumdz  (:,:) * 86400.
  call missing_bdy(ny,nz,uadvy,missv,1,1,2,1)
  call missing_bdy(ny,nz,uadvz,missv,2,2,1,0)

  ! coriolis term
  cor(:,:) = (vres(:,:)*f(:,:) - wres(:,:)*rbeta(:,:)) * 86400.
  call missing_bdy(ny,nz,cor,missv,2,2,2,1)

  ! epf, epd
  fy(:,:) = -rcos(:,:)*(rvu(:,:) - phi(:,:)*(dumdz  (:,:)+rbeta(:,:)))
  fz(:,:) = -rcos(:,:)*(rwu(:,:) + phi(:,:)*(divy_um(:,:)-f    (:,:)))
  call missing_bdy(ny,nz,fy,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  rrhocos(:,:) = rcos(:,:)*rhom(:,:)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y,missv,2,2,1,0)

  temp(:,:) = r_earth(:,:)*r_earth(:,:)
  call gradz_2nd_irr(ny,nz,fz*temp,z, temp2)
  do k=1, nz
  do j=2, ny-1
    epd_z(j,k) = temp2(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_z,missv,1,1,2,1)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,2,1)

  ! u_force 
  u_force(:,:) = cor(:,:) + uadvy(:,:) + uadvz(:,:) + epd(:,:)

  call missing_bdy(ny,nz,u_force,missv,2,2,2,1)

END subroutine tem_z

SUBROUTINE tem_z_cp(                                                     &
     nx,ny,nz,lat,z,zt,u,v,w,pt,rho,missv,                               &
     vres,wres,cor,uadvy,uadvz,epd,epd_y,epd_z,fy,fz,u_force,            &
     rueddy )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  z, zt
  real, dimension(nx,ny,nz), intent(in) ::  u, v, w, pt, rho
 
  real, dimension(ny,nz), intent(out) ::  vres, wres, fy, fz
  real, dimension(ny,nz), intent(out) ::  cor, uadvy, uadvz, epd,        &
                                          epd_y, epd_z, u_force
  real, dimension(ny,nz), intent(out) ::  rueddy
 
  real, dimension(ny,nz) ::  um, ptm, rhom, rvm, rwm
  real, dimension(ny,nz) ::  rvpt, rwpt, rvu, rwu, druwr3dz
  real, dimension(ny,nz) ::  phi, fphiy, fphiz, fphiz_t,                 &
                             dptmdy, dptmdz, divy_um, divy_um_t, dumdz
  real, dimension(ny,nz) ::  temp, temp2
  
  real, dimension(:,:,:), allocatable ::  prt, prt2, rv_prt, rw_prt
    
  integer ::  j,k
  
  call set_gridvar_z_cp(ny,nz,lat,z,zt)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), rv_prt(nx,ny,nz), rw_prt(nx,ny,nz) )
  allocate( prt2(nx,ny,nz) )

  rw_prt(:,:,:) = log(rho(:,:,:))
  prt(:,:,1:nz-1) = 0.5*(rw_prt(:,:,1:nz-1)+rw_prt(:,:,2:nz))
  prt(:,:,nz) = rw_prt(:,:,nz)
  prt(:,:,:) = exp(prt(:,:,:))

  rv_prt(:,:,:) = rho(:,:,:)*v(:,:,:)
  rw_prt(:,:,:) = prt(:,:,:)*w(:,:,:)

  um  (:,:) = sum(u     , dim=1)/float(nx)
  rvm (:,:) = sum(rv_prt, dim=1)/float(nx)
  rwm (:,:) = sum(rw_prt, dim=1)/float(nx)
  ptm (:,:) = sum(pt    , dim=1)/float(nx)
  rhom(:,:) = sum(rho   , dim=1)/float(nx)

  rv_prt(:,:,:) = rv_prt(:,:,:) - spread(rvm(:,:),1,nx)
  rw_prt(:,:,:) = rw_prt(:,:,:) - spread(rwm(:,:),1,nx)

  prt(:,:,:) = pt(:,:,:) - spread(ptm(:,:),1,nx)

  prt2(:,:,:) = rv_prt(:,:,:)
  call intpprt_m2i(nx,ny,nz, prt2)
  rvpt(:,:) = sum(prt2*prt, dim=1)/float(nx)

  rwpt(:,:) = sum(rw_prt*prt, dim=1)/float(nx)

  prt(:,:,:) = u(:,:,:) - spread(um(:,:),1,nx)

  rvu(:,:) = sum(rv_prt*prt, dim=1)/float(nx)

  prt2(:,:,:) = rw_prt(:,:,:)
  call intpprt_i2m(nx,ny,nz, prt2)
  rwu(:,:) = sum(prt2*prt, dim=1)/float(nx)

  prt2(:,:,:) = prt(:,:,:)*(rho(:,:,:)-spread(rhom(:,:),1,nx))

  ! -(rho'u')_bar / rho_bar
  rueddy(:,:) = (-1.)*sum(prt2, dim=1)/float(nx)/rhom(:,:)

  prt   (:,:,:) = prt   (:,:,:)*spread(r_earth,1,nx)  ! u*r
  rw_prt(:,:,:) = rw_prt(:,:,:)*spread(r2_t   ,1,nx)  ! rw*r2
  call gradzprt_m2i(nx,ny,nz,prt,z, prt2)
  druwr3dz(:,:) = sum(rw_prt*prt2, dim=1)
  call intp_i2m(ny,nz, druwr3dz)
  call gradz_i2m(nx,ny,nz,rw_prt,zt, prt2)
  druwr3dz(:,:) = (druwr3dz(:,:) + sum(prt*prt2, dim=1))/float(nx)

  deallocate( prt2 )
  deallocate( prt, rv_prt, rw_prt )

  ! grad, pt
  call grady_2nd(ny,nz,ptm,lat, dptmdy)

  temp(:,:) = log(ptm(:,:))
  call gradz_2nd_irr(ny,nz,temp,zt, dptmdz)
  dptmdz(:,:) = ptm(:,:)*dptmdz(:,:)

  ! grad, u
  call grady_2nd(ny,nz,um*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  divy_um_t(:,:) = divy_um(:,:)*r_earth(:,:)
  call intp_m2i(ny,nz,z,zt, divy_um_t)
  divy_um_t(:,:) = divy_um_t(:,:)/r_t(:,:)
  call gradz_m2i(ny,nz,um*r_earth,z, dumdz)
  dumdz(:,:) = dumdz(:,:)/r_t(:,:)

  ! phi
  phi(:,:) = (rvpt(:,:)*dptmdz(:,:) - rwpt(:,:)*dptmdy(:,:))/            &
             (dptmdy(:,:)*dptmdy(:,:) + dptmdz(:,:)*dptmdz(:,:))
  fphiy(:,:) = phi(:,:)*(dumdz    (:,:)+rbeta(:,:))
  fphiz(:,:) = phi(:,:)*(divy_um_t(:,:)-f    (:,:))
  fphiz_t(:,:) = fphiz(:,:)
  call intp_i2m(ny,nz, fphiy)
  call intp_i2m(ny,nz, fphiz)
!  call missing_bdy(ny,nz,phi  ,missv,1,1,1,0)
!  call missing_bdy(ny,nz,fphiy,missv,1,1,1,0)
!  call missing_bdy(ny,nz,fphiz,missv,1,1,1,0)

  ! residual mean meridional circulation
  temp(:,:) = r_t(:,:)*coslat(:,:)*phi(:,:)
  call gradz_i2m(ny,nz,temp,zt , temp2)
  do k=1, nz
  do j=2, ny-1
    vres(j,k) = (rvm(j,k) - temp2(j,k)/rcos(j,k))/rhom(j,k)
  enddo
  enddo
  call grady_2nd(ny,nz,temp,lat, temp2)
  do k=1, nz
  do j=2, ny-1
    rwm(j,k) = rwm(j,k) + temp2(j,k)/(r_t(j,k)*coslat(j,k))
  enddo
  enddo
  call missing_bdy(ny,nz,vres,missv,1,1,1,0)
 
  ! adv
  uadvy(:,:) = -vres(:,:)*divy_um(:,:) * 86400.
  uadvz(:,:) = -rwm(:,:)*dumdz(:,:)
  call intp_i2m(ny,nz, uadvz)
  uadvz(:,:) = uadvz(:,:)/rhom(:,:) * 86400.
  call missing_bdy(ny,nz,uadvy,missv,1,1,1,0)
  call missing_bdy(ny,nz,uadvz,missv,2,2,1,0)

  call intp_i2m(ny,nz, rwm)
  wres(:,:) = rwm(:,:)/rhom(:,:)
  call missing_bdy(ny,nz,wres,missv,2,2,1,0)

  ! coriolis term
  cor(:,:) = (vres(:,:)*f(:,:) - wres(:,:)*rbeta(:,:)) * 86400.
  call missing_bdy(ny,nz,cor,missv,2,2,1,0)

  ! epf, epd
  fy(:,:) = -rcos(:,:)*(rvu(:,:) - fphiy(:,:))
  fz(:,:) = -rcos(:,:)*(rwu(:,:) + fphiz(:,:))
  call missing_bdy(ny,nz,fy,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  rrhocos(:,:) = rcos(:,:)*rhom(:,:)
  r3rho  (:,:) = r3  (:,:)*rhom(:,:)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y,missv,2,2,1,0)

  call gradz_i2m(ny,nz,r3_t*fphiz_t,zt, temp)
  epd_z(:,:) = -(druwr3dz(:,:) + temp(:,:))/r3rho(:,:) * 86400.
  call missing_bdy(ny,nz,epd_z,missv,1,1,2,1)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,2,1)

  ! u_force 
  u_force(:,:) = cor(:,:) + uadvy(:,:) + uadvz(:,:) + epd(:,:)

  call missing_bdy(ny,nz,u_force,missv,2,2,2,1)

END subroutine tem_z_cp

SUBROUTINE set_gridvar_p(ny,nz,lat,p,h_scale)

  integer,             intent(in) ::  ny, nz
  real,                intent(in) ::  h_scale
  real, dimension(ny), intent(in) ::  lat
  real, dimension(nz), intent(in) ::  p

  if ( allocated(lat_pre) ) then
    if ( ny == ny_pre .and. nz == nz_pre ) then
      if ( all(lat == lat_pre) .and. all(p == ht_pre) )  RETURN
    end if
    deallocate( lat_pre, ht_pre, r_earth, zp, coslat, f, rho0 )
  end if

  allocate( r_earth(ny,nz) )
  allocate( zp(nz), coslat(ny,nz), f(ny,nz), rho0(ny,nz) )
  allocate( lat_pre(ny), ht_pre(nz) )

  r_earth(:,:) = a_earth
  zp(:) = -h_scale*(log(p(:)/1.e5))

  coslat(:,:) = spread(cos(lat(:)*deg2rad),2,nz)
  if (abs(lat(1 )) == 90.)  coslat(1 ,:) = 0.
  if (abs(lat(ny)) == 90.)  coslat(ny,:) = 0.
  f   (:,:) = spread(2.*ome_earth*sin(lat(:)*deg2rad),2,nz)
  rho0(:,:) = spread(p(:)/g/h_scale                  ,1,ny)

  ny_pre = ny  ;  nz_pre = nz
  lat_pre(:) = lat(:)
  ht_pre (:) = p  (:)

END subroutine set_gridvar_p

SUBROUTINE set_gridvar_betap(ny,nz,lat,p,lat0,h_scale)

  integer,             intent(in) ::  ny, nz
  real,                intent(in) ::  lat0, h_scale
  real, dimension(ny), intent(in) ::  lat
  real, dimension(nz), intent(in) ::  p

  if ( allocated(ht_pre) ) then
    if ( nz == nz_pre ) then
      if ( all(p == ht_pre) )  RETURN
    end if
    deallocate( ht_pre, r_earth, zp, rho0 )
  end if

  allocate( r_earth(ny,nz) )
  allocate( zp(nz), rho0(ny,nz) )
  allocate( ht_pre(nz) )

  r_earth(:,:) = a_earth
  zp(:) = -h_scale*(log(p(:)/1.e5))

  rho0(:,:) = spread(p(:)/g/h_scale,1,ny)

  nz_pre = nz
  ht_pre(:) = p(:)

END subroutine set_gridvar_betap

SUBROUTINE set_gridvar_z(ny,nz,lat,z)

  integer,             intent(in) ::  ny, nz
  real, dimension(ny), intent(in) ::  lat
  real, dimension(nz), intent(in) ::  z

  if ( allocated(lat_pre) ) then
    if ( ny == ny_pre .and. nz == nz_pre ) then
      if ( all(lat == lat_pre) .and. all(z == ht_pre) )  RETURN
    end if
    deallocate( lat_pre, ht_pre, r_earth, coslat, f, rbeta, rcos,        &
                rrhocos )
  end if

  allocate( r_earth(ny,nz), coslat(ny,nz), f(ny,nz), rbeta(ny,nz) )
  allocate( rcos(ny,nz), rrhocos(ny,nz) )
  allocate( lat_pre(ny), ht_pre(nz) )

  r_earth(:,:) = a_earth + spread(z(:),1,ny)

  coslat(:,:) = spread(cos(lat(:)*deg2rad),2,nz)
  if (abs(lat(1 )) == 90.)  coslat(1 ,:) = 0.
  if (abs(lat(ny)) == 90.)  coslat(ny,:) = 0.
  f    (:,:) = spread(2.*ome_earth*sin(lat(:)*deg2rad),2,nz)
  rbeta(:,:) = spread(2.*ome_earth*coslat(:,1),        2,nz)

  rcos(:,:) = r_earth(:,:)*coslat(:,:)

  ny_pre = ny  ;  nz_pre = nz
  lat_pre(:) = lat(:)
  ht_pre (:) = z  (:)

END subroutine set_gridvar_z

SUBROUTINE set_gridvar_z_cp(ny,nz,lat,z,zt)

  integer,             intent(in) ::  ny, nz
  real, dimension(ny), intent(in) ::  lat
  real, dimension(nz), intent(in) ::  z, zt

  if ( allocated(lat_pre) ) then
    if ( ny == ny_pre .and. nz == nz_pre ) then
      if ( all(lat == lat_pre) .and. all(z == ht_pre) )  RETURN
    end if
    deallocate( lat_pre, ht_pre, r_earth, coslat, f, rbeta, rcos,        &
                rrhocos, r3, r3rho, r_t, r2_t, r3_t )
  end if

  allocate( r_earth(ny,nz), coslat(ny,nz), f(ny,nz), rbeta(ny,nz) )
  allocate( rcos(ny,nz), rrhocos(ny,nz), r3(ny,nz), r3rho(ny,nz),        &
            r_t(ny,nz), r2_t(ny,nz), r3_t(ny,nz) )
  allocate( lat_pre(ny), ht_pre(nz) )

  r_earth(:,:) = a_earth + spread(z (:),1,ny)
  r_t    (:,:) = a_earth + spread(zt(:),1,ny)
  r3     (:,:) = r_earth(:,:)**3
  r2_t   (:,:) = r_t (:,:)*r_t(:,:)
  r3_t   (:,:) = r2_t(:,:)*r_t(:,:)

  coslat(:,:) = spread(cos(lat(:)*deg2rad),2,nz)
  if (abs(lat(1 )) == 90.)  coslat(1 ,:) = 0.
  if (abs(lat(ny)) == 90.)  coslat(ny,:) = 0.
  f    (:,:) = spread(2.*ome_earth*sin(lat(:)*deg2rad),2,nz)
  rbeta(:,:) = spread(2.*ome_earth*coslat(:,1),        2,nz)

  rcos(:,:) = r_earth(:,:)*coslat(:,:)

  ny_pre = ny  ;  nz_pre = nz
  lat_pre(:) = lat(:)
  ht_pre (:) = z  (:)

END subroutine set_gridvar_z_cp

SUBROUTINE grady_2nd(ny,nz,var,lat, grady)

  integer,                intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(ny)   , intent(in)  ::  lat
  real, dimension(ny,nz), intent(out) ::  grady

  integer ::  j,k
  real    ::  inv_2dy(ny,nz)

  do j=2, ny-1
    inv_2dy(j,1) = 1./((lat(j+1)-lat(j-1))*deg2rad)
  enddo
  inv_2dy(:,2:nz) = spread(inv_2dy(:,1),2,nz-1)
  inv_2dy(:,:) = inv_2dy(:,:)/r_earth(:,:)

  do k=1, nz
  do j=2, ny-1
    grady(j,k) = (var(j+1,k)-var(j-1,k))*inv_2dy(j,k)
  enddo
  enddo
  grady(1 ,:) = 0.
  grady(ny,:) = 0.

END subroutine grady_2nd

SUBROUTINE gradz_2nd_irr(ny,nz,var,z, gradz)

  integer,                intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(nz),    intent(in)  ::  z
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
  do j=1, ny
    gradz(j,k) = sum( coef(:,k)*var(j,k-1:k+1) )
  enddo
  enddo

  gradz(:,1 ) = (var(:,2 ) - var(:,1   ))*inv_dz_1
  gradz(:,nz) = (var(:,nz) - var(:,nz-1))*inv_dz_n

END subroutine gradz_2nd_irr

SUBROUTINE gradz_i2m_2d(ny,nz,var,zi, gradz)

  integer,                intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(nz),    intent(in)  ::  zi
  real, dimension(ny,nz), intent(out) ::  gradz

  integer ::  k
  real    ::  inv_dz(nz)

  do k=2, nz
    inv_dz(k) = 1./(zi(k) - zi(k-1))
  enddo
  do k=2, nz
    gradz(:,k) = (var(:,k) - var(:,k-1))*inv_dz(k)
  enddo
  gradz(:,1) = gradz(:,2)

END subroutine gradz_i2m_2d

SUBROUTINE gradz_i2m_3d(nx,ny,nz,var,zi, gradz)

  integer,                   intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(nz),       intent(in)  ::  zi
  real, dimension(nx,ny,nz), intent(out) ::  gradz

  integer ::  k
  real    ::  inv_dz(nz)

  do k=2, nz
    inv_dz(k) = 1./(zi(k) - zi(k-1))
  enddo
  do k=2, nz
    gradz(:,:,k) = (var(:,:,k) - var(:,:,k-1))*inv_dz(k)
  enddo
  gradz(:,:,1) = gradz(:,:,2)

END subroutine gradz_i2m_3d

SUBROUTINE gradz_m2i(ny,nz,var,z, gradz)

  integer,                intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(nz),    intent(in)  ::  z
  real, dimension(ny,nz), intent(out) ::  gradz

  integer ::  k
  real    ::  inv_dz(nz)

  do k=1, nz-1
    inv_dz(k) = 1./(z(k+1) - z(k))
  enddo
  do k=1, nz-1
    gradz(:,k) = (var(:,k+1) - var(:,k))*inv_dz(k)
  enddo
  gradz(:,nz) = gradz(:,nz-1)

END subroutine gradz_m2i

SUBROUTINE gradzwav_m2i(nk,no,ny,nz,fcr,fci,z, gradzr,gradzi)

  integer,                      intent(in)  ::  nk, no, ny, nz
  real, dimension(nk,no,ny,nz), intent(in)  ::  fcr, fci
  real, dimension(nz),          intent(in)  ::  z
  real, dimension(nk,no,ny,nz), intent(out) ::  gradzr, gradzi

  integer ::  k
  real    ::  inv_dz(nz), r_small
  complex ::  fc(nk,no,ny,nz)

  real, parameter ::  pi = 3.14159265358979323846
  real, parameter ::  twopi = 2.*3.14159265358979323846

  r_small = tiny(fcr(1,1,1,1))

  do k=1, nz-1
    inv_dz(k) = 1./(z(k+1) - z(k))
  enddo

  gradzr(:,:,:,:) = fcr(:,:,:,:)  ;  gradzi(:,:,:,:) = fci(:,:,:,:)
  call intpwav_m2i(nk,no,ny,nz, gradzr,gradzi)

  fc(:,:,:,:) = log(cmplx(fcr+r_small,fci))
  do k=1, nz-1
    fc(:,:,:,k) = fc(:,:,:,k+1) - fc(:,:,:,k)
  enddo
  where ( abs(aimag(fc(:,:,:,:nz-1))) > pi )
    fc(:,:,:,:nz-1) = fc(:,:,:,:nz-1) -                                  &
                      cmplx(0.,sign(twopi,aimag(fc(:,:,:,:nz-1))))
  end where
  fc(:,:,:,:nz-1) = fc(:,:,:,:nz-1)*                                     &
                    cmplx(gradzr(:,:,:,:nz-1),gradzi(:,:,:,:nz-1))*      &
                    spread(spread(spread(inv_dz(:nz-1),1,nk),2,no),3,ny)
  fc(:,:,:,nz) = fc(:,:,:,nz-1)

  gradzr = real (fc)
  gradzi = aimag(fc)

END subroutine gradzwav_m2i

SUBROUTINE gradzprt_m2i(nx,ny,nz,prt,z, gradz)

  use fft

  integer,                   intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  prt
  real, dimension(nz),       intent(in)  ::  z
  real, dimension(nx,ny,nz), intent(out) ::  gradz

  integer ::  j,k
  real    ::  gradz_fcr(nx,ny,1,nz), gradz_fci(nx,ny,1,nz)
  complex ::  fc(nx,ny,1,nz)

  do k=1, nz
  do j=1, ny
    call fft1d_f(nx,prt(:,j,k), fc(:,j,1,k))
  enddo
  enddo

  call gradzwav_m2i(nx,ny,1,nz,real(fc),aimag(fc),z, gradz_fcr,gradz_fci)

  fc(:,:,:,:) = cmplx(gradz_fcr,gradz_fci)
  do k=1, nz
  do j=1, ny
    call fft1d_b(nx,fc(:,j,1,k), gradz(:,j,k))
  enddo
  enddo

END subroutine gradzprt_m2i

SUBROUTINE intp_i2m(ny,nz, var)

  integer,                intent(in)    ::  ny, nz
  real, dimension(ny,nz), intent(inout) ::  var

  integer ::  k

  do k=nz, 2, -1
    var(:,k) = 0.5*(var(:,k-1) + var(:,k))
  enddo

END subroutine intp_i2m

SUBROUTINE intp_m2i(ny,nz,zm,zi, var)

  integer,                intent(in)    ::  ny, nz
  real, dimension(nz),    intent(in)    ::  zm, zi
  real, dimension(ny,nz), intent(inout) ::  var

  integer ::  k
  real    ::  coef(2,nz), dz(nz)

  dz(1:nz-1) = zm(2:nz) - zm(1:nz-1)
  do k=1, nz-1
    coef(1,k) = (zm(k+1) - zi(k))/dz(k)
    coef(2,k) = (zi(k)   - zm(k))/dz(k)
  enddo

  do k=1, nz-1
    var(:,k) = var(:,k)*coef(1,k) + var(:,k+1)*coef(2,k)
  enddo

END subroutine intp_m2i

SUBROUTINE intpwav_i2m(nk,no,ny,nz, fcr,fci)

  integer,                      intent(in)    ::  nk, no, ny, nz
  real, dimension(nk,no,ny,nz), intent(inout) ::  fcr, fci

  integer ::  k
  real    ::  r_small
  complex ::  fc(nk,no,ny,nz)

  complex, parameter ::  ipi = (0.,3.14159265358979323846)
  real,    parameter ::  pi05 = 0.5*3.14159265358979323846

  r_small = tiny(fcr(1,1,1,1))

  fc(:,:,:,:) = log(cmplx(fcr+r_small,fci))
  do k=1, nz-1
    fc(:,:,:,k) = 0.5*(fc(:,:,:,k) + fc(:,:,:,k+1))
    where ( abs(aimag(fc(:,:,:,k)-fc(:,:,:,k+1))) > pi05 )
      fc(:,:,:,k) = fc(:,:,:,k) + ipi
    end where
  enddo
  fc(:,:,:,:) = exp(fc(:,:,:,:))

  fcr(:,:,:,2:nz) = real (fc(:,:,:,1:nz-1))
  fci(:,:,:,2:nz) = aimag(fc(:,:,:,1:nz-1))

END subroutine intpwav_i2m

SUBROUTINE intpwav_m2i(nk,no,ny,nz, fcr,fci)

  integer,                      intent(in)    ::  nk, no, ny, nz
  real, dimension(nk,no,ny,nz), intent(inout) ::  fcr, fci

  real ::  fc2(nk,no,ny,2)

  fc2(:,:,:,1) = fcr(:,:,:,nz)
  fc2(:,:,:,2) = fci(:,:,:,nz)

  call intpwav_i2m(nk,no,ny,nz, fcr,fci)
  fcr(:,:,:,1:nz-1) = fcr(:,:,:,2:nz)
  fci(:,:,:,1:nz-1) = fci(:,:,:,2:nz)

  fcr(:,:,:,nz) = fc2(:,:,:,1)
  fci(:,:,:,nz) = fc2(:,:,:,2)

END subroutine intpwav_m2i

SUBROUTINE intpprt_i2m(nx,ny,nz, prt)

  use fft

  integer,                   intent(in)    ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(inout) ::  prt

  integer ::  j,k
  real    ::  r_small
  complex ::  fc(nx,ny,nz)

  complex, parameter ::  ipi = (0.,3.14159265358979323846)
  real   , parameter ::  pi05 = 0.5*3.14159265358979323846

  r_small = tiny(prt(1,1,1))

  do k=1, nz
  do j=1, ny
    call fft1d_f(nx,prt(:,j,k), fc(:,j,k))
  enddo
  enddo

  fc(:,:,:) = log(fc(:,:,:)+r_small)
  do k=1, nz-1
    fc(:,:,k) = 0.5*(fc(:,:,k) + fc(:,:,k+1))
    where ( abs(aimag(fc(:,:,k)-fc(:,:,k+1))) > pi05 )
      fc(:,:,k) = fc(:,:,k) + ipi
    end where
  enddo
  fc(:,:,:) = exp(fc(:,:,:))

  do k=2, nz
  do j=1, ny
    call fft1d_b(nx,fc(:,j,k-1), prt(:,j,k))
  enddo
  enddo

END subroutine intpprt_i2m

SUBROUTINE intpprt_m2i(nx,ny,nz, prt)

  integer,                   intent(in)    ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(inout) ::  prt

  real ::  prt1(nx,ny)

  prt1(:,:) = prt(:,:,nz)

  call intpprt_i2m(nx,ny,nz, prt)
  prt(:,:,1:nz-1) = prt(:,:,2:nz)

  prt(:,:,nz) = prt1(:,:)

END subroutine intpprt_m2i

SUBROUTINE missing_bdy(ny,nz,var,missv,nm_y1,nm_y2,nm_z1,nm_z2)

  integer, intent(in) ::  ny, nz, nm_y1, nm_y2, nm_z1, nm_z2
  real,    intent(in) ::  missv

  real, dimension(ny,nz), intent(inout) ::  var

  if (nm_y1 > 0)  var(1         :nm_y1,          :     ) = missv
  if (nm_y2 > 0)  var(ny+1-nm_y2:ny   ,          :     ) = missv
  if (nm_z1 > 0)  var(          :     ,1         :nm_z1) = missv
  if (nm_z2 > 0)  var(          :     ,nz+1-nm_z2:nz   ) = missv

END subroutine missing_bdy

END module warc3d


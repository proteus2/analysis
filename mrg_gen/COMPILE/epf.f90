MODULE epf

  implicit none

  private ::  set_gridvar_p, set_gridvar_z, set_gridvar_z_cp
  private ::  grady_2nd, gradz_2nd_irr, gradz_i2m, gradz_m2i,            &
              gradzwav_m2i, gradzprt_m2i
  private ::  intp_i2m, intp_m2i, intpwav_i2m, intpwav_m2i, intpprt_i2m, &
              intpprt_m2i
  private ::  missing_bdy

  character(len=64), dimension(8) ::  varname_epf =                      &
      (/'f_y    ','f_z    ','epd    ','epd_z  ','f_uv   ',               &
        'f_uw   ','epd_uv ','epd_uw '/)

  real, dimension(:,:),   allocatable ::  dptmdz_from_thlev
  real, dimension(:,:,:), allocatable ::  dptmdz_from_thlev3

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
    module procedure gradz_i2m_2d, gradz_i2m_3d, gradz_i2m_4d
  end interface
  interface intp_i2m
    module procedure intp_i2m_2d, intp_i2m_4d
  end interface


  CONTAINS


SUBROUTINE epf_hydro_p(                                                  &
     nx,ny,nz,lat,p,u,v,w,pt,h_scale,missv,                              &
     fy,fz,epd,epd_z,f_uv,f_uw,epd_uv,epd_uw )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  u, v, w, pt

  real, dimension(ny,nz), intent(out) ::  fy, fz, f_uv, f_uw
  real, dimension(ny,nz), intent(out) ::  epd, epd_z, epd_uv, epd_uw

  real, dimension(ny,nz) ::  um, ptm, vm, wm
  real, dimension(ny,nz) ::  phi, rvpt, rwpt, rvu, rwu
  real, dimension(ny,nz) ::  dptmdy, dptmdz, dptmd2z, divy_um, dumdz
  real, dimension(ny,nz) ::  grady, gradz, temp, epd_y

  real, dimension(:,:,:), allocatable ::  prt, v_prt, w_prt

  integer ::  j,k

  call set_gridvar_p(ny,nz,lat,p,h_scale)

  ! mean and perturbation
  allocate( prt(nx,ny,nz), v_prt(nx,ny,nz), w_prt(nx,ny,nz) )

  um  (:,:) = sum(u , dim=1)/float(nx)
  vm  (:,:) = sum(v , dim=1)/float(nx)
  wm  (:,:) = sum(w , dim=1)/float(nx)
  ptm (:,:) = sum(pt, dim=1)/float(nx)

  v_prt(:,:,:) = v(:,:,:) - spread(vm(:,:),1,nx)
  w_prt(:,:,:) = w(:,:,:) - spread(wm(:,:),1,nx)

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
  call grady_2nd(ny,nz,ptm,lat, dptmdy)
  call missing_bdy(ny,nz,dptmdy,missv,1,1,0,0)

  temp(:,:) = log(ptm(:,:))
  call gradz_2nd_irr(ny,nz,temp,zp, dptmdz)
  call gradz_2nd_irr(ny,nz,dptmdz,zp, dptmd2z)
  dptmd2z(:,:) = dptmd2z(:,:) + dptmdz(:,:)*dptmdz(:,:)
  dptmdz (:,:) = ptm(:,:)*dptmdz (:,:)
  dptmd2z(:,:) = ptm(:,:)*dptmd2z(:,:)

  call missing_bdy(ny,nz,dptmdz ,missv,0,0,1,0)
  call missing_bdy(ny,nz,dptmd2z,missv,0,0,1,0)

  ! grad, u
  call grady_2nd(ny,nz,um*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um,zp, dumdz)
  call missing_bdy(ny,nz,divy_um,missv,1,1,0,0)
  call missing_bdy(ny,nz,dumdz  ,missv,0,0,1,0)

  ! phi
  phi(:,:) = rvpt(:,:)/dptmdz(:,:)
  call missing_bdy(ny,nz,phi,missv,0,0,1,0)

  ! epf, epd
  fy  (:,:) = -a_earth*coslat(:,:)*(rvu(:,:) - phi(:,:)*dumdz(:,:))
  fz  (:,:) = -a_earth*coslat(:,:)*(rwu(:,:) + phi(:,:)*(divy_um(:,:)-f(:,:)))
  f_uv(:,:) = -a_earth*coslat(:,:)*rvu(:,:)
  f_uw(:,:) = -a_earth*coslat(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy  ,missv,0,0,1,0)
  call missing_bdy(ny,nz,fz  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uv,missv,0,0,1,0)
  call missing_bdy(ny,nz,f_uw,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = grady(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = grady(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) &
                  *86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,1,1,1,0)

  temp(:,:) = (divy_um(:,:)-f(:,:))*rvpt(:,:)
  temp(1,:) = 0.  ;  temp(ny,:) = 0.
  call gradz_2nd_irr(ny,nz,temp,zp, gradz)
  temp(:,:) = (gradz(:,:) - temp(:,:)*dptmd2z(:,:)/dptmdz(:,:))/dptmdz(:,:)
  call gradz_2nd_irr(ny,nz,rwu,zp, gradz)
  epd_z(:,:) = -(gradz(:,:)+temp(:,:))/rho0(:,:) * 86400.
  epd_uw(:,:) = -gradz(:,:)/rho0(:,:) * 86400.
  call missing_bdy(ny,nz,epd_z ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,1,0)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,1,1,1,0)

END subroutine epf_hydro_p

SUBROUTINE epf_z(                                                        &
     nx,ny,nz,lat,z,u,v,w,pt,rho,missv,                                  &
     fy,fz,epd,epd_z,f_uv,f_uw,epd_uv,epd_uw )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  z
  real, dimension(nx,ny,nz), intent(in) ::  u, v, w, pt, rho
  
  real, dimension(ny,nz), intent(out) ::  fy, fz, f_uv, f_uw
  real, dimension(ny,nz), intent(out) ::  epd, epd_z, epd_uv, epd_uw
  
  real, dimension(ny,nz) ::  um, ptm, rhom, rvm, rwm
  real, dimension(ny,nz) ::  phi, rvpt, rwpt, rvu, rwu
  real, dimension(ny,nz) ::  dptmdy, dptmdz, divy_um, dumdz
  real, dimension(ny,nz) ::  grady, gradz, temp, epd_y
  
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

  deallocate( prt, rv_prt, rw_prt )

  ! grad, pt
  call grady_2nd(ny,nz,ptm,lat, dptmdy)
  call missing_bdy(ny,nz,dptmdy,missv,1,1,0,0)

  if ( allocated(dptmdz_from_thlev) ) then
    dptmdz(:,:) = dptmdz_from_thlev(:,:)
  else
    temp(:,:) = log(ptm(:,:))
    call gradz_2nd_irr(ny,nz,temp,z, dptmdz)
    dptmdz(:,:) = ptm(:,:)*dptmdz(:,:)
  end if
  call missing_bdy(ny,nz,dptmdz,missv,0,0,1,0)

  ! grad, u
  call grady_2nd(ny,nz,um*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um*r_earth,z, dumdz)
  dumdz(:,:) = dumdz(:,:)/r_earth(:,:)
  call missing_bdy(ny,nz,divy_um,missv,1,1,0,0)
  call missing_bdy(ny,nz,dumdz  ,missv,0,0,1,0)

  ! phi
  phi(:,:) = (rvpt(:,:)*dptmdz(:,:) - rwpt(:,:)*dptmdy(:,:))/            &
             (dptmdy(:,:)*dptmdy(:,:) + dptmdz(:,:)*dptmdz(:,:))
  call missing_bdy(ny,nz,phi,missv,1,1,1,0)

  ! epf, epd
  fy  (:,:) = -rcos(:,:)*(rvu(:,:) - phi(:,:)*(dumdz  (:,:)+rbeta(:,:)))
  fz  (:,:) = -rcos(:,:)*(rwu(:,:) + phi(:,:)*(divy_um(:,:)-f    (:,:)))
  f_uv(:,:) = -rcos(:,:)*rvu(:,:)
  f_uw(:,:) = -rcos(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uv,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uw,missv,1,1,1,0)

  rrhocos(:,:) = rcos(:,:)*rhom(:,:)

  call grady_2nd(ny,nz,fy*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,2,2,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,2,2,1,0)

  temp(:,:) = r_earth(:,:)*r_earth(:,:)
  call gradz_2nd_irr(ny,nz,fz*temp,z, gradz)
  do k=1, nz
  do j=2, ny-1
    epd_z(j,k) = gradz(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,f_uw*temp,z, gradz)
  do k=1, nz
  do j=2, ny-1
    epd_uw(j,k) = gradz(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_z ,missv,1,1,2,1)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,2,1)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,2,1)

END subroutine epf_z

SUBROUTINE epf_hydro_p_fc2(                                              &
     nk,nome,ny,nz,lat,p,nt,um,ptm,h_scale,missv,                        &
     u_fc_r,pt_fc_r,v_fc_r,w_fc_r,u_fc_i,pt_fc_i,v_fc_i,w_fc_i,          &
     fy0,fz0,epd0,epd_z0,f_uv0,f_uw0,epd_uv0,epd_uw0 )

  integer,                   intent(in) ::  nk, nome, ny, nz, nt
  real,                      intent(in) ::  h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(ny,nz,nt), intent(in) ::  um, ptm

  real, dimension(nk,nome,ny,nz), intent(in) ::  u_fc_r, pt_fc_r,        &
                         v_fc_r, w_fc_r, u_fc_i, pt_fc_i, v_fc_i, w_fc_i

  real, dimension(nk,nome,ny,nz), intent(out) ::  fy0, fz0, f_uv0, f_uw0
  real, dimension(nk,nome,ny,nz), intent(out) ::  epd0, epd_z0, epd_uv0, &
                                                  epd_uw0

  real, dimension(ny,nz) ::  fy, fz, f_uv, f_uw
  real, dimension(ny,nz) ::  epd, epd_z, epd_uv, epd_uw 
  real, dimension(ny,nz) ::  rvpt, rwpt, rvu, rwu, fphiy, fphiz
  real, dimension(ny,nz) ::  phiy, phiz
  real, dimension(ny,nz) ::  dptmdy, dptmdz, dptmd2z, divy_um, dumdz
  real, dimension(ny,nz) ::  grady, gradz, temp, epd_y

  integer ::  j,k, i,n

  call set_gridvar_p(ny,nz,lat,p,h_scale)

  phiy(:,:) = 0.  ;  phiz(:,:) = 0.

  L_BTIME:  DO n=1, nt

  ! grad, pt
  call grady_2nd(ny,nz,ptm(:,:,n),lat, dptmdy)
  call missing_bdy(ny,nz,dptmdy,missv,1,1,0,0)

  temp(:,:) = log(ptm(:,:,n))
  call gradz_2nd_irr(ny,nz,temp,zp, dptmdz)
  call gradz_2nd_irr(ny,nz,dptmdz,zp, dptmd2z)
  dptmd2z(:,:) = dptmd2z(:,:) + dptmdz(:,:)*dptmdz(:,:)
  dptmdz (:,:) = ptm(:,:,n)*dptmdz (:,:)
  dptmd2z(:,:) = ptm(:,:,n)*dptmd2z(:,:)

  call missing_bdy(ny,nz,dptmdz ,missv,0,0,1,0)
  call missing_bdy(ny,nz,dptmd2z,missv,0,0,1,0)

  ! grad, u
  call grady_2nd(ny,nz,um(:,:,n)*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um(:,:,n),zp, dumdz)
  call missing_bdy(ny,nz,divy_um,missv,1,1,0,0)
  call missing_bdy(ny,nz,dumdz  ,missv,0,0,1,0)

  divy_um(:,:) = divy_um(:,:) - f(:,:)

  phiy(:,:) = phiy(:,:) + dumdz  (:,:)/dptmdz(:,:)/float(nt)
  phiz(:,:) = phiz(:,:) + divy_um(:,:)/dptmdz(:,:)/float(nt)

  ENDDO  L_BTIME

  call missing_bdy(ny,nz,phiy,missv,0,0,1,0)
  call missing_bdy(ny,nz,phiz,missv,1,1,1,0)

  fy0     = 0.  ;  fz0     = 0.  ;  epd0    = 0.  ;  epd_z0  = 0.
  f_uv0   = 0.  ;  f_uw0   = 0.  ;  epd_uv0 = 0.  ;  epd_uw0 = 0.

  L_OME:  DO n=1, nome
  L_KWN:  DO i=2, nk

  rvpt(:,:) = v_fc_r(i,n,:,:)*pt_fc_r(i,n,:,:) + &
              v_fc_i(i,n,:,:)*pt_fc_i(i,n,:,:)
  rwpt(:,:) = w_fc_r(i,n,:,:)*pt_fc_r(i,n,:,:) + &
              w_fc_i(i,n,:,:)*pt_fc_i(i,n,:,:)
  rvu (:,:) = v_fc_r(i,n,:,:)*u_fc_r (i,n,:,:) + &
              v_fc_i(i,n,:,:)*u_fc_i (i,n,:,:)
  rwu (:,:) = w_fc_r(i,n,:,:)*u_fc_r (i,n,:,:) + &
              w_fc_i(i,n,:,:)*u_fc_i (i,n,:,:)

  rvpt(:,:) = rho0(:,:)*rvpt(:,:)
  rwpt(:,:) = rho0(:,:)*rwpt(:,:)
  rvu (:,:) = rho0(:,:)*rvu (:,:)
  rwu (:,:) = rho0(:,:)*rwu (:,:)

  if ( all(rvpt == rwpt) .and. all(rvu == rwu) )  CYCLE
  ! all-zero

  ! phi
  fphiy(:,:) = rvpt(:,:)*phiy(:,:)
  fphiz(:,:) = rvpt(:,:)*phiz(:,:)
  call missing_bdy(ny,nz,fphiy,missv,0,0,1,0)
  call missing_bdy(ny,nz,fphiz,missv,1,1,1,0)

  ! epf, epd
  fy  (:,:) = -a_earth*coslat(:,:)*(rvu(:,:) - fphiy(:,:))
  fz  (:,:) = -a_earth*coslat(:,:)*(rwu(:,:) + fphiz(:,:))
  f_uv(:,:) = -a_earth*coslat(:,:)*rvu(:,:)
  f_uw(:,:) = -a_earth*coslat(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy  ,missv,0,0,1,0)
  call missing_bdy(ny,nz,fz  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uv,missv,0,0,1,0)
  call missing_bdy(ny,nz,f_uw,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = grady(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = grady(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) &
                  *86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,1,1,1,0)

  temp(:,:) = (divy_um(:,:)-f(:,:))*rvpt(:,:)
  temp(1,:) = 0.  ;  temp(ny,:) = 0.
  call gradz_2nd_irr(ny,nz,temp,zp, gradz)
  temp(:,:) = (gradz(:,:) - temp(:,:)*dptmd2z(:,:)/dptmdz(:,:))/dptmdz(:,:)
  call gradz_2nd_irr(ny,nz,rwu,zp, gradz)
  epd_z(:,:) = -(gradz(:,:)+temp(:,:))/rho0(:,:) * 86400.
  epd_uw(:,:) = -gradz(:,:)/rho0(:,:) * 86400.
  call missing_bdy(ny,nz,epd_z ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,1,0)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,1,1,1,0)

  fy0    (i,n,:,:) = fy    (:,:)
  fz0    (i,n,:,:) = fz    (:,:)
  epd0   (i,n,:,:) = epd   (:,:)
  epd_z0 (i,n,:,:) = epd_z (:,:)
  f_uv0  (i,n,:,:) = f_uv  (:,:)
  f_uw0  (i,n,:,:) = f_uw  (:,:)
  epd_uv0(i,n,:,:) = epd_uv(:,:)
  epd_uw0(i,n,:,:) = epd_uw(:,:)

  ENDDO  L_KWN
  ENDDO  L_OME

END subroutine epf_hydro_p_fc2

SUBROUTINE epf_z_fc2(                                                    &
     nk,nome,ny,nz,lat,z,nt,um,ptm,rhom,missv,                           &
     u_fc_r,pt_fc_r,rv_fc_r,rw_fc_r,u_fc_i,pt_fc_i,rv_fc_i,rw_fc_i,      &
     fy0,fz0,epd0,epd_z0,f_uv0,f_uw0,epd_uv0,epd_uw0 )

  integer,                        intent(in) ::  nk, nome, ny, nz, nt
  real,                           intent(in) ::  missv
  real, dimension(ny),            intent(in) ::  lat
  real, dimension(nz),            intent(in) ::  z
  real, dimension(ny,nz,nt),      intent(in) ::  um, ptm, rhom
  real, dimension(nk,nome,ny,nz), intent(in) ::  u_fc_r, pt_fc_r,        &
                     rv_fc_r, rw_fc_r, u_fc_i, pt_fc_i, rv_fc_i, rw_fc_i
  
  real, dimension(nk,nome,ny,nz), intent(out) ::  fy0, fz0, f_uv0, f_uw0
  real, dimension(nk,nome,ny,nz), intent(out) ::  epd0, epd_z0, epd_uv0, &
                                                  epd_uw0
 
  real, dimension(ny,nz) ::  fy, fz, f_uv, f_uw
  real, dimension(ny,nz) ::  epd, epd_z, epd_uv, epd_uw 
  real, dimension(ny,nz) ::  rvpt, rwpt, rvu, rwu, fphiy, fphiz
  real, dimension(ny,nz) ::  phi1, phi2, phi1y, phi1z, phi2y, phi2z
  real, dimension(ny,nz) ::  dptmdy, dptmdz, divy_um, dumdz, rhomm
  real, dimension(ny,nz) ::  grady, gradz, temp, epd_y
  
  integer ::  j,k, i,n
 
  call set_gridvar_z(ny,nz,lat,z)

  rhomm(:,:) = sum(rhom(:,:,1:nt), dim=3)/float(nt)

  phi1y(:,:) = 0.  ;  phi1z(:,:) = 0.
  phi2y(:,:) = 0.  ;  phi2z(:,:) = 0.

  L_BTIME:  DO n=1, nt

  ! grad, pt
  call grady_2nd(ny,nz,ptm(:,:,n),lat, dptmdy)
  call missing_bdy(ny,nz,dptmdy,missv,1,1,0,0)

  if ( allocated(dptmdz_from_thlev3) ) then
    dptmdz(:,:) = dptmdz_from_thlev3(:,:,n)
  else
    temp(:,:) = log(ptm(:,:,n))
    call gradz_2nd_irr(ny,nz,temp,z, dptmdz)
    dptmdz(:,:) = ptm(:,:,n)*dptmdz(:,:)
  end if
  call missing_bdy(ny,nz,dptmdz,missv,0,0,1,0)

  ! grad, u
  call grady_2nd(ny,nz,um(:,:,n)*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um(:,:,n)*r_earth,z, dumdz)
  dumdz(:,:) = dumdz(:,:)/r_earth(:,:)
  call missing_bdy(ny,nz,divy_um,missv,1,1,0,0)
  call missing_bdy(ny,nz,dumdz  ,missv,0,0,1,0)

  temp(:,:) = dptmdy(:,:)*dptmdy(:,:) + dptmdz(:,:)*dptmdz(:,:)
  phi1(:,:) = dptmdz(:,:)/temp(:,:)
  phi2(:,:) = dptmdy(:,:)/temp(:,:)

  dumdz  (:,:) = dumdz  (:,:)+rbeta(:,:)
  divy_um(:,:) = divy_um(:,:)-f    (:,:)

  phi1y(:,:) = phi1y(:,:) + phi1(:,:)*dumdz  (:,:)/float(nt)
  phi2y(:,:) = phi2y(:,:) + phi2(:,:)*dumdz  (:,:)/float(nt)
  phi1z(:,:) = phi1z(:,:) + phi1(:,:)*divy_um(:,:)/float(nt)
  phi2z(:,:) = phi2z(:,:) + phi2(:,:)*divy_um(:,:)/float(nt)

  ENDDO  L_BTIME

  call missing_bdy(ny,nz,phi1y,missv,1,1,1,0)
  call missing_bdy(ny,nz,phi1z,missv,1,1,1,0)
  call missing_bdy(ny,nz,phi2y,missv,1,1,1,0)
  call missing_bdy(ny,nz,phi2z,missv,1,1,1,0)

  rrhocos(:,:) = rcos(:,:)*rhomm(:,:)

  fy0     = 0.  ;  fz0     = 0.  ;  epd0    = 0.  ;  epd_z0  = 0.
  f_uv0   = 0.  ;  f_uw0   = 0.  ;  epd_uv0 = 0.  ;  epd_uw0 = 0.

  L_OME:  DO n=1, nome
  L_KWN:  DO i=2, nk

  rvpt(:,:) = rv_fc_r(i,n,:,:)*pt_fc_r(i,n,:,:) + &
              rv_fc_i(i,n,:,:)*pt_fc_i(i,n,:,:)
  rwpt(:,:) = rw_fc_r(i,n,:,:)*pt_fc_r(i,n,:,:) + &
              rw_fc_i(i,n,:,:)*pt_fc_i(i,n,:,:)
  rvu (:,:) = rv_fc_r(i,n,:,:)*u_fc_r (i,n,:,:) + &
              rv_fc_i(i,n,:,:)*u_fc_i (i,n,:,:)
  rwu (:,:) = rw_fc_r(i,n,:,:)*u_fc_r (i,n,:,:) + &
              rw_fc_i(i,n,:,:)*u_fc_i (i,n,:,:)

  if ( all(rvpt == rwpt) .and. all(rvu == rwu) )  CYCLE
  ! all-zero

  ! phi
  fphiy(:,:) = rvpt(:,:)*phi1y(:,:) - rwpt(:,:)*phi2y(:,:)
  fphiz(:,:) = rvpt(:,:)*phi1z(:,:) - rwpt(:,:)*phi2z(:,:)
  call missing_bdy(ny,nz,fphiy,missv,1,1,1,0)
  call missing_bdy(ny,nz,fphiz,missv,1,1,1,0)

  ! epf, epd
  fy  (:,:) = -rcos(:,:)*(rvu(:,:) - fphiy(:,:))
  fz  (:,:) = -rcos(:,:)*(rwu(:,:) + fphiz(:,:))
  f_uv(:,:) = -rcos(:,:)*rvu(:,:)
  f_uw(:,:) = -rcos(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uv,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uw,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,2,2,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,2,2,1,0)

  temp(:,:) = r_earth(:,:)*r_earth(:,:)
  call gradz_2nd_irr(ny,nz,fz*temp,z, gradz)
  do k=1, nz
  do j=2, ny-1
    epd_z(j,k) = gradz(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,f_uw*temp,z, gradz)
  do k=1, nz
  do j=2, ny-1
    epd_uw(j,k) = gradz(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_z ,missv,1,1,2,1)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,2,1)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,2,1)

  fy0    (i,n,:,:) = fy    (:,:)
  fz0    (i,n,:,:) = fz    (:,:)
  epd0   (i,n,:,:) = epd   (:,:)
  epd_z0 (i,n,:,:) = epd_z (:,:)
  f_uv0  (i,n,:,:) = f_uv  (:,:)
  f_uw0  (i,n,:,:) = f_uw  (:,:)
  epd_uv0(i,n,:,:) = epd_uv(:,:)
  epd_uw0(i,n,:,:) = epd_uw(:,:)

  ENDDO  L_KWN
  ENDDO  L_OME

END subroutine epf_z_fc2

SUBROUTINE epf_z_cp(                                                     &
     nx,ny,nz,lat,z,zt,u,v,w,pt,rho,missv,                               &
     fy,fz,epd,epd_z,f_uv,f_uw,epd_uv,epd_uw )

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  z, zt
  real, dimension(nx,ny,nz), intent(in) ::  u, v, w, pt, rho
  
  real, dimension(ny,nz), intent(out) ::  fy, fz, f_uv, f_uw
  real, dimension(ny,nz), intent(out) ::  epd, epd_z, epd_uv, epd_uw
  
  real, dimension(ny,nz) ::  um, ptm, rhom, rvm, rwm
  real, dimension(ny,nz) ::  rvpt, rwpt, rvu, rwu, druwr3dz
  real, dimension(ny,nz) ::  phi, fphiy, fphiz, fphiz_t,                 &
                             dptmdy, dptmdz, divy_um_t, dumdz
  real, dimension(ny,nz) ::  grady, gradz, temp, epd_y
  
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
  call missing_bdy(ny,nz,dptmdy,missv,1,1,0,0)

  temp(:,:) = log(ptm(:,:))
  call gradz_2nd_irr(ny,nz,temp,zt, dptmdz)
  dptmdz(:,:) = ptm(:,:)*dptmdz(:,:)
  call missing_bdy(ny,nz,dptmdz,missv,0,0,1,0)

  ! grad, u
  call grady_2nd(ny,nz,um*coslat,lat, divy_um_t)
  do k=1, nz
  do j=2, ny-1
    divy_um_t(j,k) = divy_um_t(j,k)/coslat(j,k)
  enddo
  enddo
  divy_um_t(:,:) = divy_um_t(:,:)*r_earth(:,:)
  call intp_m2i(ny,nz,z,zt, divy_um_t)
  divy_um_t(:,:) = divy_um_t(:,:)/r_t(:,:)
  call gradz_m2i(ny,nz,um*r_earth,z, dumdz)
  dumdz(:,:) = dumdz(:,:)/r_t(:,:)

  call missing_bdy(ny,nz,divy_um_t,missv,1,1,0,0)
  call missing_bdy(ny,nz,dumdz    ,missv,0,0,0,0)

  ! phi
  phi(:,:) = (rvpt(:,:)*dptmdz(:,:) - rwpt(:,:)*dptmdy(:,:))/            &
             (dptmdy(:,:)*dptmdy(:,:) + dptmdz(:,:)*dptmdz(:,:))
  fphiy(:,:) = phi(:,:)*(dumdz    (:,:)+rbeta(:,:))
  fphiz(:,:) = phi(:,:)*(divy_um_t(:,:)-f    (:,:))
  fphiz_t(:,:) = fphiz(:,:)
  call intp_i2m(ny,nz, fphiy)
  call intp_i2m(ny,nz, fphiz)
  call missing_bdy(ny,nz,fphiy,missv,1,1,1,0)
  call missing_bdy(ny,nz,fphiz,missv,1,1,1,0)

  ! epf, epd
  fy  (:,:) = -rcos(:,:)*(rvu(:,:) - fphiy(:,:))
  fz  (:,:) = -rcos(:,:)*(rwu(:,:) + fphiz(:,:))
  f_uv(:,:) = -rcos(:,:)*rvu(:,:)
  f_uw(:,:) = -rcos(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uv,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uw,missv,1,1,1,0)

  rrhocos(:,:) = rcos(:,:)*rhom(:,:)
  r3rho  (:,:) = r3  (:,:)*rhom(:,:)

  call grady_2nd(ny,nz,fy*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,2,2,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,2,2,1,0)

  epd_uw(:,:) = -druwr3dz(:,:)/r3rho(:,:) * 86400.
  call gradz_i2m(ny,nz,r3_t*fphiz_t,zt, gradz)
  epd_z(:,:) = epd_uw(:,:) - gradz(:,:)/r3rho(:,:)*86400.
  call missing_bdy(ny,nz,epd_z ,missv,1,1,2,1)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,2,1)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,2,1)

END subroutine epf_z_cp

SUBROUTINE epf_z_cp_fc2(                                                 &
     nk,nome,ny,nz,lat,z,zt,nt,um,ptm,rhom,missv,                        &
     u_fc_r,pt_fc_r,rv_fc_r,rw_fc_r,u_fc_i,pt_fc_i,rv_fc_i,rw_fc_i,      &
     fy0,fz0,epd0,epd_z0,f_uv0,f_uw0,epd_uv0,epd_uw0 )

  integer,                        intent(in) ::  nk, nome, ny, nz, nt
  real,                           intent(in) ::  missv
  real, dimension(ny),            intent(in) ::  lat
  real, dimension(nz),            intent(in) ::  z, zt
  real, dimension(ny,nz,nt),      intent(in) ::  um, ptm, rhom
  real, dimension(nk,nome,ny,nz), intent(inout) ::  u_fc_r, pt_fc_r,     &
                     rv_fc_r, rw_fc_r, u_fc_i, pt_fc_i, rv_fc_i, rw_fc_i
  
  real, dimension(nk,nome,ny,nz), intent(out) ::  fy0, fz0, f_uv0, f_uw0
  real, dimension(nk,nome,ny,nz), intent(out) ::  epd0, epd_z0, epd_uv0, &
                                                  epd_uw0
 
  real, dimension(ny,nz) ::  fy, fz, f_uv, f_uw
  real, dimension(ny,nz) ::  epd, epd_z, epd_uv, epd_uw 
  real, dimension(ny,nz) ::  rvpt, rwpt, rvu, rwu, druwr3dz
  real, dimension(ny,nz) ::  phi1, phi2, phi1y, phi1z, phi2y, phi2z,     &
                             fphiy, fphiz, fphiz_t,                      &
                             dptmdy, dptmdz, divy_um_t, dumdz, rhomm
  real, dimension(ny,nz) ::  grady, gradz, temp, epd_y
  
  real, dimension(nk,nome,ny,nz) ::  tmp4d, gradz_fc_r, gradz_fc_i

  integer ::  j,k, i,n
 
  call set_gridvar_z_cp(ny,nz,lat,z,zt)

  rhomm(:,:) = sum(rhom(:,:,1:nt), dim=3)/float(nt)

  phi1y(:,:) = 0.  ;  phi1z(:,:) = 0.
  phi2y(:,:) = 0.  ;  phi2z(:,:) = 0.

  L_BTIME:  DO n=1, nt

  ! grad, pt
  call grady_2nd(ny,nz,ptm(:,:,n),lat, dptmdy)
  call missing_bdy(ny,nz,dptmdy,missv,1,1,0,0)

  temp(:,:) = log(ptm(:,:,n))
  call gradz_2nd_irr(ny,nz,temp,zt, dptmdz)
  dptmdz(:,:) = ptm(:,:,n)*dptmdz(:,:)
  call missing_bdy(ny,nz,dptmdz,missv,0,0,1,0)

  ! grad, u
  call grady_2nd(ny,nz,um(:,:,n)*coslat,lat, divy_um_t)
  do k=1, nz
  do j=2, ny-1
    divy_um_t(j,k) = divy_um_t(j,k)/coslat(j,k)
  enddo
  enddo
  divy_um_t(:,:) = divy_um_t(:,:)*r_earth(:,:)
  call intp_m2i(ny,nz,z,zt, divy_um_t)
  divy_um_t(:,:) = divy_um_t(:,:)/r_t(:,:)
  call gradz_m2i(ny,nz,um(:,:,n)*r_earth,z, dumdz)
  dumdz(:,:) = dumdz(:,:)/r_t(:,:)
  call missing_bdy(ny,nz,divy_um_t,missv,1,1,0,0)
  call missing_bdy(ny,nz,dumdz  ,missv,0,0,0,0)

  temp(:,:) = dptmdy(:,:)*dptmdy(:,:) + dptmdz(:,:)*dptmdz(:,:)
  phi1(:,:) = dptmdz(:,:)/temp(:,:)
  phi2(:,:) = dptmdy(:,:)/temp(:,:)

  dumdz    (:,:) = dumdz    (:,:)+rbeta(:,:)
  divy_um_t(:,:) = divy_um_t(:,:)-f    (:,:)

  phi1y(:,:) = phi1y(:,:) + phi1(:,:)*dumdz    (:,:)/float(nt)
  phi2y(:,:) = phi2y(:,:) + phi2(:,:)*dumdz    (:,:)/float(nt)
  phi1z(:,:) = phi1z(:,:) + phi1(:,:)*divy_um_t(:,:)/float(nt)
  phi2z(:,:) = phi2z(:,:) + phi2(:,:)*divy_um_t(:,:)/float(nt)

  ENDDO  L_BTIME

  call missing_bdy(ny,nz,phi1y,missv,1,1,1,0)
  call missing_bdy(ny,nz,phi1z,missv,1,1,1,0)
  call missing_bdy(ny,nz,phi2y,missv,1,1,1,0)
  call missing_bdy(ny,nz,phi2z,missv,1,1,1,0)

  rrhocos(:,:) = rcos(:,:)*rhomm(:,:)
  r3rho  (:,:) = r3  (:,:)*rhomm(:,:)

  fy0     = 0.  ;  fz0     = 0.  ;  epd0    = 0.  ;  epd_z0  = 0.
  f_uv0   = 0.  ;  f_uw0   = 0.  ;  epd_uv0 = 0.  ;  epd_uw0 = 0.

  tmp4d  (:,:,:,:) = rv_fc_r(:,:,:,:)*u_fc_r (:,:,:,:) + &
                     rv_fc_i(:,:,:,:)*u_fc_i (:,:,:,:)
  call intpwav_m2i(nk,nome,ny,nz, rv_fc_r,rv_fc_i)
  rv_fc_r(:,:,:,:) = rv_fc_r(:,:,:,:)*pt_fc_r(:,:,:,:) + &
                     rv_fc_i(:,:,:,:)*pt_fc_i(:,:,:,:)
  rv_fc_i(:,:,:,:) = tmp4d  (:,:,:,:)
  pt_fc_r(:,:,:,:) = rw_fc_r(:,:,:,:)*pt_fc_r(:,:,:,:) + &
                     rw_fc_i(:,:,:,:)*pt_fc_i(:,:,:,:)
  tmp4d  (:,:,:,:) = rw_fc_r(:,:,:,:)
  pt_fc_i(:,:,:,:) = rw_fc_i(:,:,:,:)
  call intpwav_i2m(nk,nome,ny,nz, tmp4d,pt_fc_i)
  pt_fc_i(:,:,:,:) = tmp4d  (:,:,:,:)*u_fc_r (:,:,:,:) + &
                     pt_fc_i(:,:,:,:)*u_fc_i (:,:,:,:)

  u_fc_r (:,:,:,:) = u_fc_r (:,:,:,:)*spread(spread(r_earth,1,nk),2,nome)
  u_fc_i (:,:,:,:) = u_fc_i (:,:,:,:)*spread(spread(r_earth,1,nk),2,nome)
  rw_fc_r(:,:,:,:) = rw_fc_r(:,:,:,:)*spread(spread(r2_t   ,1,nk),2,nome)
  rw_fc_i(:,:,:,:) = rw_fc_i(:,:,:,:)*spread(spread(r2_t   ,1,nk),2,nome)
  call gradzwav_m2i(nk,nome,ny,nz,u_fc_r,u_fc_i,z, gradz_fc_r,gradz_fc_i)
  tmp4d(:,:,:,:) = rw_fc_r(:,:,:,:)*gradz_fc_r(:,:,:,:) +                &
                   rw_fc_i(:,:,:,:)*gradz_fc_i(:,:,:,:)
  call intp_i2m(nk,nome,ny,nz, tmp4d)
  call gradz_i2m(nk,nome,ny,nz,rw_fc_r,zt, gradz_fc_r)
  call gradz_i2m(nk,nome,ny,nz,rw_fc_i,zt, gradz_fc_i)
  rw_fc_r(:,:,:,:) = tmp4d(:,:,:,:) +                                    &
                     ( u_fc_r(:,:,:,:)*gradz_fc_r(:,:,:,:) +             &
                       u_fc_i(:,:,:,:)*gradz_fc_i(:,:,:,:) )

  L_OME:  DO n=1, nome
  L_KWN:  DO i=2, nk

  rvpt    (:,:) = rv_fc_r(i,n,:,:)
  rwpt    (:,:) = pt_fc_r(i,n,:,:)
  rvu     (:,:) = rv_fc_i(i,n,:,:)
  rwu     (:,:) = pt_fc_i(i,n,:,:)
  druwr3dz(:,:) = rw_fc_r(i,n,:,:)

  if ( all(rvpt == rwpt) .and. all(rvu == rwu) )  CYCLE
  ! all-zero

  ! phi
  fphiy(:,:) = rvpt(:,:)*phi1y(:,:) - rwpt(:,:)*phi2y(:,:)
  call intp_i2m(ny,nz, fphiy)
  fphiz(:,:) = rvpt(:,:)*phi1z(:,:) - rwpt(:,:)*phi2z(:,:)
  fphiz_t(:,:) = fphiz(:,:)
  call intp_i2m(ny,nz, fphiz)
  call missing_bdy(ny,nz,fphiy,missv,1,1,1,0)
  call missing_bdy(ny,nz,fphiz,missv,1,1,1,0)

  ! epf, epd
  fy  (:,:) = -rcos(:,:)*(rvu(:,:) - fphiy(:,:))
  fz  (:,:) = -rcos(:,:)*(rwu(:,:) + fphiz(:,:))
  f_uv(:,:) = -rcos(:,:)*rvu(:,:)
  f_uw(:,:) = -rcos(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz  ,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uv,missv,1,1,1,0)
  call missing_bdy(ny,nz,f_uw,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, grady)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = grady(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,2,2,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,2,2,1,0)

  epd_uw(:,:) = -druwr3dz(:,:)/r3rho(:,:) * 86400.
  call gradz_i2m(ny,nz,r3_t*fphiz_t,zt, gradz)
  epd_z(:,:) = epd_uw(:,:) - gradz(:,:)/r3rho(:,:)*86400.
  call missing_bdy(ny,nz,epd_z ,missv,1,1,2,1)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,2,1)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,2,1)

  fy0    (i,n,:,:) = fy    (:,:)
  fz0    (i,n,:,:) = fz    (:,:)
  epd0   (i,n,:,:) = epd   (:,:)
  epd_z0 (i,n,:,:) = epd_z (:,:)
  f_uv0  (i,n,:,:) = f_uv  (:,:)
  f_uw0  (i,n,:,:) = f_uw  (:,:)
  epd_uv0(i,n,:,:) = epd_uv(:,:)
  epd_uw0(i,n,:,:) = epd_uw(:,:)

  ENDDO  L_KWN
  ENDDO  L_OME

END subroutine epf_z_cp_fc2

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

SUBROUTINE gradz_i2m_4d(nk,no,ny,nz,var,zi, gradz)

  integer,                      intent(in)  ::  nk, no, ny, nz
  real, dimension(nk,no,ny,nz), intent(in)  ::  var
  real, dimension(nz),          intent(in)  ::  zi
  real, dimension(nk,no,ny,nz), intent(out) ::  gradz

  integer ::  k
  real    ::  inv_dz(nz)

  do k=2, nz
    inv_dz(k) = 1./(zi(k) - zi(k-1))
  enddo
  do k=2, nz
    gradz(:,:,:,k) = (var(:,:,:,k) - var(:,:,:,k-1))*inv_dz(k)
  enddo
  gradz(:,:,:,1) = gradz(:,:,:,2)

END subroutine gradz_i2m_4d

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

SUBROUTINE intp_i2m_2d(ny,nz, var)

  integer,                intent(in)    ::  ny, nz
  real, dimension(ny,nz), intent(inout) ::  var

  call intp_i2m_4d(ny,1,1,nz, var)

END subroutine intp_i2m_2d

SUBROUTINE intp_i2m_4d(nk,no,ny,nz, var)

  integer,                      intent(in)    ::  nk, no, ny, nz
  real, dimension(nk,no,ny,nz), intent(inout) ::  var

  integer ::  k

  do k=nz, 2, -1
    var(:,:,:,k) = 0.5*(var(:,:,:,k-1) + var(:,:,:,k))
  enddo

END subroutine intp_i2m_4d

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

END module epf


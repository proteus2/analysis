MODULE epf

  use util,  only: grady_2nd, gradz_2nd_irr, missing_bdy
  use const_glob,  only:  g, kappa, a_earth, ome_earth, deg2rad

  implicit none

  private ::  set_gridvar_p, set_gridvar_z
  private ::  g, kappa, a_earth, ome_earth, deg2rad

  character(len=64), dimension(8) ::  varname_epf =                      &
      (/'f_y    ','f_z    ','epd    ','epd_z  ','f_uv   ',               &
        'f_uw   ','epd_uv ','epd_uw '/)

  real, dimension(:,:),   allocatable ::  dptmdz_from_thlev
  real, dimension(:,:,:), allocatable ::  dptmdz_from_thlev3

  integer, dimension(4,8) ::  jk_bnd

  integer,                           private ::  ny_pre, nz_pre
  real, dimension(:),   allocatable, private ::  lat_pre, ht_pre, zp
  real, dimension(:,:), allocatable, private ::  coslat, f, rho0, rbeta, &
                                                 r_earth, rcos, rrhocos, &
                                                 r3, r3rho, r_t, r2_t,   &
                                                 r3_t


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
  real, dimension(ny,nz) ::  dptmdy, dptmdz, divy_um, dumdz
  real, dimension(ny,nz) ::  temp, epd_y

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

  ! epf, epd
  fy  (:,:) = -a_earth*coslat(:,:)*(rvu(:,:) - phi(:,:)*dumdz(:,:))
  fz  (:,:) = -a_earth*coslat(:,:)*(rwu(:,:) + phi(:,:)*(divy_um(:,:)-f(:,:)))
  f_uv(:,:) = -a_earth*coslat(:,:)*rvu(:,:)
  f_uw(:,:) = -a_earth*coslat(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy,missv,0,0,1,0)
!mv  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, epd_uv)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = epd_uv(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) &
                  *86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,1,1,1,0)

  call gradz_2nd_irr(ny,nz,fz,zp, epd_z)
  epd_z(:,:) = epd_z(:,:)/(rho0(:,:)*a_earth*coslat(:,:)) * 86400.

  call gradz_2nd_irr(ny,nz,f_uw,zp, epd_uw)
  epd_uw(:,:) = epd_uw(:,:)/(rho0(:,:)*a_earth*coslat(:,:)) * 86400.
  call missing_bdy(ny,nz,epd_z ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,1,0)

  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,1,1,1,0)

  jk_bnd(:,1) = (/0,0,1,0/)  ! fy
  jk_bnd(:,2) = (/1,1,1,0/)  ! fz
  jk_bnd(:,3) = (/1,1,1,0/)  ! epd
  jk_bnd(:,4) = (/1,1,1,0/)  ! epd_z
  jk_bnd(:,5) = (/0,0,0,0/)  ! f_uv
  jk_bnd(:,6) = (/0,0,0,0/)  ! f_uw
  jk_bnd(:,7) = (/1,1,1,0/)  ! epd_uv
  jk_bnd(:,8) = (/1,1,1,0/)  ! epd_uw

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
  real, dimension(ny,nz) ::  temp, epd_y
  
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

  ! epf, epd
  fy  (:,:) = -rcos(:,:)*(rvu(:,:) - phi(:,:)*(dumdz  (:,:)+rbeta(:,:)))
  fz  (:,:) = -rcos(:,:)*(rwu(:,:) + phi(:,:)*(divy_um(:,:)-f    (:,:)))
  f_uv(:,:) = -rcos(:,:)*rvu(:,:)
  f_uw(:,:) = -rcos(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  rrhocos(:,:) = rcos(:,:)*rhom(:,:)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, epd_uv)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = epd_uv(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,2,2,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,2,2,1,0)

  temp(:,:) = r_earth(:,:)*r_earth(:,:)
  call gradz_2nd_irr(ny,nz,fz*temp,z, epd_z)
  do k=1, nz
  do j=2, ny-1
    epd_z(j,k) = epd_z(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,f_uw*temp,z, epd_uw)
  do k=1, nz
  do j=2, ny-1
    epd_uw(j,k) = epd_uw(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_z ,missv,1,1,2,1)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,2,1)

  epd(:,:) = epd_y(:,:) + epd_z(:,:)
  call missing_bdy(ny,nz,epd,missv,2,2,2,1)

  jk_bnd(:,1) = (/1,1,1,0/)  ! fy
  jk_bnd(:,2) = (/1,1,1,0/)  ! fz
  jk_bnd(:,3) = (/2,2,2,1/)  ! epd
  jk_bnd(:,4) = (/1,1,2,1/)  ! epd_z
  jk_bnd(:,5) = (/0,0,0,0/)  ! f_uv
  jk_bnd(:,6) = (/0,0,0,0/)  ! f_uw
  jk_bnd(:,7) = (/2,2,1,0/)  ! epd_uv
  jk_bnd(:,8) = (/1,1,2,1/)  ! epd_uw

END subroutine epf_z

SUBROUTINE epf_hydro_p_fc2(                                              &
     nk,nome,ny,nz,lat,p,nt,um,ptm,h_scale,vpt,vu,wpt,wu,missv,          &
     fy0,fz0,epd0,epd_z0,f_uv0,f_uw0,epd_uv0,epd_uw0 )

  integer,                   intent(in) ::  nk, nome, ny, nz, nt
  real,                      intent(in) ::  h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(ny,nz,nt), intent(in) ::  um, ptm
  real, dimension(nk,nome,ny,nz), intent(in) ::  vpt, vu, wpt, wu

  real, dimension(nk,nome,ny,nz), intent(out) ::  fy0, fz0, f_uv0, f_uw0
  real, dimension(nk,nome,ny,nz), intent(out) ::  epd0, epd_z0, epd_uv0, &
                                                  epd_uw0

  real, dimension(ny,nz) ::  fy, fz, f_uv, f_uw
  real, dimension(ny,nz) ::  epd, epd_z, epd_uv, epd_uw 
  real, dimension(ny,nz) ::  rvpt, rwpt, rvu, rwu, fphiy, fphiz
  real, dimension(ny,nz) ::  phiy, phiz
  real, dimension(ny,nz) ::  dptmdy, dptmdz, divy_um, dumdz
  real, dimension(ny,nz) ::  temp, epd_y

  integer ::  j,k, i,n

  call set_gridvar_p(ny,nz,lat,p,h_scale)

  phiy(:,:) = 0.  ;  phiz(:,:) = 0.

  L_BTIME:  DO n=1, nt

  ! grad, pt
!  call grady_2nd(ny,nz,ptm(:,:,n),lat, dptmdy)

  temp(:,:) = log(ptm(:,:,n))
  call gradz_2nd_irr(ny,nz,temp,zp, dptmdz)
  dptmdz(:,:) = ptm(:,:,n)*dptmdz(:,:)

  ! grad, u
  call grady_2nd(ny,nz,um(:,:,n)*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um(:,:,n),zp, dumdz)

  divy_um(:,:) = divy_um(:,:) - f(:,:)

  phiy(:,:) = phiy(:,:) + dumdz  (:,:)/dptmdz(:,:)/float(nt)
  phiz(:,:) = phiz(:,:) + divy_um(:,:)/dptmdz(:,:)/float(nt)

  ENDDO  L_BTIME

!!  call missing_bdy(ny,nz,phiy,missv,0,0,1,0)
!!  call missing_bdy(ny,nz,phiz,missv,1,1,1,0)

  fy0     = 0.  ;  fz0     = 0.  ;  epd0    = 0.  ;  epd_z0  = 0.
  f_uv0   = 0.  ;  f_uw0   = 0.  ;  epd_uv0 = 0.  ;  epd_uw0 = 0.

  L_OME:  DO n=1, nome
  L_KWN:  DO i=2, nk

  rvpt(:,:) = rho0(:,:)*vpt(i,n,:,:)
  rwpt(:,:) = rho0(:,:)*wpt(i,n,:,:)
  rvu (:,:) = rho0(:,:)*vu (i,n,:,:)
  rwu (:,:) = rho0(:,:)*wu (i,n,:,:)

  if (sum(abs(rvpt)) == 0.)  CYCLE

  ! phi
  fphiy(:,:) = rvpt(:,:)*phiy(:,:)
  fphiz(:,:) = rvpt(:,:)*phiz(:,:)
!!  call missing_bdy(ny,nz,fphiy,missv,0,0,1,0)
!!  call missing_bdy(ny,nz,fphiz,missv,1,1,1,0)

  ! epf, epd
  fy  (:,:) = -a_earth*coslat(:,:)*(rvu(:,:) - fphiy(:,:))
  fz  (:,:) = -a_earth*coslat(:,:)*(rwu(:,:) + fphiz(:,:))
  f_uv(:,:) = -a_earth*coslat(:,:)*rvu(:,:)
  f_uw(:,:) = -a_earth*coslat(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy,missv,0,0,1,0)
!mv  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, epd_uv)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = epd_uv(j,k)/a_earth/(coslat(j,k)*coslat(j,k))/rho0(j,k) &
                  *86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,1,1,1,0)

  call gradz_2nd_irr(ny,nz,fz,zp, epd_z)
  epd_z(:,:) = epd_z(:,:)/(rho0(:,:)*a_earth*coslat(:,:)) * 86400.
  call gradz_2nd_irr(ny,nz,f_uw,zp, epd_uw)
  epd_uw(:,:) = epd_uw(:,:)/(rho0(:,:)*a_earth*coslat(:,:)) * 86400.
  call missing_bdy(ny,nz,epd_z ,missv,1,1,1,0)
  call missing_bdy(ny,nz,epd_uw,missv,1,1,1,0)

  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

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

  jk_bnd(:,1) = (/0,0,1,0/)  ! fy
  jk_bnd(:,2) = (/1,1,1,0/)  ! fz
  jk_bnd(:,3) = (/1,1,1,0/)  ! epd
  jk_bnd(:,4) = (/1,1,1,0/)  ! epd_z
  jk_bnd(:,5) = (/0,0,0,0/)  ! f_uv
  jk_bnd(:,6) = (/0,0,0,0/)  ! f_uw
  jk_bnd(:,7) = (/1,1,1,0/)  ! epd_uv
  jk_bnd(:,8) = (/1,1,1,0/)  ! epd_uw

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
  real, dimension(ny,nz) ::  temp, epd_y
  
  integer ::  j,k, i,n
 
  call set_gridvar_z(ny,nz,lat,z)

  rhomm(:,:) = sum(rhom(:,:,1:nt), dim=3)/float(nt)

  phi1y(:,:) = 0.  ;  phi1z(:,:) = 0.
  phi2y(:,:) = 0.  ;  phi2z(:,:) = 0.

  L_BTIME:  DO n=1, nt

  ! grad, pt
  call grady_2nd(ny,nz,ptm(:,:,n),lat, dptmdy)

  if ( allocated(dptmdz_from_thlev3) ) then
    dptmdz(:,:) = dptmdz_from_thlev3(:,:,n)
  else
    temp(:,:) = log(ptm(:,:,n))
    call gradz_2nd_irr(ny,nz,temp,z, dptmdz)
    dptmdz(:,:) = ptm(:,:,n)*dptmdz(:,:)
  end if

  ! grad, u
  call grady_2nd(ny,nz,um(:,:,n)*coslat,lat, divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/coslat(j,k)
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,um(:,:,n)*r_earth,z, dumdz)
  dumdz(:,:) = dumdz(:,:)/r_earth(:,:)

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

!  call missing_bdy(ny,nz,phi1y,missv,1,1,1,0)
!  call missing_bdy(ny,nz,phi1z,missv,1,1,1,0)
!  call missing_bdy(ny,nz,phi2y,missv,1,1,1,0)
!  call missing_bdy(ny,nz,phi2z,missv,1,1,1,0)

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

  if (sum(abs(rvpt)) == 0.)  CYCLE

  ! phi
  fphiy(:,:) = rvpt(:,:)*phi1y(:,:) - rwpt(:,:)*phi2y(:,:)
  fphiz(:,:) = rvpt(:,:)*phi1z(:,:) - rwpt(:,:)*phi2z(:,:)
!  call missing_bdy(ny,nz,fphiy,missv,1,1,1,0)
!  call missing_bdy(ny,nz,fphiz,missv,1,1,1,0)

  ! epf, epd
  fy  (:,:) = -rcos(:,:)*(rvu(:,:) - fphiy(:,:))
  fz  (:,:) = -rcos(:,:)*(rwu(:,:) + fphiz(:,:))
  f_uv(:,:) = -rcos(:,:)*rvu(:,:)
  f_uw(:,:) = -rcos(:,:)*rwu(:,:)
  call missing_bdy(ny,nz,fy,missv,1,1,1,0)
  call missing_bdy(ny,nz,fz,missv,1,1,1,0)

  call grady_2nd(ny,nz,fy*coslat,lat, epd_y)
  do k=1, nz
  do j=2, ny-1
    epd_y(j,k) = epd_y(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call grady_2nd(ny,nz,f_uv*coslat,lat, epd_uv)
  do k=1, nz
  do j=2, ny-1
    epd_uv(j,k) = epd_uv(j,k)/coslat(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call missing_bdy(ny,nz,epd_y ,missv,2,2,1,0)
  call missing_bdy(ny,nz,epd_uv,missv,2,2,1,0)

  temp(:,:) = r_earth(:,:)*r_earth(:,:)
  call gradz_2nd_irr(ny,nz,fz*temp,z, epd_z)
  do k=1, nz
  do j=2, ny-1
    epd_z(j,k) = epd_z(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
  enddo
  enddo
  call gradz_2nd_irr(ny,nz,f_uw*temp,z, epd_uw)
  do k=1, nz
  do j=2, ny-1
    epd_uw(j,k) = epd_uw(j,k)/temp(j,k)/rrhocos(j,k) * 86400.
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

  jk_bnd(:,1) = (/1,1,1,0/)  ! fy
  jk_bnd(:,2) = (/1,1,1,0/)  ! fz
  jk_bnd(:,3) = (/2,2,2,1/)  ! epd
  jk_bnd(:,4) = (/1,1,2,1/)  ! epd_z
  jk_bnd(:,5) = (/0,0,0,0/)  ! f_uv
  jk_bnd(:,6) = (/0,0,0,0/)  ! f_uw
  jk_bnd(:,7) = (/2,2,1,0/)  ! epd_uv
  jk_bnd(:,8) = (/1,1,2,1/)  ! epd_uw

END subroutine epf_z_fc2

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

END module epf


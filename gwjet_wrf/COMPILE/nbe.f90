MODULE nbe

  implicit none

  logical, save ::  l_gp_discont_lon = .False.

  character(len=64), dimension(5) ::  varname_nbe_p =                    &
      (/'div    ','nbe    ','term2  ','curv   ','adv    '/)
  character(len=64), dimension(8) ::  varname_nbe_z =                    &
      (/'div    ','nbe_maj','ageo   ','two_jac','nl_min ','curv   ','adv    ','fvor   '/)
!      (/'div    ','nbe    ','term2  ','curv   ','adv    ','nonl   '/)

  private ::  set_gridvar
  private ::  gradx_2nd, grady_2nd, grady_cos_2nd, missing_bdy

  integer,                           private ::  nx_pre, ny_pre
  real, dimension(:),   allocatable, private ::  lat_pre, inv_2dx,       &
                                                 inv_2dy
  real, dimension(:,:), allocatable, private ::  coslat, sinlat, f, beta

  real, parameter, private ::  a_earth = 6371229.
  real, parameter, private ::  pi = 3.14159265358979323846
  real, parameter, private ::  ome_earth = 7.292116e-5
  real, parameter, private ::  deg2rad = pi/180.


  CONTAINS


SUBROUTINE nbe_p(                                                        &
     nx,ny,lat,u,v,gp,omega,us,vs,missv,                                 &
     div,nbe,term2,curv,adv )

  integer,                intent(in) ::  nx, ny
  real,                   intent(in) ::  missv
  real, dimension(ny),    intent(in) ::  lat
  real, dimension(nx,ny), intent(in) ::  u, v, gp, omega, us, vs

  real, dimension(nx,ny), intent(out) ::  div, nbe, term2, curv, adv

  call set_gridvar(nx,ny,lat)

  div(:,:) = gradx_2nd(u) + grady_cos_2nd(v)

  nbe(:,:) = (gradx_2nd(f*v) - grady_cos_2nd(f*u))                       &
           + 2.*(gradx_2nd(u)*grady_2nd(v) - gradx_2nd(v)*grady_2nd(u))  &
           - ( gradx2_2nd(gp) + grady_cos_2nd(grady_2nd(gp)) )

  term2(:,:) = (-1.)*div(:,:)*div(:,:) - gradx_2nd(omega)*us(:,:)        &
             - grady_2nd(omega)*vs(:,:)

  curv(:,:) = (-1.)*grady_2nd((u(:,:)*u(:,:) + v(:,:)*v(:,:))*           &
              sinlat(:,:))/(a_earth*coslat(:,:))

  adv(:,:) = (-1.)*( u(:,:)*gradx_2nd(div) + v(:,:)*grady_2nd(div)       &
           + omega(:,:)*(gradx_2nd(us) + grady_cos_2nd(vs)) )

  call missing_bdy(nx,ny,div  ,missv,0,0,1,1)
  call missing_bdy(nx,ny,nbe  ,missv,0,0,2,2)
  call missing_bdy(nx,ny,term2,missv,0,0,1,1)
  call missing_bdy(nx,ny,curv ,missv,0,0,1,1)
  call missing_bdy(nx,ny,adv  ,missv,0,0,2,2)

END subroutine nbe_p

SUBROUTINE nbe_z(                                                        &
     nx,ny,lat,u,v,p,inv_rho,w,us,vs,missv,n_smth9,                      &
     div,nbe,adv,nbe_maj,curv,nl_min,ageo,two_jac,fvor)

  integer,                intent(in) ::  nx, ny, n_smth9
  real,                   intent(in) ::  missv
  real, dimension(ny),    intent(in) ::  lat
  real, dimension(nx,ny), intent(in) ::  u, v, p, inv_rho, w, us, vs

  real, dimension(nx,ny), intent(out) ::  div, nbe, adv, nbe_maj, curv,  &
                                          nl_min, ageo, two_jac, fvor

  integer                ::  j, is
  real, dimension(nx,ny) ::  div_pgf, tmp

  call set_gridvar(nx,ny,lat)

  div(:,:) = gradx_2nd(u) + grady_cos_2nd(v)

  div_pgf(:,:) = (-1.)*( gradx2c_2nd(p,inv_rho) +                        &
                         grady2c_2nd(p,coslat*inv_rho)/coslat )

  ! 9-pt local smoothing (1-2-1 smoothing)
  do is=1, n_smth9
    tmp = div_pgf
    div_pgf(2:nx-1,:) = 0.25*(2.*tmp(2:nx-1,:) + tmp(1:nx-2,:) + tmp(3:nx,:))
    div_pgf(1 ,:) = 0.25*(2.*tmp(1 ,:) + tmp(nx  ,:) + tmp(2,:))
    div_pgf(nx,:) = 0.25*(2.*tmp(nx,:) + tmp(nx-1,:) + tmp(1,:))
    tmp = div_pgf
    do j=2, ny-1
      div_pgf(:,j) = 0.25*(2.*tmp(:,j) + tmp(:,j-1) + tmp(:,j+1))
    enddo
  enddo

  fvor(:,:) = f*(gradx_2nd(v) - grady_cos_2nd(u)) - beta*u
!!  fvor(:,:) = f*(gradx_2nd(v) - grady_cos_2nd(u)) - beta*(u + a_earth*gradx_2nd(w))

  ageo(:,:) = div_pgf + fvor

  two_jac(:,:) = 2.*(gradx_2nd(u)*grady_2nd(v) - gradx_2nd(v)*grady_2nd(u))

  nbe_maj(:,:) = ageo(:,:) + two_jac(:,:)

  nl_min(:,:) = (-1.)*div(:,:)*div(:,:)                                  &
              - (gradx_2nd(w)*us(:,:) + grady_2nd(w)*vs(:,:))
!!  nl_min(:,:) = (-1.)*div(:,:)*(div(:,:)+w(:,:)/a_earth)                 &
!!              - (gradx_2nd(w)*(us(:,:)+u(:,:)/a_earth) +                 &
!!                 grady_2nd(w)*(vs(:,:)+v(:,:)/a_earth))

  curv(:,:) = (-1.)*grady_2nd((u(:,:)*u(:,:) + v(:,:)*v(:,:))*           &
              sinlat(:,:))/(a_earth*coslat(:,:))

  nbe(:,:) = nbe_maj(:,:) + curv(:,:) + nl_min(:,:)

  adv(:,:) = (-1.)*( u(:,:)*gradx_2nd(div) + v(:,:)*grady_2nd(div)       &
           + w(:,:)*(gradx_2nd(us) + grady_cos_2nd(vs)) )

  call missing_bdy(nx,ny,div    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,nbe    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,nbe_maj,missv,0,0,1,1)
  call missing_bdy(nx,ny,ageo   ,missv,0,0,1,1)
  call missing_bdy(nx,ny,two_jac,missv,0,0,1,1)
  call missing_bdy(nx,ny,nl_min ,missv,0,0,1,1)
  call missing_bdy(nx,ny,curv   ,missv,0,0,1,1)
  call missing_bdy(nx,ny,adv    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,fvor   ,missv,0,0,1,1)

END subroutine nbe_z

SUBROUTINE set_gridvar(nx,ny,lat)

  integer,             intent(in) ::  nx, ny
  real, dimension(ny), intent(in) ::  lat

  if ( allocated(lat_pre) ) then
    if ( nx == nx_pre .and. ny == ny_pre ) then
      if ( all(lat == lat_pre) )  RETURN
    end if
    deallocate( lat_pre, coslat, sinlat, f, beta, inv_2dx, inv_2dy )
  end if

  allocate( coslat(nx,ny), sinlat(nx,ny), f(nx,ny), beta(nx,ny),         &
            inv_2dx(ny), inv_2dy(ny) )
  allocate( lat_pre(ny) )

  coslat(:,:) = spread(cos(lat(:)*deg2rad),1,nx)
  sinlat(:,:) = spread(sin(lat(:)*deg2rad),1,nx)
  if (abs(lat(1 )) == 90.)  coslat(:,1 ) = 0.
  if (abs(lat(ny)) == 90.)  coslat(:,ny) = 0.
  f   (:,:) = spread(2.*ome_earth*sinlat(1,:)        ,1,nx)
  beta(:,:) = spread(2.*ome_earth*coslat(1,:)/a_earth,1,nx)

  inv_2dx(2:ny-1) = 0.5/(360./float(nx)*deg2rad*a_earth*coslat(1,2:ny-1))
  inv_2dy(2:ny-1) = 1./((lat(3:ny)-lat(1:ny-2))*deg2rad*a_earth)

  nx_pre = nx  ;  ny_pre = ny
  lat_pre(:) = lat(:)

END subroutine set_gridvar

FUNCTION gradx_2nd(var)

  real, dimension(:,:), intent(in)  ::  var

  real, dimension(size(var,1),size(var,2)) ::  gradx_2nd

  integer ::  i, nxi, nyi
  real    ::  inv_2dxp(size(var,2))

  nxi = size(var,1)  ;  nyi = size(var,2)

  inv_2dxp(:) = inv_2dx(:)
  inv_2dxp(1) = 1.e20  ;  inv_2dxp(nyi) = 1.e20
  do i=2, nxi-1
    gradx_2nd(i,:) = (var(i+1,:)-var(i-1,:))*inv_2dxp(:)
  enddo
  gradx_2nd(1  ,:) = (var(2,:)-var(nxi  ,:))*inv_2dxp(:)
  gradx_2nd(nxi,:) = (var(1,:)-var(nxi-1,:))*inv_2dxp(:)
  gradx_2nd(:,1  ) = 0.
  gradx_2nd(:,nyi) = 0.

END function gradx_2nd

FUNCTION grady_2nd(var)

  real, dimension(:,:), intent(in)  ::  var

  real, dimension(size(var,1),size(var,2)) ::  grady_2nd

  integer ::  j, nyi

  nyi = size(var,2)

  do j=2, nyi-1
    grady_2nd(:,j) = (var(:,j+1)-var(:,j-1))*inv_2dy(j)
  enddo
  grady_2nd(:,1  ) = 0.
  grady_2nd(:,nyi) = 0.

END function grady_2nd

FUNCTION grady_cos_2nd(var)

  real, dimension(:,:), intent(in)  ::  var

  real, dimension(size(var,1),size(var,2)) ::  grady_cos_2nd

  integer ::  nyi

  nyi = size(var,2)

  grady_cos_2nd(:,:) = grady_2nd(var*coslat)

  grady_cos_2nd(:,2:nyi-1) = grady_cos_2nd(:,2:nyi-1)/coslat(:,2:nyi-1)

END function grady_cos_2nd

FUNCTION gradx2_2nd(var)  ! correction needed (for 2nd order)

  real, dimension(:,:), intent(in)  ::  var

  real, dimension(size(var,1),size(var,2)) ::  gradx2_2nd

  integer ::  i, nxi, nyi
  real    ::  inv_2dxp(size(var,2)), gradx(size(var,1),size(var,2)),     &
              tmp(size(var,1),size(var,2),2)

  nxi = size(var,1)  ;  nyi = size(var,2)

  inv_2dxp(:) = inv_2dx(:)
  inv_2dxp(1) = 1.e20  ;  inv_2dxp(nyi) = 1.e20

  do i=2, nxi-1
    gradx(i,:) = (var(i+1,:)-var(i-1,:))
  enddo
  gradx(1  ,:) = (var(2,:)-var(nxi  ,:))
  gradx(nxi,:) = (var(1,:)-var(nxi-1,:))
  gradx(:,1  ) = 0.
  gradx(:,nyi) = 0.

  if ( l_gp_discont_lon ) then  ! e.g., MERRA GP
    if ( var(nxi/2+1,nyi/2+1) > 1.e3 ) then  ! var: gp
      tmp(:,:,1) = abs(gradx(:,:))
      where ( tmp(:,:,1) < 1.e-10 )
        tmp(:,:,1) = 1.e20
      end where
      tmp(:,:,2) = spread( minval(tmp(:,:,1), dim=1), 1,nxi )
      tmp(:,:,1) = abs(gradx(:,:))
      where ( tmp(:,:,1) < 10.*tmp(:,:,2) )
        gradx(:,:) = 0.
      end where
    end if
  end if

  gradx(:,:) = gradx(:,:)*spread(inv_2dxp(:),1,nxi)

  gradx2_2nd(:,:) = gradx_2nd(gradx)

END function gradx2_2nd

FUNCTION gradx2c_2nd(var,c)

  real, dimension(:,:), intent(in)  ::  var, c

  real, dimension(size(var,1),size(var,2)) ::  gradx2c_2nd

  integer ::  j, nxi, nyi
  real    ::  inv_dx2p(size(var,2)), gradx(0:size(var,1))

  nxi = size(var,1)  ;  nyi = size(var,2)

  inv_dx2p(:) = (inv_2dx(:)*2.)**2
!  inv_dxp(1) = 1.e20  ;  inv_dxp(nyi) = 1.e20
  do j=2, nyi-1
    gradx(1:nxi-1) = (var(2:nxi,j)-var(1:nxi-1,j))*    &
                     (c  (2:nxi,j)+c  (1:nxi-1,j))*0.5
    gradx(nxi) = (var(1,j)-var(nxi,j))*    &
                 (c  (1,j)+c  (nxi,j))*0.5
    gradx(0) = gradx(nxi)
    gradx2c_2nd(1:nxi,j) = (gradx(1:nxi)-gradx(0:nxi-1))*inv_dx2p(j)
  enddo
  gradx2c_2nd(:,1) = 0.  ;  gradx2c_2nd(:,nyi) = 0.

END function gradx2c_2nd

FUNCTION grady2c_2nd(var,c)

  real, dimension(:,:), intent(in)  ::  var, c

  real, dimension(size(var,1),size(var,2)) ::  grady2c_2nd

  integer ::  j, nxi, nyi
  real    ::  inv_dyp(size(var,2)), inv_dyp_i(size(var,2)), &
              grady(size(var,1),size(var,2))

  nxi = size(var,1)  ;  nyi = size(var,2)

  inv_dyp_i(1:nyi-1) = 1./((lat_pre(2:nyi)-lat_pre(1:nyi-1))*deg2rad*a_earth)
  inv_dyp(:) = inv_2dy(:)*2.
  do j=1, nyi-1
    grady(:,j) = (var(:,j+1)-var(:,j))*inv_dyp_i(j)* &
                 (c  (:,j+1)+c  (:,j))*0.5
  enddo
  do j=2, nyi-1
    grady2c_2nd(:,j) = (grady(:,j)-grady(:,j-1))*inv_dyp(j)
    ! This centered-method calculation is accurate only for a uniform 
    ! y-grid size, say dy(j) = dy(j+1)
  enddo
  grady2c_2nd(:,1) = 0.  ;  grady2c_2nd(:,nyi) = 0.

END function grady2c_2nd

SUBROUTINE missing_bdy(nx,ny,var,missv,nm_x1,nm_x2,nm_y1,nm_y2)

  integer, intent(in) ::  nx, ny, nm_x1, nm_x2, nm_y1, nm_y2
  real,    intent(in) ::  missv

  real, dimension(nx,ny), intent(inout) ::  var

  if (nm_x1 > 0)  var(1         :nm_x1,          :     ) = missv
  if (nm_x2 > 0)  var(nx+1-nm_x2:nx   ,          :     ) = missv
  if (nm_y1 > 0)  var(          :     ,1         :nm_y1) = missv
  if (nm_y2 > 0)  var(          :     ,ny+1-nm_y2:ny   ) = missv

END subroutine missing_bdy

END module nbe


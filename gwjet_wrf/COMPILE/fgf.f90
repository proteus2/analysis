MODULE fgf
! adiabatic

  implicit none

  character(len=64), dimension(6) ::  varname_fgf_2d =                   &
      (/'pt_grad','fgf    ','div    ','def    ','fgf_div','fgf_def'/)

  private ::  set_gridvar
  private ::  gradx_2nd, grady_2nd, grady_cos_2nd, missing_bdy

  integer,                           private ::  nx_pre, ny_pre
  real, dimension(:),   allocatable, private ::  lat_pre, inv_2dx,       &
                                                 inv_2dy
  real, dimension(:,:), allocatable, private ::  coslat, tanlat

  real, parameter, private ::  a_earth = 6371229.
  real, parameter, private ::  deg2rad = 3.14159265358979323846/180.


  CONTAINS


SUBROUTINE fgf_2d(                                                       &
     nx,ny,lat,u,v,pt,missv,                                             &
     pt_grad,fgf,div,def,fgf_div,fgf_def )

  integer,                intent(in) ::  nx, ny
  real,                   intent(in) ::  missv
  real, dimension(ny),    intent(in) ::  lat
  real, dimension(nx,ny), intent(in) ::  u, v, pt

  real, dimension(nx,ny), intent(out) ::  pt_grad, fgf, div, def,        &
                                          fgf_div, fgf_def

  real, dimension(nx,ny) ::  pt_x, pt_y, pt_x_2, pt_y_2, def_st, def_sh

  call set_gridvar(nx,ny,lat)

  pt_x(:,:) = gradx_2nd(pt)
  pt_y(:,:) = grady_2nd(pt)
  pt_x_2(:,:) = pt_x(:,:)*pt_x(:,:)
  pt_y_2(:,:) = pt_y(:,:)*pt_y(:,:)

  div(:,:) = gradx_2nd(u) + grady_cos_2nd(v)

  pt_grad(:,:) = pt_x_2(:,:) + pt_y_2(:,:)
  fgf_div(:,:) = -0.5*div(:,:)*pt_grad(:,:)
  pt_grad(:,:) = sqrt(pt_grad)

!!===== CALCULATE FGF USING DIV and DEF ==================================

  ! v_y = 0.5*(div - def_st)
  ! u_y = 0.5*(def_sh - vor)
  def_st(:,:) = div(:,:) - 2.*grady_2nd(v)
  def_sh(:,:) = gradx_2nd(v) - grady_cos_2nd(u) + 2.*grady_2nd(u)
!  def_st(:,:) = gradx_2nd(u) - v(:,:)*tanlat/a_earth - grady_2nd(v)
!  def_sh(:,:) = gradx_2nd(v) + grady_2nd(u) + u(:,:)*tanlat/a_earth

  def(:,:) = sqrt(def_st(:,:)*def_st(:,:) + def_sh(:,:)*def_sh(:,:))

  fgf_def(:,:) = -0.5*(def_st(:,:)*(pt_x_2(:,:) - pt_y_2(:,:))) - &
                 def_sh(:,:)*pt_x(:,:)*pt_y(:,:)

  fgf(:,:) = fgf_div(:,:) + fgf_def(:,:)

!!========================================================================

!!===== CALCULATE FGF WITHOUT DIV and DEF ================================
!
!  fgf(:,:) = (-1.)*( pt_x_2(:,:)*(gradx_2nd(u) -                         &
!                                  v(:,:)*tanlat/a_earth) +               &
!                     pt_y_2(:,:)*grady_2nd(v) +                          &
!                     pt_x(:,:)*pt_y(:,:)*(gradx_2nd(v) + grady_2nd(u) +  &
!                                          u(:,:)*tanlat/a_earth) )
!
!  fgf_def(:,:) = fgf(:,:) - fgf_div(:,:)
!
!!========================================================================

  pt_grad(:,:) = pt_grad(:,:)*1.e5  ! [ K / (100 km) ]
  fgf    (:,:) = fgf    (:,:)*(1.e10*3600.)  ! [ K^2 / (100 km)^2 / hr ]
  fgf_div(:,:) = fgf_div(:,:)*(1.e10*3600.)  ! [ K^2 / (100 km)^2 / hr ]
  fgf_def(:,:) = fgf_def(:,:)*(1.e10*3600.)  ! [ K^2 / (100 km)^2 / hr ]

  call missing_bdy(nx,ny,pt_grad,missv,0,0,1,1)
  call missing_bdy(nx,ny,fgf    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,div    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,def    ,missv,0,0,1,1)
  call missing_bdy(nx,ny,fgf_div,missv,0,0,1,1)
  call missing_bdy(nx,ny,fgf_def,missv,0,0,1,1)

END subroutine fgf_2d

SUBROUTINE set_gridvar(nx,ny,lat)

  integer,             intent(in) ::  nx, ny
  real, dimension(ny), intent(in) ::  lat

  if ( allocated(lat_pre) ) then
    if ( nx == nx_pre .and. ny == ny_pre ) then
      if ( all(lat == lat_pre) )  RETURN
    end if
    deallocate( lat_pre, coslat, tanlat, inv_2dx, inv_2dy )
  end if

  allocate( coslat(nx,ny), tanlat(nx,ny), inv_2dx(ny), inv_2dy(ny) )
  allocate( lat_pre(ny) )

  coslat(:,:) = spread(cos(lat(:)*deg2rad),1,nx)
  tanlat(:,:) = spread(tan(lat(:)*deg2rad),1,nx)
  if (abs(lat(1)) == 90.) then
    coslat(:,1) = 0.  ;  tanlat(:,1) = sign(1.e20,lat(1))
  end if
  if (abs(lat(ny)) == 90.) then
    coslat(:,ny) = 0.  ;  tanlat(:,ny) = sign(1.e20,lat(ny))
  end if

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

SUBROUTINE missing_bdy(nx,ny,var,missv,nm_x1,nm_x2,nm_y1,nm_y2)

  integer, intent(in) ::  nx, ny, nm_x1, nm_x2, nm_y1, nm_y2
  real,    intent(in) ::  missv

  real, dimension(nx,ny), intent(inout) ::  var

  if (nm_x1 > 0)  var(1         :nm_x1,          :     ) = missv
  if (nm_x2 > 0)  var(nx+1-nm_x2:nx   ,          :     ) = missv
  if (nm_y1 > 0)  var(          :     ,1         :nm_y1) = missv
  if (nm_y2 > 0)  var(          :     ,ny+1-nm_y2:ny   ) = missv

END subroutine missing_bdy

END module fgf


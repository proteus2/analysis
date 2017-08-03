MODULE util

  use const_glob,  only:  a_earth, deg2rad

  implicit none

  public ::  gradx_2nd, grady_2nd, gradz_2nd_irr, grady_wgt_2nd,         &
             gradxx_4nd, gradxx_2nd
  public ::  get_ind_bnd, missing_bdy
  public ::  lowpass_k, filter121_y

  private ::  a_earth, deg2rad

  interface grady_2nd
    module procedure grady_2nd_3d, grady_2nd_2d
  end interface
  interface gradz_2nd_irr
    module procedure gradz_2nd_irr_3d, gradz_2nd_irr_2d, gradz_2nd_irr_1d
  end interface
  interface grady_wgt_2nd
    module procedure grady_wgt_2nd_3d, grady_wgt_2nd_2d
  end interface
  interface gradz_wgt_2nd_irr
    module procedure gradz_wgt_2nd_irr_3d, gradz_wgt_2nd_irr_2d,         &
                     gradz_wgt_2nd_irr_1d
  end interface
  interface missing_bdy
    module procedure missing_bdy_3d, missing_bdy_2d
  end interface
  interface filter121_y
    module procedure filter121_y_3d, filter121_y_yz
  end interface
 

  CONTAINS


SUBROUTINE gradx_2nd(nx,ny,nz,var,lat,dlon, gradx)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real                     , intent(in)  ::  dlon
  real, dimension(nx,ny,nz), intent(out) ::  gradx

  real    ::  inv_2dx(ny)
  integer ::  j,k

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

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real, dimension(nx,ny,nz), intent(out) ::  grady

  real    ::  inv_2dy(ny)
  integer ::  j,k

  inv_2dy(2:ny-1) = 1./((lat(3:ny) - lat(1:ny-2))*deg2rad*a_earth)

  do k=1, nz
  do j=2, ny-1
    grady(:,j,k) = (var(:,j+1,k) - var(:,j-1,k))*inv_2dy(j)
  enddo
  enddo
  grady(:,1 ,:) = 0.
  grady(:,ny,:) = 0.

END subroutine grady_2nd_3d

SUBROUTINE gradz_2nd_irr_3d(nx,ny,nz,var,z, gradz)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(nz)      , intent(in)  ::  z
  real, dimension(nx,ny,nz), intent(out) ::  gradz

  real    ::  coef(3,nz), inv_dz_1, inv_dz_n
  integer ::  k

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coef(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coef(:,k))
  enddo
  inv_dz_1 = 1./(z(2 ) - z(1   ))
  inv_dz_n = 1./(z(nz) - z(nz-1))

  do k=2, nz-1
    gradz(:,:,k) = coef(1,k)*var(:,:,k-1) +                              &
                   coef(2,k)*var(:,:,k  ) +                              &
                   coef(3,k)*var(:,:,k+1)
  enddo
  gradz(:,:,1 ) = (var(:,:,2 ) - var(:,:,1   ))*inv_dz_1
  gradz(:,:,nz) = (var(:,:,nz) - var(:,:,nz-1))*inv_dz_n

END subroutine gradz_2nd_irr_3d

SUBROUTINE grady_wgt_2nd_3d(nx,ny,nz,var,lat,coslat, grady)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat, coslat
  real, dimension(nx,ny,nz), intent(out) ::  grady

  real    ::  inv_2dyw(ny), tmp(nx,ny)
  integer ::  j,k

  inv_2dyw(2:ny-1) = 1./((lat(3:ny) - lat(1:ny-2))*deg2rad*a_earth)/     &
                     coslat(2:ny-1)

  do k=1, nz
    do j=1, ny
      tmp(:,j) = var(:,j,k)*coslat(j)
    enddo
    do j=2, ny-1
      grady(:,j,k) = (tmp(:,j+1) - tmp(:,j-1))*inv_2dyw(j)
    enddo
  enddo
  grady(:,1 ,:) = 0.
  grady(:,ny,:) = 0.

END subroutine grady_wgt_2nd_3d

SUBROUTINE gradz_wgt_2nd_irr_3d(nx,ny,nz,var,z,rho, gradz)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(nz)      , intent(in)  ::  z, rho
  real, dimension(nx,ny,nz), intent(out) ::  gradz

  real    ::  coefw(3,nz), inv_dzw_1, inv_dzw_n, tmp(nx,ny,nz)
  integer ::  k

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coefw(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coefw(:,k))
    coefw(:,k) = coefw(:,k)/rho(k)
  enddo
  inv_dzw_1 = 1./(z(2 ) - z(1   ))/rho(1 )
  inv_dzw_n = 1./(z(nz) - z(nz-1))/rho(nz)

  do k=1, nz
    tmp(:,:,k) = var(:,:,k)*rho(k)
  enddo

  do k=2, nz-1
    gradz(:,:,k) = coefw(1,k)*tmp(:,:,k-1) +                             &
                   coefw(2,k)*tmp(:,:,k  ) +                             &
                   coefw(3,k)*tmp(:,:,k+1)
  enddo
  gradz(:,:,1 ) = (tmp(:,:,2 ) - tmp(:,:,1   ))*inv_dzw_1
  gradz(:,:,nz) = (tmp(:,:,nz) - tmp(:,:,nz-1))*inv_dzw_n

END subroutine gradz_wgt_2nd_irr_3d

SUBROUTINE grady_2nd_2d(ny,nz,var,lat, grady)

  integer               , intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(ny)   , intent(in)  ::  lat
  real, dimension(ny,nz), intent(out) ::  grady

  real    ::  inv_2dy(ny)
  integer ::  j,k

  inv_2dy(2:ny-1) = 1./((lat(3:ny) - lat(1:ny-2))*deg2rad*a_earth)

  do k=1, nz
  do j=2, ny-1
    grady(j,k) = (var(j+1,k) - var(j-1,k))*inv_2dy(j)
  enddo
  enddo
  grady(1 ,:) = 0.
  grady(ny,:) = 0.

END subroutine grady_2nd_2d

SUBROUTINE gradz_2nd_irr_2d(ny,nz,var,z, gradz)

  integer               , intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(nz)   , intent(in)  ::  z
  real, dimension(ny,nz), intent(out) ::  gradz

  real    ::  coef(3,nz), inv_dz_1, inv_dz_n
  integer ::  k

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coef(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coef(:,k))
  enddo
  inv_dz_1 = 1./(z(2 ) - z(1   ))
  inv_dz_n = 1./(z(nz) - z(nz-1))

  do k=2, nz-1
    gradz(:,k) = coef(1,k)*var(:,k-1) +                                  &
                 coef(2,k)*var(:,k  ) +                                  &
                 coef(3,k)*var(:,k+1)
  enddo
  gradz(:,1 ) = (var(:,2 ) - var(:,1   ))*inv_dz_1
  gradz(:,nz) = (var(:,nz) - var(:,nz-1))*inv_dz_n

END subroutine gradz_2nd_irr_2d

SUBROUTINE gradz_wgt_2nd_irr_2d(ny,nz,var,z,rho, gradz)

  integer               , intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(nz)   , intent(in)  ::  z, rho
  real, dimension(ny,nz), intent(out) ::  gradz

  real    ::  coefw(3,nz), inv_dzw_1, inv_dzw_n, tmp(ny,nz)
  integer ::  k

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coefw(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coefw(:,k))
    coefw(:,k) = coefw(:,k)/rho(k)
  enddo
  inv_dzw_1 = 1./(z(2 ) - z(1   ))/rho(1 )
  inv_dzw_n = 1./(z(nz) - z(nz-1))/rho(nz)

  do k=1, nz
    tmp(:,k) = var(:,k)*rho(k)
  enddo

  do k=2, nz-1
    gradz(:,k) = coefw(1,k)*tmp(:,k-1) +                                 &
                 coefw(2,k)*tmp(:,k  ) +                                 &
                 coefw(3,k)*tmp(:,k+1)
  enddo
  gradz(:,1 ) = (tmp(:,2 ) - tmp(:,1   ))*inv_dzw_1
  gradz(:,nz) = (tmp(:,nz) - tmp(:,nz-1))*inv_dzw_n

END subroutine gradz_wgt_2nd_irr_2d

SUBROUTINE grady_wgt_2nd_2d(ny,nz,var,lat,coslat, grady)

  integer               , intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(ny)   , intent(in)  ::  lat, coslat
  real, dimension(ny,nz), intent(out) ::  grady

  real    ::  inv_2dyw(ny), tmp(ny,nz)
  integer ::  j,k

  inv_2dyw(2:ny-1) = 1./((lat(3:ny) - lat(1:ny-2))*deg2rad*a_earth)/     &
                     coslat(2:ny-1)

  do k=1, nz
    tmp(:,k) = var(:,k)*coslat(:)
  enddo

  do k=1, nz
  do j=2, ny-1
    grady(j,k) = (tmp(j+1,k) - tmp(j-1,k))*inv_2dyw(j)
  enddo
  enddo
  grady(1 ,:) = 0.
  grady(ny,:) = 0.

END subroutine grady_wgt_2nd_2d

SUBROUTINE gradz_2nd_irr_1d(nz,var,z, gradz)

  integer            , intent(in)  ::  nz
  real, dimension(nz), intent(in)  ::  var
  real, dimension(nz), intent(in)  ::  z
  real, dimension(nz), intent(out) ::  gradz

  real    ::  coef(3,nz), inv_dz_1, inv_dz_n
  integer ::  k

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coef(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coef(:,k))
  enddo
  inv_dz_1 = 1./(z(2 ) - z(1   ))
  inv_dz_n = 1./(z(nz) - z(nz-1))

  do k=2, nz-1
    gradz(k) = coef(1,k)*var(k-1) +                                      &
               coef(2,k)*var(k  ) +                                      &
               coef(3,k)*var(k+1)
  enddo
  gradz(1 ) = (var(2 ) - var(1   ))*inv_dz_1
  gradz(nz) = (var(nz) - var(nz-1))*inv_dz_n

END subroutine gradz_2nd_irr_1d

SUBROUTINE gradz_wgt_2nd_irr_1d(nz,var,z,rho, gradz)

  integer            , intent(in)  ::  nz
  real, dimension(nz), intent(in)  ::  var
  real, dimension(nz), intent(in)  ::  z, rho
  real, dimension(nz), intent(out) ::  gradz

  real    ::  coefw(3,nz), inv_dzw_1, inv_dzw_n, tmp(nz)
  integer ::  k

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coefw(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coefw(:,k))
    coefw(:,k) = coefw(:,k)/rho(k)
  enddo
  inv_dzw_1 = 1./(z(2 ) - z(1   ))/rho(1 )
  inv_dzw_n = 1./(z(nz) - z(nz-1))/rho(nz)

  tmp(:) = var(:)*rho(:)

  do k=2, nz-1
    gradz(k) = coefw(1,k)*tmp(k-1) +                                     &
               coefw(2,k)*tmp(k  ) +                                     &
               coefw(3,k)*tmp(k+1)
  enddo
  gradz(1 ) = (tmp(2 ) - tmp(1   ))*inv_dzw_1
  gradz(nz) = (tmp(nz) - tmp(nz-1))*inv_dzw_n

END subroutine gradz_wgt_2nd_irr_1d

SUBROUTINE gradxx_2nd(nx,ny,nz,var,lat,dlon, gradxx)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real                     , intent(in)  ::  dlon
  real, dimension(nx,ny,nz), intent(out) ::  gradxx

  real    ::  inv_dx2(ny)
  integer ::  j,k

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

SUBROUTINE gradxx_4nd(nx,ny,nz,var,lat,dlon, gradxx)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real                     , intent(in)  ::  dlon
  real, dimension(nx,ny,nz), intent(out) ::  gradxx

  real, dimension(-1:nx+2,ny,nz) ::  var_b

  real    ::  inv_dx2(ny)
  integer ::  j,k

  real, parameter ::  c0 = -2.5, c1 = 4./3., c2 = -1./12.

  inv_dx2(:) = 1.0/((dlon*deg2rad*a_earth)*cos(lat(:)*deg2rad))**2
  if (abs(lat(1 )) == 90.)  inv_dx2(1 ) = 0.
  if (abs(lat(ny)) == 90.)  inv_dx2(ny) = 0.

  var_b(1:nx,:,:) = var(:,:,:)
  var_b(-1:0,:,:) = var(nx-1:nx,:,:)
  var_b(nx+1:nx+2,:,:) = var(1:2,:,:)

  gradxx(1:nx,:,:) = c0*var_b(1:nx,:,:) +                                &
                     c1*(var_b(0 :nx-1,:,:) + var_b(2:nx+1,:,:)) +       &
                     c2*(var_b(-1:nx-2,:,:) + var_b(3:nx+2,:,:))

  do k=1, nz
  do j=1, ny
    gradxx(:,j,k) = gradxx(:,j,k)*inv_dx2(j)
  enddo
  enddo

END subroutine gradxx_4nd

SUBROUTINE get_ind_bnd(nx,x,bnd, i1b,i2b)

  integer            , intent(in)  ::  nx
  real, dimension(nx), intent(in)  ::  x
  real               , intent(in)  ::  bnd(2)
  integer            , intent(out) ::  i1b, i2b

  real    ::  minbnd, maxbnd
  integer ::  i

  i1b = nx+1  ;  i2b = 0

  minbnd = min(bnd(1),bnd(2))  ;  maxbnd = max(bnd(1),bnd(2))
  do i=1, nx
    if ( x(i) >= minbnd .and. x(i) <= maxbnd ) then
      i1b = i  ;  EXIT
    end if
  enddo
  do i=nx, 1, -1
    if ( x(i) >= minbnd .and. x(i) <= maxbnd ) then
      i2b = i  ;  EXIT
    end if
  enddo

END subroutine get_ind_bnd

SUBROUTINE missing_bdy_3d(nx,ny,nz,var,missv,nm_y1,nm_y2,nm_z1,nm_z2)

  integer, intent(in) ::  nx, ny, nz, nm_y1, nm_y2, nm_z1, nm_z2
  real   , intent(in) ::  missv

  real, dimension(nx,ny,nz), intent(inout) ::  var

  if (nm_y1 > 0)  var(:,1         :nm_y1,          :     ) = missv
  if (nm_y2 > 0)  var(:,ny+1-nm_y2:ny   ,          :     ) = missv
  if (nm_z1 > 0)  var(:,          :     ,1         :nm_z1) = missv
  if (nm_z2 > 0)  var(:,          :     ,nz+1-nm_z2:nz   ) = missv

END subroutine missing_bdy_3d

SUBROUTINE missing_bdy_2d(ny,nz,var,missv,nm_y1,nm_y2,nm_z1,nm_z2)

  integer, intent(in) ::  ny, nz, nm_y1, nm_y2, nm_z1, nm_z2
  real   , intent(in) ::  missv

  real, dimension(ny,nz), intent(inout) ::  var

  if (nm_y1 > 0)  var(1         :nm_y1,          :     ) = missv
  if (nm_y2 > 0)  var(ny+1-nm_y2:ny   ,          :     ) = missv
  if (nm_z1 > 0)  var(          :     ,1         :nm_z1) = missv
  if (nm_z2 > 0)  var(          :     ,nz+1-nm_z2:nz   ) = missv

END subroutine missing_bdy_2d

SUBROUTINE lowpass_k(var,k_high)

  integer, intent(in) ::  k_high

  real, dimension(:,:,:), intent(inout) ::  var

  real, dimension(size(var,2),size(var,3))        ::  vm
  real, dimension(size(var,1),k_high)             ::  coskx, sinkx
  real, dimension(size(var,2),size(var,3),k_high) ::  coef_c, coef_s

  integer ::  nx, ny, nz
  integer ::  i,j,k, ik

  double precision, parameter ::  twopi = 2.d0*3.14159265358979323846d0

  nx = size(var,1)  ;  ny = size(var,2)  ;  nz = size(var,3)

  vm(:,:) = sum(var, dim=1)/float(nx)

  var(:,:,:) = var(:,:,:) - spread(vm(:,:),1,nx)

  do ik=1, k_high
  do i=1, nx
    coskx(i,ik) = real(cos(dble(ik*(i-1))/dble(nx)*twopi))
    sinkx(i,ik) = real(sin(dble(ik*(i-1))/dble(nx)*twopi))
  enddo
  enddo
 
  do ik=1, k_high
  do k=1, nz
  do j=1, ny
    coef_c(j,k,ik) = sum(var(:,j,k)*coskx(:,ik))
    coef_s(j,k,ik) = sum(var(:,j,k)*sinkx(:,ik))
  enddo
  enddo
  enddo
  coef_c(:,:,:) = coef_c(:,:,:)*(2./float(nx))
  coef_s(:,:,:) = coef_s(:,:,:)*(2./float(nx))

  var(:,:,:) = 0.
  do ik=1, k_high
  do k=1, nz
  do j=1, ny
    var(:,j,k) = var(:,j,k) +                                            &
        (coef_c(j,k,ik)*coskx(:,ik) + coef_s(j,k,ik)*sinkx(:,ik))
  enddo
  enddo
  enddo
  var(:,:,:) = var(:,:,:) + spread(vm(:,:),1,nx)

END subroutine lowpass_k

SUBROUTINE filter121_y_3d(var,vtype)

  real, dimension(:,:,:), intent(inout) ::  var

  character(len=*), intent(in), optional ::  vtype

  real, dimension(size(var,1),size(var,2),size(var,3)) ::  tmp
  real, dimension(size(var,1)) ::  cosx, sinx
  real, dimension(size(var,3)) ::  ns90, sn90

  real    ::  coef_c1, coef_s1, coef_c2, coef_s2
  integer ::  nx, ny, nz
  integer ::  i,k

  double precision, parameter ::  twopi = 2.d0*3.14159265358979323846d0

  nx = size(var,1)  ;  ny = size(var,2)  ;  nz = size(var,3)

  tmp(:,:,:) = var(:,:,:)

  var(:,2:ny-1,:) = 0.5*tmp(:,2:ny-1,:) +                                &
                    0.25*( tmp(:,1:ny-2,:) + tmp(:,3:ny,:) )

  if ( present(vtype) ) then

    if ( trim(vtype) == 'u' .or. trim(vtype) == 'U' .or.                 &
         trim(vtype) == 'v' .or. trim(vtype) == 'V' ) then
      do i=1, nx
        cosx(i) = real(cos(dble(i-1)/dble(nx)*twopi))
        sinx(i) = real(sin(dble(i-1)/dble(nx)*twopi))
      enddo
      do k=1, nz
        coef_c1 = sum(tmp(:,2   ,k)*cosx(:))*(2./float(nx))
        coef_s1 = sum(tmp(:,2   ,k)*sinx(:))*(2./float(nx))
        coef_c2 = sum(tmp(:,ny-1,k)*cosx(:))*(2./float(nx))
        coef_s2 = sum(tmp(:,ny-1,k)*sinx(:))*(2./float(nx))
        var(:,1 ,k) = 0.5*( tmp(:,1 ,k) +                                &
                            (coef_c1*cosx(:) + coef_s1*sinx(:)) )
        var(:,ny,k) = 0.5*( tmp(:,ny,k) +                                &
                            (coef_c2*cosx(:) + coef_s2*sinx(:)) )
      enddo
    else  ! scalar, same values at all x grids
      ns90(:) = 0.5*( tmp(1,1 ,:) + sum(tmp(:,2   ,:),dim=1)/float(nx) )
      sn90(:) = 0.5*( tmp(1,ny,:) + sum(tmp(:,ny-1,:),dim=1)/float(nx) )
      var(:,1 ,:) = spread(ns90(:),1,nx)
      var(:,ny,:) = spread(sn90(:),1,nx)
    end if

  else

    var(:,1 ,:) = 0.5*( tmp(:,1 ,:) + tmp(:,2   ,:) )
    var(:,ny,:) = 0.5*( tmp(:,ny,:) + tmp(:,ny-1,:) )

  end if

END subroutine filter121_y_3d

SUBROUTINE filter121_y_yz(var,vtype)

  real, dimension(:,:), intent(inout) ::  var

  character(len=*), intent(in), optional ::  vtype

  real, dimension(size(var,1),size(var,2)) ::  tmp

  integer ::  ny, nz
  integer ::  k

  ny = size(var,1)  ;  nz = size(var,2)

  tmp(:,:) = var(:,:)

  var(2:ny-1,:) = 0.5*tmp(2:ny-1,:) +                                    &
                  0.25*( tmp(1:ny-2,:) + tmp(3:ny,:) )

  if ( present(vtype) ) then

    if ( trim(vtype) == 'u' .or. trim(vtype) == 'U' .or.                 &
         trim(vtype) == 'v' .or. trim(vtype) == 'V' ) then
      var(1 ,k) = 0.
      var(ny,k) = 0.
    else
      var(1 ,:) = 0.5*( tmp(1 ,:) + tmp(2   ,:) )
      var(ny,:) = 0.5*( tmp(ny,:) + tmp(ny-1,:) )
    end if

  else

    var(1 ,:) = 0.5*( tmp(1 ,:) + tmp(2   ,:) )
    var(ny,:) = 0.5*( tmp(ny,:) + tmp(ny-1,:) )

  end if

END subroutine filter121_y_yz

END module util


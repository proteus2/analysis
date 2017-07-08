MODULE util

  use const_glob,  only:  a_earth, deg2rad

  implicit none

  public ::  gradx_2nd, grady_2nd, gradz_2nd_irr, gradxx_2nd
  public ::  missing_bdy

  private ::  a_earth, deg2rad

  interface grady_2nd
    module procedure grady_2nd_3d, grady_2nd_2d
  end interface
  interface gradz_2nd_irr
    module procedure gradz_2nd_irr_3d, gradz_2nd_irr_2d, gradz_2nd_irr_1d
  end interface
  interface missing_bdy
    module procedure missing_bdy_3d, missing_bdy_2d
  end interface


  CONTAINS


SUBROUTINE gradx_2nd(nx,ny,nz,var,lat,dlon, gradx)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real                     , intent(in)  ::  dlon
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

  integer                  , intent(in)  ::  nx, ny, nz
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

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(nz)      , intent(in)  ::  z
  real, dimension(nx,ny,nz), intent(out) ::  gradz

  integer ::  k
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

  integer               , intent(in)  ::  ny, nz
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

  integer               , intent(in)  ::  ny, nz
  real, dimension(ny,nz), intent(in)  ::  var
  real, dimension(nz)   , intent(in)  ::  z
  real, dimension(ny,nz), intent(out) ::  gradz

  integer ::  k
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

SUBROUTINE gradz_2nd_irr_1d(nz,var,z, gradz)

  integer            , intent(in)  ::  nz
  real, dimension(nz), intent(in)  ::  var
  real, dimension(nz), intent(in)  ::  z
  real, dimension(nz), intent(out) ::  gradz

  integer ::  k
  real    ::  coef(3,nz), inv_dz_1, inv_dz_n

  do k=2, nz-1
!    call fdcoef(1,2,z(k),z(k-1:k+1), coef(:,k))
    call fdcoef_1d2o(z(k-1:k+1), coef(:,k))
  enddo
  inv_dz_1 = 1./(z(2 ) - z(1   ))
  inv_dz_n = 1./(z(nz) - z(nz-1))

  do k=2, nz-1
    gradz(k) = coef(1,k)*var(k-1) +  &
               coef(2,k)*var(k  ) +  &
               coef(3,k)*var(k+1)
  enddo
  gradz(1 ) = (var(2 ) - var(1   ))*inv_dz_1
  gradz(nz) = (var(nz) - var(nz-1))*inv_dz_n

END subroutine gradz_2nd_irr_1d

SUBROUTINE gradxx_2nd(nx,ny,nz,var,lat,dlon, gradxx)

  integer                  , intent(in)  ::  nx, ny, nz
  real, dimension(nx,ny,nz), intent(in)  ::  var
  real, dimension(ny)      , intent(in)  ::  lat
  real                     , intent(in)  ::  dlon
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

END module util


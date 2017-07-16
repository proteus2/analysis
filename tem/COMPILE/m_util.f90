MODULE util

  use const_glob,  only:  a_earth, deg2rad

  implicit none

  public ::  gradx_2nd, grady_2nd, gradz_2nd_irr, gradxx_4nd, gradxx_2nd
  public ::  get_ind_bnd, missing_bdy
  public ::  lowpass_k

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

  real    ::  coef(3,nz), inv_dz_1, inv_dz_n
  integer ::  k

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

  real    ::  inv_2dy(ny)
  integer ::  j,k

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

  real    ::  coef(3,nz), inv_dz_1, inv_dz_n
  integer ::  k

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

  real    ::  coef(3,nz), inv_dz_1, inv_dz_n
  integer ::  k

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

END module util


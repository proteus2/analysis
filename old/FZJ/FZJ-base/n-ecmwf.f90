PROGRAM n_ecmwf

  use netcdfio

  implicit none

  integer, parameter ::  i_out3d = 1

  integer, parameter ::  ix = 480, nx = 84  ! 119.75 - 140.5
  integer, parameter ::  iy = 243, ny = 77  ! 29.5 - 10.5
  integer, parameter ::  iz = 31 , nz = 101 ! 15 - 65

  real, dimension(nx,ny,nz) ::  temp, nb
  real, dimension(nx,ny)    ::  area
  real, dimension(nz)       ::  t0, n0

  real ::  lon(nx), lat(ny), z(nz)


  integer ::  i,j,k, ncid
  real    ::  cst1, cst2

  real, dimension(:),     allocatable ::  dat1d
  real, dimension(:,:,:), allocatable ::  dat3d

  real, parameter ::  g = 9.81, cp = 1004.

  allocate( dat1d(ny) )
  allocate( dat3d(nz,ny,nx) )

  call opennc('ECMWF/ecmwfr_ana_ml_T_06070700.nc',ncid)
  call geta1d_yh(ncid,'lon'   ,ix,nx, lon  )
  call geta1d_yh(ncid,'lat'   ,iy,ny, dat1d)
  call geta1d_yh(ncid,'height',iz,nz, z    )
  call geta4d_yh(ncid,'TEMP',iz,nz,iy,ny,ix,nx,1,1, dat3d)
  call closenc(ncid)

  do j=1, ny
    lat(j) = dat1d(ny+1-j)
  enddo

  do k=1, nz
  do j=1, ny
    temp(:,j,k) = dat3d(k,ny+1-j,:)
  enddo
  enddo

  deallocate( dat3d )  ;  deallocate( dat1d )

  cst1 = z(3) - z(1)
  cst2 = g/cp
  do k=2, nz-1
    nb(:,:,k) = g/temp(:,:,k)*((temp(:,:,k+1)-temp(:,:,k-1))/cst1+cst2)
  enddo
  nb(:,:,1 ) = nb(:,:,2   )
  nb(:,:,nz) = nb(:,:,nz-1)

  do j=1, ny
    area(:,j) = cos(lat(j)*3.1415926/180.)
  enddo

  do k=1, nz
    t0(k) = sum(temp(:,:,k)*area(:,:)) / sum(area)
    n0(k) = sum(nb  (:,:,k)*area(:,:)) / sum(area)
  enddo

  nb = sqrt(nb)
  n0 = sqrt(n0)

  call out1d_yh('ECMWF/ec-npro.nc',2,(/'T0','N0'/),(/t0,n0/), &
                'z',nz,z,'averaged profiles')

  if (i_out3d == 1)  &
     call out3d_yh('ECMWF/ec-temp.nc',1,(/'T'/),temp, &
                   'lon',nx,lon,'lat',ny,lat,'z',nz,z,'temperture')


  STOP

END program n_ecmwf

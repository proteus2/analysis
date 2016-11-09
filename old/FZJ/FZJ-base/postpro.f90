PROGRAM post_prog

  use netcdfio

  implicit none

  integer, parameter ::  date = 7, time = 0

  integer ::  nx, ny, nz
  integer ::  i, ncid
  character(len=256) ::  fhead, fname

  real, parameter ::  g = 9.81, kappa = 0.286

  real, dimension(:),     allocatable ::  x, y, z
  real, dimension(:,:,:), allocatable ::  t3d, p3d, pb3d

  fhead = '/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/wrfout_d01_2006-07-'
  fhead = '/data1/ksy/WRFV3/WRFV3/run/wrf_real_input_em.d01.2006-07-'
  write(fname,'(a,i2.2,a,i2.2,a)') trim(fhead),date,'_',time,':00:00'

  call opennc(fname,ncid)
  call dilen(ncid,'west_east'  ,nx)
  call dilen(ncid,'south_north',ny)
  call dilen(ncid,'bottom_top' ,nz)
  allocate( t3d(nx,ny,nz), p3d(nx,ny,nz), pb3d(nx,ny,nz) )
  call geta4d_yh(ncid,'T' ,1,nx,1,ny,1,nz,1,1, t3d )
  call geta4d_yh(ncid,'P' ,1,nx,1,ny,1,nz,1,1, p3d )
  call geta4d_yh(ncid,'PB',1,nx,1,ny,1,nz,1,1, pb3d)
  p3d = p3d + pb3d
  deallocate( pb3d )
  t3d = (t3d+300.) * (p3d/1.e5)**kappa
  call closenc(ncid)

  allocate( x(nx), y(ny), z(nz) )
  do i=1, nx
    x(i) = float(i-1)
  enddo
  do i=1, ny
    y(i) = float(i-1)
  enddo
  do i=1, nz
    z(i) = float(i-1)
  enddo

  write(fname,'(a,i2.2,a,i2.2,a)') 'wrfinpp_d01_2006-07-',date,'_',time,':00:00.nc'
!  call out3d_yh(fname,2,(/'T','P'/),(/t3d,p3d/), &
!                'x',nx,x,'y',ny,y,'z',nz,z,'pp')
call out1d('zz.nc',2,(/'T','P'/),(/t3d(85,136,:),p3d(85,136,:)/),'z',nz,z,'pp')


  STOP

END program post_prog

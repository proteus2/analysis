    program MF 

    use netcdfio
    use fft

    implicit none

    integer :: i, j, k, t, file, tt, kk, t1
    integer :: ncid
    integer :: day, hour
    integer :: dim1, dim2, dim3, dim4

    character(len=2), parameter ::  speriod='p1'
    integer, parameter :: start_day=7, start_hour=1
    integer, parameter :: nfile=29
    integer, parameter :: nxt=186, nyt=186, nz=51, nt1=6, iz1 = 35
    integer, parameter :: nt=nfile*nt1+1
    integer, parameter :: ngrid=21

    integer, parameter :: sx=32
    integer, parameter :: ex=116
    integer, parameter :: sy=32
    integer, parameter :: ey=116 

    integer, parameter :: sx1=sx-10 ! for background
    integer, parameter :: ex1=ex-10
    integer, parameter :: sy1=sy-10
    integer, parameter :: ey1=ey-10

    integer, parameter :: nx1=166, ny1=166 ! background

    integer, parameter :: nx=ex-sx+1
    integer, parameter :: ny=ey-sy+1
 
    real, parameter :: dx=27000., dy=27000., dz=500., dt=600.
    real, parameter :: g=9.81

    character(len=100) :: ncfilename
    character(len=50) :: varname
    character(len=50) :: d1name,d2name, d3name, d4name
    character(len=50) :: title

    real, dimension(nxt) :: x_org
    real, dimension(nyt) :: y_org
    real, dimension(nx) :: x
    real, dimension(ny) :: y
    real, dimension(0:nz) :: z1
    real, dimension(nz) :: z, tavg, tsd
    real, dimension(nt) :: h
    real, dimension(nz,nt) :: tavg2d, tano

    real, dimension(nx,ny,nz,nt) :: tb
 
    real, allocatable, dimension(:,:,:,:) :: th_org, th4d
 
!== READ WRF OUTPUT =====================================================
 
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',start_day,'_',start_hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/w/'//ncfilename, ncid)
      call get1d(ncid, 'west_east', nxt, x_org)
      do i=sx, ex, 1
        x(i-sx+1)=(x_org(i)-1)*27.
      end do
      call get1d(ncid, 'south_north', nyt, y_org)
      do j=sy, ey, 1
        y(j-sy+1)=(y_org(j)-1)*27.
      end do
      call geta1d_yh(ncid, 'height', iz1-1, nz+1, z1)
      z(:) = z1(1:)
    call closenc(ncid)

!!!

    allocate(th4d(nx,ny,1:nz,nt))

    day=start_day
    hour=start_hour

    t=0
    do tt=1, nfile  ! Loop: file 

    allocate(th_org(nxt,nyt,1:nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    varname='temperature'
    call geta4d_yh(ncid, varname, 1,nxt, 1,nyt, iz1,nz, 1,nt1, th_org)
    call closenc(ncid)

    do t1=1, nt1
      t=t+1
      do k=1, nz
      do j=1, ny
      do i=1, nx
        th4d(i,j,k,t)=sum(th_org(i+sx-1-ngrid/2:i+sx-1+ngrid/2, &
                                 j+sy-1-ngrid/2:j+sy-1+ngrid/2, &
                                 k,t1))
      end do 
      end do 
      end do 
      th4d(:,:,:,t) = th4d(:,:,:,t)/float(ngrid*ngrid)
    end do

    deallocate(th_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file
    allocate(th_org(nxt,nyt,1:nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    varname='temperature'
    call geta4d_yh(ncid, varname, 1,nxt, 1,nyt, iz1,nz, 1,1, th_org)
    call closenc(ncid)

    t=t+1
    do k=1, nz
    do j=1, ny
    do i=1, nx
      th4d(i,j,k,t)=sum(th_org(i+sx-1-ngrid/2:i+sx-1+ngrid/2, &
                               j+sy-1-ngrid/2:j+sy-1+ngrid/2, &
                               k,1))
    end do
    end do
    end do
    th4d(:,:,:,t) = th4d(:,:,:,t)/float(ngrid*ngrid)

    deallocate(th_org)

    tb(:,:,:,:) = th4d(:,:,:,:)

    deallocate(th4d)

!!!

  do t=1, nt
  do k=1, nz
    tavg2d(k,t) = sum(tb(:,:,k,t))/float(nx*ny)
  end do
  end do

  do k=1, nz
    tavg(k) = sum(tavg2d(k,:))/float(nt)
    tano(k,:) = tavg2d(k,:) - tavg(k)
  enddo

  do k=1, nz
    tsd(k) = sqrt( sum((tb(:,:,k,:)-tavg(k))**2)/float(nx*ny*nt) )
  enddo


!== 3D FFT with respect to x, y, and t ================================

    do t=1, nt
      h(t) = 1. + float(t-1)*dt/3600.
    enddo

    dim1=nx
    dim2=ny
    dim3=nz
    dim4=nt
    ncfilename='res/t_background_'//speriod//'_region.nc'
    d1name='x_km'
    d2name='y_km'
    d3name='height_m'
    d4name='time_min'
    title='background'

    call out4d_yh(ncfilename, 1, (/'tb'/), (/tb/), &
               d1name, dim1, x, &
               d2name, dim2, y, &
               d3name, dim3, z, &
               d4name, dim4, h, &
               title)

    ncfilename='res/t_background2_'//speriod//'_region.nc'
    call out1d_yh(ncfilename, 2, (/'tbavg','tbsd'/), (/tavg,tsd/), &
               d3name, dim3, z, &
               title)

    ncfilename='res/t_background3_'//speriod//'_region.nc'
    call out2d_yh(ncfilename, 2, (/'tbavg','tbano'/), (/tavg2d,tano/), &
               d3name, dim3, z, &
               d4name, dim4, h, &
               title)

    end program

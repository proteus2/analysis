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
    real, dimension(nz) :: z, uavg, vavg, navg, usd, vsd, nsd
    real, dimension(nt) :: h

    real, dimension(nx,ny,nz,nt) :: ub, vb, nb
 
    real, allocatable, dimension(:,:,:,:) :: ub_org, vb_org, th_org, th4d
 
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

    allocate(ub_org(nx1,ny1,nz,nt))
    call opennc('/data1/ksy/ewiniar-MF/0700-0806/running-mean-back/u_background.nc', ncid)
    call geta4d_yh(ncid, 'ub', 1,nx1, 1,ny1, iz1,nz, 1,nt, ub_org)
    call closenc(ncid)

    allocate(vb_org(nx1,ny1,nz,nt))
    call opennc('/data1/ksy/ewiniar-MF/0700-0806/running-mean-back/v_background.nc', ncid)
    call geta4d_yh(ncid, 'vb', 1,nx1, 1,ny1, iz1,nz, 1,nt, vb_org)
    call closenc(ncid)

!!!

    allocate(th4d(nx,ny,0:nz,nt))

    day=start_day
    hour=start_hour

    t=0
    do tt=1, nfile  ! Loop: file 

    allocate(th_org(nxt,nyt,0:nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_thp_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/perturb_th/'//ncfilename, ncid)
    varname='perturbation_theta'
    call geta4d_yh(ncid, varname, 1,nxt, 1,nyt, iz1-1,nz+1, 1,nt1, th_org)
    call closenc(ncid)

    do t1=1, nt1
      t=t+1
      do k=0, nz
      do j=1, ny
      do i=1, nx
        th4d(i,j,k,t)=sum(th_org(i+sx-1-ngrid/2:i+sx-1+ngrid/2, &
                                 j+sy-1-ngrid/2:j+sy-1+ngrid/2, &
                                 k,t1))
      end do 
      end do 
      end do 
      th4d(:,:,:,t) = th4d(:,:,:,t)/float(ngrid*ngrid) + 300.
    end do

    deallocate(th_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file
    allocate(th_org(nxt,nyt,0:nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_thp_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/perturb_th/'//ncfilename, ncid)
    varname='perturbation_theta'
    call geta4d_yh(ncid, varname, 1,nxt, 1,nyt, iz1-1,nz+1, 1,1, th_org)
    call closenc(ncid)

    t=t+1
    do k=0, nz
    do j=1, ny
    do i=1, nx
      th4d(i,j,k,t)=sum(th_org(i+sx-1-ngrid/2:i+sx-1+ngrid/2, &
                               j+sy-1-ngrid/2:j+sy-1+ngrid/2, &
                               k,1))
    end do
    end do
    end do
    th4d(:,:,:,t) = th4d(:,:,:,t)/float(ngrid*ngrid) + 300.

    deallocate(th_org)

    do k=1, nz-1
      nb(:,:,k,:) = g*log(th4d(:,:,k+1,:)/th4d(:,:,k-1,:))/(z1(k+1)-z1(k-1))
    enddo
    nb(:,:,nz,:) = g*log(th4d(:,:,nz,:)/th4d(:,:,nz-1,:))/(z1(nz)-z1(nz-1))
    nb = sqrt(nb)

    deallocate(th4d)

!!!

    do t=1, nt
    do k=1, nz
    do j=sy1, ey1, 1
    do i=sx1, ex1, 1
      ub(i-sx1+1,j-sy1+1,k,t)=ub_org(i,j,k,t) 
      vb(i-sx1+1,j-sy1+1,k,t)=vb_org(i,j,k,t)
    end do
    end do
    end do
    end do

  do k=1, nz
    uavg(k) = sum(ub(:,:,k,:))/float(nx*ny*nt)
    vavg(k) = sum(vb(:,:,k,:))/float(nx*ny*nt)
    navg(k) = sum(nb(:,:,k,:))/float(nx*ny*nt)
  enddo

  do k=1, nz
    usd(k) = sqrt( sum((ub(:,:,k,:)-uavg(k))**2)/float(nx*ny*nt) )
    vsd(k) = sqrt( sum((vb(:,:,k,:)-vavg(k))**2)/float(nx*ny*nt) )
    nsd(k) = sqrt( sum((nb(:,:,k,:)-navg(k))**2)/float(nx*ny*nt) )
  enddo


!== 3D FFT with respect to x, y, and t ================================

    do t=1, nt
      h(t) = 1. + float(t-1)*dt/3600.
    enddo

    dim1=nx
    dim2=ny
    dim3=nz
    dim4=nt
    ncfilename='res/uvn_background_'//speriod//'_region.nc'
    d1name='x_km'
    d2name='y_km'
    d3name='height_m'
    d4name='time_min'
    title='background'

    call out4d_yh(ncfilename, 3, (/'ub','vb','nb'/), (/ub,vb,nb/), &
               d1name, dim1, x, &
               d2name, dim2, y, &
               d3name, dim3, z, &
               d4name, dim4, h, &
               title)

    ncfilename='res/uvn_background2_'//speriod//'_region.nc'
    call out1d_yh(ncfilename, 6, (/'ubavg','vbavg','nbavg','ubsd','vbsd', &
                       'nbsd'/), (/uavg,vavg,navg,usd,vsd,nsd/), &
               d3name, dim3, z, &
               title)

    end program  

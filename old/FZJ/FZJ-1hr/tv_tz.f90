    program MF 

    use netcdfio
    use fft

    implicit none

    integer :: i, j, k, t, file, tt, kk, l1, l2, iz
    integer :: ncid
    integer :: day, hour
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: start_day=7, start_hour=0
    integer, parameter :: nfile=30
    integer, parameter :: nxt=186, nyt=186, nz=85, nt1=1
    integer, parameter :: nt=nfile*nt1+1
    integer, parameter :: ngrid=21
    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants
    real, parameter :: g=9.81
    character(len=2), parameter ::  speriod='p1'

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
 
    real, parameter :: dx=27000., dy=27000., dz=500., dt=3600.

    character(len=100) :: ncfilename, ifdir
    character(len=50) :: varname
    character(len=50) :: d1name,d2name, d3name, d4name
    character(len=50) :: title

    real, dimension(nxt) :: x_org
    real, dimension(nyt) :: y_org
    real, dimension(nx) :: x
    real, dimension(ny) :: y
    real, dimension(nz) :: z
    real, dimension(nt) :: h

    real, dimension(nx,ny,nz,nt) :: tb4d
    real, dimension(nt,nz) :: tvpro
 
    real, allocatable, dimension(:,:) :: t2d
    real, allocatable, dimension(:,:,:,:) :: t_org, t4d, tb_org

!== READ WRF OUTPUT =====================================================
 
    day=start_day
    hour=start_hour
    t=0
   
    allocate(t4d(nx,ny,nz,nt))
 
    do file=1, nfile  ! Loop: file 

    allocate(t_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    call geta4d_yh(ncid, 'temperature', 1,nxt, 1,nyt, 1,nz, 1,nt1, t_org)
    if(t.eq.0) then
      call get1d(ncid, 'west_east', nxt, x_org)
      do i=sx, ex, 1
        x(i-sx+1)=(x_org(i)-1-4)*27.
      end do
      call get1d(ncid, 'south_north', nyt, y_org)
      do j=sy, ey, 1
        y(j-sy+1)=(y_org(j)-1-4)*27.
      end do
      call get1d(ncid, 'height', nz, z)
    end if
    call closenc(ncid)

    do tt=1, nt1
      t=t+1
      do k=1, nz
      do j=sy, ey
      do i=sx, ex
        t4d (i-sx+1,j-sy+1,k,t)=t_org(i,j,k,tt)
        tb4d(i-sx+1,j-sy+1,k,t)=sum(t_org(i-ngrid/2:i+ngrid/2, &
                                          j-ngrid/2:j+ngrid/2, &
                                          k,tt))
      end do
      end do
      end do
      tb4d(:,:,:,t) = tb4d(:,:,:,t)/float(ngrid*ngrid)
    end do 

    deallocate(t_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file

    allocate(t_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    call get4d(ncid, 'temperature', nxt, nyt, nz, 1, t_org)
    call closenc(ncid)

    t=t+1
    do k=1, nz
    do j=sy, ey, 1
    do i=sx, ex, 1
      t4d (i-sx+1,j-sy+1,k,t)=t_org(i,j,k,1)
      tb4d(i-sx+1,j-sy+1,k,t)=sum(t_org(i-ngrid/2:i+ngrid/2, &
                                        j-ngrid/2:j+ngrid/2, &
                                        k,1))
    end do
    end do
    end do
    tb4d(:,:,:,t) = tb4d(:,:,:,t)/float(ngrid*ngrid)

    deallocate(t_org)

!== 3D FFT with respect to x, y, and t ================================

    allocate(t2d(nx,ny))

    do k=1, nz
    do t=1, nt

      t2d(:,:)=t4d(:,:,k,t)-tb4d(:,:,k,t)

      tvpro(t,k) = sum(t2d(:,:)**2)/float(nx*ny)

    end do
    end do

!== Write Output ========================================================

    do t=1, nt
      h(t) = 0. + float(t-1)*dt/3600.
    enddo

    ncfilename='res/tv_tz_'//speriod//'.nc'
    call out2d_yh(ncfilename, 1, (/'tv'/), tvpro, &
                  't', nt, h, &
                  'z', nz, z, &
                  'T variance')


    end program  

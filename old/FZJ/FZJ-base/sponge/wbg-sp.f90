    program background

    use netcdfio

    implicit none

    integer :: i, j, k, t, file, tt, kk 
    integer :: ncid
    integer :: day, hour
    integer :: sx, ex, sy, ey
    integer :: nx, ny
    integer :: l1, l2

    character(len=2), parameter ::  speriod='p1'
    integer, parameter :: start_day=7, start_hour=1
    integer, parameter :: nfile=29
    integer, parameter :: nxt=186, nyt=186, nz=101, nt1=6
    integer, parameter :: nt=nfile*nt1+1
    integer, parameter :: ngrid=21

    real, parameter :: dt=600., dz=500., dx=27000., dy=27000.

    character(len=100) :: ncfilename

    real, dimension(nz) :: z
    real, dimension(nt) :: time
    
    real, allocatable, dimension(:) :: x, y
    real, allocatable, dimension(:,:,:,:) :: u_org, ub4d

    sx=1+ngrid/2
    ex=nxt-ngrid/2
    sy=1+ngrid/2
    ey=nyt-ngrid/2

    nx=ex-sx+1
    ny=ey-sy+1

    print *, nx, ny
    print *, sx, ex, sy, ey
 
!== READ WRF OUTPUT =====================================================

    day=start_day
    hour=start_hour
    t=0
   
    allocate(ub4d(nx,ny,nz,nt)) 
 
    ub4d(:,:,:,:)=0.

    do file=1, nfile  ! Loop: file 

    allocate(u_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/data14/kyh/ewiniar/0700-0806/new_data/org/w/'//ncfilename, ncid)
    call get4d(ncid, 'vertical_wind', nxt, nyt, nz, nt1, u_org)
    call get1d(ncid, 'height', nz, z)
    call closenc(ncid)

    do tt=1, nt1
      t=t+1
      do l2=-ngrid/2, ngrid/2
      do l1=-ngrid/2, ngrid/2
        ub4d(:,:,:,t) = ub4d(:,:,:,t) + u_org(sx+l1:ex+l1,sy+l2:ey+l2,:,tt)
      end do
      end do
      ub4d(:,:,:,t) = ub4d(:,:,:,t)/float(ngrid*ngrid)
    end do 

    deallocate(u_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file

    allocate(u_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/data14/kyh/ewiniar/0700-0806/new_data/org/w/'//ncfilename, ncid)
    call get4d(ncid, 'vertical_wind', nxt, nyt, nz, 1, u_org)
    call closenc(ncid)

    t=t+1

    do l2=-ngrid/2, ngrid/2
    do l1=-ngrid/2, ngrid/2
      ub4d(:,:,:,t) = ub4d(:,:,:,t) + u_org(sx+l1:ex+l1,sy+l2:ey+l2,:,1)
    end do
    end do
    ub4d(:,:,:,t)=ub4d(:,:,:,t)/float(ngrid*ngrid)

    deallocate(u_org)

!== Write Output ========================================================


    allocate(x(nx))
    allocate(y(ny))
    do i=1, nx
      x(i)=(sx-1+i-1)*dx/1000.
    end do
    do j=1, ny
      y(j)=(sy-1+j-1)*dy/1000.
    end do
    do t=1, nt
      time(t)=float(t-1)*int(dt/60.)+60.
    end do  

    ncfilename='/data14/kyh/ewiniar/w_background_'//speriod//'.nc'
    call out4d(ncfilename, 'wb', ub4d, &
               'x_km', nx, x, &
               'y_km', ny, y, &
               'height_m', nz, z, &
               'time_min', nt, time, &
               'wb')

    deallocate(ub4d)

    end program  

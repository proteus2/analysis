    program background

    use netcdfio
    use parafit2d

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
    integer, parameter :: nxt=186, nyt=186, nz=85, nt1=6
    integer, parameter :: nt=nfile*nt1+1
    integer, parameter :: ngrid=21

    real, parameter :: dt=600., dz=500., dx=27000., dy=27000.

    character(len=100) :: ncfilename

    real                :: a(5), bfit(5:nxt-4,5:nyt-4)
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

    allocate( x(nxt-8), y(nyt-8) )
    do i=1, nxt-8
      x(i)=float(i)/float(nxt-8)*2.
    end do
    do j=1, nyt-8
      y(j)=float(j)/float(nyt-8)*2.
    end do
    x(:) = x(:) - sum(x)/float(nxt-8)
    y(:) = y(:) - sum(y)/float(nyt-8)
 
!== READ WRF OUTPUT =====================================================

    day=start_day
    hour=start_hour
    t=0
   
    allocate(ub4d(nx,ny,nz,nt)) 
 
    ub4d(:,:,:,:)=0.

    do file=1, nfile  ! Loop: file 

    allocate(u_org(5:nxt-4,5:nyt-4,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_v_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/v_wind/'//ncfilename, ncid)
    call geta4d_yh(ncid, 'meridional_wind', 5,nxt-8, 5,nyt-8, 1,nz, 1,nt1, u_org)
    call get1d(ncid, 'height', nz, z)
    call closenc(ncid)

    do tt=1, nt1
      t=t+1
      do k=1, nz
        call para_fit_2d(nxt-8,nyt-8,u_org(:,:,k,tt),a,bfit)
        ub4d(:,:,k,t) = bfit(sx:ex,sy:ey)
      enddo
    end do

    deallocate(u_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file

    allocate(u_org(5:nxt-4,5:nyt-4,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_v_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/v_wind/'//ncfilename, ncid)
    call geta4d_yh(ncid, 'meridional_wind', 5,nxt-8, 5,nyt-8, 1,nz, 1,1, u_org)
    call closenc(ncid)

    t=t+1
    do k=1, nz
      call para_fit_2d(nxt-8,nyt-8,u_org(:,:,k,1),a,bfit)
      ub4d(:,:,k,t) = bfit(sx:ex,sy:ey)
    enddo

    deallocate(u_org)

!== Write Output ========================================================

    deallocate( x, y )

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

    ncfilename='res/v_background_'//speriod//'.nc'
    call out4d(ncfilename, 'vb', ub4d, &
               'x_km', nx, x, &
               'y_km', ny, y, &
               'height_m', nz, z, &
               'time_min', nt, time, &
               'vb')

    deallocate(ub4d)

    end program  

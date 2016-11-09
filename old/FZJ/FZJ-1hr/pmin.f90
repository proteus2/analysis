    program MF 

    use netcdfio
    use fft

    implicit none

    integer :: i, j, k, t, file, tt, kk, l1, l2, iz
    integer :: ncid
    integer :: day, hour
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: start_day=9, start_hour=12
    integer, parameter :: nfile=30
    integer, parameter :: nxt=186, nyt=186, nz=85, nt1=1
    integer, parameter :: nt=nfile*nt1+1
    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants
    character(len=2), parameter ::  speriod='p3'

    integer, parameter :: sx=43
    integer, parameter :: ex=129
    integer, parameter :: sy=72
    integer, parameter :: ey=158 

    integer, parameter :: sx1=sx-10 ! for background
    integer, parameter :: ex1=ex-10
    integer, parameter :: sy1=sy-10
    integer, parameter :: ey1=ey-10

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
    real, dimension(nt)    ::  h, pmin
    integer, dimension(nt,2)  ::  minl

    real, allocatable, dimension(:,:,:) :: p3d, p_org

!== READ WRF OUTPUT =====================================================
 
    day=start_day
    hour=start_hour
    t=0
   
    allocate(p3d(nx,ny,nt))
 
    do file=1, nfile  ! Loop: file 

    allocate(p_org(nxt,nyt,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'pslv_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0912-1018/new_data/new/pslv/'//ncfilename, ncid)
    call geta3d_yh(ncid, 'pslv', 1,nxt, 1,nyt, 1,nt1, p_org)
    if(t.eq.0) then
      call get1d(ncid, 'west_east', nxt, x_org)
      do i=sx, ex, 1
        x(i-sx+1)=(x_org(i)-1-4)*27.
      end do
      call get1d(ncid, 'south_north', nyt, y_org)
      do j=sy, ey, 1
        y(j-sy+1)=(y_org(j)-1-4)*27.
      end do
    end if
    call closenc(ncid)

    do tt=1, nt1
      t=t+1
      do j=sy, ey, 1
      do i=sx, ex, 1
        p3d(i-sx+1,j-sy+1,t)=p_org(i,j,tt)
      end do
      end do
    end do

    deallocate(p_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file

    allocate(p_org(nxt,nyt,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'pslv_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0912-1018/new_data/new/pslv/'//ncfilename, ncid)
    call get3d(ncid, 'pslv', nxt, nyt, 1, p_org)
    call closenc(ncid)

    t=t+1
    do j=sy, ey, 1
    do i=sx, ex, 1
      p3d(i-sx+1,j-sy+1,t)=p_org(i,j,1)
    end do
    end do

    deallocate(p_org)

    minl(1,:) = minloc(p3d(:,:,1))
    do t=2, nt
      minl(t,:) = minloc(p3d(minl(t-1,1)-10:minl(t-1,1)+10, &
                             minl(t-1,2)-10:minl(t-1,2)+10,t))
      minl(t,:) = minl(t,:) + minl(t-1,:)-10 - 1
    enddo

    do t=1, nt
      pmin(t) = p3d(minl(t,1),minl(t,2),t)
    enddo

!== Write Output ========================================================

    do t=1, nt
      h(t) = 0. + float(t-1)*dt/3600.
    enddo

    ncfilename='res/slp-min_'//speriod//'.nc'
    call out1d_yh(ncfilename, 3, (/'pmin','x','y'/), &
               (/pmin,x(minl(:,1)),y(minl(:,2))/), &
               'time',nt,h, &
               'sea-level pressure')


    end program  

    program MF 

    use netcdfio
    use fft

    implicit none

    integer, parameter ::  cf = 1

    integer :: i, j, k, t, file, tt, kk 
    integer :: ncid
    integer :: day, hour
    integer :: dim1, dim2, dim3, dim4

    character(len=2), parameter ::  speriod='p1'
    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants
    real, parameter :: lat = 20.
    integer, parameter :: start_day=7, start_hour=1
    integer, parameter :: nfile=29
    integer, parameter :: nxt=186, nyt=186, nz=101, nt1=6
    integer, parameter :: nt=nfile*nt1+1

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

    character(len=100) :: ncfilename, ifdir
    character(len=50) :: varname
    character(len=50) :: d1name,d2name, d3name, d4name
    character(len=50) :: title

    real, dimension(nxt) :: x_org
    real, dimension(nyt) :: y_org
    real, dimension(nx) :: x
    real, dimension(ny) :: y
    real, dimension(nz) :: z, u0, v0, nb
    real, dimension(nz) :: rho_tmean, w2

    real, dimension(nx,ny,nz,nt) :: ub, vb, wb
 
    real, dimension(nx,ny,nt) :: w3d, ww3d

    real, dimension(nx) :: freq_x
    real, dimension(ny) :: freq_y
    real, dimension(nt) :: freq_t
    double complex, dimension(nx,ny,nt) ::  coefw

    real, allocatable, dimension(:,:,:,:) :: w4d, w_org 
    real, allocatable, dimension(:,:,:,:) :: wb_org

    real    ::  f0, cc, omeh

    f0 = 2.*omega*sin(lat*pi/180.)

!== READ WRF OUTPUT =====================================================
 
    day=start_day
    hour=start_hour
    t=0
   
    allocate(w4d(nx,ny,nz,nt))
 
    do file=1, nfile  ! Loop: file 

    allocate(w_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/data14/kyh/ewiniar/0700-0806/new_data/org/w/'//ncfilename, ncid)
    call get4d(ncid, 'vertical_wind', nxt, nyt, nz, nt1, w_org)
    if(t.eq.0) then
      call get1d(ncid, 'west_east', nxt, x_org)
      do i=sx, ex, 1
        x(i-sx+1)=(x_org(i)-1)*27.
      end do
      call get1d(ncid, 'south_north', nyt, y_org)
      do j=sy, ey, 1
        y(j-sy+1)=(y_org(j)-1)*27.
      end do
      call get1d(ncid, 'height', nz, z)
    end if
    call closenc(ncid)


    do tt=1, nt1
      t=t+1
      do k=1, nz
      do j=sy, ey, 1
      do i=sx, ex, 1
        w4d(i-sx+1,j-sy+1,k,t)=w_org(i,j,k,tt) 
      end do
      end do
      end do
    end do 

    deallocate(w_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file

    allocate(w_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/data14/kyh/ewiniar/0700-0806/new_data/org/w/'//ncfilename, ncid)
    call get4d(ncid, 'vertical_wind', nxt, nyt, nz, 1, w_org)
    call closenc(ncid)

    t=t+1
    do k=1, nz
    do j=sy, ey, 1
    do i=sx, ex, 1
      w4d(i-sx+1,j-sy+1,k,t)=w_org(i,j,k,1) 
    end do
    end do
    end do

    deallocate(w_org)

!    ncfilename='/data4/atmosdata/T_SPECTRUM/'//speriod//'/rho_ave_profile_'//speriod//'.nc'
!    call opennc(ncfilename, ncid)
!    call get1d(ncid,'rho', nz, rho_tmean)
!    call closenc(ncid)
rho_tmean(:) = exp(-(z(:)-20.e3)/7.e3)

!== Background wind (21x21 running mean) ==============================

    allocate(wb_org(nx1,ny1,nz,nt))
    call opennc('/data14/kyh/ewiniar/w_background_'//speriod//'.nc', ncid)
    call get4d(ncid, 'wb', nx1, ny1, nz, nt, wb_org)
    call closenc(ncid)

    do t=1, nt
    do k=1, nz
    do j=sy1, ey1, 1
    do i=sx1, ex1, 1
      wb(i-sx1+1,j-sy1+1,k,t)=wb_org(i,j,k,t) !0.
    end do
    end do
    end do
    end do

    ifdir='/data4/atmosdata/T_SPECTRUM/'//speriod

!    ncfilename=trim(ifdir)//'/u_profile_ave_'//speriod//'.nc'
!    call opennc(ncfilename, ncid)
!    call get1d(ncid,'u_ave', nz, u0)
!    call closenc(ncid)
!
!    ncfilename=trim(ifdir)//'/v_profile_ave_'//speriod//'.nc'
!    call opennc(ncfilename, ncid)
!    call get1d(ncid,'v_ave', nz, v0)
!    call closenc(ncid)
!
!    ncfilename=trim(ifdir)//'/n_profile_ave_'//speriod//'.nc'
!    call opennc(ncfilename, ncid)
!    call get1d(ncid,'n_ave', nz, nb)
!    call closenc(ncid)

!== 3D FFT with respect to x, y, and t ================================

    w4d = w4d - wb

    w2 = 0.

    do k=1, nz

      w3d(:,:,:) = w4d(:,:,k,:)

      call fft3df(nx, ny, nt, dx, dy, dt, w3d, freq_x, freq_y, freq_t, coefw)

      do t=1, nt
      do j=1, ny
      do i=1, nx
!        omeh = 2.*pi*( -freq_t(t) - u0(k)*freq_x(i) - v0(k)*freq_y(j) )
!        cc = omeh*omeh
!        if ( cc <= f0*f0*float(cf) .or. cc >= nb(k)*nb(k)/float(cf) ) then
!          coefw(i,j,t) = 0.
!        end if
      end do
      end do
      end do
      coefw(1,1,:) = 0.

      call fft3db(nx, ny, nt, coefw, w3d)

      ww3d(:,:,:) = w3d(:,:,:)*w3d(:,:,:)

      do t=1, nt
      do j=1, ny
      do i=1, nx
        w2(k) = w2(k) + ww3d(i,j,t)
      end do
      end do
      end do

      w2(k) = rho_tmean(k) * w2(k) / float(nx*ny*nt)

    end do


!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') '/data14/kyh/ewiniar/w2_'//speriod//'_f',cf,'.nc'
    call out1d_yh(ncfilename, 1, (/'w2'/), &
               (/w2/), &
               'height_m',nz,z,  &
               'mfx')

    end program  

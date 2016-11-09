    program MF 

    use netcdfio
    use fft

    implicit none

    integer, parameter ::  cf = 10

    integer :: i, j, k, t, file, tt, kk, l1, l2, iz
    integer :: ncid
    integer :: day, hour
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: start_day=7, start_hour=1
    integer, parameter :: nfile=29
    integer, parameter :: nxt=186, nyt=186, nz=85, nt1=6
    integer, parameter :: nt=nfile*nt1+1
    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants
    real, parameter :: g=9.81
    real, parameter :: lat = 20.
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
 
    real, parameter :: dx=27000., dy=27000., dz=500., dt=600.

    character(len=100) :: ncfilename, ifdir
    character(len=50) :: varname
    character(len=50) :: d1name,d2name, d3name, d4name
    character(len=50) :: title

    real, dimension(nxt) :: x_org
    real, dimension(nyt) :: y_org
    real, dimension(nx) :: x
    real, dimension(ny) :: y
    real, dimension(nz) :: z, t0, u0, v0, nb
    real, dimension(nt) :: h

    real, dimension(nx) :: freq_x
    real, dimension(ny) :: freq_y
    real, dimension(nt) :: freq_t

    real, dimension(nx, ny, nz, nt, 4) :: var
    real, dimension(nx,ny,nz,nt) :: tb, ub, vb 
 
    real, allocatable, dimension(:,:,:) :: t3d, u3d, v3d
    real, allocatable, dimension(:,:,:,:) :: t_org, u_org, v_org
    real, allocatable, dimension(:,:,:,:) :: t4d, u4d, v4d
    real, allocatable, dimension(:,:,:,:) :: tb_org, ub_org, vb_org

    double complex, allocatable, dimension(:,:,:) :: coeft, coefu, coefv
    double complex, allocatable, dimension(:,:,:) :: coef_up, coef_dn, coef_up2, coef_dn2

    real    ::  f0, cc, omeh
    complex ::  dd

    f0 = 2.*omega*sin(lat*pi/180.)

!== READ WRF OUTPUT =====================================================
 
    day=start_day
    hour=start_hour
    t=0
   
    allocate(t4d(nx,ny,nz,nt), u4d(nx,ny,nz,nt), v4d(nx,ny,nz,nt))
 
    do file=1, nfile  ! Loop: file 

    allocate(t_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    call get4d(ncid, 'temperature', nxt, nyt, nz, nt1, t_org)
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

    allocate(u_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_u_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/u_wind/'//ncfilename, ncid)
    call get4d(ncid, 'zonal_wind', nxt, nyt, nz, nt1, u_org)
    call closenc(ncid)

    allocate(v_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_v_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/v_wind/'//ncfilename, ncid)
    call get4d(ncid, 'meridional_wind', nxt, nyt, nz, nt1, v_org)
    call closenc(ncid)

    do tt=1, nt1
      t=t+1
      do k=1, nz
      do j=sy, ey, 1
      do i=sx, ex, 1
        t4d(i-sx+1,j-sy+1,k,t)=t_org(i,j,k,tt)
        u4d(i-sx+1,j-sy+1,k,t)=u_org(i,j,k,tt)
        v4d(i-sx+1,j-sy+1,k,t)=v_org(i,j,k,tt)
      end do
      end do
      end do
    end do 

    deallocate(t_org, u_org, v_org)

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

    allocate(u_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_u_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/u_wind/'//ncfilename, ncid)
    call get4d(ncid, 'zonal_wind', nxt, nyt, nz, 1, u_org)
    call closenc(ncid)

    allocate(v_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_v_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/v_wind/'//ncfilename, ncid)
    call get4d(ncid, 'meridional_wind', nxt, nyt, nz, 1, v_org)
    call closenc(ncid)

    t=t+1
    do k=1, nz
    do j=sy, ey, 1
    do i=sx, ex, 1
      t4d(i-sx+1,j-sy+1,k,t)=t_org(i,j,k,1) 
      u4d(i-sx+1,j-sy+1,k,t)=u_org(i,j,k,1)
      v4d(i-sx+1,j-sy+1,k,t)=v_org(i,j,k,1)
    end do
    end do
    end do

    deallocate(t_org, u_org, v_org)

!== Background wind (21x21 running mean) ==============================

    allocate(tb_org(nx1,ny1,nz,nt))
    call opennc('/data1/ksy/ewiniar-MF/0700-0806/running-mean-back/t_background.nc', ncid)
    call get4d(ncid, 'tb', nx1, ny1, nz, nt, tb_org)
    call closenc(ncid)

    allocate(ub_org(nx1,ny1,nz,nt))
    call opennc('/data1/ksy/ewiniar-MF/0700-0806/running-mean-back/u_background.nc', ncid)
    call get4d(ncid, 'ub', nx1, ny1, nz, nt, ub_org)
    call closenc(ncid)

    allocate(vb_org(nx1,ny1,nz,nt))
    call opennc('/data1/ksy/ewiniar-MF/0700-0806/running-mean-back/v_background.nc', ncid)
    call get4d(ncid, 'vb', nx1, ny1, nz, nt, vb_org)
    call closenc(ncid)

    do t=1, nt
    do k=1, nz
    do j=sy1, ey1, 1
    do i=sx1, ex1, 1
      tb(i-sx1+1,j-sy1+1,k,t)=tb_org(i,j,k,t)
      ub(i-sx1+1,j-sy1+1,k,t)=ub_org(i,j,k,t)
      vb(i-sx1+1,j-sy1+1,k,t)=vb_org(i,j,k,t)
    end do
    end do
    end do
    end do

    ifdir='/data4/atmosdata/T_SPECTRUM/'//speriod

    ncfilename=trim(ifdir)//'/u_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'u_ave', nz, u0)
    call closenc(ncid)

    ncfilename=trim(ifdir)//'/v_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'v_ave', nz, v0)
    call closenc(ncid)

    ncfilename=trim(ifdir)//'/n_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'n_ave', nz, nb)
    call closenc(ncid)

    ncfilename='res/t_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'t_ave', nz, t0)
    call closenc(ncid)

!== 3D FFT with respect to x, y, and t ================================

    allocate(t3d(nx,ny,nt), u3d(nx,ny,nt), v3d(nx,ny,nt))
    allocate(coeft(nx,ny,nt), coefu(nx,ny,nt), coefv(nx,ny,nt))

    do k=1, nz

      do j=1, ny
      do t=1, nt
      do i=1, nx
        t3d(i,j,t)=t4d(i,j,k,t)-tb(i,j,k,t)
        u3d(i,j,t)=u4d(i,j,k,t)-ub(i,j,k,t)
        v3d(i,j,t)=v4d(i,j,k,t)-vb(i,j,k,t)
      end do
      end do
      end do

      call fft3df(nx, ny, nt, dx, dy, dt, t3d, freq_x, freq_y, freq_t, coeft)
      call fft3df(nx, ny, nt, dx, dy, dt, u3d, freq_x, freq_y, freq_t, coefu)
      call fft3df(nx, ny, nt, dx, dy, dt, v3d, freq_x, freq_y, freq_t, coefv)

      allocate(coef_up(nx,ny,nt), coef_dn(nx,ny,nt))
      allocate(coef_up2(nx,ny,nt), coef_dn2(nx,ny,nt))
      coef_up(:,:,:) = 0.
      coef_dn(:,:,:) = 0.
      coef_up2(:,:,:) = 0.
      coef_dn2(:,:,:) = 0.

      do t=1, nt
      do j=1, ny
      do i=1, nx/2+1
        omeh = 2.*pi*( -freq_t(t) - u0(k)*freq_x(i) - v0(k)*freq_y(j) )
        cc = omeh*omeh
if (omeh > 0. .and. freq_y(j) == 0.) then
        if ( cc > f0*f0*float(cf) .and. cc < nb(k)*nb(k)/float(cf) ) then
          dd = (coefu(i,j,t)*freq_x(i)+coefv(i,j,t)*freq_y(j)) / &
               sqrt(freq_x(i)*freq_x(i)+freq_y(j)*freq_y(j)) *   &
               (0.,-1.)*nb(k)*t0(k)/g
          coef_up(i,j,t) = 0.5*(coeft(i,j,t) + dd)
          coef_dn(i,j,t) = 0.5*(coeft(i,j,t) - dd)

          coef_up(i,j,t) = coef_up(i,j,t)*(0.,1.)*g/nb(k)/t0(k)
          coef_dn(i,j,t) = coef_dn(i,j,t)*(0.,1.)*g/nb(k)/t0(k)
          coef_up2(i,j,t) = coef_up(i,j,t)*omeh/nb(k)
          coef_dn2(i,j,t) = coef_dn(i,j,t)*omeh/nb(k)
          coef_dn(i,j,t) = -coef_dn(i,j,t)
        end if
end if
      end do
      end do
      end do
      coef_up(1,:,:) = coef_up(1,:,:)*0.5
      coef_dn(1,:,:) = coef_dn(1,:,:)*0.5
      coef_up(1,1,:) = 0.
      coef_dn(1,1,:) = 0.
      coef_up2(1,:,:) = coef_up2(1,:,:)*0.5
      coef_dn2(1,:,:) = coef_dn2(1,:,:)*0.5
      coef_up2(1,1,:) = 0.
      coef_dn2(1,1,:) = 0.

      call fft3db(nx, ny, nt, coef_up, var(:,:,k,:,1))
      call fft3db(nx, ny, nt, coef_dn, var(:,:,k,:,2))
      call fft3db(nx, ny, nt, coef_up2, var(:,:,k,:,3))
      call fft3db(nx, ny, nt, coef_dn2, var(:,:,k,:,4))

    end do

    deallocate(coeft, coefu, coefv)

!== Write Output ========================================================

    do t=1, nt
      h(t) = 1. + float(t-1)*dt/3600.
    enddo

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/field/wp_'//speriod//'_f',cf,'-e1.nc'
    call out4d_yh(ncfilename, 4, (/'u_up','u_dn','w_up','w_dn'/), &
               var, &
               'west_east' , nx, x, &
               'sout_north', ny, y, &
               'height'    , nz, z, &
               'time'      , nt, h, &
               'T field')

    end program  

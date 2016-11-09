    program MF 

    use netcdfio
    use fft

    implicit none

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

    real, dimension(nx) :: freq_x
    real, dimension(ny) :: freq_y
    real, dimension(nt) :: freq_t
    real, dimension(nx/2+1,ny,nt) :: psd, tup, tdn

    real, dimension(nx/2+1, (ny/2)*2+1, (nt/2)*2+1, nz, 3) :: var
    real, dimension(nx,ny,nz,nt) :: tb, ub, vb 
 
    real, allocatable, dimension(:) :: axis1, axis2, axis3, axis4
    real, allocatable, dimension(:,:,:) :: t3d, u3d, v3d
    real, allocatable, dimension(:,:,:,:) :: t_org, u_org, v_org
    real, allocatable, dimension(:,:,:,:) :: t4d, u4d, v4d
    real, allocatable, dimension(:,:,:,:) :: tb_org, ub_org, vb_org

    double complex, allocatable, dimension(:,:,:) :: coeft, coefu, coefv

    real    ::  cc, f0
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

      do t=1, nt
      do j=1, ny
      do i=1, nx/2+1
        psd(i,j,t)=real(cdabs(coeft(i,j,t)))**2
        cc = (2.*pi*( -freq_t(t) - u0(k)*freq_x(i) - v0(k)*freq_y(j) ))**2
        if ( cc > f0*f0 .and. cc < nb(k)*nb(k) ) then
          dd = (coefu(i,j,t)*freq_x(i)+coefv(i,j,t)*freq_y(j)) / &
               sqrt(freq_x(i)*freq_x(i)+freq_y(j)*freq_y(j)) *   &
               (0.,-1.)*nb(k)*t0(k)/g
          tup(i,j,t) = 0.25*real(cdabs(coeft(i,j,t) + dd))**2
          tdn(i,j,t) = 0.25*real(cdabs(coeft(i,j,t) - dd))**2
        else
          tup(i,j,t) = 0.
          tdn(i,j,t) = 0.
        end if
      end do
      end do
      end do
      psd(:,:,:) = psd(:,:,:)*2*dx*dy*dt/float(nx)/float(ny)/float(nt)
      tup(:,:,:) = tup(:,:,:)*2*dx*dy*dt/float(nx)/float(ny)/float(nt)
      tdn(:,:,:) = tdn(:,:,:)*2*dx*dy*dt/float(nx)/float(ny)/float(nt)
      psd(1,:,:) = psd(1,:,:)*0.5
      tup(1,:,:) = tup(1,:,:)*0.5
      tdn(1,:,:) = tdn(1,:,:)*0.5
      psd(1,1,:) = 0.
      tup(1,1,:) = 0.
      tdn(1,1,:) = 0.


    do i=1, nx/2+1
      do j=1, ny/2
        do t=1, nt/2
          var(i,j,t,k,1)=psd(i,ny/2+j+mod(ny,2),nt/2+t+mod(nt,2))
          var(i,j,t,k,2)=tup(i,ny/2+j+mod(ny,2),nt/2+t+mod(nt,2))
          var(i,j,t,k,3)=tdn(i,ny/2+j+mod(ny,2),nt/2+t+mod(nt,2))
        end do
        do t=nt/2+1, (nt/2)*2+1, 1
          var(i,j,t,k,1)=psd(i,ny/2+j+mod(ny,2),t-nt/2)
          var(i,j,t,k,2)=tup(i,ny/2+j+mod(ny,2),t-nt/2)
          var(i,j,t,k,3)=tdn(i,ny/2+j+mod(ny,2),t-nt/2)
        end do
      end do
      do j=ny/2+1, (ny/2)*2+1, 1
        do t=1, nt/2
          var(i,j,t,k,1)=psd(i,j-ny/2,nt/2+t+mod(nt,2))
          var(i,j,t,k,2)=tup(i,j-ny/2,nt/2+t+mod(nt,2))
          var(i,j,t,k,3)=tdn(i,j-ny/2,nt/2+t+mod(nt,2))
        end do
        do t=nt/2+1, (nt/2)*2+1, 1
          var(i,j,t,k,1)=psd(i,j-ny/2,t-nt/2)
          var(i,j,t,k,2)=tup(i,j-ny/2,t-nt/2)
          var(i,j,t,k,3)=tdn(i,j-ny/2,t-nt/2)
        end do
      end do
    end do

    end do

    deallocate(coeft, coefu, coefv)

!== Write Output ========================================================

    allocate(axis1(nx/2+1))
    allocate(axis2((ny/2)*2+1))
    allocate(axis3((nt/2)*2+1))
    allocate(axis4(nz))

    do i=1, nx/2+1
      axis1(i)=freq_x(i)
    end do

    axis2(1)=-1*freq_y(ny/2+1)
    do j=2, ny/2
      axis2(j)=freq_y(ny/2+j+mod(ny,2))
    end do
    do j=ny/2+1, (ny/2)*2+1, 1
      axis2(j)=freq_y(j-ny/2)
    end do

    axis3(1)=-1*freq_t(nt/2+1)
    do t=2, nt/2
      axis3(t)=freq_t(nt/2+t+mod(nt,2))
    end do
    do t=nt/2+1, (nt/2)*2+1, 1
      axis3(t)=freq_t(t-nt/2)
    end do

    do k=1, nz
      axis4(k)=z(k)
    end do

    dim1=nx/2+1
    dim2=(ny/2)*2+1
    dim3=(nt/2)*2+1
    dim4=nz
    d1name='k'
    d2name='l'
    d3name='f'
    d4name='height_m'
    title='psd3d_tp_k_l_f_z'

    ncfilename='res/psd_tp_k_l_f_'//speriod//'.nc'
    call out4d_yh(ncfilename, 3, (/'psd ','t2up','t2dn'/), &
               var, &
               d1name, dim1, axis1, &
               d2name, dim2, axis2, &
               d3name, dim3, axis3, &
               d4name, dim4, axis4, &
               title)


    end program  

    program MF 

    use netcdfio
    use fft

    implicit none

    integer :: i, j, k, t, file, tt, kk 
    integer :: ncid
    integer :: day, hour
    integer :: dim1, dim2, dim3, dim4

    character(len=2), parameter ::  speriod='p1'
    integer, parameter :: start_day=7, start_hour=1
    integer, parameter :: nfile=29
    integer, parameter :: nxt=186, nyt=186, nz=85, nt1=1
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
 
    real, parameter :: dx=27000., dy=27000., dz=500., dt=3600.

    character(len=100) :: ncfilename
    character(len=50) :: varname
    character(len=50) :: d1name,d2name, d3name, d4name
    character(len=50) :: title

    real, dimension(nxt) :: x_org
    real, dimension(nyt) :: y_org
    real, dimension(nx) :: x
    real, dimension(ny) :: y
    real, dimension(nz) :: z
    real, dimension(nz) :: rho_tmean, wu, wv, wu_p, wv_p, wu_n, wv_n

    real, dimension(nx,ny,nz,nt) :: ub, vb, wb

    real, dimension(nx,ny,nt) :: wu3d, wv3d
    real, allocatable, dimension(:,:,:,:) :: w4d, w_org
    real, allocatable, dimension(:,:,:,:) :: u4d, u_org
    real, allocatable, dimension(:,:,:,:) :: v4d, v_org
    real, allocatable, dimension(:,:,:,:) :: ub_org, vb_org, wb_org

!== READ WRF OUTPUT =====================================================
 
    day=start_day
    hour=start_hour
    t=0
   
    allocate(w4d(nx,ny,nz,nt))
    allocate(u4d(nx,ny,nz,nt))
    allocate(v4d(nx,ny,nz,nt))
 
    do file=1, nfile  ! Loop: file 

    allocate(w_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/w/'//ncfilename, ncid)
    call geta4d_yh(ncid, 'vertical_wind', 1,nxt, 1,nyt, 1,nz, 1,nt1, w_org)
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


    allocate(u_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_u_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/u_wind/'//ncfilename, ncid)
    call geta4d_yh(ncid, 'zonal_wind', 1,nxt, 1,nyt, 1,nz, 1,nt1, u_org)
    call closenc(ncid)



    allocate(v_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_v_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/v_wind/'//ncfilename, ncid)
    call geta4d_yh(ncid, 'meridional_wind', 1,nxt, 1,nyt, 1,nz, 1,nt1, v_org)
    call closenc(ncid)

    do tt=1, nt1
      t=t+1
      do k=1, nz
      do j=sy, ey, 1
      do i=sx, ex, 1
        w4d(i-sx+1,j-sy+1,k,t)=w_org(i,j,k,tt) 
        u4d(i-sx+1,j-sy+1,k,t)=u_org(i,j,k,tt) 
        v4d(i-sx+1,j-sy+1,k,t)=v_org(i,j,k,tt) 
      end do
      end do
      end do
    end do 

    deallocate(w_org)
    deallocate(u_org)
    deallocate(v_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file

    allocate(w_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_w_d01_2006-07-',day,'_',hour,':00:00'
    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/w/'//ncfilename, ncid)
    call get4d(ncid, 'vertical_wind', nxt, nyt, nz, 1, w_org)
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
      w4d(i-sx+1,j-sy+1,k,t)=w_org(i,j,k,1) 
      u4d(i-sx+1,j-sy+1,k,t)=u_org(i,j,k,1) 
      v4d(i-sx+1,j-sy+1,k,t)=v_org(i,j,k,1) 
    end do
    end do
    end do

    deallocate(w_org)
    deallocate(u_org)
    deallocate(v_org)

    ncfilename='/data4/atmosdata/T_SPECTRUM/'//speriod//'/rho_ave_profile_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'rho', nz, rho_tmean)
    call closenc(ncid)

!== Background wind (21x21 running mean) ==============================

    allocate(ub_org(nx1,ny1,nz,175))
    call opennc('/data1/ksy/ewiniar-MF/0700-0806/running-mean-back/u_background.nc', ncid)
    call get4d(ncid, 'ub', nx1, ny1, nz, 175, ub_org)
    call closenc(ncid)

    allocate(vb_org(nx1,ny1,nz,175))
    call opennc('/data1/ksy/ewiniar-MF/0700-0806/running-mean-back/v_background.nc', ncid)
    call get4d(ncid, 'vb', nx1, ny1, nz, 175, vb_org)
    call closenc(ncid)

    allocate(wb_org(nx1,ny1,nz,175))
    call opennc('../res/w_background_'//speriod//'.nc', ncid)
    call get4d(ncid, 'wb', nx1, ny1, nz, 175, wb_org)
    call closenc(ncid)

    do t=1, nt
    do k=1, nz
    do j=sy1, ey1, 1
    do i=sx1, ex1, 1
      wb(i-sx1+1,j-sy1+1,k,t)=wb_org(i,j,k,t*6-5) !0.
      ub(i-sx1+1,j-sy1+1,k,t)=ub_org(i,j,k,t*6-5) 
      vb(i-sx1+1,j-sy1+1,k,t)=vb_org(i,j,k,t*6-5) 
    end do
    end do
    end do
    end do

!== 3D FFT with respect to x, y, and t ================================

    w4d = w4d - wb
    u4d = u4d - ub
    v4d = v4d - vb

    wu_p = 0.
    wu_n = 0.
    wv_p = 0.
    wv_n = 0.

    do k=1, nz

      wu3d(:,:,:) = w4d(:,:,k,:)*u4d(:,:,k,:)
      wv3d(:,:,:) = w4d(:,:,k,:)*v4d(:,:,k,:)

      do t=1, nt
      do j=1, ny
      do i=1, nx
        if (wu3d(i,j,t) > 0.) then
          wu_p(k) = wu_p(k) + wu3d(i,j,t)
        else
          wu_n(k) = wu_n(k) + wu3d(i,j,t)
        end if
        if (wv3d(i,j,t) > 0.) then
          wv_p(k) = wv_p(k) + wv3d(i,j,t)
        else
          wv_n(k) = wv_n(k) + wv3d(i,j,t)
        end if
      end do
      end do
      end do

      wu_p(k) = rho_tmean(k) * wu_p(k) / float(nx*ny*nt)
      wu_n(k) = rho_tmean(k) * wu_n(k) / float(nx*ny*nt)
      wv_p(k) = rho_tmean(k) * wv_p(k) / float(nx*ny*nt)
      wv_n(k) = rho_tmean(k) * wv_n(k) / float(nx*ny*nt)

    end do

    wu = wu_p + wu_n
    wv = wv_p + wv_n


!== Write Output ========================================================

    ncfilename='res/wind/mean-mf/mean-mfx_'//speriod//'.nc'
    call out1d_yh(ncfilename, 3, (/'mfx  ','mfx_p','mfx_n'/), &
               (/wu,wu_p,wu_n/), &
               'height_m',nz,z,  &
               'mfx')

    ncfilename='res/wind/mean-mf/mean-mfy_'//speriod//'.nc'
    call out1d_yh(ncfilename, 3, (/'mfy  ','mfy_p','mfy_n'/), &
               (/wv,wv_p,wv_n/), &
               'height_m',nz,z,  &
               'mfy')


    end program  

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
    integer, parameter :: nxt=186, nyt=186, nz=85, nt1=6
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

    character(len=100) :: ncfilename
    character(len=50) :: varname
    character(len=50) :: d1name,d2name, d3name, d4name
    character(len=50) :: title

    real, dimension(nxt) :: x_org
    real, dimension(nyt) :: y_org
    real, dimension(nx) :: x
    real, dimension(ny) :: y
    real, dimension(nz) :: z
    real, dimension(nz) :: rho_tmean

    real, dimension(nx) :: freq_x
    real, dimension(ny) :: freq_y
    real, dimension(nt) :: freq_t
    real, dimension(nx/2+1,ny,nz,nt) :: copsdx, copsdy
    real, dimension(nx/2+1, (ny/2)*2+1, (nt/2)*2+1, nz) :: var
    real, dimension(nx,ny,nz,nt) :: ub, vb, wb
 
    real, allocatable, dimension(:) :: axis1, axis2, axis3, axis4
    real, allocatable, dimension(:,:,:) :: u3d, v3d, w3d
    real, allocatable, dimension(:,:,:,:) :: w4d, w_org 
    real, allocatable, dimension(:,:,:,:) :: u4d, u_org 
    real, allocatable, dimension(:,:,:,:) :: v4d, v_org 
    real, allocatable, dimension(:,:,:,:) :: ub_org, vb_org, wb_org
 
    double complex, allocatable, dimension(:,:,:) :: coefu, coefv, coefw

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

    allocate(ub_org(nx1,ny1,nz,nt))
    call opennc('res/u_background_'//speriod//'.nc', ncid)
    call get4d(ncid, 'ub', nx1, ny1, nz, nt, ub_org)
    call closenc(ncid)

    allocate(vb_org(nx1,ny1,nz,nt))
    call opennc('res/v_background_'//speriod//'.nc', ncid)
    call get4d(ncid, 'vb', nx1, ny1, nz, nt, vb_org)
    call closenc(ncid)

    allocate(wb_org(nx1,ny1,nz,nt))
    call opennc('res/w_background_'//speriod//'.nc', ncid)
    call get4d(ncid, 'wb', nx1, ny1, nz, nt, wb_org)
    call closenc(ncid)

    do t=1, nt
    do k=1, nz
    do j=sy1, ey1, 1
    do i=sx1, ex1, 1
      wb(i-sx1+1,j-sy1+1,k,t)=wb_org(i,j,k,t) !0.
      ub(i-sx1+1,j-sy1+1,k,t)=ub_org(i,j,k,t) 
      vb(i-sx1+1,j-sy1+1,k,t)=vb_org(i,j,k,t) 
    end do
    end do
    end do
    end do

!== 3D FFT with respect to x, y, and t ================================

    allocate(u3d(nx,ny,nt))
    allocate(v3d(nx,ny,nt))
    allocate(w3d(nx,ny,nt))
    allocate(coefu(nx,ny,nt))
    allocate(coefv(nx,ny,nt))
    allocate(coefw(nx,ny,nt))

    do k=1, nz
    do j=1, ny
    do i=1, nx/2+1
    do t=1, nt
      copsdx(i,j,k,t)=0.
      copsdy(i,j,k,t)=0.
    end do
    end do
    end do
    end do

    do k=1, nz

      do j=1, ny
      do t=1, nt
      do i=1, nx
        w3d(i,j,t)=w4d(i,j,k,t)-wb(i,j,k,t)
        u3d(i,j,t)=u4d(i,j,k,t)-ub(i,j,k,t)
        v3d(i,j,t)=v4d(i,j,k,t)-vb(i,j,k,t)
      end do
      end do
      end do

      call fft3df(nx, ny, nt, dx, dy, dt, w3d, freq_x, freq_y, freq_t, coefw)
      call fft3df(nx, ny, nt, dx, dy, dt, u3d, freq_x, freq_y, freq_t, coefu)
      call fft3df(nx, ny, nt, dx, dy, dt, v3d, freq_x, freq_y, freq_t, coefv)

      do t=1, nt
      do j=1, ny
      do i=1, nx/2+1
        copsdx(i,j,k,t)=rho_tmean(k)*real(coefu(i,j,t)*conjg(coefw(i,j,t)))
        copsdy(i,j,k,t)=rho_tmean(k)*real(coefv(i,j,t)*conjg(coefw(i,j,t)))
      end do
      end do
      end do

    end do

    copsdx(:,:,:,:) = copsdx(:,:,:,:)*2*dx*dy*dt/float(nx)/float(ny)/float(nt)
    copsdy(:,:,:,:) = copsdy(:,:,:,:)*2*dx*dy*dt/float(nx)/float(ny)/float(nt)
    copsdx(1,:,:,:) = copsdx(1,:,:,:)*0.5
    copsdy(1,:,:,:) = copsdy(1,:,:,:)*0.5
    copsdx(1,1,:,:) = 0.
    copsdy(1,1,:,:) = 0.

    deallocate(coefu)
    deallocate(coefv)
    deallocate(coefw)


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

    do k=1, nz
    do i=1, nx/2+1
      do j=1, ny/2
        do t=1, nt/2
          var(i,j,t,k)=copsdx(i,ny/2+j+mod(ny,2),k,nt/2+t+mod(nt,2))
        end do
        do t=nt/2+1, (nt/2)*2+1, 1
          var(i,j,t,k)=copsdx(i,ny/2+j+mod(ny,2),k,t-nt/2)
        end do
      end do
      do j=ny/2+1, (ny/2)*2+1, 1
        do t=1, nt/2
          var(i,j,t,k)=copsdx(i,j-ny/2,k,nt/2+t+mod(nt,2))
        end do
        do t=nt/2+1, (nt/2)*2+1, 1
          var(i,j,t,k)=copsdx(i,j-ny/2,k,t-nt/2)
        end do
      end do
    end do
    end do
      
    dim1=nx/2+1
    dim2=(ny/2)*2+1
    dim3=(nt/2)*2+1
    dim4=nz
    ncfilename='res/wind/mf_k_l_f/mfx_k_l_f_'//speriod//'.nc'
    varname='mfx'
    d1name='k'
    d2name='l'
    d3name='f'
    d4name='height_m'
    title='mf3d_k_l_f_z'

    call out4d(ncfilename, varname, var, &
               d1name, dim1, axis1, &
               d2name, dim2, axis2, &
               d3name, dim3, axis3, &
               d4name, dim4, axis4, &
               title)

    do k=1, nz
    do i=1, nx/2+1
      do j=1, ny/2
        do t=1, nt/2
          var(i,j,t,k)=copsdy(i,ny/2+j+mod(ny,2),k,nt/2+t+mod(nt,2))
        end do
        do t=nt/2+1, (nt/2)*2+1, 1
          var(i,j,t,k)=copsdy(i,ny/2+j+mod(ny,2),k,t-nt/2)
        end do
      end do
      do j=ny/2+1, (ny/2)*2+1, 1
        do t=1, nt/2
          var(i,j,t,k)=copsdy(i,j-ny/2,k,nt/2+t+mod(nt,2))
        end do
        do t=nt/2+1, (nt/2)*2+1, 1
          var(i,j,t,k)=copsdy(i,j-ny/2,k,t-nt/2)
        end do
      end do
    end do
    end do

    ncfilename='res/wind/mf_k_l_f/mfy_k_l_f_'//speriod//'.nc'
    varname='mfy'

    call out4d(ncfilename, varname, var, &
               d1name, dim1, axis1, &
               d2name, dim2, axis2, &
               d3name, dim3, axis3, &
               d4name, dim4, axis4, &
               title)

    end program  

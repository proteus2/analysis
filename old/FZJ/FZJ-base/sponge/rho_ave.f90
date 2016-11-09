    program MF 

    use netcdfio

    implicit none

    integer :: i, j, k, t, file, tt, kk 
    integer :: ncid
    integer :: day, hour
    integer :: dim1, dim2

    integer, parameter :: start_day=7, start_hour=1
    integer, parameter :: nfile=29
    integer, parameter :: nxt=186, nyt=186, nz=101, nt1=6
    integer, parameter :: nt=nfile*nt1+1

    integer, parameter :: sx=32
    integer, parameter :: ex=116
    integer, parameter :: sy=32
    integer, parameter :: ey=116 

    integer, parameter :: nx=ex-sx+1
    integer, parameter :: ny=ey-sy+1
 
    real, parameter :: r_d=287.

    character(len=100) :: ncfilename
    character(len=50) :: varname
    character(len=50) :: d1name,d2name, d3name, d4name
    character(len=50) :: title

    real, dimension(nz) :: z
    
    real, allocatable, dimension(:) :: rho_tmean
    real, allocatable, dimension(:,:) :: rhomean

    real, allocatable, dimension(:,:,:,:) :: t_org 
    real, allocatable, dimension(:,:,:,:) :: pb_org 
    real, allocatable, dimension(:,:,:,:) :: pp_org 
    real, allocatable, dimension(:,:,:,:) :: rho4d
    real, allocatable, dimension(:) :: axis1, axis2, x, y
 
!== READ WRF OUTPUT =====================================================
 
    day=start_day
    hour=start_hour
    t=0
   
    allocate(rho4d(nx,ny,nz,nt))
 
    do file=1, nfile  ! Loop: file 

    allocate(t_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
!    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    varname='temperature'
    call get4d(ncid, varname, nxt, nyt, nz, nt1, t_org)
    if(t.eq.0) then
      varname='height'
      call get1d(ncid, varname, nz, z)
    end if
    call closenc(ncid)

    allocate(pb_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_pb_d01_2006-07-',day,'_',hour,':00:00'
!    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    varname='basic_pressure'
    call get4d(ncid, varname, nxt, nyt, nz, nt1, pb_org)
    call closenc(ncid)

    allocate(pp_org(nxt,nyt,nz,nt1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_pp_d01_2006-07-',day,'_',hour,':00:00'
!    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    varname='perturbation_pressure'
    call get4d(ncid, varname, nxt, nyt, nz, nt1, pp_org)
    call closenc(ncid)
   
    do tt=1, nt1
      t=t+1
      do k=1, nz
      do j=sy, ey, 1
      do i=sx, ex, 1
        rho4d(i-sx+1,j-sy+1,k,t)=(pb_org(i,j,k,tt)+pp_org(i,j,k,tt))/r_d/t_org(i,j,k,tt) 
      end do
      end do
      end do
    end do 

    deallocate(t_org)
    deallocate(pb_org)
    deallocate(pp_org)

    hour=hour+1
    if(hour.eq.24) then
      day=day+1
      hour=0
    end if

    end do  ! End of Loop: file

!-- For last file

    allocate(t_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_t_d01_2006-07-',day,'_',hour,':00:00'
!    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/new/temperature/'//ncfilename, ncid)
    varname='temperature'
    call get4d(ncid, varname, nxt, nyt, nz, 1, t_org)
    call closenc(ncid)

    allocate(pb_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_pb_d01_2006-07-',day,'_',hour,':00:00'
!    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    varname='basic_pressure'
    call get4d(ncid, varname, nxt, nyt, nz, 1, pb_org)
    call closenc(ncid)

    allocate(pp_org(nxt,nyt,nz,1))
    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'new_pp_d01_2006-07-',day,'_',hour,':00:00'
!    call opennc('/export30/ksy/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    call opennc('/data2/ksy/exp30-back/ewiniar/WRF-0.1mb/0700-0806/new_data/org/pressure/'//ncfilename, ncid)
    varname='perturbation_pressure'
    call get4d(ncid, varname, nxt, nyt, nz, 1, pp_org)
    call closenc(ncid)



    t=t+1
    do k=1, nz
    do j=sy, ey, 1
    do i=sx, ex, 1
      rho4d(i-sx+1,j-sy+1,k,t)=(pb_org(i,j,k,1)+pp_org(i,j,k,1))/r_d/t_org(i,j,k,1) 
    end do
    end do
    end do

    deallocate(t_org)
    deallocate(pb_org)
    deallocate(pp_org)


!== Calculate mean-state (z,t) ========================================
 
    allocate(rhomean(nz,nt))  
    allocate(rho_tmean(nz))

    do t=1, nt
    do k=1, nz
      rhomean(k,t)=0.
      do i=1, nx
      do j=1, ny
        rhomean(k,t)=rhomean(k,t)+rho4d(i,j,k,t)
      end do
      end do
      rhomean(k,t)=rhomean(k,t)/nx/ny
    end do
    end do

    do k=1, nz
      rho_tmean(k)=0.
      do t=1, nt
        rho_tmean(k)=rho_tmean(k)+rhomean(k,t)
      end do
      rho_tmean(k)=rho_tmean(k)/nt
    end do


!== Write Output ========================================================


    allocate(axis2(nt))

    do t=1, nt
      axis2(t)=float(t)
    end do  

    dim1=nz
    dim2=nt
    varname='rho'
    d1name='height_m'
    d2name='time'
    title='rho_profile_domain-mean'

    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'rho_profile_1-29.nc'
    call out2d(ncfilename, varname, rhomean, &
               d1name, dim1, z, &
               d2name, dim2, axis2, &
               title)


!-- PSD

    dim1=nz
    ncfilename='rho_ave_profile_1-29.nc'
    varname='rho'
    d1name='height_m'
    title='rho_profile_time-domain-mean'

    call out1d(ncfilename, varname, rho_tmean, &
               d1name, dim1, z, &
               title)

    allocate(x(nx))
    allocate(y(ny))
    do i=1, nx
      x(i)=float(i-1)*27.
    end do
    do j=1, ny
      y(j)=float(j-1)*27.
    end do

    write(ncfilename,'(a,i2.2,a,i2.2,a)') 'rho4d_1-29.nc'
    call out4d(ncfilename, 'rho', rho4d, &
               'x', nx, x, &
               'y', ny, y, &
               'height_m', nz, z, &
               'time', nt, axis2, &
               'rho')

    deallocate(axis2) 

    end program  

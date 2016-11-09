    program  power_spectra 

    use netcdfio

    implicit none

    integer, parameter ::  cf = 1

    integer :: nx, ny, nt
    integer :: i, j, k, t, l1, l2, kk
    integer :: ncid
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: nz=85

    integer, parameter :: nc=31, max_lz=30

    real, parameter :: dz=500.
    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants
    real, parameter :: r_d=287.

    real, parameter :: dc=1.0  ! interval of phase speed
 
    real :: f0, lat, absomeh
    real :: lz
    real :: dx, dy, dt

    real, dimension(nz) :: z
    real, dimension(nz) :: uave, vave, nave
    real, dimension(nc) :: lzgrd
    real, dimension(nc,nz) :: mfx2, mfy2

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:,:,:,:) :: mfx, mfy, m2
 
    character(len=100) :: ncfilename, ifdir
    character(len=50) :: varname
    character(len=2)  ::  speriod

    call getarg(1,speriod)

    nx = 43  ;  ny = 85  ;  nt = 31
    if (speriod == 'p1')  lat = 20.
    if (speriod == 'p2')  lat = 27.
    if (speriod == 'p3') then
      nx = 44  ;  ny = 87  ;  nt = 31
      lat = 34.
    end if

    allocate( freq_x(nx), freq_y(ny), freq_t(nt) )

!== Calculate average f0 at this domain

    f0=2.*omega*sin(lat*pi/180.) ! f at center of domain (24.8N)
    print *, 'f0 at this domain', f0

!== READ 3D MF SPECTRUM =================================================
 
    allocate(mfx(nx,ny,nt,nz))

    write(ncfilename,'(a,i0.0,a)') 'res/wind/mf_k_l_f/mfx_k_l_f_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfx'
    call get4d(ncid, varname, nx, ny, nt, nz, mfx)
    call get1d(ncid, 'k', nx, freq_x)
    call get1d(ncid, 'l', ny, freq_y)
    call get1d(ncid, 'f', nt, freq_t)
    call get1d(ncid, 'height_m', nz, z)
    call closenc(ncid)

    allocate(mfy(nx,ny,nt,nz))

    write(ncfilename,'(a,i0.0,a)') 'res/wind/mf_k_l_f/mfy_k_l_f_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfy'
    call get4d(ncid, varname, nx, ny, nt, nz, mfy)
    call closenc(ncid)

    allocate(m2(nx,ny,nt,nz))
    call opennc('res/m2_'//speriod//'.nc', ncid)
    call get4d(ncid, 'm2', nx, ny, nt, nz, m2)
    call closenc(ncid)

!== Divide cospectra to bins (phase speed, propagation direction)

    dx=freq_x(2)-freq_x(1)
    dy=freq_y(2)-freq_y(1)
    dt=freq_t(2)-freq_t(1)

    print *, 'dk=',dx
    print *, 'dl=',dy
    print *, 'df=',dt
  
    do l1=1, nc
      lzgrd(l1)=(l1-1)*dc
    end do
 
    mfx2(:,:)=0.
    mfy2(:,:)=0.

    do k=1, nz
      do t=1, nt
      do j=1, ny
      do i=1, nx
        if ( i == 1 .and. j == ny/2+1 )  CYCLE
        if ( m2(i,j,t,k) < 0.)  CYCLE
        lz=2.*pi/sqrt(m2(i,j,t,k))/1.e3
        do l1=1, nc
          if(lz.gt.lzgrd(l1)-dc/2. .and. lz.le.lzgrd(l1)+dc/2.) then
            mfx2(l1,k)=mfx2(l1,k)+mfx(i,j,t,k)*dx*dy*dt/dc
            mfy2(l1,k)=mfy2(l1,k)+mfy(i,j,t,k)*dx*dy*dt/dc
          end if
        end do
      end do
      end do
      end do
    end do 

    deallocate(mfx)
    deallocate(mfy)

!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/wind/mf_lz/mfx_lz_'//speriod//'_f',cf,'.nc'
    call out2d_yh(ncfilename, 1, (/'mfx'/), mfx2, &
               'lz', nc, lzgrd, &
               'height', nz, z, &
               'mfx')

    write(ncfilename,'(a,i0.0,a)') 'res/wind/mf_lz/mfy_lz_'//speriod//'_f',cf,'.nc'
    call out2d_yh(ncfilename, 3, (/'mfy'/), mfy2, &
               'lz', nc, lzgrd, &
               'height', nz, z, &
               'mfy')

    end program  

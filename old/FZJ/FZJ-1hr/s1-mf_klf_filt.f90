    program  power_spectra 

    use netcdfio

    implicit none

    integer, parameter ::  cf = 1

    integer :: nx, ny, nt
    integer :: i, j, k, t, l1, l2, kk
    integer :: ncid
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: nz=85

    integer, parameter :: nc=26, max_cp=100
    integer, parameter :: nphi=36

    real, parameter :: dz=500.
    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants
    real, parameter :: r_d=287.

    real, parameter :: dc=4.  ! interval of phase speed
    real, parameter :: dphi=10.  ! deg, interval of propagation direction (phi)
 
    real :: f0, lat, omeh2
    real :: cp, phi
    real :: dx, dy, dt
    real :: lo_phi, up_phi

    real, dimension(nz) :: z
    real, dimension(nz) :: uave, vave, nave
    real, dimension(nc) :: cpgrd
    real, dimension(nphi) :: phigrd
    real, dimension(nphi,nc,nz) :: mfx2, mfy2

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:,:,:,:) :: mfx, mfy
 
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

!== READ Background Profile

    ifdir='/data4/atmosdata/T_SPECTRUM/'//speriod

    ncfilename=trim(ifdir)//'/u_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'u_ave', nz, uave)
    call closenc(ncid)

    ncfilename=trim(ifdir)//'/v_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'v_ave', nz, vave)
    call closenc(ncid)

    ncfilename=trim(ifdir)//'/n_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'n_ave', nz, nave)
    call closenc(ncid)


!== READ 3D MF SPECTRUM =================================================
 
    allocate(mfx(nx,ny,nt,nz)) 

    ncfilename='res/wind/mf_k_l_f/mfx_k_l_f_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    varname='mfx'
    call get4d(ncid, varname, nx, ny, nt, nz, mfx)
    call get1d(ncid, 'k', nx, freq_x)
    call get1d(ncid, 'l', ny, freq_y)
    call get1d(ncid, 'f', nt, freq_t)
    call get1d(ncid, 'height_m', nz, z)
    call closenc(ncid)

    allocate(mfy(nx,ny,nt,nz)) 

    ncfilename='res/wind/mf_k_l_f/mfy_k_l_f_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    varname='mfy'
    call get4d(ncid, varname, nx, ny, nt, nz, mfy)
    call closenc(ncid)


!== Filter out the components that don't satisfy vertical propagation condition of IGWs


    do t=1, nt
    do j=1, ny
    do i=1, nx
    do k=1, nz
      omeh2 = (2.*pi*( -freq_t(t)-uave(k)*freq_x(i)-vave(k)*freq_y(j) ))**2
      if( omeh2 <= f0*f0*float(cf) .or. omeh2 >= nave(k)*nave(k)/float(cf) ) then
        mfx(i,j,t,k)=0.
        mfy(i,j,t,k)=0.
      end if
    enddo
    enddo
    enddo
    enddo
    mfx(1,ny/2+1,:,:) = 0.
    mfy(1,ny/2+1,:,:) = 0.

!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/wind/mf_k_l_f/mfx_k_l_f_'//speriod//'_f',cf,'.nc'
    varname='mfx'

    call out4d(ncfilename, varname, mfx, &
               'k', nx, freq_x, &
               'l', ny, freq_y, &
               'f', nt, freq_t, &
               'height_m', nz, z,'mfx_filt')

    write(ncfilename,'(a,i0.0,a)') 'res/wind/mf_k_l_f/mfy_k_l_f_'//speriod//'_f',cf,'.nc'
    varname='mfy'

    call out4d(ncfilename, varname, mfy, &
               'k', nx, freq_x, &
               'l', ny, freq_y, &
               'f', nt, freq_t, &
               'height_m', nz, z,'mfy_filt')

    end program  

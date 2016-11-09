    program  power_spectra 

    use netcdfio
    use regress

    implicit none

    integer :: nx, ny, nt
    integer :: i, j, k, t, l1, l2, kk
    integer :: ncid
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: nz=85

    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants

    real :: f0, lat
    real :: dx, dy, dt

    real, dimension(nz) :: z
    real, dimension(nz) :: uave, vave, nave

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:)       :: kwn, lwn, ome, ome_h, kvec, temp
    real, allocatable, dimension(:,:,:,:) :: m2
 
    character(len=100) :: ncfilename, ifdir
    character(len=50)  :: varname
    character(len=2)   :: speriod

    call getarg(1,speriod)

    nx = 43  ;  ny = 85  ;  nt = 175
    if (speriod == 'p1')  lat = 20.
    if (speriod == 'p2')  lat = 27.
    if (speriod == 'p3') then
      nx = 44  ;  ny = 87  ;  nt = 175
      lat = 34.
    end if

    allocate( freq_x(nx), freq_y(ny), freq_t(nt) )
    allocate( kwn   (nx), lwn   (ny), ome   (nt) )

    allocate( ome_h(nx), kvec(nx) )

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
 
    ncfilename='res/psd_tp_k_l_f_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid, 'k', nx, freq_x)
    call get1d(ncid, 'l', ny, freq_y)
    call get1d(ncid, 'f', nt, freq_t)
    call get1d(ncid, 'height_m', nz, z)
    call closenc(ncid)


    kwn(:) = freq_x(:)*2.*pi
    lwn(:) = freq_y(:)*2.*pi
    ome(:) = freq_t(:)*2.*pi * (-1.)

    allocate(m2(nx,ny,nt,nz))
    do k=1, nz
    do t=1, nt
    do j=1, ny
      ome_h(:) = ome(t)-uave(k)*kwn(:)-vave(k)*lwn(j)
      kvec (:) = sqrt( kwn(:)*kwn(:) + lwn(j)*lwn(j) )
      m2(:,j,t,k) = kvec(:)*kvec(:)*(nave(k)*nave(k)-ome_h(:)*ome_h(:)) &
                    / (ome_h(:)*ome_h(:)) - 5.e-9
    enddo
    enddo
    enddo
    m2(1,ny/2+1,:,:) = 1.e32

print*, minval(m2(:,:,:,35:))

!== Write Output ========================================================

    ncfilename='res/m2_'//speriod//'.nc'

    call out4d_yh(ncfilename, 1, (/'m2  '/), m2, &
               'k', nx, freq_x, &
               'l', ny, freq_y, &
               'f', nt, freq_t, &
               'height_m', nz, z, &
               'm2')

    end program  

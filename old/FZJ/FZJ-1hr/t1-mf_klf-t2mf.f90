    program  power_spectra 

    use netcdfio
    use regress

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
    real, parameter :: g=9.81
    real, parameter :: r_d=287.

    real, parameter :: dc=4.  ! interval of phase speed
    real, parameter :: dphi=10.  ! deg, interval of propagation direction (phi)
 
    real :: cp, phi, omeh2, f0, lat
    real :: dx, dy, dt
    real :: lo_phi, up_phi

    real, dimension(nz) :: z
    real, dimension(nz) :: uave, vave, nave, rave, tave
    real, dimension(nc) :: cpgrd
    real, dimension(nphi) :: phigrd

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:)       :: kwn, lwn, ome, ome_h, kvec, temp
    real, allocatable, dimension(:,:,:,:) :: mfx, mfy, tup, tdn, tupx, tdnx, tupy, tdny, psd, tb
 
    character(len=100) :: ncfilename, ifdir
    character(len=50) ::  varname
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

    ncfilename=trim(ifdir)//'/rho_ave_profile_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'rho', nz, rave)
    call closenc(ncid)

    ncfilename='../res/t_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'t_ave', nz, tave)
    call closenc(ncid)

!== READ 3D MF SPECTRUM =================================================
 
    allocate(mfx(nx,ny,nt,nz)) 
    allocate(mfy(nx,ny,nt,nz))
    allocate(tup(nx,ny,nt,nz),tdn(nx,ny,nt,nz))
    allocate(tupx(nx,ny,nt,nz),tdnx(nx,ny,nt,nz),tupy(nx,ny,nt,nz),tdny(nx,ny,nt,nz))

    ncfilename='res/psd_tp_k_l_f_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get4d(ncid, 't2up', nx, ny, nt, nz, tup)
    call get4d(ncid, 't2dn', nx, ny, nt, nz, tdn)
    call get1d(ncid, 'k', nx, freq_x)
    call get1d(ncid, 'l', ny, freq_y)
    call get1d(ncid, 'f', nt, freq_t)
    call get1d(ncid, 'height_m', nz, z)
    call closenc(ncid)

!== FILTER

    do t=1, nt
    do j=1, ny
    do i=1, nx
    do k=1, nz
      omeh2 = (2.*pi*( -freq_t(t)-uave(k)*freq_x(i)-vave(k)*freq_y(j) ))**2
      if( omeh2 <= f0*f0*float(cf) .or. omeh2 >= nave(k)*nave(k)/float(cf) ) then
        tup(i,j,t,k)=0.
        tdn(i,j,t,k)=0.
      end if
    enddo
    enddo
    enddo
    enddo
    tup(1,ny/2+1,:,:) = 0.
    tdn(1,ny/2+1,:,:) = 0.

!== Deduce horizontal MF spectrum from T spectrum =======================

    kwn(:) = freq_x(:)*2.*pi
    lwn(:) = freq_y(:)*2.*pi
    ome(:) = freq_t(:)*2.*pi * (-1.)
    allocate(temp(nx))

    do k=1, nz
    do t=1, nt
    do j=1, ny
      ome_h(:) = ome(t)-uave(k)*kwn(:)-vave(k)*lwn(j)
      kvec (:) = sqrt( kwn(:)*kwn(:) + lwn(j)*lwn(j) )
!!!
      temp(:) = rave(k)*g*g/nave(k)**3 * ome_h(:)/kvec(:)    &
                     /(tave(k)*tave(k))
!!!
      tup(:,j,t,k) = tup(:,j,t,k)*temp(:)
      tdn(:,j,t,k) = -tdn(:,j,t,k)*temp(:)
    enddo
    enddo
    enddo
    tup(1,ny/2+1,:,:) = 0.
    tdn(1,ny/2+1,:,:) = 0.

    do j=1, ny
      tupy(:,j,:,:) = tup(:,j,:,:)*lwn(j)
      tdny(:,j,:,:) = tdn(:,j,:,:)*lwn(j)
    enddo
    mfy = tupy + tdny

    do i=1, nx
      tupx(i,:,:,:) = tup(i,:,:,:)*kwn(i)
      tdnx(i,:,:,:) = tdn(i,:,:,:)*kwn(i)
    enddo
    mfx = tupx + tdnx


!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_k_l_f/mfx_k_l_f_'//speriod//'_f',cf,'.nc'
    varname='mfx'

    call out4d_yh(ncfilename, 3, (/'mfx   ','mfx_up','mfx_dn'/), (/mfx,tupx,tdnx/), &
               'k', nx, freq_x, &
               'l', ny, freq_y, &
               'f', nt, freq_t, &
               'height_m', nz, z, &
               'mfx')

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_k_l_f/mfy_k_l_f_'//speriod//'_f',cf,'.nc'
    varname='mfy'

    call out4d_yh(ncfilename, 3, (/'mfy   ','mfy_up','mfy_dn'/), (/mfy,tupy,tdny/), &
               'k', nx, freq_x, &
               'l', ny, freq_y, &
               'f', nt, freq_t, &
               'height_m', nz, z, &
               'mfy')

    end program  

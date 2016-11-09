    program  power_spectra 

    use netcdfio

    implicit none

    integer, parameter ::  cf = 1

    integer :: nx, ny, nt
    integer :: i, j, k, t, l1, l2, kk
    integer :: ncid
    integer :: dim1, dim2, dim3, dim4

    integer, parameter :: nz=85

    integer, parameter :: nc=21, max_kh=10.e-06
    integer, parameter :: nphi=36

    real, parameter :: dz=500.
    real, parameter :: omega=7.2921E-5 ! same as in WRF constants
    real, parameter :: pi=3.1415926  ! same as in WRF constants
    real, parameter :: r_d=287.

    real, parameter :: dc=0.5e-06  ! interval of phase speed
    real, parameter :: dphi=10.  ! deg, interval of propagation direction (phi)
 
    real :: f0, lat, absomeh
    real :: kh, phi
    real :: dx, dy, dt
    real :: lo_phi, up_phi

    real, dimension(nz) :: z
    real, dimension(nz) :: uave, vave, nave
    real, dimension(nc) :: khgrd
    real, dimension(nphi) :: phigrd
    real, dimension(nphi,nc,nz) :: mfx2_up, mfx2_dn, mfy2_up, mfy2_dn

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:,:,:,:) :: mfx_up, mfx_dn, mfy_up, mfy_dn
 
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
 
    allocate(mfx_up(nx,ny,nt,nz), mfx_dn(nx,ny,nt,nz)) 

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_k_l_f/mfx_k_l_f_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfx_up'
    call get4d(ncid, varname, nx, ny, nt, nz, mfx_up)
    varname='mfx_dn'
    call get4d(ncid, varname, nx, ny, nt, nz, mfx_dn)
    call get1d(ncid, 'k', nx, freq_x)
    call get1d(ncid, 'l', ny, freq_y)
    call get1d(ncid, 'f', nt, freq_t)
    call get1d(ncid, 'height_m', nz, z)
    call closenc(ncid)

    allocate(mfy_up(nx,ny,nt,nz), mfy_dn(nx,ny,nt,nz))

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_k_l_f/mfy_k_l_f_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)
    varname='mfy_up'
    call get4d(ncid, varname, nx, ny, nt, nz, mfy_up)
    varname='mfy_dn'
    call get4d(ncid, varname, nx, ny, nt, nz, mfy_dn)
    call closenc(ncid)

!== Divide cospectra to bins (phase speed, propagation direction)

    dx=freq_x(2)-freq_x(1)
    dy=freq_y(2)-freq_y(1)
    dt=freq_t(2)-freq_t(1)

    print *, 'dk=',dx
    print *, 'dl=',dy
    print *, 'df=',dt
  
    do l1=1, nc
      khgrd(l1)=(l1-1)*dc
    end do
    do l2=1, nphi
      phigrd(l2)=(l2-1)*dphi
    end do
 
    mfx2_up(:,:,:)=0.  ;  mfx2_dn(:,:,:)=0.
    mfy2_up(:,:,:)=0.  ;  mfy2_dn(:,:,:)=0.

    do k=1, nz
      do t=1, nt
      do j=1, ny
      do i=1, nx
        if ( i == 1 .and. j == ny/2+1 )  CYCLE
        do l2=1, nphi
          phi=sqrt(freq_x(i)*freq_x(i)+freq_y(j)*freq_y(j))
          kh=phi
          phi=acos(freq_x(i)/phi)
          if(-1*freq_t(t).ge.0. .and. freq_y(j).ge.0.) phi=phi
          if(-1*freq_t(t).lt.0. .and. freq_y(j).lt.0.) phi=pi-phi
          if(-1*freq_t(t).lt.0. .and. freq_y(j).ge.0.) phi=phi+pi
          if(-1*freq_t(t).ge.0. .and. freq_y(j).lt.0.) phi=2.*pi-phi
          phi=phi*180./pi
          lo_phi=phigrd(l2)-dphi/2.
          up_phi=phigrd(l2)+dphi/2.
          if( ( (lo_phi < 0.) .and. (phi <= dphi/2. .or. phi > lo_phi+360.) ) &
              .or. (phi > lo_phi .and. phi <= up_phi) ) then
            do l1=1, nc
              if(kh.gt.khgrd(l1)-dc/2. .and. kh.le.khgrd(l1)+dc/2.) then
                mfx2_up(l2,l1,k)=mfx2_up(l2,l1,k)+mfx_up(i,j,t,k)
                mfy2_up(l2,l1,k)=mfy2_up(l2,l1,k)+mfy_up(i,j,t,k)
                mfx2_dn(l2,l1,k)=mfx2_dn(l2,l1,k)+mfx_dn(i,j,t,k)
                mfy2_dn(l2,l1,k)=mfy2_dn(l2,l1,k)+mfy_dn(i,j,t,k)
              end if
            end do
          end if
        end do

      end do
      end do
      end do
    end do 

    mfx2_up(:,:,:) = mfx2_up(:,:,:)*dx*dy*dt/dc/dphi
    mfy2_up(:,:,:) = mfy2_up(:,:,:)*dx*dy*dt/dc/dphi
    mfx2_dn(:,:,:) = mfx2_dn(:,:,:)*dx*dy*dt/dc/dphi
    mfy2_dn(:,:,:) = mfy2_dn(:,:,:)*dx*dy*dt/dc/dphi

    deallocate(mfx_up, mfx_dn)
    deallocate(mfy_up, mfy_dn)

!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_kh_phi/mfx_kh_phi_'//speriod//'_f',cf,'.nc'
    call out3d_yh(ncfilename, 3, (/'mfx   ','mfx_up','mfx_dn'/), &
               (/mfx2_up+mfx2_dn,mfx2_up,mfx2_dn/), &
               'direction', nphi, phigrd, &
               'kh', nc, khgrd, &
               'height', nz, z, &
               'mfx')

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_kh_phi/mfy_kh_phi_'//speriod//'_f',cf,'.nc'
    call out3d_yh(ncfilename, 3, (/'mfy   ','mfy_up','mfy_dn'/), &
               (/mfy2_up+mfy2_dn,mfy2_up,mfy2_dn/), &
               'direction', nphi, phigrd, &
               'kh', nc, khgrd, &
               'height', nz, z, &
               'mfy')

    end program  

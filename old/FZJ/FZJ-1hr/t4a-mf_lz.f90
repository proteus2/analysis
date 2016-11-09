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
    real :: lz, omeh
    real :: dx, dy, dt, delta

    real, dimension(nz) :: z
    real, dimension(nz) :: uave, vave, nave
    real, dimension(nc) :: lzgrd
    real, dimension(nc,nz) :: mfx2, mfy2, mfx2_up, mfx2_dn, mfy2_up, mfy2_dn
    real, dimension(nc,nz) :: mfx2_e, mfy2_n, mfx2_e_up, mfx2_e_dn, mfy2_n_up, mfy2_n_dn
    real, dimension(nc,nz) :: mfx2_w, mfy2_s, mfx2_w_up, mfx2_w_dn, mfy2_s_up, mfy2_s_dn

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:,:,:,:) :: mfx_up, mfx_dn, mfy_up, mfy_dn, m2
 
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
 
    mfx2_up(:,:)=0.  ;  mfx2_dn(:,:)=0.
    mfy2_up(:,:)=0.  ;  mfy2_dn(:,:)=0.
    mfx2_e_up(:,:)=0.  ;  mfx2_e_dn(:,:)=0.
    mfy2_n_up(:,:)=0.  ;  mfy2_n_dn(:,:)=0.
    mfx2_w_up(:,:)=0.  ;  mfx2_w_dn(:,:)=0.
    mfy2_s_up(:,:)=0.  ;  mfy2_s_dn(:,:)=0.

    do k=1, nz
      do t=1, nt
      do j=1, ny
      do i=1, nx
        if ( i == 1 .and. j == ny/2+1 )  CYCLE
        if ( m2(i,j,t,k) < 0.)  CYCLE
        lz=2.*pi/sqrt(m2(i,j,t,k))/1.e3
        omeh = 2.*pi*(-freq_t(t)-uave(k)*freq_x(i)-vave(k)*freq_y(j))
        do l1=1, nc
          if(lz.gt.lzgrd(l1)-dc/2. .and. lz.le.lzgrd(l1)+dc/2.) then

            mfx2_up(l1,k)=mfx2_up(l1,k)+mfx_up(i,j,t,k)
            mfy2_up(l1,k)=mfy2_up(l1,k)+mfy_up(i,j,t,k)
            mfx2_dn(l1,k)=mfx2_dn(l1,k)+mfx_dn(i,j,t,k)
            mfy2_dn(l1,k)=mfy2_dn(l1,k)+mfy_dn(i,j,t,k)

            if ( omeh*freq_x(i) > 0. ) then
              mfx2_e_up(l1,k) = mfx2_e_up(l1,k) + mfx_up(i,j,t,k)
              mfx2_e_dn(l1,k) = mfx2_e_dn(l1,k) + mfx_dn(i,j,t,k)
            end if
            if ( omeh*freq_x(i) < 0. ) then
              mfx2_w_up(l1,k) = mfx2_w_up(l1,k) + mfx_up(i,j,t,k)
              mfx2_w_dn(l1,k) = mfx2_w_dn(l1,k) + mfx_dn(i,j,t,k)
            end if
            if ( omeh*freq_y(i) > 0. ) then
              mfy2_n_up(l1,k) = mfy2_n_up(l1,k) + mfy_up(i,j,t,k)
              mfy2_n_dn(l1,k) = mfy2_n_dn(l1,k) + mfy_dn(i,j,t,k)
            end if
            if ( omeh*freq_y(i) < 0. ) then
              mfy2_s_up(l1,k) = mfy2_s_up(l1,k) + mfy_up(i,j,t,k)
              mfy2_s_dn(l1,k) = mfy2_s_dn(l1,k) + mfy_dn(i,j,t,k)
            end if
 
          end if
        end do
 
      end do
      end do
      end do
    end do 

    deallocate(mfx_up, mfx_dn)
    deallocate(mfy_up, mfy_dn)

    delta = dx*dy*dt/dc

    mfx2_up  =mfx2_up  *delta
    mfy2_up  =mfy2_up  *delta
    mfx2_dn  =mfx2_dn  *delta
    mfy2_dn  =mfy2_dn  *delta
    mfx2_e_up=mfx2_e_up*delta
    mfy2_n_up=mfy2_n_up*delta
    mfx2_e_dn=mfx2_e_dn*delta
    mfy2_n_dn=mfy2_n_dn*delta
    mfx2_w_up=mfx2_w_up*delta
    mfy2_s_up=mfy2_s_up*delta
    mfx2_w_dn=mfx2_w_dn*delta
    mfy2_s_dn=mfy2_s_dn*delta


    mfx2   = mfx2_up   + mfx2_dn
    mfy2   = mfy2_up   + mfy2_dn
    mfx2_e = mfx2_e_up + mfx2_e_dn
    mfx2_w = mfx2_w_up + mfx2_w_dn
    mfy2_n = mfy2_n_up + mfy2_n_dn
    mfy2_s = mfy2_s_up + mfy2_s_dn

!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_lz/mfx_lz_'//speriod//'_f',cf,'.nc'
    call out2d_yh(ncfilename, 9, &
               (/'mfx     ','mfx_e   ','mfx_w   ',   &
                 'mfx_up  ','mfx_up_e','mfx_up_w',   &
                 'mfx_dn  ','mfx_dn_e','mfx_dn_w'/), &
               (/mfx2   ,mfx2_e   ,mfx2_w   ,   &
                 mfx2_up,mfx2_e_up,mfx2_w_up,   &
                 mfx2_dn,mfx2_e_dn,mfx2_w_dn/), &
               'lz', nc, lzgrd, &
               'height', nz, z, &
               'mfx')

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_lz/mfy_lz_'//speriod//'_f',cf,'.nc'
    call out2d_yh(ncfilename, 9, &
               (/'mfy     ','mfy_n   ','mfy_s   ',   &
                 'mfy_up  ','mfy_up_n','mfy_up_s',   &
                 'mfy_dn  ','mfy_dn_n','mfy_dn_s'/), &
               (/mfy2   ,mfy2_n   ,mfy2_s   ,   &
                 mfy2_up,mfy2_n_up,mfy2_s_up,   &
                 mfy2_dn,mfy2_n_dn,mfy2_s_dn/), &
               'lz', nc, lzgrd, &
               'height', nz, z, &
               'mfy')

    end program  

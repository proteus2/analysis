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
 
    real :: cp, phi, omeh
    real :: dx, dy, dt, delta
    real :: lo_phi, up_phi

    real, dimension(nz) :: z, uave, vave
    real, dimension(nz) :: mfxa_p, mfxa_m, mfya_p, mfya_m, mfxa, mfya
    real, dimension(nz) :: mfxa_e, mfxa_w, mfya_n, mfya_s

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:,:,:,:) :: mfx, mfy
 
    character(len=100) :: ncfilename, ifdir
    character(len=50) :: varname
    character(len=2)  ::  speriod

    call getarg(1,speriod)

    nx = 43  ;  ny = 85  ;  nt = 175
    if (speriod == 'p3') then
      nx = 44  ;  ny = 87  ;  nt = 175
    end if

    allocate( freq_x(nx), freq_y(ny), freq_t(nt) )

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

    ifdir='/data4/atmosdata/T_SPECTRUM/'//speriod

    ncfilename=trim(ifdir)//'/u_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'u_ave', nz, uave)
    call closenc(ncid)
    
    ncfilename=trim(ifdir)//'/v_profile_ave_'//speriod//'.nc'
    call opennc(ncfilename, ncid)
    call get1d(ncid,'v_ave', nz, vave)
    call closenc(ncid)

!== Divide cospectra to bins (phase speed, propagation direction)

    mfxa_p(:) = 0.  ;  mfxa_m(:) = 0.
    mfya_p(:) = 0.  ;  mfya_m(:) = 0.

    mfxa_e(:) = 0.  ;  mfxa_w(:) = 0.
    mfya_n(:) = 0.  ;  mfya_s(:) = 0.

    do k=1, nz
    do t=1, nt
    do j=1, ny
    do i=1, nx

      if ( mfx(i,j,t,k) > 0. ) then
        mfxa_p(k) = mfxa_p(k) + mfx(i,j,t,k)
      else
        mfxa_m(k) = mfxa_m(k) + mfx(i,j,t,k)
      end if
      if ( mfy(i,j,t,k) > 0. ) then
        mfya_p(k) = mfya_p(k) + mfy(i,j,t,k)
      else
        mfya_m(k) = mfya_m(k) + mfy(i,j,t,k)
      end if

      omeh = 2.*pi*(-freq_t(t)-uave(k)*freq_x(i)-vave(k)*freq_y(j))
      if ( omeh*freq_x(i) > 0. ) then
        mfxa_e(k) = mfxa_e(k) + mfx(i,j,t,k)
      end if
      if ( omeh*freq_x(i) < 0. ) then
        mfxa_w(k) = mfxa_w(k) + mfx(i,j,t,k)
      end if
      if ( omeh*freq_y(i) > 0. ) then
        mfya_n(k) = mfya_n(k) + mfy(i,j,t,k)
      end if
      if ( omeh*freq_y(i) < 0. ) then
        mfya_s(k) = mfya_s(k) + mfy(i,j,t,k)
      end if

    enddo
    enddo
    enddo
    enddo

    delta = (freq_x(2)-freq_x(1))*(freq_y(2)-freq_y(1))* &
            (freq_t(2)-freq_t(1))

    mfxa_p(:) = mfxa_p(:) * delta
    mfxa_m(:) = mfxa_m(:) * delta
    mfya_p(:) = mfya_p(:) * delta
    mfya_m(:) = mfya_m(:) * delta

    mfxa_e(:) = mfxa_e(:) * delta
    mfxa_w(:) = mfxa_w(:) * delta
    mfya_n(:) = mfya_n(:) * delta
    mfya_s(:) = mfya_s(:) * delta

    mfxa(:) = mfxa_p(:) + mfxa_m(:)
    mfya(:) = mfya_p(:) + mfya_m(:)

!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/wind/mfa/mfa-klf_'//speriod//'_f',cf,'.nc'

    call out1d(ncfilename, 10,                                      &
               (/'mfx     ','mfx_p   ','mfx_m   ','mfx_e   ','mfx_w   ',   &
                 'mfy     ','mfy_p   ','mfy_m   ','mfy_n   ','mfy_s   '/), &
               (/mfxa   ,mfxa_p   ,mfxa_m   ,mfxa_e   ,mfxa_w   ,   &
                 mfya   ,mfya_p   ,mfya_m   ,mfya_n   ,mfya_s   /), &
               'height', nz, z,                                     &
               'mf profile')


    end program  

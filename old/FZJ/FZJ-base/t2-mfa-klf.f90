    program  power_spectra 

    use netcdfio

    implicit none

    integer, parameter ::  cf = 10

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
    real, dimension(nz) :: mfxa_p_up, mfxa_m_up, mfya_p_up, mfya_m_up, mfxa_up, mfya_up
    real, dimension(nz) :: mfxa_p_dn, mfxa_m_dn, mfya_p_dn, mfya_m_dn, mfxa_dn, mfya_dn
    real, dimension(nz) :: mfxa_e_up, mfxa_w_up, mfya_n_up, mfya_s_up
    real, dimension(nz) :: mfxa_e_dn, mfxa_w_dn, mfya_n_dn, mfya_s_dn
    real, dimension(nz) :: mfxa_e, mfxa_w, mfya_n, mfya_s

    real, allocatable, dimension(:)       :: freq_x, freq_y, freq_t
    real, allocatable, dimension(:,:,:,:) :: mfx, mfy, mfx_up, mfy_up, mfx_dn, mfy_dn
 
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
    allocate(mfx_up(nx,ny,nt,nz))
    allocate(mfx_dn(nx,ny,nt,nz))

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_k_l_f/mfx_k_l_f_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)

    varname='mfx'
    call get4d(ncid, varname, nx, ny, nt, nz, mfx)
    varname='mfx_up'
    call get4d(ncid, varname, nx, ny, nt, nz, mfx_up)
    varname='mfx_dn'
    call get4d(ncid, varname, nx, ny, nt, nz, mfx_dn)

    call get1d(ncid, 'k', nx, freq_x)
    call get1d(ncid, 'l', ny, freq_y)
    call get1d(ncid, 'f', nt, freq_t)
    call get1d(ncid, 'height_m', nz, z)
    call closenc(ncid)

    allocate(mfy(nx,ny,nt,nz)) 
    allocate(mfy_up(nx,ny,nt,nz))
    allocate(mfy_dn(nx,ny,nt,nz))

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mf_k_l_f/mfy_k_l_f_'//speriod//'_f',cf,'.nc'
    call opennc(ncfilename, ncid)

    varname='mfy'
    call get4d(ncid, varname, nx, ny, nt, nz, mfy)
    varname='mfy_up'
    call get4d(ncid, varname, nx, ny, nt, nz, mfy_up)
    varname='mfy_dn'
    call get4d(ncid, varname, nx, ny, nt, nz, mfy_dn)

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

    mfxa_p   (:) = 0.  ;  mfxa_m   (:) = 0.
    mfya_p   (:) = 0.  ;  mfya_m   (:) = 0.
    mfxa_p_up(:) = 0.  ;  mfxa_m_up(:) = 0.
    mfya_p_up(:) = 0.  ;  mfya_m_up(:) = 0.
    mfxa_p_dn(:) = 0.  ;  mfxa_m_dn(:) = 0.
    mfya_p_dn(:) = 0.  ;  mfya_m_dn(:) = 0.

    mfxa_e   (:) = 0.  ;  mfxa_w   (:) = 0.
    mfya_n   (:) = 0.  ;  mfya_s   (:) = 0.
    mfxa_e_up(:) = 0.  ;  mfxa_w_up(:) = 0.
    mfya_n_up(:) = 0.  ;  mfya_s_up(:) = 0.
    mfxa_e_dn(:) = 0.  ;  mfxa_w_dn(:) = 0.
    mfya_n_dn(:) = 0.  ;  mfya_s_dn(:) = 0.

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
      if ( mfx_up(i,j,t,k) > 0. ) then
        mfxa_p_up(k) = mfxa_p_up(k) + mfx_up(i,j,t,k)
      else
        mfxa_m_up(k) = mfxa_m_up(k) + mfx_up(i,j,t,k)
      end if
      if ( mfy_up(i,j,t,k) > 0. ) then
        mfya_p_up(k) = mfya_p_up(k) + mfy_up(i,j,t,k)
      else
        mfya_m_up(k) = mfya_m_up(k) + mfy_up(i,j,t,k)
      end if
      if ( mfx_dn(i,j,t,k) > 0. ) then
        mfxa_p_dn(k) = mfxa_p_dn(k) + mfx_dn(i,j,t,k)
      else
        mfxa_m_dn(k) = mfxa_m_dn(k) + mfx_dn(i,j,t,k)
      end if
      if ( mfy_dn(i,j,t,k) > 0. ) then
        mfya_p_dn(k) = mfya_p_dn(k) + mfy_dn(i,j,t,k)
      else
        mfya_m_dn(k) = mfya_m_dn(k) + mfy_dn(i,j,t,k)
      end if

      omeh = 2.*pi*(-freq_t(t)-uave(k)*freq_x(i)-vave(k)*freq_y(j))
      if ( omeh*freq_x(i) > 0. ) then
        mfxa_e   (k) = mfxa_e   (k) + mfx   (i,j,t,k)
        mfxa_e_up(k) = mfxa_e_up(k) + mfx_up(i,j,t,k)
        mfxa_e_dn(k) = mfxa_e_dn(k) + mfx_dn(i,j,t,k)
      end if
      if ( omeh*freq_x(i) < 0. ) then
        mfxa_w   (k) = mfxa_w   (k) + mfx   (i,j,t,k)
        mfxa_w_up(k) = mfxa_w_up(k) + mfx_up(i,j,t,k)
        mfxa_w_dn(k) = mfxa_w_dn(k) + mfx_dn(i,j,t,k)
      end if
      if ( omeh*freq_y(i) > 0. ) then
        mfya_n   (k) = mfya_n   (k) + mfy   (i,j,t,k)
        mfya_n_up(k) = mfya_n_up(k) + mfy_up(i,j,t,k)
        mfya_n_dn(k) = mfya_n_dn(k) + mfy_dn(i,j,t,k)
      end if
      if ( omeh*freq_y(i) < 0. ) then
        mfya_s   (k) = mfya_s   (k) + mfy   (i,j,t,k)
        mfya_s_up(k) = mfya_s_up(k) + mfy_up(i,j,t,k)
        mfya_s_dn(k) = mfya_s_dn(k) + mfy_dn(i,j,t,k)
      end if

    enddo
    enddo
    enddo
    enddo

    delta = (freq_x(2)-freq_x(1))*(freq_y(2)-freq_y(1))* &
            (freq_t(2)-freq_t(1))

    mfxa_p   (:) = mfxa_p   (:) * delta
    mfxa_m   (:) = mfxa_m   (:) * delta
    mfya_p   (:) = mfya_p   (:) * delta
    mfya_m   (:) = mfya_m   (:) * delta
    mfxa_p_up(:) = mfxa_p_up(:) * delta
    mfxa_m_up(:) = mfxa_m_up(:) * delta
    mfya_p_up(:) = mfya_p_up(:) * delta
    mfya_m_up(:) = mfya_m_up(:) * delta
    mfxa_p_dn(:) = mfxa_p_dn(:) * delta
    mfxa_m_dn(:) = mfxa_m_dn(:) * delta
    mfya_p_dn(:) = mfya_p_dn(:) * delta
    mfya_m_dn(:) = mfya_m_dn(:) * delta

    mfxa_e   (:) = mfxa_e   (:) * delta
    mfxa_w   (:) = mfxa_w   (:) * delta
    mfya_n   (:) = mfya_n   (:) * delta
    mfya_s   (:) = mfya_s   (:) * delta
    mfxa_e_up(:) = mfxa_e_up(:) * delta
    mfxa_w_up(:) = mfxa_w_up(:) * delta
    mfya_n_up(:) = mfya_n_up(:) * delta
    mfya_s_up(:) = mfya_s_up(:) * delta
    mfxa_e_dn(:) = mfxa_e_dn(:) * delta
    mfxa_w_dn(:) = mfxa_w_dn(:) * delta
    mfya_n_dn(:) = mfya_n_dn(:) * delta
    mfya_s_dn(:) = mfya_s_dn(:) * delta

    mfxa   (:) = mfxa_p   (:) + mfxa_m   (:)
    mfya   (:) = mfya_p   (:) + mfya_m   (:)
    mfxa_up(:) = mfxa_p_up(:) + mfxa_m_up(:)
    mfya_up(:) = mfya_p_up(:) + mfya_m_up(:)
    mfxa_dn(:) = mfxa_p_dn(:) + mfxa_m_dn(:)
    mfya_dn(:) = mfya_p_dn(:) + mfya_m_dn(:)

!== Write Output ========================================================

    write(ncfilename,'(a,i0.0,a)') 'res/temp2/mfa/mfa-klf_'//speriod//'_f',cf,'.nc'

    call out1d(ncfilename, 30,                                      &
               (/'mfx     ','mfx_p   ','mfx_m   ','mfx_e   ','mfx_w   ',   &
                 'mfx_up  ','mfx_up_p','mfx_up_m','mfx_up_e','mfx_up_w',   &
                 'mfx_dn  ','mfx_dn_p','mfx_dn_m','mfx_dn_e','mfx_dn_w',   &
                 'mfy     ','mfy_p   ','mfy_m   ','mfy_n   ','mfy_s   ',   &
                 'mfy_up  ','mfy_up_p','mfy_up_m','mfy_up_n','mfy_up_s',   &
                 'mfy_dn  ','mfy_dn_p','mfy_dn_m','mfy_dn_n','mfy_dn_s'/), &
               (/mfxa   ,mfxa_p   ,mfxa_m   ,mfxa_e   ,mfxa_w   ,   &
                 mfxa_up,mfxa_p_up,mfxa_m_up,mfxa_e_up,mfxa_w_up,   &
                 mfxa_dn,mfxa_p_dn,mfxa_m_dn,mfxa_e_dn,mfxa_w_dn,   &
                 mfya   ,mfya_p   ,mfya_m   ,mfya_n   ,mfya_s   ,   &
                 mfya_up,mfya_p_up,mfya_m_up,mfya_n_up,mfya_s_up,   &
                 mfya_dn,mfya_p_dn,mfya_m_dn,mfya_n_dn,mfya_s_dn/), &
               'height', nz, z,                                     &
               'mf profile')


    end program  

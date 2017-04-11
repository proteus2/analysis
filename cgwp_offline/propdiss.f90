SUBROUTINE propdiss(  &
    nx, ny, nz, ncol,                  &
    nphi, nc, phi_dir, c_phase, dc,     &
    igwp, jgwp, u_tlev, v_tlev, nbv_tlev, rho, rho_tlev, lat, r_tlev,    &
    kcta, mfs_ct,                         &
! diagnostics
    l_drag_u_on, l_drag_v_on, drag_u, drag_v,                            &
    l_mflx_u_on, l_mflx_v_on,                                            &
    mflx_east, mflx_west, mflx_north, mflx_south,                        &
    l_mflx_u_ctop_on, l_mflx_v_ctop_on,                                  &
    mflx_e_ctop, mflx_w_ctop, mflx_n_ctop, mflx_s_ctop,                  &
    l_spec_on, diag_spec)

!
! PURPOSE:  To calculate the convective gravity wave drag
!
! METHOD:
!
!   1) Interpolate data and gather at the GWDC columns into 2-D array.
!   2) Call <gw_ctop> to calculate the cloud-top momentum flux spectrum
!   3) Obtain momentum flux profiles (Warner and McIntyre, 1999, EPS; 
!      Song and Chun, 2006, JKMS)
!   4) Calculate convective gravity wave drag
!
! HISTORY:
!

  USE subr_common
  USE param_gwp

  implicit none

! SUBROUTINE ARGUMENTS

! spatial dimensions
 
  integer, intent(in) ::  nx, ny, nz, ncol

! spectrum parameters

  integer                , intent(in) ::  nphi, nc
  real, dimension(nphi)  , intent(in) ::  phi_dir
  real, dimension(-nc:nc), intent(in) ::  c_phase
  real                   , intent(in) ::  dc

! data arrays

  integer, dimension(ncol), intent(in) ::  igwp, jgwp, kcta
  real, dimension(ncol,nz), intent(in) ::  u_tlev, v_tlev, nbv_tlev,     &
                                           rho, rho_tlev, r_tlev
  real, dimension(ncol)   , intent(in) ::  lat

  real, dimension(-nc:nc,ncol,nphi), intent(in) ::  mfs_ct

! switches

  logical, intent(in) ::  l_drag_u_on, l_drag_v_on,                      &
                          l_mflx_u_on, l_mflx_v_on,                      &
                          l_mflx_u_ctop_on, l_mflx_v_ctop_on, l_spec_on

! output variables

  real, dimension(nx,ny,nz), intent(out) ::  drag_u, drag_v,             &
                                             mflx_east, mflx_west,       &
                                             mflx_north, mflx_south 
  real, dimension(nx,ny)   , intent(out) ::  mflx_e_ctop, mflx_w_ctop,   &
                                             mflx_n_ctop, mflx_s_ctop

  real, dimension(nx,ny,nz,-nc:nc,nphi*2), intent(out) ::  diag_spec

! LOCAL VARIABLES

! data arrays

  real                          ::  mfct_pos, mfct_neg
  real, dimension(ncol)         ::  f2
  real, dimension(nphi)         ::  phi_rad, cosphi, sinphi
  real, dimension(-nc:nc)       ::  c_intr, mfsp, mfsp_s
  real, dimension(nz)           ::  mf_pos, mf_neg
  real, dimension(ncol,nz)      ::  mf_x, mf_y, fact_s
  real, dimension(ncol,nz,nphi) ::  ub_tlev

! work arrays

  real    ::  tmp, b0, w0, ome_min, pm1, c2mp
  integer ::  tmpi, ipos, ineg
  integer ::  i,j,k,l,ic,iphi   ! loop counters

! parameters and constants

!  real, parameter ::  g = 9.80665
  real, parameter ::  two_omega = 2.*7.292116E-5
  real, parameter ::  beta_eq = 2.3e-11

  include 'c_math.inc'

  if (ncol < 1)  RETURN

!------------------------------------------------------------------
!  1-1. INTERPOLATE WINDS TO THETA-GRID AND GATHER AT GWDC ARRAY.
!------------------------------------------------------------------

  call get_wm_hg2cgwp

  mflx_east (:,:,:) = 0.
  mflx_west (:,:,:) = 0.
  mflx_north(:,:,:) = 0.
  mflx_south(:,:,:) = 0.

  diag_spec(:,:,:,:,:) = 0.

!------------------------------------------------------------------
!  1-3. PREPARE SOME VARIABLES USED IN GWDC CALCULATIONS.
!------------------------------------------------------------------

  phi_rad(:) = phi_dir(:)*deg2rad
  cosphi(:) = cos(phi_rad(:))
  sinphi(:) = sin(phi_rad(:))

  f2(:) = (two_omega*two_omega)*sin(lat(:)*deg2rad)**2

  call basic_u_phi(u_tlev,v_tlev,phi_dir, ub_tlev)

!------------------------------------------------------------------
!  3. OPTAIN MOMENTUM FLUX PROFILES
!------------------------------------------------------------------

  ! for calculating saturated spectrum
  tmp  = beta_wm/(sqrt(2.0)*pi)
  pm1  = p_wm - 1.0
  c2mp = 2.0 - p_wm
  do k=1, nz
  do l=1, ncol
    ome_min = sqrt(max(f2(l), nbv_tlev(l,k)*beta_eq/mstar_wm))
    b0 = pm1*ome_min**pm1/(1.0-(ome_min/nbv_tlev(l,k))**pm1)
!    if (p_wm == 1.0)  b0 = 1.0/log(nbv_tlev(l,k)/ome_min)
    w0 = (nbv_tlev(l,k)**c2mp-ome_min**c2mp)/c2mp
    fact_s(l,k) = abs(tmp*b0*w0*rho_tlev(l,k))
  enddo
  enddo

  mf_x(:,:) = 0.
  mf_y(:,:) = 0.


  L_PHI:  DO iphi=1, nphi
  L_COL:  DO l=1, ncol

  i = igwp(l)
  j = jgwp(l)

  mfsp(:) = mfs_ct(:,l,iphi)*dc
  c_intr(:) = c_phase(:) - ub_tlev(l,kcta(l),iphi)

  tmpi = minloc(abs(c_intr),1) - (nc+1)

  ineg = tmpi - 1
  ipos = tmpi + 1

  ! diagnostics - unfiltered cloud-top momentum flux
  ! saved at mflx_xxxx(:,:,1)

  mfct_neg = sum(mfsp(-nc:ineg))
  mfct_pos = sum(mfsp(ipos:nc))

  if ( l_mflx_u_on ) then
    if (cosphi(iphi) > 0.) then
      mflx_east(i,j,1) = mflx_east(i,j,1) + mfct_pos*cosphi(iphi)
      mflx_west(i,j,1) = mflx_west(i,j,1) + mfct_neg*cosphi(iphi)
    else
      mflx_east(i,j,1) = mflx_east(i,j,1) + mfct_neg*cosphi(iphi)
      mflx_west(i,j,1) = mflx_west(i,j,1) + mfct_pos*cosphi(iphi)
    end if
  end if
  if ( l_mflx_v_on ) then
    mflx_north(i,j,1) = mflx_north(i,j,1) + mfct_pos*sinphi(iphi)
    mflx_south(i,j,1) = mflx_south(i,j,1) + mfct_neg*sinphi(iphi)
  end if

!yhspec+ :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ! diagnostics - unfiltered cloud-top momentum flux spectrum
  ! saved at diag_spec(:,:,1,:,:)
  if ( l_spec_on ) then
    diag_spec(i,j,1,ipos:nc,iphi) = mfsp(ipos:nc)
    diag_spec(i,j,1,-ineg:nc,nphi+iphi) = -mfsp(ineg:-nc:-1)
  end if
!yhspec- :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! momentum flux profiles
  mf_neg(:) = 0.
  mf_pos(:) = 0.

  do k=kcta(l), nz

    ! calculate the saturation spectrum using Warner and
    ! McIntyre's method but as a function of phase speed,
    ! and apply the saturation condition

    ! for components having negative MF
    mfsp_s(-nc:ineg) = fact_s(l,k)*c_intr(-nc:ineg)/nbv_tlev(l,k)*dc
    do ic=-nc, ineg
      mfsp(ic) = max(mfsp_s(ic), mfsp(ic))
    enddo
    mf_neg(k) = sum(mfsp(-nc:ineg))

    ! for components having positive MF
    mfsp_s(ipos:nc) = fact_s(l,k)*c_intr(ipos:nc)/nbv_tlev(l,k)*dc
    do ic=ipos, nc
      mfsp(ic) = min(mfsp_s(ic), mfsp(ic))
    enddo
    mf_pos(k) = sum(mfsp(ipos:nc))

    ! calculate zonal and meridional MF
    tmp = mf_pos(k) + mf_neg(k)
    mf_x(l,k) = mf_x(l,k) + tmp*cosphi(iphi)
    mf_y(l,k) = mf_y(l,k) + tmp*sinphi(iphi)

!yhspec+ :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if ( l_spec_on ) then
      if (k == kcta(l)) then
        diag_spec(i,j,2,ipos:nc,iphi) = mfsp(ipos:nc)
        diag_spec(i,j,2,-ineg:nc,nphi+iphi) = -mfsp(ineg:-nc:-1)
      end if
      diag_spec(i,j,k,ipos:nc,iphi) = mfsp(ipos:nc)
      diag_spec(i,j,k,-ineg:nc,nphi+iphi) = -mfsp(ineg:-nc:-1)
    end if
!yhspec- :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! check the critical levels between k and k+1 levels,
    ! and prepare the next-level calculations
    if (k /= nz) then

      ! update c_intr for the next level
      c_intr(:) = c_phase(:) - ub_tlev(l,k+1,iphi)

      ! update mfsp, ipos, ineg
      tmpi = ineg
      do ic=tmpi, -nc, -1
        if (c_intr(ic) < 0.)  EXIT
        mfsp(ic) = 0.
        ineg = ic - 1
      enddo
      tmpi = ipos
      do ic=tmpi, nc
        if (c_intr(ic) > 0.)  EXIT
        mfsp(ic) = 0.
        ipos = ic + 1
      enddo

    end if  ! k /= nz

  enddo  ! k

  ! diagnostics
  if ( l_mflx_u_on ) then
    if (cosphi(iphi) > 0.) then
      mflx_east(i,j,3:) = mflx_east(i,j,3:) + mf_pos(3:)*cosphi(iphi)
      mflx_west(i,j,3:) = mflx_west(i,j,3:) + mf_neg(3:)*cosphi(iphi)
    else
      mflx_east(i,j,3:) = mflx_east(i,j,3:) + mf_neg(3:)*cosphi(iphi)
      mflx_west(i,j,3:) = mflx_west(i,j,3:) + mf_pos(3:)*cosphi(iphi)
    end if
  end if
  if ( l_mflx_v_on ) then
    mflx_north(i,j,3:) = mflx_north(i,j,3:) + mf_pos(3:)*sinphi(iphi)
    mflx_south(i,j,3:) = mflx_south(i,j,3:) + mf_neg(3:)*sinphi(iphi)
  end if

  if ( l_mflx_u_ctop_on .or. l_mflx_v_ctop_on ) then
    k = kcta(l)
    if ( l_mflx_u_ctop_on ) then
      if (cosphi(iphi) > 0.) then
        mflx_e_ctop(i,j) = mflx_e_ctop(i,j) + mf_pos(k)*cosphi(iphi)
        mflx_w_ctop(i,j) = mflx_w_ctop(i,j) + mf_neg(k)*cosphi(iphi)
      else
        mflx_e_ctop(i,j) = mflx_e_ctop(i,j) + mf_neg(k)*cosphi(iphi)
        mflx_w_ctop(i,j) = mflx_w_ctop(i,j) + mf_pos(k)*cosphi(iphi)
      end if
    end if
    if ( l_mflx_v_ctop_on ) then
      mflx_n_ctop(i,j) = mflx_n_ctop(i,j) + mf_pos(k)*sinphi(iphi)
      mflx_s_ctop(i,j) = mflx_s_ctop(i,j) + mf_neg(k)*sinphi(iphi)
    end if
  end if

  ENDDO  L_COL
  ENDDO  L_PHI

!------------------------------------------------------------------
!  4. CALCULATE CONVECTIVE GRAVITY WAVE DRAG
!------------------------------------------------------------------

  if ( l_drag_u_on .or. l_drag_v_on ) then

    drag_u(:,:,:) = 0.
    drag_v(:,:,:) = 0.

    do l=1, ncol
      i = igwp(l)
      j = jgwp(l)
      do k=kcta(l)+1, nz
        tmp = 1.0/rho(l,k)/(r_tlev(l,k) - r_tlev(l,k-1))
        drag_u(i,j,k) = (mf_x(l,k-1) - mf_x(l,k))*tmp
        drag_v(i,j,k) = (mf_y(l,k-1) - mf_y(l,k))*tmp
      enddo
      k = kcta(l)
      tmp = 1.0/rho(l,k)/(r_tlev(l,k) - r_tlev(l,k-1))
      drag_u(i,j,k) = (mf_x(l,nz) - mf_x(l,k))*tmp
      drag_v(i,j,k) = (mf_y(l,nz) - mf_y(l,k))*tmp
      ! as the outward flux of momentum at the model top is allowed,
      ! drags below the cloud top must be reduced for the momentum
      ! budget of the model (i.e., mf(nz)*tmp)
    enddo

  end if

  RETURN

END SUBROUTINE propdiss


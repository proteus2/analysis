MODULE prop_diss

  implicit none

  integer, parameter ::  kfmin = 3

CONTAINS


SUBROUTINE propdiss(  &
    nx, ny, nz, ncol,                  &
    nphi, nc, phi_dir, c_phase, dc,     &
    igwp, jgwp, u_tlev, v_tlev, nbv_tlev, rho_tlev, lat,                 &
    kcta, mfs_ct,                         &
! diagnostics
    l_mflx_u_on, l_mflx_v_on,                                            &
    mflx_east, mflx_west, mflx_north, mflx_south,                        &
    l_mflx_u_ctop_on, l_mflx_v_ctop_on,                                  &
    mflx_e_ctop, mflx_w_ctop, mflx_n_ctop, mflx_s_ctop,                  &
    l_spec_on, diag_spec )

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

  integer, intent(in) ::  nx, ny, nz, ncol

  ! spectrum parameters
  integer                , intent(in) ::  nphi, nc
  real, dimension(nphi)  , intent(in) ::  phi_dir
  real, dimension(-nc:nc), intent(in) ::  c_phase
  real                   , intent(in) ::  dc

  ! data arrays
  integer, dimension(ncol), intent(in) ::  igwp, jgwp, kcta
  real, dimension(ncol,nz), intent(in) ::  u_tlev, v_tlev, nbv_tlev,     &
                                           rho_tlev
  real, dimension(ncol)   , intent(in) ::  lat
  real, dimension(-nc:nc,ncol,nphi), intent(in) ::  mfs_ct

  ! switches
  logical, intent(in) ::  l_mflx_u_on, l_mflx_v_on,                      &
                          l_mflx_u_ctop_on, l_mflx_v_ctop_on, l_spec_on

  ! output variables
  real, dimension(ncol,nz,nphi), intent(out)  ::  mf_pos, mf_neg
  real, dimension(ncol,nz), intent(out) ::  mflx_east, mflx_west,        &
                                            mflx_north, mflx_south 
  real, dimension(nx,ny)   , intent(out) ::  mflx_e_ctop, mflx_w_ctop,   &
                                             mflx_n_ctop, mflx_s_ctop

  real, dimension(nx,ny,nz,-nc:nc,nphi*2), intent(out) ::  diag_spec

! LOCAL VARIABLES

! data arrays

  real                          ::  mfct_pos, mfct_neg
  real, dimension(ncol)         ::  f2
  real, dimension(nphi)         ::  phi_rad, cosphi, sinphi
  real, dimension(-nc:nc)       ::  c_intr
  real, dimension(ncol,nz,nphi) ::  ub_tlev

! work arrays

  real    ::  tmp, b0, w0, ome_min, pm1, c2mp
  integer ::  tmpi, ipos, ineg
  integer ::  i,j,k,l,ic,iphi   ! loop counters

! parameters and constants

  include 'c_math.inc'

  if (ncol < 1)  RETURN

  do iphi=1, nphi
  do l=1, ncol
    c_intr(:) = c_phase(:) - ub_tlev(l,kcta(l),iphi)
    tmpi = minloc(abs(c_intr),1) - (nc+1)

    ineg = tmpi - 1
    ipos = tmpi + 1

    diag_spec(l,ipos:nc,iphi) = mfs_ct(ipos:nc,l,iphi)
    diag_spec(l,-ineg:nc,nphi+iphi) = -mfs_ct(ineg:-nc:-1,l,iphi)
  enddo
  enddo

!------------------------------------------------------------------
!  1-1. INTERPOLATE WINDS TO THETA-GRID AND GATHER AT GWDC ARRAY.
!------------------------------------------------------------------

  call get_wm_hg2cgwp

  mflx_east (:,:) = 0.
  mflx_west (:,:) = 0.
  mflx_north(:,:) = 0.
  mflx_south(:,:) = 0.

  diag_spec(:,:,:,:,:) = 0.

!------------------------------------------------------------------
!  1-3. PREPARE SOME VARIABLES USED IN GWDC CALCULATIONS.
!------------------------------------------------------------------

  phi_rad(:) = phi_dir(:)*deg2rad
  cosphi(:) = cos(phi_rad(:))
  sinphi(:) = sin(phi_rad(:))

  call basic_u_phi(u_tlev,v_tlev,phi_dir, ub_tlev)

  ! momentum flux profiles
  mf_neg(:,:,:) = 0.
  mf_pos(:,:,:) = 0.


  L_PHI:  DO iphi=1, nphi
  L_COL:  DO l=1, ncol

  c_intr(:) = c_phase(:) - ub_tlev(l,kcta(l),iphi)

  tmpi = minloc(abs(c_intr),1) - (nc+1)

  ineg = tmpi - 1
  ipos = tmpi + 1

  mfct_neg = sum(mfs_ct(-nc:ineg,l,iphi))
  mfct_pos = sum(mfs_ct(ipos:nc ,l,iphi))

  ! diagnostics - unfiltered cloud-top momentum flux
  if ( l_mflx_u_on ) then
    if (cosphi(iphi) > 0.) then
      mflx_east(l) = mflx_east(l) + mfct_pos*cosphi(iphi)
      mflx_west(l) = mflx_west(l) + mfct_neg*cosphi(iphi)
    else
      mflx_east(l) = mflx_east(l) + mfct_neg*cosphi(iphi)
      mflx_west(l) = mflx_west(l) + mfct_pos*cosphi(iphi)
    end if
  end if
  if ( l_mflx_v_on ) then
    mflx_north(l) = mflx_north(l) + mfct_pos*sinphi(iphi)
    mflx_south(l) = mflx_south(l) + mfct_neg*sinphi(iphi)
  end if

  ENDDO  L_COL
  ENDDO  L_PHI

  RETURN

END SUBROUTINE propdiss

END module prop_diss


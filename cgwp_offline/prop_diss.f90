MODULE prop_diss

  implicit none

  logical ::  l_spec_o, l_spec_ctop_o

  real, dimension(:,:,:,:), allocatable ::  diag_spec_col

  real, dimension(:,:,:), allocatable ::  diag_spec_ct_col

CONTAINS


SUBROUTINE propdiss(  &
    ncol, nz,      &
    u_flev, v_flev, nbv_flev, rho_flev, lat,                 &
    kcta, mfs_ct,                         &
    mf_pos, mf_neg )
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

  USE param_gwp
  USE subr_common

  implicit none

! SUBROUTINE ARGUMENTS

  integer, intent(in) ::  ncol, nz

  ! data arrays
  integer, dimension(ncol), intent(in) ::  kcta
  real, dimension(ncol,nz), intent(in) ::  u_flev, v_flev, nbv_flev,     &
                                           rho_flev
  real, dimension(ncol)   , intent(in) ::  lat
  real, dimension(-nc:nc,ncol,nphi), intent(in) ::  mfs_ct

  ! output variables
  real, dimension(ncol,nz,nphi), intent(out)  ::  mf_pos, mf_neg

! LOCAL VARIABLES

! data arrays

  real, dimension(ncol)         ::  f2
  real, dimension(nphi)         ::  phi_rad, cosphi, sinphi
  real, dimension(-nc:nc)       ::  c_intr, mfsp, mfsp_s
  real, dimension(ncol,nz)      ::  fact_s
  real, dimension(ncol,nz,nphi) ::  ub_flev

! work arrays

  real    ::  tmp, b0, w0, ome_min, pm1, c2mp
  integer ::  tmpi, ipos, ineg
  integer ::  k,l,ic,iphi   ! loop counters

! parameters and constants

!  real, parameter ::  g = 9.80665
  real, parameter ::  two_omega = 2.*7.292116E-5
  real, parameter ::  beta_eq = 2.3e-11

  include 'c_math.inc'

  if (ncol < 1)  RETURN

  if ( l_spec_o ) then
    if ( allocated(diag_spec_col) )  deallocate( diag_spec_col )
    allocate( diag_spec_col(ncol,nz,-nc:nc,nphi*2) )
    diag_spec_col(:,:,:,:) = 0.
  end if
  if ( l_spec_ctop_o ) then
    if ( allocated(diag_spec_ct_col) )  deallocate( diag_spec_ct_col )
    allocate( diag_spec_ct_col(ncol,-nc:nc,nphi*2) )
    diag_spec_ct_col(:,:,:) = 0.
  end if

  call get_wm_hg2cgwp

  f2(:) = (two_omega*two_omega)*sin(lat(:)*deg2rad)**2

  call basic_u_phi(u_flev,v_flev,phi_deg, ub_flev)

  ! for calculating saturation spectrum
  tmp  = beta_wm/(sqrt(2.0)*pi)
  pm1  = p_wm - 1.0
  c2mp = 2.0 - p_wm
  do k=1, nz
  do l=1, ncol
    ome_min = sqrt(max(f2(l), nbv_flev(l,k)*beta_eq/mstar_wm))
    b0 = pm1*ome_min**pm1/(1.0-(ome_min/nbv_flev(l,k))**pm1)
!    if (p_wm == 1.0)  b0 = 1.0/log(nbv_flev(l,k)/ome_min)
    w0 = (nbv_flev(l,k)**c2mp-ome_min**c2mp)/c2mp
    fact_s(l,k) = abs(tmp*b0*w0*rho_flev(l,k))
  enddo
  enddo

  ! momentum flux profiles
  mf_pos(:,:,:) = 0.  ;  mf_neg(:,:,:) = 0.


  L_PHI:  DO iphi=1, nphi
  L_COL:  DO l=1, ncol

  mfsp(:) = mfs_ct(:,l,iphi)
  c_intr(:) = c_phase(:) - ub_flev(l,kcta(l),iphi)

  tmpi = minloc(abs(c_intr),1) - (nc+1)

  ineg = tmpi - 1
  ipos = tmpi + 1

  do k=kcta(l), nz

    ! calculate the saturation spectrum using Warner and
    ! McIntyre's method but as a function of phase speed,
    ! and apply the saturation condition

    ! for components having positive MF
    mfsp_s(ipos:nc) = fact_s(l,k)*c_intr(ipos:nc)/nbv_flev(l,k)
    do ic=ipos, nc
      mfsp(ic) = min(mfsp_s(ic), mfsp(ic))
    enddo
    mf_pos(l,k,iphi) = sum(mfsp(ipos:nc))*dc

    ! for components having negative MF
    mfsp_s(-nc:ineg) = fact_s(l,k)*c_intr(-nc:ineg)/nbv_flev(l,k)
    do ic=-nc, ineg
      mfsp(ic) = max(mfsp_s(ic), mfsp(ic))
    enddo
    mf_neg(l,k,iphi) = sum(mfsp(-nc:ineg))*dc

    if ( l_spec_o ) then
      diag_spec_col(l,k,ipos:nc,iphi) = mfsp(ipos:nc)
      diag_spec_col(l,k,-ineg:nc,nphi+iphi) = -mfsp(ineg:-nc:-1)
    end if
    if ( l_spec_ctop_o .and. k == kcta(l) ) then
      diag_spec_ct_col(l,ipos:nc,iphi) = mfsp(ipos:nc)
      diag_spec_ct_col(l,-ineg:nc,nphi+iphi) = -mfsp(ineg:-nc:-1)
    end if

    ! check the critical levels between k and k+1 levels,
    ! and prepare the next-level calculations
    if (k /= nz) then

      ! update c_intr for the next level
      c_intr(:) = c_phase(:) - ub_flev(l,k+1,iphi)

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

  ENDDO  L_COL
  ENDDO  L_PHI

  RETURN

END SUBROUTINE propdiss

SUBROUTINE calc_drag(                                                    &
    ncol, nz, z_flev, mf_pos, mf_neg, rho, i_conserve,                   &
    drag_u, drag_v, kcta )

  USE param_gwp,  ONLY: nphi, phi_deg

  ! i_conserve |  option for the column momentum conservation
  !            |  0: not conserved (e.g., offline calculation)
  !            |  1: conserved by vanishing the flux at the lid
  !            |  2: conserved by reducing the parameterized
  !            |     counteractive drag below the cloud top in a
  !            |     practical sense (not physical but inducing only
  !            |     little change in wind tendency)

  integer                      , intent(in) ::  ncol, nz
  real, dimension(ncol,nz)     , intent(in) ::  z_flev, rho
  real, dimension(ncol,nz,nphi), intent(in) ::  mf_pos, mf_neg
  integer                      , intent(in) ::  i_conserve

  integer, dimension(ncol), intent(in), optional ::  kcta

  real, dimension(ncol,nz), intent(out) ::  drag_u, drag_v

  real, dimension(nphi)    ::  cosphi, sinphi
  real, dimension(ncol,nz) ::  tmp, mf_x, mf_y
  integer                  ::  k,l,iphi

  include 'c_math.inc'

  mf_x(:,:) = 0.  ;  mf_y(:,:) = 0.

  cosphi(:) = cos(phi_deg(:)*deg2rad)
  sinphi(:) = sin(phi_deg(:)*deg2rad)

  do iphi=1, nphi
    tmp(:,:) = mf_pos(:,:,iphi) + mf_neg(:,:,iphi)
    mf_x(:,:) = mf_x(:,:) + tmp(:,:)*cosphi(iphi)
    mf_y(:,:) = mf_y(:,:) + tmp(:,:)*sinphi(iphi)
  enddo

  if (i_conserve == 1) then
    mf_x(:,nz) = 0.  ;  mf_y(:,nz) = 0.
  else
    if ( i_conserve == 2 .and. ( .not. present(kcta) ) ) then
      write(6,*) 'ERROR: calc_drag'
      write(6,*) '  kcta not provided when i_conserve = 2.'
      STOP
    end if
  end if

  drag_u(:,1) = 0.  ;  drag_v(:,1) = 0.

  do k=2, nz
    tmp(:,k) = 1.0/(rho(:,k)*(z_flev(:,k) - z_flev(:,k-1)))
    drag_u(:,k) = (mf_x(:,k-1) - mf_x(:,k))*tmp(:,k)
    drag_v(:,k) = (mf_y(:,k-1) - mf_y(:,k))*tmp(:,k)
  enddo

  if (i_conserve < 2)  RETURN

  ! If the outward flux of momentum at the model lid is allowed,
  ! drags below the cloud top must be reduced somehow in order to
  ! conserve the momentum budget of the model
 
  if (i_conserve == 2) then
    do l=1, ncol
      k = kcta(l)
      tmp(l,1) = 1.0/rho(l,k)/(z_flev(l,k) - z_flev(l,k-1))
      drag_u(l,k) = (mf_x(l,nz) - mf_x(l,k))*tmp(l,1)
      drag_v(l,k) = (mf_y(l,nz) - mf_y(l,k))*tmp(l,1)
    enddo
  end if

  RETURN

END subroutine calc_drag

END module prop_diss


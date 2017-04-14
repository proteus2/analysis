MODULE prop_diss

  implicit none

  real, dimension(:,:,:), allocatable ::  mf_pos, mf_neg

  real, dimension(:,:,:,:), allocatable ::  diag_spec

  real, dimension(:,:,:), allocatable ::  diag_spec_ctop

  real, dimension(:,:), allocatable ::  mflx_east, mflx_west,            &
                                        mflx_north, mflx_south 

  real, dimension(:), allocatable ::  mflx_ct_east, mflx_ct_west,        &
                                      mflx_ct_north, mflx_ct_south

  real, dimension(:,:), allocatable ::  drag_u, drag_v


CONTAINS


SUBROUTINE propdiss(  &
    ncol, nz,      &
    u_flev, v_flev, nbv_flev, rho_flev, lat,                 &
    kcta, mfs_ct )
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
  USE switch_dump,  ONLY: l_spec_o, l_spec_ctop_o
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

! LOCAL VARIABLES

! data arrays

  real, dimension(ncol)    ::  f2
  real, dimension(nz)      ::  ub_flev
  real, dimension(-nc:nc)  ::  c_intr, mfsp, mfsp_s
  real, dimension(ncol,nz) ::  fact_s

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

  if ( allocated(mf_pos) )  deallocate( mf_pos, mf_neg )
  allocate( mf_pos(ncol,nz,nphi), mf_neg(ncol,nz,nphi) )
  mf_pos(:,:,:) = 0.  ;  mf_neg(:,:,:) = 0.

  if ( l_spec_o ) then
    if ( allocated(diag_spec) )  deallocate( diag_spec )
    allocate( diag_spec(ncol,nz,-nc:nc,nphi*2) )
    diag_spec(:,:,:,:) = 0.
  end if
  if ( l_spec_ctop_o ) then
    if ( allocated(diag_spec_ctop) )  deallocate( diag_spec_ctop )
    allocate( diag_spec_ctop(ncol,-nc:nc,nphi*2) )
    diag_spec_ctop(:,:,:) = 0.
  end if

  call get_wm_hg2cgwp

  f2(:) = (two_omega*two_omega)*sin(lat(:)*deg2rad)**2

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


  L_PHI:  DO iphi=1, nphi
  L_COL:  DO l=1, ncol

  mfsp(:) = mfs_ct(:,l,iphi)

  ub_flev(:) = u_flev(l,:)*cosphi(iphi) + v_flev(l,:)*sinphi(iphi)

  c_intr(:) = c_phase(:) - ub_flev(kcta(l))

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
      diag_spec(l,k,ipos:nc,iphi) = mfsp(ipos:nc)
      diag_spec(l,k,-ineg:nc,nphi+iphi) = -mfsp(ineg:-nc:-1)
    end if
    if ( l_spec_ctop_o .and. k == kcta(l) ) then
      diag_spec_ctop(l,ipos:nc,iphi) = mfsp(ipos:nc)
      diag_spec_ctop(l,-ineg:nc,nphi+iphi) = -mfsp(ineg:-nc:-1)
    end if

    ! check the critical levels between k and k+1 levels,
    ! and prepare the next-level calculations
    if (k /= nz) then

      ! update c_intr for the next level
      c_intr(:) = c_phase(:) - ub_flev(k+1)

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
    ncol, nz, z_flev, rho, i_conserve,                                   &
    kcta )

  USE param_gwp  ,  ONLY: nphi, cosphi, sinphi
  USE switch_dump,  ONLY: l_drag_u_o, l_drag_v_o

  ! i_conserve |  option for the column momentum conservation
  !            |  0: not conserved (e.g., offline calculation)
  !            |  1: conserved by vanishing the flux at the lid
  !            |  2: conserved by reducing the parameterized
  !            |     counteractive drag below the cloud top in a
  !            |     practical sense (not physical but inducing only
  !            |     little change in wind tendency)

  integer                      , intent(in) ::  ncol, nz
  real, dimension(ncol,nz)     , intent(in) ::  z_flev, rho
  integer                      , intent(in) ::  i_conserve

  integer, dimension(ncol), intent(in), optional ::  kcta

  real, dimension(ncol,nz) ::  tmp, mf_x, mf_y
  integer                  ::  k,l,iphi

  include 'c_math.inc'

  if ( allocated(drag_u) )  deallocate( drag_u, drag_v )
  allocate( drag_u(ncol,nz), drag_v(ncol,nz) )
  drag_u(:,:) = 0.  ;  drag_v(:,:) = 0.

  mf_x(:,:) = 0.  ;  mf_y(:,:) = 0.

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

SUBROUTINE mflux_ewns(ncol,nz)

  USE param_gwp  ,  ONLY: nphi, cosphi, sinphi
  USE switch_dump,  ONLY: l_mflx_u_o, l_mflx_v_o

  implicit none

  integer, intent(in) ::  ncol, nz

  integer ::  iphi

  if (ncol < 1)  RETURN

  if ( l_mflx_u_o ) then

    if ( allocated(mflx_east) )  deallocate( mflx_east, mflx_west )
    allocate( mflx_east(ncol,nz), mflx_west(ncol,nz) )
    mflx_east(:,:) = 0.  ;  mflx_west(:,:) = 0.

    do iphi=1, nphi
      if (cosphi(iphi) > 0.) then
        mflx_east(:,:) = mflx_east(:,:) + mf_pos(:,:,iphi)*cosphi(iphi)
        mflx_west(:,:) = mflx_west(:,:) + mf_neg(:,:,iphi)*cosphi(iphi)
      else
        mflx_east(:,:) = mflx_east(:,:) + mf_neg(:,:,iphi)*cosphi(iphi)
        mflx_west(:,:) = mflx_west(:,:) + mf_pos(:,:,iphi)*cosphi(iphi)
      end if
    enddo 
    mflx_west(:,:) = mflx_west(:,:)*(-1.)

  end if

  if ( l_mflx_v_o ) then

    if ( allocated(mflx_north) )  deallocate( mflx_north, mflx_south )
    allocate( mflx_north(ncol,nz), mflx_south(ncol,nz) )
    mflx_north(:,:) = 0.  ;  mflx_south(:,:) = 0.

    do iphi=1, nphi
      mflx_north(:,:) = mflx_north(:,:) + mf_pos(:,:,iphi)*sinphi(iphi)
      mflx_south(:,:) = mflx_south(:,:) + mf_neg(:,:,iphi)*sinphi(iphi)
    enddo
    mflx_south(:,:) = mflx_south(:,:)*(-1.)

  end if

  RETURN

END SUBROUTINE mflux_ewns

SUBROUTINE mflux_ewns_ctop(ncol,nz,kcta)

  USE param_gwp  ,  ONLY: nphi, cosphi, sinphi
  USE switch_dump,  ONLY: l_mflx_u_ctop_o, l_mflx_v_ctop_o

  implicit none

  integer                 , intent(in) ::  ncol, nz
  integer, dimension(ncol), intent(in) ::  kcta

  integer ::  k,l,iphi

  if (ncol < 1)  RETURN

  if ( l_mflx_u_ctop_o ) then

    if ( allocated(mflx_ct_east) )  deallocate( mflx_ct_east,            &
                                                mflx_ct_west )
    allocate( mflx_ct_east(ncol), mflx_ct_west(ncol) )
    mflx_ct_east(:) = 0.  ;  mflx_ct_west(:) = 0.

    do iphi=1, nphi
    do l=1, ncol
      k = kcta(l)
      if (cosphi(iphi) > 0.) then
        mflx_ct_east(l) = mflx_ct_east(l) + mf_pos(l,k,iphi)*cosphi(iphi)
        mflx_ct_west(l) = mflx_ct_west(l) + mf_neg(l,k,iphi)*cosphi(iphi)
      else
        mflx_ct_east(l) = mflx_ct_east(l) + mf_neg(l,k,iphi)*cosphi(iphi)
        mflx_ct_west(l) = mflx_ct_west(l) + mf_pos(l,k,iphi)*cosphi(iphi)
      end if
    enddo
    enddo
    mflx_ct_west(:) = mflx_ct_west(:)*(-1.)

  end if

  if ( l_mflx_v_ctop_o ) then

    if ( allocated(mflx_ct_north) )  deallocate( mflx_ct_north,          &
                                                 mflx_ct_south )
    allocate( mflx_ct_north(ncol), mflx_ct_south(ncol) )
    mflx_ct_north(:) = 0.  ;  mflx_ct_south(:) = 0.

    do l=1, ncol
      k = kcta(l)
      mflx_ct_north(l) = sum(mf_pos(l,k,:)*sinphi(:))
      mflx_ct_south(l) = sum(mf_neg(l,k,:)*sinphi(:))
    enddo
    mflx_ct_south(:) = mflx_ct_south(:)*(-1.)

  end if

  RETURN

END SUBROUTINE mflux_ewns_ctop

END module prop_diss


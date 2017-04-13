MODULE diag_prop_diss

  implicit none

  logical ::  l_mflx_u_o, l_mflx_v_o
  logical ::  l_mflx_u_ctop_o, l_mflx_v_ctop_o

  real, dimension(:,:), allocatable ::  mflx_east_col, mflx_west_col,    &
                                        mflx_north_col, mflx_south_col 

  real, dimension(:), allocatable ::  mflx_ct_e_col, mflx_ct_w_col,      &
                                      mflx_ct_n_col, mflx_ct_s_col

CONTAINS


SUBROUTINE mflux_ewns(ncol,nz,mf_pos,mf_neg)

  USE param_gwp,  ONLY: nphi, phi_deg

  implicit none

  integer                      , intent(in) ::  ncol, nz
  real, dimension(ncol,nz,nphi), intent(in) ::  mf_pos, mf_neg

  real, dimension(nphi) ::  cosphi, sinphi
  integer               ::  iphi

  include 'c_math.inc'

  if (ncol < 1)  RETURN

  if ( l_mflx_u_o ) then

    if ( allocated(mflx_east_col) )  deallocate( mflx_east_col )
    if ( allocated(mflx_west_col) )  deallocate( mflx_west_col )
    allocate( mflx_east_col(ncol,nz), mflx_west_col(ncol,nz) )
    mflx_east_col(:,:) = 0.  ;  mflx_west_col(:,:) = 0.

    cosphi(:) = cos(phi_deg(:)*deg2rad)

    do iphi=1, nphi
      if (cosphi(iphi) > 0.) then
        mflx_east_col(:,:) = mflx_east_col(:,:) + mf_pos(:,:,iphi)*cosphi(iphi)
        mflx_west_col(:,:) = mflx_west_col(:,:) + mf_neg(:,:,iphi)*cosphi(iphi)
      else
        mflx_east_col(:,:) = mflx_east_col(:,:) + mf_neg(:,:,iphi)*cosphi(iphi)
        mflx_west_col(:,:) = mflx_west_col(:,:) + mf_pos(:,:,iphi)*cosphi(iphi)
      end if
    enddo 
    mflx_west_col(:,:) = mflx_west_col(:,:)*(-1.)

  end if

  if ( l_mflx_v_o ) then

    if ( allocated(mflx_north_col) )  deallocate( mflx_north_col )
    if ( allocated(mflx_south_col) )  deallocate( mflx_south_col )
    allocate( mflx_north_col(ncol,nz), mflx_south_col(ncol,nz) )
    mflx_north_col(:,:) = 0.  ;  mflx_south_col(:,:) = 0.

    sinphi(:) = sin(phi_deg(:)*deg2rad)

    do iphi=1, nphi
      mflx_north_col(:,:) = mflx_north_col(:,:) + mf_pos(:,:,iphi)*sinphi(iphi)
      mflx_south_col(:,:) = mflx_south_col(:,:) + mf_neg(:,:,iphi)*sinphi(iphi)
    enddo
    mflx_south_col(:,:) = mflx_south_col(:,:)*(-1.)

  end if

  RETURN

END SUBROUTINE mflux_ewns

SUBROUTINE mflux_ewns_ctop(ncol,nz,mf_pos,mf_neg,kcta)

  USE param_gwp,  ONLY: nphi, phi_deg

  implicit none

  integer                      , intent(in) ::  ncol, nz
  real, dimension(ncol,nz,nphi), intent(in) ::  mf_pos, mf_neg
  integer, dimension(ncol)     , intent(in) ::  kcta

  real, dimension(nphi) ::  cosphi, sinphi
  integer               ::  k,l,iphi

  include 'c_math.inc'

  if (ncol < 1)  RETURN

  if ( l_mflx_u_ctop_o ) then

    if ( allocated(mflx_ct_e_col) )  deallocate( mflx_ct_e_col )
    if ( allocated(mflx_ct_w_col) )  deallocate( mflx_ct_w_col )
    allocate( mflx_ct_e_col(ncol), mflx_ct_w_col(ncol) )
    mflx_ct_e_col(:) = 0.  ;  mflx_ct_w_col(:) = 0.

    cosphi(:) = cos(phi_deg(:)*deg2rad)

    do iphi=1, nphi
    do l=1, ncol
      k = kcta(l)
      if (cosphi(iphi) > 0.) then
        mflx_ct_e_col(l) = mflx_ct_e_col(l) + mf_pos(l,k,iphi)*cosphi(iphi)
        mflx_ct_w_col(l) = mflx_ct_w_col(l) - mf_neg(l,k,iphi)*cosphi(iphi)
      else
        mflx_ct_e_col(l) = mflx_ct_e_col(l) + mf_neg(l,k,iphi)*cosphi(iphi)
        mflx_ct_w_col(l) = mflx_ct_w_col(l) - mf_pos(l,k,iphi)*cosphi(iphi)
      end if
    enddo
    enddo
    mflx_ct_w_col(:) = mflx_ct_w_col(:)*(-1.)

  end if

  if ( l_mflx_v_ctop_o ) then

    if ( allocated(mflx_ct_n_col) )  deallocate( mflx_ct_n_col )
    if ( allocated(mflx_ct_s_col) )  deallocate( mflx_ct_s_col )
    allocate( mflx_ct_n_col(ncol), mflx_ct_s_col(ncol) )

    sinphi(:) = sin(phi_deg(:)*deg2rad)

    do l=1, ncol
      k = kcta(l)
      mflx_ct_n_col(l) = sum(mf_pos(l,k,:)*sinphi(:))
      mflx_ct_s_col(l) = sum(mf_neg(l,k,:)*sinphi(:))*(-1.)
    enddo
    mflx_ct_s_col(:) = mflx_ct_s_col(:)*(-1.)

  end if

  RETURN

END SUBROUTINE mflux_ewns_ctop

END module diag_prop_diss


MODULE zonal_average

  USE switch_dump
  USE prop_diss,  ONLY: mflx_ct_east, mflx_ct_west,                      &
                        mflx_ct_north, mflx_ct_south,                    &
                        mflx_east, mflx_west, mflx_north, mflx_south,    &
                        drag_u, drag_v, diag_spec_ctop, diag_spec

  implicit none

  real, dimension(:)      , allocatable ::  zm_mflx_ct_east,             &
                                            zm_mflx_ct_west,             &
                                            zm_mflx_ct_north,            &
                                            zm_mflx_ct_south
  real, dimension(:,:)    , allocatable ::  zm_mflx_east,                &
                                            zm_mflx_west,                &
                                            zm_mflx_north,               &
                                            zm_mflx_south,               &
                                            zm_drag_u,                   &
                                            zm_drag_v
  real, dimension(:,:,:)  , allocatable ::  zm_diag_spec_ctop
  real, dimension(:,:,:,:), allocatable ::  zm_diag_spec

CONTAINS


SUBROUTINE alloc_zonal_avg(ny,nz,nc,nphi)

  integer, intent(in) ::  ny, nz, nc, nphi

  if ( l_mflx_u_ctop_o ) then
    allocate( zm_mflx_ct_east(ny), zm_mflx_ct_west(ny) )
    zm_mflx_ct_east(:) = 0.  ;  zm_mflx_ct_west(:) = 0.
  end if
  if ( l_mflx_v_ctop_o ) then
    allocate( zm_mflx_ct_north(ny), zm_mflx_ct_south(ny) )
    zm_mflx_ct_north(:) = 0.  ;  zm_mflx_ct_south(:) = 0.
  end if

  if ( l_mflx_u_o ) then
    allocate( zm_mflx_east(ny,nz), zm_mflx_west(ny,nz) )
    zm_mflx_east(:,:) = 0.  ;  zm_mflx_west(:,:) = 0.
  end if
  if ( l_mflx_v_o ) then
    allocate( zm_mflx_north(ny,nz), zm_mflx_south(ny,nz) )
    zm_mflx_north(:,:) = 0.  ;  zm_mflx_south(:,:) = 0.
  end if

  if ( l_drag_u_o .or. l_drag_v_o ) then
    allocate( zm_drag_u(ny,nz), zm_drag_v(ny,nz) )
    zm_drag_u(:,:) = 0.  ;  zm_drag_v(:,:) = 0.
  end if

  if ( l_spec_ctop_o ) then
    allocate( zm_diag_spec_ctop(nc*2+1,nphi*2,ny) )
    zm_diag_spec_ctop(:,:,:) = 0.
  end if

  if ( l_spec_o ) then
    allocate( zm_diag_spec(nc*2+1,nphi*2,nz,ny) )
    zm_diag_spec(:,:,:,:) = 0.
  end if

END subroutine alloc_zonal_avg

SUBROUTINE dealloc_zonal_avg

  if ( l_mflx_u_ctop_o )  deallocate( zm_mflx_ct_east, zm_mflx_ct_west )
  if ( l_mflx_v_ctop_o )  deallocate( zm_mflx_ct_north, zm_mflx_ct_south )

  if ( l_mflx_u_o )  deallocate( zm_mflx_east, zm_mflx_west )
  if ( l_mflx_v_o )  deallocate( zm_mflx_north, zm_mflx_south )

  if ( l_drag_u_o .or. l_drag_v_o )  deallocate( zm_drag_u, zm_drag_v )

  if ( l_spec_ctop_o )  deallocate( zm_diag_spec_ctop )

  if ( l_spec_o )  deallocate( zm_diag_spec )

END subroutine dealloc_zonal_avg

SUBROUTINE zonal_avg(j_col,nx, w_acc)

  integer, dimension(:), intent(in) ::  j_col
  integer              , intent(in) ::  nx
  real                 , intent(in), optional ::  w_acc

  real    ::  wgt
  integer ::  ncol
  logical ::  l_initialize
  integer ::  j,l

  ncol = size(j_col)

  wgt = 1./float(nx)
  l_initialize = .True.
  if ( present(w_acc) ) then
    wgt = wgt*w_acc
    l_initialize = .False.
  end if

  if ( l_mflx_u_ctop_o ) then
    if ( l_initialize ) then
      zm_mflx_ct_east(:) = 0.  ;  zm_mflx_ct_west(:) = 0.
    end if
    do l=1, ncol
      j = j_col(l)
      zm_mflx_ct_east(j) = zm_mflx_ct_east(j) + mflx_ct_east(l)*wgt
      zm_mflx_ct_west(j) = zm_mflx_ct_west(j) + mflx_ct_west(l)*wgt
    enddo
  end if
  if ( l_mflx_v_ctop_o ) then
    if ( l_initialize ) then
      zm_mflx_ct_north(:) = 0.  ;  zm_mflx_ct_south(:) = 0.
    end if
    do l=1, ncol
      j = j_col(l)
      zm_mflx_ct_north(j) = zm_mflx_ct_north(j) + mflx_ct_north(l)*wgt
      zm_mflx_ct_south(j) = zm_mflx_ct_south(j) + mflx_ct_south(l)*wgt
    enddo
  end if
  if ( l_mflx_u_o ) then
    if ( l_initialize ) then
      zm_mflx_east(:,:) = 0.  ;  zm_mflx_west(:,:) = 0.
    end if
    do l=1, ncol
      j = j_col(l)
      zm_mflx_east(j,:) = zm_mflx_east(j,:) + mflx_east(l,:)*wgt
      zm_mflx_west(j,:) = zm_mflx_west(j,:) + mflx_west(l,:)*wgt
    enddo
  end if
  if ( l_mflx_v_o ) then
    if ( l_initialize ) then
      zm_mflx_north(:,:) = 0.  ;  zm_mflx_south(:,:) = 0.
    end if
    do l=1, ncol
      j = j_col(l)
      zm_mflx_north(j,:) = zm_mflx_north(j,:) + mflx_north(l,:)*wgt
      zm_mflx_south(j,:) = zm_mflx_south(j,:) + mflx_south(l,:)*wgt
    enddo
  end if
  if ( l_drag_u_o .or. l_drag_v_o ) then
    if ( l_initialize ) then
      zm_drag_u(:,:) = 0.  ;  zm_drag_v(:,:) = 0.
    end if
    do l=1, ncol
      j = j_col(l)
      zm_drag_u(j,:) = zm_drag_u(j,:) + drag_u(l,:)*wgt
      zm_drag_v(j,:) = zm_drag_v(j,:) + drag_v(l,:)*wgt
    enddo
  end if
  if ( l_spec_ctop_o ) then
    if ( l_initialize )  zm_diag_spec_ctop(:,:,:) = 0.
    do l=1, ncol
      j = j_col(l)
      zm_diag_spec_ctop(:,:,j) = zm_diag_spec_ctop(:,:,j) +              &
          diag_spec_ctop(:,:,l)*wgt
    enddo
  end if
  if ( l_spec_o ) then
    if ( l_initialize )  zm_diag_spec(:,:,:,:) = 0.
    do l=1, ncol
      j = j_col(l)
      zm_diag_spec(:,:,:,j) = zm_diag_spec(:,:,:,j) +                    &
          diag_spec(:,:,:,l)*wgt
    enddo
  end if

END subroutine zonal_avg

END module zonal_average


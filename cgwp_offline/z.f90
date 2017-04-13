PROGRAM cgwp

  USE param_gwp
  USE prop_diss
  USE diag_prop_diss
  USE subr_common

  implicit none

  allocate( phi_deg(nphi) )
  nc = 30
  dc = 2.
  nphi = 2
  phi_deg = (/45.,135./)
  cfactor = 125.
 
  allocate( c_phase(-nc:nc) )
  do ic=-nc, nc
    c_phase(ic) = float(ic)*dc
  enddo

!  call src

  l_mflx_u_o = .True.  ;  l_mflx_v_o = .True.
  l_mflx_u_ctop_o = .True.  ;  l_mflx_v_ctop_o = .True.
  l_mflx_regrid = .True.
  l_spec_o = .True.  ;  l_spec_ctop_o = .True.

  call propdiss(ncol,nz,     &
    u_flev, v_flev, nbv_flev, rho_flev, lat,                 &
    kcta, mfs_ct,                         &
    mf_pos, mf_neg )
  ! diag_spec_col ; diag_spec_ct_col

  call mflux_ewns_ctop(ncol,nz,mf_pos,mf_neg,kcta)
  ! mflx_XXXX_col
  call mflux_ewns     (ncol,nz,mf_pos,mf_neg)
  ! mflx_ct_X_col
 
  call calc_drag(ncol,nz,z_flev_col,mf_x,mf_y,rho_col,0,                 &
                 drag_u_col,drag_v_col)

  if ( l_mflx_regrid ) then

  if ( l_spec_o ) then
    call column_to_grid3d_sp(icol,jcol,diag_spec_col,1, diag_spec_grd)
  end if
  if ( l_spec_ctop_o ) then
    call column_to_grid2d_sp(icol,jcol,diag_spec_ct_col,1, diag_spec_ct_grd)
  end if
  call column_to_grid3d(icol,jcol,mflx_east_col,1, mflx_east_grd)
  call column_to_grid3d(icol,jcol,drag_u_col,1, drag_u_grd)
  call column_to_grid3d(icol,jcol,drag_v_col,0, drag_v_grd)

  end if

  STOP

END program cgwp


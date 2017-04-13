PROGRAM cgwp

  USE param_gwp
  USE switch_dump
  USE prop_diss
  USE subr_common

  implicit none

  integer ::  nx = 3, ny = 2, nz = 30
  integer ::  ncol

  real, dimension(:,:), allocatable ::  u_flev, v_flev, nbv_flev
  real, dimension(:,:), allocatable ::  rho_flev, z_flev_col

  real, dimension(:,:), allocatable ::  rho_col
  real, dimension(:)  , allocatable ::  lat

  real, dimension(:,:,:), allocatable ::  mfs_ct
  integer, dimension(:) , allocatable ::  kcta, icol, jcol

  nc = 30
  dc = 2.
  nphi = 2
  allocate( phi_deg(nphi) )
  phi_deg = (/45.,135./)
  cfactor = 125.
 
  call set_spec_param

  ncol = 4

  call switch_defaults
  l_mflx_regrid = .False.
!  l_drag_u_o = .True.  ;  l_drag_v_o = .True.
!  l_mflx_u_o = .True.  ;  l_mflx_v_o = .True.
!  l_mflx_u_ctop_o = .True.  ;  l_mflx_v_ctop_o = .True.
!  l_spec_o = .True.  ;  l_spec_ctop_o = .True.

  allocate( u_flev(ncol,nz), v_flev(ncol,nz) )
  allocate( nbv_flev(ncol,nz) )
  allocate( rho_flev(ncol,nz), z_flev_col(ncol,nz) )
  allocate( lat(ncol) )

  allocate( rho_col(ncol,nz) )

  allocate( icol(ncol), jcol(ncol) )

  if ( l_mflx_regrid )  call alloc_grd_output

!  call src

  allocate( kcta(ncol) )
  allocate( mfs_ct(-nc:nc,ncol,nphi) )

  call propdiss(ncol,nz,     &
    u_flev, v_flev, nbv_flev, rho_flev, lat,                 &
    kcta, mfs_ct )
  ! mf_pos ; mf_neg ; diag_spec ; diag_spec_ctop

  if ( l_mflx_u_ctop_o .or. l_mflx_v_ctop_o )                            &
     call mflux_ewns_ctop(ncol,nz,kcta)
     ! mflx_XXXX
  if ( l_mflx_u_o .or. l_mflx_v_o )  call mflux_ewns(ncol,nz)
     ! mflx_ct_XXXX
 
  if ( l_drag_u_o .or. l_drag_v_o )                                      &
     call calc_drag(ncol,nz,z_flev_col,rho_col,0)


  if ( l_mflx_regrid ) then

!  if ( l_spec_o ) then
!    call column_to_grid3d_sp(icol,jcol,diag_spec,1, diag_spec_grd)
!  end if
!  if ( l_spec_ctop_o ) then
!    call column_to_grid2d_sp(icol,jcol,diag_spec_ctop,1, diag_spec_ct_grd)
!  end if

!  call column_to_grid3d(icol,jcol,mflx_east,1, mflx_east_grd)
!
!  if ( l_drag_u_o ) then
!    call column_to_grid3d(icol,jcol,drag_u,1, drag_u_grd)
!  end if
!  if ( l_drag_v_o ) then
!    call column_to_grid3d(icol,jcol,drag_v,0, drag_v_grd)
!  end if

  end if

  STOP

END program cgwp

SUBROUTINE alloc_grd_output

!  if ( l_drag_u_o )  allocate( drag_u_grd(nx,ny,nz) )
!  if ( l_drag_v_o )  allocate( drag_v_grd(nx,ny,nz) )
!  if ( l_mflx_u_o )  allocate( mflx_east_grd(nx,ny,nz),                  &
!                               mflx_west_grd(nx,ny,nz) )
!  if ( l_mflx_v_o )  allocate( mflx_north_grd(nx,ny,nz),                 &
!                               mflx_south_grd(nx,ny,nz) )
!  if ( l_mflx_u_ctop_o )  allocate( mflx_ct_e_grd(nx,ny),                &
!                                    mflx_ct_w_grd(nx,ny) )
!  if ( l_mflx_v_ctop_o )  allocate( mflx_ct_n_grd(nx,ny),                &
!                                    mflx_ct_s_grd(nx,ny) )

END subroutine alloc_grd_output


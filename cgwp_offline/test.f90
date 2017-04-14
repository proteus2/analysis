PROGRAM cgwp

  USE param_gwp
  USE switch_dump
  USE prop_diss
  USE subr_common
  USE netio

  implicit none

  integer ::  nx = 3, ny = 2, nz = 30
  integer ::  ncol, nvo

  real, dimension(:,:), allocatable ::  u_flev, v_flev, nbv_flev
  real, dimension(:,:), allocatable ::  rho_flev, z_flev_col

  real, dimension(:,:), allocatable ::  rho_col
  real, dimension(:)  , allocatable ::  lat

  real, dimension(:,:,:), allocatable ::  mfs_ct
  integer, dimension(:) , allocatable ::  kcta, icol, jcol

  character(len=128) ::  file_o
  integer ::  k, l, iv

  type(vset), dimension(:), allocatable ::  set

  nc = 30
  dc = 2.
  nphi = 2
  allocate( phi_deg(nphi) )
  phi_deg = (/45.,135./)
  cfactor = 125.
 
  call set_spec_param

  call switch_defaults
!  l_spec_o = .False.
  l_mflx_regrid = .False.

  ncol = 4

  allocate( u_flev(ncol,nz), v_flev(ncol,nz) )
  allocate( nbv_flev(ncol,nz) )
  allocate( rho_flev(ncol,nz), z_flev_col(ncol,nz) )
  allocate( lat(ncol) )

  allocate( rho_col(ncol,nz) )

  allocate( icol(ncol), jcol(ncol) )

!  call src

  allocate( kcta(ncol) )
  allocate( mfs_ct(-nc:nc,ncol,nphi) )

  u_flev = 0.  ;  v_flev = 0.
  nbv_flev(:,:) = 2.6e-2
  do k=1, nz
    z_flev_col(:,k) = 800.*float(k)
    rho_flev(:,k) = exp(-z_flev_col(:,k)/7.e3)
    rho_col (:,k) = exp(-(z_flev_col(:,k)-400.)/7.e3)
  enddo
  lat = (/-20.,-20.,0.,0./)
  icol = (/1,2,1,3/)
  jcol = (/1,1,1,2/)

  kcta = (/6,9,12,9/)
  mfs_ct(:,:,:) = 0.
  mfs_ct(-nc:-1,:,:) = -1.
  mfs_ct(1:nc,:,:) = 1.

file_o = '/Users/kyh/analy/cgwp_offline/zzz.nc'

  call propdiss(ncol,nz,     &
    u_flev, v_flev, nbv_flev, rho_flev, lat,                 &
    kcta, mfs_ct )
  ! mf_pos ; mf_neg ; diag_spec ; diag_spec_ctop

  if ( l_mflx_u_ctop_o .or. l_mflx_v_ctop_o )                            &
     call mflux_ewns_ctop(ncol,nz,kcta)
     ! mflx_ct_XXXX
  if ( l_mflx_u_o .or. l_mflx_v_o )  call mflux_ewns(ncol,nz)
     ! mflx_XXXX
 
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


  call get_nv_output

  call put_vars_set

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nvo,set,'CGWP offline calculation')

! END

  call finalize

  STOP


CONTAINS


SUBROUTINE get_nv_output

  nvo = 0
  if ( l_mflx_u_ctop_o )  nvo = nvo + 2
  if ( l_mflx_v_ctop_o )  nvo = nvo + 2
  if ( l_mflx_u_o      )  nvo = nvo + 2
  if ( l_mflx_v_o      )  nvo = nvo + 2
  if ( l_drag_u_o      )  nvo = nvo + 1
  if ( l_drag_v_o      )  nvo = nvo + 1

END subroutine get_nv_output

SUBROUTINE put_vars_set

  character(len=32) ::  axis(4)
  integer           ::  ndim(4)

  real, dimension(ncol) ::  coln
  real, dimension(nz)   ::  z

  do l=1, ncol
    coln(l) = float(l)
  enddo
  z(:) = z_flev_col(1,:)

  allocate( set(nvo) )

  iv = 0

  axis = (/'case','    ','    ','    '/)
  ndim = (/ncol,1,1,1/)
  if ( l_mflx_u_ctop_o ) then
    call xxxx1d(iv,'mflx_ct_east',mflx_ct_east,axis,ndim, coln)
    call xxxx1d(iv,'mflx_ct_west',mflx_ct_west,axis,ndim, coln)
  end if
  if ( l_mflx_v_ctop_o ) then
    call xxxx1d(iv,'mflx_ct_north',mflx_ct_north,axis,ndim, coln)
    call xxxx1d(iv,'mflx_ct_south',mflx_ct_south,axis,ndim, coln)
  end if

  axis = (/'z   ','case','    ','    '/)
  ndim = (/nz,ncol,1,1/)
  if ( l_mflx_u_o ) then
    call xxxx2d(iv,'mflx_east',transpose(mflx_east),axis,ndim, z,coln)
    call xxxx2d(iv,'mflx_west',transpose(mflx_west),axis,ndim, z,coln)
  end if
  if ( l_mflx_v_o ) then
    call xxxx2d(iv,'mflx_north',transpose(mflx_north),axis,ndim, z,coln)
    call xxxx2d(iv,'mflx_south',transpose(mflx_south),axis,ndim, z,coln)
  end if
  if ( l_drag_u_o ) then
    call xxxx2d(iv,'drag_u',transpose(drag_u),axis,ndim, z,coln)
  end if
  if ( l_drag_v_o ) then
    call xxxx2d(iv,'drag_v',transpose(drag_v),axis,ndim, z,coln)
  end if

END subroutine put_vars_set

SUBROUTINE xxxx1d(ivl,varname,out1d,axis,ndim,axis1)

  character(len=*)     , intent(in) ::  varname, axis(4)
  real   , dimension(:), intent(in) ::  out1d
  integer, dimension(4), intent(in) ::  ndim
  real   , dimension(:), intent(in) ::  axis1

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axis,ndim)

  set(ivl)%axis1(:) = axis1(:)
 
  allocate( set(ivl)%var_out(ndim(1),ndim(2),ndim(3),ndim(4)) )
  set(ivl)%var_out(:,1,1,1) = out1d(:)

END subroutine xxxx1d

SUBROUTINE xxxx2d(ivl,varname,out2d,axis,ndim,axis1,axis2)

  character(len=*)       , intent(in) ::  varname, axis(4)
  real   , dimension(:,:), intent(in) ::  out2d
  integer, dimension(4)  , intent(in) ::  ndim
  real   , dimension(:)  , intent(in) ::  axis1, axis2

  integer, intent(inout) ::  ivl

  call setdim0(ivl,varname,axis,ndim)

  set(ivl)%axis1(:) = axis1(:)
  set(ivl)%axis2(:) = axis2(:)
 
  allocate( set(ivl)%var_out(ndim(1),ndim(2),ndim(3),ndim(4)) )
  set(ivl)%var_out(:,:,1,1) = out2d(:,:)

END subroutine xxxx2d

SUBROUTINE setdim0(ivl,varname,axis,ndim)

  character(len=*)     , intent(in) ::  varname, axis(4)
  integer, dimension(4), intent(in) ::  ndim

  integer, intent(inout) ::  ivl

  ivl = ivl + 1

  set(ivl)%vname = trim(varname)

  set(ivl)%axis(:) = axis(:)
  set(ivl)%nd  (:) = ndim(:)
  allocate( set(ivl)%axis1(ndim(1)) )
  allocate( set(ivl)%axis2(ndim(2)) )
  allocate( set(ivl)%axis3(ndim(3)) )
  allocate( set(ivl)%axis4(ndim(4)) )
  set(ivl)%axis1(:) = -999.
  set(ivl)%axis2(:) = -999.
  set(ivl)%axis3(:) = -999.
  set(ivl)%axis4(:) = -999.

END subroutine setdim0

SUBROUTINE finalize

!  deallocate( var4d, var3d )
!  deallocate( lon, lat, p, t, t2pt )
  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

END subroutine finalize


END program cgwp


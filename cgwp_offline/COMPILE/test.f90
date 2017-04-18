PROGRAM cgwp

  USE param_gwp
  USE switch_dump
  USE mflx_cldtop
  USE prop_diss
  USE subr_common
  USE netio

  implicit none

  integer ::  nx = 3, ny = 2, nz = 30
  integer ::  ncol, nvo

  real, dimension(:,:), allocatable ::  u_flev, v_flev, nbv_flev
  real, dimension(:,:), allocatable ::  rho_flev, z_flev

  real, dimension(:,:), allocatable ::  rho_uvlev
  real, dimension(:)  , allocatable ::  lat

  character(len=128) ::  file_o
  integer ::  k, l, iv

  nc = 30
  dc = 2.
  nphi = 2
  allocate( phi_deg(nphi) )
  phi_deg = (/45.,135./)
  cfactor = 125.
 
  call set_spec_param

  call switch_defaults

!  l_spec_on = .False.

  call get_nv_output(nvo)

  ncol = 4

  allocate( u_flev(ncol,nz), v_flev(ncol,nz) )
  allocate( nbv_flev(ncol,nz) )
  allocate( rho_flev(ncol,nz), z_flev(ncol,nz) )
  allocate( lat(ncol) )

  allocate( rho_uvlev(ncol,nz) )

  u_flev = 0.  ;  v_flev = 0.
  nbv_flev(:,:) = 2.6e-2
  do k=1, nz
    z_flev(:,k) = 800.*float(k)
    rho_flev(:,k) = exp(-z_flev(:,k)/7.e3)
    rho_uvlev(:,k) = exp(-(z_flev(:,k)-400.)/7.e3)
  enddo
  lat = (/-20.,-20.,0.,0./)

file_o = '/Users/kyh/analy/cgwp_offline/zzz.nc'

!  call gw_ctop(ncol,nz,  &
!     &   eta_theta_levels,                                              &
!     &   u_flev,v_flev,rho_flev,z_flev,nbv_flev,         &
!     &   heat,heatmax,u_sfc,v_sfc,kcb,kct,                            &
!     &   t_flev                                      &
!     &   )
  ! kcta ; mfs_ct ; diag_znwcq

  allocate( kcta(ncol) )
  allocate( mfs_ct(-nc:nc,ncol,nphi) )
  kcta = (/6,9,12,9/)
  mfs_ct(:,:,:) = 0.
  mfs_ct(-nc:-1,:,:) = -1.
  mfs_ct(1:nc,:,:) = 1.

  call propdiss(ncol,nz,     &
    u_flev, v_flev, nbv_flev, rho_flev, lat,                 &
    kcta, mfs_ct )
  ! mf_pos ; mf_neg ; diag_spec_ctop ; diag_spec

  ! mflx_ct_XXXX
  if ( l_mflx_u_ctop_o .or. l_mflx_v_ctop_o )                            &
     call mflux_ewns_ctop(ncol,nz,kcta)

  ! mflx_XXXX
  if ( l_mflx_u_o .or. l_mflx_v_o )  call mflux_ewns(ncol,nz)
 
  ! drag
  if ( l_drag_u_o .or. l_drag_v_o )                                      &
     call calc_drag(ncol,nz,z_flev,rho_uvlev,0)


  call put_vars_set

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nvo,set,'CGWP offline calculation')

! END

  call finalize

  STOP


CONTAINS


SUBROUTINE put_vars_set

  character(len=32) ::  axis(4)
  integer           ::  ndim(4)

  real, dimension(ncol) ::  coln
  real, dimension(nz)   ::  z

  do l=1, ncol
    coln(l) = float(l)
  enddo
  z(:) = z_flev(1,:)

  allocate( set(nvo) )

  iv = 0

  axis = (/'case','    ','    ','    '/)
  ndim = (/ncol,1,1,1/)
  if ( l_mflx_u_ctop_o ) then
    call defset(iv,'mflx_ct_east',mflx_ct_east,axis,ndim, coln)
    call defset(iv,'mflx_ct_west',mflx_ct_west,axis,ndim, coln)
  end if
  if ( l_mflx_v_ctop_o ) then
    call defset(iv,'mflx_ct_north',mflx_ct_north,axis,ndim, coln)
    call defset(iv,'mflx_ct_south',mflx_ct_south,axis,ndim, coln)
  end if

  axis = (/'z   ','case','    ','    '/)
  ndim = (/nz,ncol,1,1/)
  if ( l_mflx_u_o ) then
    call defset(iv,'mflx_east',transpose(mflx_east),axis,ndim, z,coln)
    call defset(iv,'mflx_west',transpose(mflx_west),axis,ndim, z,coln)
  end if
  if ( l_mflx_v_o ) then
    call defset(iv,'mflx_north',transpose(mflx_north),axis,ndim, z,coln)
    call defset(iv,'mflx_south',transpose(mflx_south),axis,ndim, z,coln)
  end if
  if ( l_drag_u_o ) then
    call defset(iv,'drag_u',transpose(drag_u),axis,ndim, z,coln)
  end if
  if ( l_drag_v_o ) then
    call defset(iv,'drag_v',transpose(drag_v),axis,ndim, z,coln)
  end if

  axis = (/'c_ph','dir ','case','    '/)
  ndim = (/nc*2+1,nphi*2,ncol,1/)
  if ( l_spec_ctop_o ) then
    call defset(iv,'mflx_ct_spec',diag_spec_ctop,axis,ndim,              &
                c_phase,phi_deg2,coln)
  end if
 
  axis = (/'c_ph','dir ','z   ','case'/)
  ndim = (/nc*2+1,nphi*2,nz,ncol/)
  if ( l_spec_o ) then
    call defset(iv,'mflx_spec',diag_spec,axis,ndim,                      &
                c_phase,phi_deg2,z,coln)
  end if
 
END subroutine put_vars_set

SUBROUTINE finalize

!  deallocate( var4d, var3d )
!  deallocate( lon, lat, p, t, t2pt )
  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo
  deallocate( set )

END subroutine finalize


END program cgwp


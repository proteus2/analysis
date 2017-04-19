PROGRAM cgwp

  USE param_gwp
  USE switch_dump
  USE mflx_ctop_sc05
  USE prop_diss
  USE subr_common
  USE netio

  implicit none

  integer ::  nx = 3, ny = 2, nz = 30
  integer ::  ncol, nvo

  real, dimension(:,:), allocatable ::  u_flev, v_flev, nbv_flev,        &
                                        rho_flev, z_flev

  ! only for SC05
  real   , dimension(:,:), allocatable ::  heat_flev, t_flev
  integer, dimension(:)  , allocatable ::  kcb, kct

  ! only for propdiss
  real, dimension(:)  , allocatable ::  f_cor

  ! only for drag
  real, dimension(:,:), allocatable ::  rho_dlev

  character(len=128) ::  file_o
  integer ::  k,l,iv

  real, parameter ::  two_omega = 2.*7.292116e-5

  include 'c_math.inc'   ! deg2rad

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

  allocate( u_flev(ncol,nz), v_flev(ncol,nz), u_sfc(ncol), v_sfc(ncol) )
  allocate( nbv_flev(ncol,nz), rho_flev(ncol,nz), z_flev(ncol,nz) )

  ! for SC05
  allocate( heat_flev(ncol,nz), t_flev(ncol,nz), kcb(ncol), kct(ncol) )
  ! for propdiss
  allocate( f_cor(ncol) )
  ! for drag
  allocate( rho_dlev(ncol,nz) )

  u_flev = 0.  ;  v_flev = 0.  ;  u_sfc = 0.  ; v_sfc = 0.
  t_flev = 270.
  nbv_flev(:,:) = 2.6e-2
  do k=1, nz
    z_flev(:,k) = 800.*float(k)
    rho_flev(:,k) = exp(-z_flev(:,k)/7.e3)
    rho_dlev(:,k) = exp(-(z_flev(:,k)-400.)/7.e3)
    heat_flev(:,k) = (1./3600.)*( 1. - ((z_flev(:,k)-4.e3)/2.e3)**2 )
u_flev(:,k) = 30.*sin(z_flev(:,k)/10.e3 * 2.*3.141592)
  enddo
  kcb(:) = 1  ;  kct(:) = 20
  f_cor = (/-10.,-10.,0.,0./)
  f_cor(:) = two_omega*sin(f_cor(:)*deg2rad)

file_o = '/Users/kyh/analy/cgwp_offline/zzz.nc'

  ! u_sfc, v_sfc: allocate and specify, if possible


  allocate( kcta(ncol) )
  kcta = (/6,9,12,9/)
 
  call args_sc05(ncol,nz,u_flev,v_flev,t_flev,nbv_flev,rho_flev,z_flev,  &
                 heat_flev,kcb,kct)  ! opt: z_ref
  ! kcta ; diag_znwcq

  ! u_ct, v_ct, u_cb, v_cb, t_ct
  ! zcta, zcba, n_q, n_ct, rho_ct, cqx, cqy, heatmax
  ! kcta


  allocate( mfs_ct(-nc:nc,ncol,nphi) )
  mfs_ct(:,:,:) = 0.
  mfs_ct(-nc:-1,:,:) = -1.
  mfs_ct(1:nc,:,:) = 1.

  call calc_sc05(ncol,nz)  ! opt: shear_ct
  ! mfs_ct

  deallocate( u_sfc, v_sfc, t_flev, kcb, kct )


  call propdiss(ncol,nz,     &
    u_flev, v_flev, nbv_flev, rho_flev, f_cor,                 &
    kcta, mfs_ct )
  ! mf_pos  ;  mf_neg
  ! diag_spec_ctop  ;  diag_spec
  ! mflx_ct_XXXX  ;  mflx_XXXX
 
  deallocate( f_cor )


  ! drag
  if ( l_drag_u_o .or. l_drag_v_o )                                      &
     call calc_drag(ncol,nz,z_flev,rho_dlev,0)


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

  deallocate( u_flev, v_flev, nbv_flev, rho_flev, z_flev )
  deallocate( rho_dlev )
  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo
  deallocate( set )

END subroutine finalize


END program cgwp


PROGRAM cgwp

  USE param_gwp
  USE switch_dump
  USE mflx_ctop_sc05
  USE prop_diss
  USE subr_common
  USE netio

  implicit none

  integer, parameter ::  nz = 60
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

  real, dimension(nz) ::  z_dlev

  character(len=128) ::  file_i, file_o
  integer ::  k,l,iv, ncid, tmpi

  real, parameter ::  two_omega = 2.*7.292116e-5

  include 'c_math.inc'   ! deg2rad

  nc = 160
  dc = 0.5
  nphi = 2
  allocate( phi_deg(nphi) )
  phi_deg = (/45.,135./)
!  phi_deg = (/0.,90./)
  cfactor = 125.
 
  call set_spec_param

  call switch_defaults

!  l_spec_on = .False.

  call get_nv_output(nvo)

  ncol = 21

  allocate( heatmax(ncol) )
  allocate( u_ct(ncol), v_ct(ncol), u_cb(ncol), v_cb(ncol), t_ct(ncol) )
  allocate( u_sfc(ncol), v_sfc(ncol), cqx(ncol), cqy(ncol), rho_ct(ncol) )
  allocate( n_q(ncol), n_ct(ncol), zcta(ncol), zcba(ncol) )
  allocate( kcta(ncol) )

!  file_i = '/data18/kyh/dat/L60CGW/dchm_pdf/'//  &
!           'uanuj.dchm-midlev_pdf.1979-2006.01-12.nc'
  file_i = '/data18/kyh/dat/L60CGW/dchm_pdf/'//  &
           'uanuj.dchm-nonmidlev_pdf.1979-2006.01-12.nc'

  call opennc(file_i,ncid)
  call geta2d(ncid,'dchmax',32,ncol,1,1,heatmax)
  heatmax(:) = heatmax(:)/3600.
  call geta2d(ncid,'u_ct'  ,32,ncol,1,1,u_ct   )
  call geta2d(ncid,'v_ct'  ,32,ncol,1,1,v_ct   )
  call geta2d(ncid,'u_cb'  ,32,ncol,1,1,u_cb   )
  call geta2d(ncid,'v_cb'  ,32,ncol,1,1,v_cb   )
  call geta2d(ncid,'t_ct'  ,32,ncol,1,1,t_ct   )
  call geta2d(ncid,'u_sfc' ,32,ncol,1,1,u_sfc  )
  call geta2d(ncid,'v_sfc' ,32,ncol,1,1,v_sfc  )
  call geta2d(ncid,'cq_x'  ,32,ncol,1,1,cqx    )
  call geta2d(ncid,'cq_y'  ,32,ncol,1,1,cqy    )
  call geta2d(ncid,'rho_ct',32,ncol,1,1,rho_ct )
  call geta2d(ncid,'n_q'   ,32,ncol,1,1,n_q    )
  call geta2d(ncid,'n_ct'  ,32,ncol,1,1,n_ct   )
  call geta2d(ncid,'zcta'  ,32,ncol,1,1,zcta   )
  call geta2d(ncid,'zcba'  ,32,ncol,1,1,zcba   )
  call closenc(ncid)

  allocate( z_flev(ncol,nz) )

!  do k=1, nz
!    z_flev(:,k) = 500.*float(k)
!  enddo
  z_dlev = (/10.00335, 49.99991, 130.0014, 249.9995, 410.0026, 610.0023, &
    849.9984, 1130, 1449.997, 1810, 2209.999, 2650.004, 3129.996, 3650.002, &
    4210.004, 4810.003, 5449.999, 6129.999, 6849.996, 7609.998, 8409.996, &
    9250, 10130, 11050, 12010, 13010, 14050, 15130, 16250, 17410, 18590, &
    19770.05, 20950.33, 22131.32, 23313.88, 24499.48, 25690.25, 26889.18, &
    28100.23, 29328.48, 30580.27, 31863.35, 33186.98, 34562.13, 36001.55, &
    37520, 39134.29, 40863.49, 42729.05, 44754.91, 46967.73, 49396.89, &
    52074.75, 55036.74, 58321.52, 61971.07, 66030.91, 70550.15, 75581.72, &
    81182.44/)
  do l=1, ncol
    z_flev(l,:) = (/19.99828, 80.00153, 180.0014, 319.9977, 499.9991, 719.997, &
    979.9999, 1279.999, 1620.004, 1999.996, 2420.002, 2879.996, 3380.004, &
    3920, 4500, 5119.998, 5780, 6479.998, 7220.002, 8000.002, 8819.999, &
    9680.001, 10580, 11520, 12500, 13520, 14580, 15680, 16820, 18000, &
    19180.01, 20360.1, 21540.57, 22722.06, 23905.7, 25093.26, 26287.25, &
    27491.12, 28709.35, 29947.62, 31212.93, 32513.76, 33860.2, 35264.05, &
    36739.05, 38300.94, 39967.63, 41759.35, 43698.74, 45811.09, 48124.36, &
    50669.41, 53480.09, 56593.4, 60049.64, 63892.5, 68169.3, 72930.99, &
    78232.44, 84132.44/)
  enddo
 

  do l=1, ncol
    tmpi = minloc(abs(z_flev(l,:) - zcba(l)),1)
    kcta(l) = tmpi + minloc(abs(z_flev(l,tmpi+1:) - zcta(l)),1)
  enddo

  allocate( u_flev(ncol,nz), v_flev(ncol,nz) )
  allocate( nbv_flev(ncol,nz), rho_flev(ncol,nz) )

  ! for SC05
  allocate( heat_flev(ncol,nz), t_flev(ncol,nz), kcb(ncol), kct(ncol) )
  ! for propdiss
  allocate( f_cor(ncol) )
  ! for drag
  allocate( rho_dlev(ncol,nz) )

  u_flev = 0.  ;  v_flev = 0.  ;  u_sfc = 0.  ; v_sfc = 0.
!  t_flev = 270.
!  nbv_flev(:,:) = 2.6e-2
  do k=1, nz
  do l=1, ncol
    u_flev(l,k) = u_ct(l)
    v_flev(l,k) = v_ct(l)
    nbv_flev(l,k) = n_ct(l)
    rho_flev(l,k) = rho_ct(l)*exp(-(z_flev(l,k)-zcta(l))/7.e3)
    rho_dlev(l,k) = rho_ct(l)*exp(-(z_dlev(k)-zcta(l))/7.e3)
!    heat_flev(l,k) = (1./3600.)*( 1. - ((z_flev(l,k)-4.e3)/2.e3)**2 )
  enddo
  enddo
!  kcb(:) = 1  ;  kct(:) = 20
  f_cor(:) = 3.
  f_cor(:) = two_omega*sin(f_cor(:)*deg2rad)

file_o = './zzz.nc'

  ! u_sfc, v_sfc: allocate and specify, if possible


!  call args_sc05(ncol,nz,u_flev,v_flev,t_flev,nbv_flev,rho_flev,z_flev,  &
!                 heat_flev,kcb,kct)  ! opt: z_ref
!  ! kcta ; diag_znwcq

  ! u_ct, v_ct, u_cb, v_cb, t_ct
  ! zcta, zcba, n_q, n_ct, rho_ct, cqx, cqy, heatmax
  ! kcta


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


SUBROUTINE input_prof

  use netio

  implicit none

  character(len=128) ::  fdir

  fdir = '/data11/data-arch/ERA-I-nr/2005/11'

  call opennc(trim(fdir)//'/era-int_f.u.anal.00.ml.200511.nc',ncid)

  call closenc(ncid)

END subroutine input_prof

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


PROGRAM cgwp

  USE param_gwp
  USE switch_dump
  USE mflx_ctop_sc05
  USE prop_diss
  USE subr_common
  USE netio

  implicit none

  logical, parameter ::  l_sc05var_specified = .True.
  integer, parameter ::  nz = 60
  integer, parameter ::  nz_src = 30
  integer, parameter ::  nx_ei = 360
  integer, parameter ::  nx = 192, nk_fc = nx_ei/4
 
  integer ::  ncol, nvo

  real, dimension(nz) ::  p_dlev_ref, zp_dlev_ref, zp_flev_ref,          &
                          lnp_flev_ref

  real, dimension(nx,nz) ::  u_grd_dl, v_grd_dl, t_grd_dl
  real, dimension(nx,nz) ::  heat_grd_dl
 
  real, dimension(:,:), allocatable ::  u_fl_b, v_fl_b, nbv_fl_b,        &
                                        rho_fl_b, t_fl_b
  real, dimension(:,:), allocatable ::  heat_fl

  ! only for args_sc05
  real, dimension(:,:) , allocatable ::  z_fl_s, u_fl_s, v_fl_s, t_fl_s, &
                                         nbv_fl_s, rho_fl_s, heat_fl_s
  integer, dimension(:), allocatable ::  kcb, kct

  ! only for propdiss
  real, dimension(:), allocatable ::  f_cor

  ! only for calc_drag
  real, dimension(:,:), allocatable ::  p_dlev, lnp_flev

  character(len=128) ::  file_i, file_o
  integer            ::  i,k,l,iv, tmpi
  real               ::  tmp

  real, parameter ::  two_omega = 2.*7.292116e-5
  real, parameter ::  g = 9.80665
  real, parameter ::  rd = 287.05
  real, parameter ::  kappa = rd/1005.0
  real, parameter ::  n2bv_min = 1.e-6

  include 'c_math.inc'   ! deg2rad

  nc = 160  ! 30
  dc = 0.5  ! 2.
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

!-----------------------------------------------------------------------
!  READ BACKGROUND VARIABLES AND CONVECTIVE HEATING
!-----------------------------------------------------------------------

  call read_erai

  ! define vertical coordinates
  zp_flev_ref(1:nz-1) = 0.5*(zp_dlev_ref(1:nz-1) + zp_dlev_ref(2:nz))
  zp_flev_ref(nz) = 2.*zp_dlev_ref(nz) - zp_flev_ref(nz-1)
  lnp_flev_ref(:) = (-zp_flev_ref(:)/7.e3) + log(1.e5)

  allocate( u_fl_b(ncol,nz), v_fl_b(ncol,nz), nbv_fl_b(ncol,nz),         &
            rho_fl_b(ncol,nz), t_fl_b(ncol,nz) )

  do l=1, ncol
!    i = (l) ??
i=180
    u_fl_b(l,1:nz-1) = 0.5*(u_grd_dl(i,1:nz-1) + u_grd_dl(i,2:nz))
    v_fl_b(l,1:nz-1) = 0.5*(v_grd_dl(i,1:nz-1) + v_grd_dl(i,2:nz))
    t_fl_b(l,1:nz-1) = 0.5*(t_grd_dl(i,1:nz-1) + t_grd_dl(i,2:nz))
    u_fl_b(l,nz) = u_grd_dl(i,nz)
    v_fl_b(l,nz) = v_grd_dl(i,nz)
    t_fl_b(l,nz) = t_grd_dl(i,nz)
!p_coord+
    nbv_fl_b(l,1:nz-1) = 7.e3*(g*g)/(rd*t_fl_b(l,1:nz-1))*               &
        log( (t_grd_dl(i,2:nz)/t_grd_dl(i,1:nz-1))*                      &
             (p_dlev_ref(1:nz-1)/p_dlev_ref(2:nz))**kappa )/             &
        (zp_dlev_ref(2:nz) - zp_dlev_ref(1:nz-1))
!p_coord-
  enddo
  nbv_fl_b(:,nz) = nbv_fl_b(:,nz-1)
  do k=1, nz
  do l=1, ncol
    nbv_fl_b(l,k) = sqrt( max(n2bv_min,nbv_fl_b(l,k)) )
  enddo
  enddo

  do k=1, nz
    rho_fl_b(:,k) = exp(lnp_flev_ref(k))/(rd*t_fl_b(:,k))
  enddo

  ! u_sfc, v_sfc: allocate and specify, if possible

  allocate( heat_fl(ncol,nz) )
  heat_fl(:,:) = 0.
!  do l=1, ncol
!!    i = (l) ??
!i=180
!    heat_fl(l,1:nz-1) = 0.5*(q_grd_dl(i,1:nz-1) + q_grd_dl(i,2:nz))
!    heat_fl(l,nz) = q_grd_dl(i,nz)
!  enddo

!-----------------------------------------------------------------------
!  EXTRACT THE VARIABLES USED FOR CALCULATION OF SC05
!-----------------------------------------------------------------------

  SPEC_SC05:  IF ( .not. l_sc05var_specified ) then

  allocate( kcb(ncol), kct(ncol) )
  kcb = 1  ;  kct = nz_src

  allocate( z_fl_s(ncol,nz_src) )
  allocate( u_fl_s   (ncol,nz_src), v_fl_s   (ncol,nz_src),              &
            t_fl_s   (ncol,nz_src), nbv_fl_s (ncol,nz_src),              &
            rho_fl_s (ncol,nz_src), heat_fl_s(ncol,nz_src) )

!  z_fl_s   (:,:) = zp_flev_ref(1:nz_src)  ! z-coord
!  z_fl_s   (:,:) = z_fl_b   (:,1:nz_src)  ! p-coord
  u_fl_s   (:,:) = u_fl_b   (:,1:nz_src)
  v_fl_s   (:,:) = v_fl_b   (:,1:nz_src)
  t_fl_s   (:,:) = t_fl_b   (:,1:nz_src)
  nbv_fl_s (:,:) = nbv_fl_b (:,1:nz_src)
  rho_fl_s (:,:) = rho_fl_b (:,1:nz_src)
  heat_fl_s(:,:) = heat_fl  (:,1:nz_src)

  call args_sc05( ncol, nz_src, z_fl_s, u_fl_s, v_fl_s, t_fl_s,          &
                  nbv_fl_s, rho_fl_s, heat_fl_s, kcb, kct )  ! opt: z_ref
  ! OUT |  u_ct ; v_ct ; u_cb ; v_cb ; t_ct ; n_q ; n_ct ; rho_ct ;
  !     |  zcta ; zcba ; cqx ; cqy ; heatmax ; kcta
  !     |  (u_sfc ; v_sfc, if not allocated before)
  ! OUT |  diag_znwcq

  deallocate( z_fl_s )
  deallocate( u_fl_s, v_fl_s, t_fl_s, nbv_fl_s, rho_fl_s, heat_fl_s )
  deallocate( kcb, kct )

  ELSE

  call sc05vars

  if ( allocated(kcta) )  deallocate( kcta )
  allocate( kcta(ncol) )

  do l=1, ncol
!z_coord+
!    kcta(l) = minloc(abs(z_flev_ref(:) - zcta(l)),1)
!!UM    tmpi = minloc(abs(z_flev_ref(:) - zcba(l)),1)
!!UM    kcta(l) = tmpi + minloc(abs(z_flev_ref(tmpi+1:) - zcta(l)),1)
!z_coord-
!p_coord+
    tmp = 7.e3*log(1.e5/(rho_ct(l)*rd*t_ct(l)))
!p_coord-
    kcta(l) = minloc(abs(zp_flev_ref(:) - tmp),1)
    kcta(l) = max(3,kcta(l))
  enddo

  END IF  SPEC_SC05

  deallocate( heat_fl )

!-----------------------------------------------------------------------
!  SC05 CALCULATION
!-----------------------------------------------------------------------

  call calc_sc05(ncol)  ! opt: shear_ct
  ! IN  |  output variables from args_sc05:
  !     |  u_ct ; v_ct ; u_cb ; v_cb ; t_ct ; n_q ; n_ct ; rho_ct ;
  !     |  zcta ; zcba ; cqx ; cqy ; heatmax ;
  !     |  u_sfc ; v_sfc
  ! OUT |  mfs_ct

  ! u_sfc and v_sfc must be deallocated here to reset their values, if
  ! in a loop
  deallocate( u_sfc, v_sfc )

  deallocate( u_ct, v_ct, u_cb, v_cb, t_ct, n_q, n_ct, rho_ct,           &
              zcta, zcba, cqx, cqy, heatmax )

!-----------------------------------------------------------------------
!  OBTAIN MOMENTUM FLUX PROFILE ABOVE THE LAUNCH LEVEL
!-----------------------------------------------------------------------
 
  allocate( f_cor(ncol) )
  f_cor(:) = 3.
  f_cor(:) = two_omega*sin(f_cor(:)*deg2rad)

  call propdiss( ncol, nz, u_fl_b, v_fl_b, nbv_fl_b, rho_fl_b, f_cor,    &
                 kcta, mfs_ct )
  ! OUT |  mflx_ct_XXXX ; mflx_XXXX ; mf_pos ; mf_neg
  ! OUT |  diag_spec_ctop ; diag_spec
 
  deallocate( f_cor )
  deallocate( mfs_ct )
  deallocate( u_fl_b, v_fl_b, nbv_fl_b, rho_fl_b, t_fl_b )

!-----------------------------------------------------------------------
!  CALCULATE GRAVITY WAVE DRAG
!-----------------------------------------------------------------------
 
  if ( l_drag_u_o .or. l_drag_v_o ) then

    allocate( p_dlev(ncol,nz), lnp_flev(ncol,nz) )
    p_dlev  (:,:) = spread(p_dlev_ref  ,1,ncol)
    lnp_flev(:,:) = spread(lnp_flev_ref,1,ncol)

    call calc_drag_lnp(ncol,nz,lnp_flev,p_dlev,0)
    ! IN  |  mf_pos ; mf_neg
    ! OUT |  drag_u ; drag_v
 
!    deallocate( rho_dl_b )
    deallocate( p_dlev, lnp_flev )

  end if

  deallocate( kcta )

!-----------------------------------------------------------------------
!  
!-----------------------------------------------------------------------
 
  call put_vars_set

! DUMP

  file_o = './zzz.nc'

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nvo,set,'CGWP offline calculation')

! END

  call finalize

  STOP


CONTAINS


SUBROUTINE sc05vars

  integer ::  ncid

  if ( allocated(u_ct) )  deallocate( u_ct, v_ct, u_cb, v_cb, t_ct,      &
                                      n_q, n_ct, rho_ct, zcta, zcba,     &
                                      cqx, cqy, heatmax )
  allocate( u_ct(ncol), v_ct(ncol), u_cb(ncol), v_cb(ncol), t_ct(ncol),  &
            n_q(ncol), n_ct(ncol), rho_ct(ncol), zcta(ncol), zcba(ncol), &
            cqx(ncol), cqy(ncol), heatmax(ncol) )

  allocate( u_sfc(ncol), v_sfc(ncol) )

!  file_i = '../../../dat/L60CGW/dchm_pdf/'//  &
!           'uanuj.dchm-midlev_pdf.1979-2006.01-12.nc'
  file_i = '../../../dat/L60CGW/dchm_pdf/'//  &
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

END subroutine sc05vars

SUBROUTINE read_erai

  implicit none

  character(len=128) ::  fdir, fname_ei(3)

  real, dimension(nx,nz,3)    ::  var_grd
  real, dimension(nx)         ::  ix_grd
  real, dimension(nx_ei,nz,3) ::  var_ei
  real, dimension(nx_ei)      ::  ix_ei, lon_ei
  real, dimension(nx_ei,nz)   ::  tmp_ei
  real, dimension(0:nk_fc)    ::  cc, cs, ifc

  integer ::  ii, ik, ivi

  real, parameter ::  twopi = 6.283185 !3

  fdir = '../ERA-I_eq'

  fname_ei(1) = trim(fdir)//'/era-interim_eq_2005-11-01_u.dat'
  fname_ei(2) = trim(fdir)//'/era-interim_eq_2005-11-01_v.dat'
  fname_ei(3) = trim(fdir)//'/era-interim_eq_2005-11-01_t.dat'
  do ivi=1, 3
    call read_erai_lon_data(10,fname_ei(ivi),lon_ei,zp_dlev_ref,         &
                            var_ei(:,:,ivi),nx_ei,nz)
  enddo
  if (lon_ei(1) /= 0.) then
    if (lon_ei(1) == -180.) then
      do ivi=1, 3
        tmp_ei(:,:) = var_ei(:,:,ivi)
        var_ei(1:nx_ei/2,:,ivi) = tmp_ei(nx_ei/2+1:,:)
        var_ei(nx_ei/2+1:,:,ivi) = tmp_ei(1:nx_ei/2,:)
      enddo
      tmp_ei(:,1) = lon_ei(:)
      lon_ei(1:nx_ei/2) = tmp_ei(nx_ei/2+1:,1)
      lon_ei(nx_ei/2+1:) = tmp_ei(1:nx_ei/2,1)
    else
      print*, 'Check longitudes in ERA-I.'  ;  STOP
    end if
  end if

  do ii=1, nx_ei
    ix_ei(ii) = float(ii-1)
  enddo
  do ii=1, nx
    ix_grd(ii) = float(ii-1)
  enddo
  do ik=0, nk_fc
    ifc(ik) = float(ik)
  enddo

  do ivi=1, 3
  do k=1, nz
    do ik=0, nk_fc
      cc(ik) = sum(var_ei(:,k,ivi)*cos(ix_ei(:)*ifc(ik)/float(nx_ei)*twopi))
      cs(ik) = sum(var_ei(:,k,ivi)*sin(ix_ei(:)*ifc(ik)/float(nx_ei)*twopi))
    enddo
    cc(:) = cc(:)/float(nx_ei)
    cs(:) = cs(:)/float(nx_ei)
    do ii=1, nx
      var_grd(ii,k,ivi) = cc(0) + 2.*(                                   &
          sum(cc(1:)*cos(ix_grd(ii)*ifc(1:)/float(nx)*twopi)) +          &
          sum(cs(1:)*sin(ix_grd(ii)*ifc(1:)/float(nx)*twopi)) )
    enddo
  enddo
  enddo
 
  u_grd_dl(:,:) = var_grd(:,:,1)
  v_grd_dl(:,:) = var_grd(:,:,2)
  t_grd_dl(:,:) = var_grd(:,:,3)

  zp_dlev_ref(:) = zp_dlev_ref(:)*1.e3  ! [m]

  p_dlev_ref(:) = 1.e5*exp(-zp_dlev_ref(:)/7.e3)

END subroutine read_erai

SUBROUTINE put_vars_set

  character(len=32) ::  axis(4)
  integer           ::  ndim(4)

  real, dimension(ncol) ::  coln
  real, dimension(nz)   ::  z

  do l=1, ncol
    coln(l) = float(l)
  enddo
  z(:) = zp_dlev_ref(:)

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
 
  z(:) = zp_flev_ref(:)
  axis = (/'c_ph','dir ','z   ','case'/)
  ndim = (/nc*2+1,nphi*2,nz,ncol/)
  if ( l_spec_o ) then
    call defset(iv,'mflx_spec',diag_spec,axis,ndim,                      &
                c_phase,phi_deg2,z,coln)
  end if
 
END subroutine put_vars_set

SUBROUTINE finalize

  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo
  deallocate( set )

END subroutine finalize


END program cgwp


      subroutine read_erai_lon_data (nf,fname,lon,zkm,data,nlon,nlev) 
      character*80 fname,label1,label2,label3
      real lon(nlon),zkm(nlev),data(nlon,nlev)  
      if (nlon.ne.360) then 
         write(6,*) 'nlon must be 360' 
         stop
      endif 
      if (nlev.ne.60) then 
         write(6,*) 'nlev must be 60' 
         stop
      endif 
      open (nf,file=fname,status='old')
      read (nf,5000) label1
      read (nf,5010) lon
      read (nf,5000) label2
      read (nf,5020) zkm
      read (nf,5000) label3
      do l=1,nlev 
         read (nf,5020) (data(i,l),i=1,nlon)
      enddo
      close (nf)         
!      write(6,*) label1
!      write(6,*) lon
!      write(6,*) label2
!      write(6,*) zkm
!      write(6,*) label3
!      do l=1,nlev 
!         write(6,*) (data(i,l),i=1,nlon)
!      enddo
 5000 format(a80)
 5010 format(360f8.2)
 5020 format(360f13.8)
      return
      end subroutine


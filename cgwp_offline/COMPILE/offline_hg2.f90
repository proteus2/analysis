PROGRAM cgwp

  USE param_gwp
  USE switch_dump
  USE mflx_ctop_sc05
  USE prop_diss
  USE subr_common
  USE hadgem
  USE netio
  USE zonal_average

  implicit none

  real    ::  nc_dc(2)
  integer ::  phi0_dphi(2), nz_src, nx_b0, nz_b, nk_b
 
  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM_S/ NC_DC, PHI0_DPHI, CFACTOR, NZ_SRC, LAT_RNG
  namelist /PARAM_B/ BETA_WM, NX_B0, NZ_B, NK_B
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, VAR_I_NAME, FILE_I_HEAD2, FILE_I_FORM2, &
                    FILE_I_XXXX2, VAR_I_NAME2, FILE_ALT, VAR_ALT, FILE_O

  integer, parameter ::  nzu = 33

  integer ::  imon, ihour, i_time, iy2(2), ny2, ncol, ncolm, day1_0
  integer ::  nvo

  real, dimension(:), allocatable ::  p_dlev_ref, zp_dlev_ref,           &
                                      zp_flev_ref, lnp_flev_ref

  real, dimension(:)  , allocatable ::  p_grd, zp_grd, lat_grd
  real, dimension(:,:), allocatable ::  u_grd, v_grd, t_grd
 
  real, dimension(:,:,:), allocatable ::  var_sc, var_u, var_v, ind_mc

  real   , dimension(:,:), allocatable ::  u_fl_b, v_fl_b, nbv_fl_b,     &
                                           rho_fl_b, t_fl_b
  integer, dimension(:,:), allocatable ::  ij_col

  ! only for propdiss
  real, dimension(:), allocatable ::  f_cor

  ! only for calc_drag
  real, dimension(:,:), allocatable ::  p_dlev, lnp_flev

  integer ::  i,j,k,l, tmpi
  real    ::  tmp, wgt_t, wgt_acc

  real, parameter ::  two_omega = 2.*7.292116e-5
  real, parameter ::  g = 9.80665
  real, parameter ::  rd = 287.05
  real, parameter ::  cp = 1005.0
  real, parameter ::  n2bv_min = 1.e-6

  include 'c_math.inc'   ! deg2rad

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM_S)  ;  read(10, PARAM_B)
  read(10, FILEIO)
  close(10)

!=======================================================================
!  SET PARAMETERS, GET AXES, AND INITIALIZE ARRAYS
!=======================================================================

  year = yyyy
  mon  = mm(1)

  call initialize

  i_time = 0

  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------
  if (opt_30d == 0)  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
  i_time = i_time + 1

  day_from_ref = get_dayfromref(year,mon,date,hour)

  ! problem found in the data for 00 UTC 1 each month
  if ( date == 1 .and. hour == 0 ) then
    hour = hour + 24/nhour
    i_time = i_time - 1
    CYCLE
  end if

!=======================================================================
!  READ BACKGROUND FLOW VARIABLES (AND CONVECTIVE HEATING)
!=======================================================================

  call read_erai

!=======================================================================
!  EXTRACT VARIABLES USED FOR CALCULATION OF SC05
!=======================================================================
 
!-------------------------------------------------
!
! ! Opt 1) extract from 3-D fields
!
!  call args_sc05( ncol, nz_src, z_fl_s, u_fl_s, v_fl_s, t_fl_s,          &
!                  nbv_fl_s, rho_fl_s, heat_fl_s, kcb, kct )  ! opt: z_ref
!  ! OUT |  u_ct ; v_ct ; u_cb ; v_cb ; t_ct ; n_q ; n_ct ; rho_ct ;
!  !     |  zcta ; zcba ; cqx ; cqy ; heatmax ; kcta
!  !     |  (u_sfc ; v_sfc, if not allocated before)
!  ! OUT |  diag_znwcq
!
!-------------------------------------------------

 ! Opt 2) get from UM results

  call switch_para_in

  call sc05vars_um

  call switch_para_in

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

!-------------------------------------------------

!=======================================================================
!  SC05 CALCULATION
!=======================================================================

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

!=======================================================================
!  RE-GRID BACKGROUND FLOW PROFILES
!=======================================================================

  allocate( u_fl_b(ncol,nz_b), v_fl_b(ncol,nz_b), nbv_fl_b(ncol,nz_b),   &
            rho_fl_b(ncol,nz_b), t_fl_b(ncol,nz_b) )
  allocate( f_cor(ncol) )

  do l=1, ncol
    i = ij_col(l,1)
    j = ij_col(l,2)
    u_fl_b(l,1:nz_b-1) = 0.5*(u_grd(i,1:nz_b-1) + u_grd(i,2:nz_b))
    v_fl_b(l,1:nz_b-1) = 0.5*(v_grd(i,1:nz_b-1) + v_grd(i,2:nz_b))
    t_fl_b(l,1:nz_b-1) = 0.5*(t_grd(i,1:nz_b-1) + t_grd(i,2:nz_b))
    u_fl_b(l,nz_b) = u_grd(i,nz_b)
    v_fl_b(l,nz_b) = v_grd(i,nz_b)
    t_fl_b(l,nz_b) = t_grd(i,nz_b)
!p_coord+
    do k=1, nz_b-1
      nbv_fl_b(l,k) = (g*g)/t_fl_b(l,k)*( 1./cp + 7.e3/(rd*t_fl_b(l,k))* &
          ( (t_grd(i,k+1) - t_grd(i,k))/  &
            (zp_dlev_ref(k+1) - zp_dlev_ref(k)) ) )
    enddo
!p_coord-
    f_cor(l) = two_omega*sin(lat_grd(j)*deg2rad)
  enddo
  nbv_fl_b(:,nz_b) = nbv_fl_b(:,nz_b-1)
  do k=1, nz_b
  do l=1, ncol
    nbv_fl_b(l,k) = sqrt( max(n2bv_min,nbv_fl_b(l,k)) )
  enddo
  enddo

  do k=1, nz_b
    rho_fl_b(:,k) = exp(lnp_flev_ref(k))/(rd*t_fl_b(:,k))
  enddo

!=======================================================================
!  OBTAIN MOMENTUM FLUX PROFILE ABOVE THE LAUNCH LEVEL
!=======================================================================
 
  call propdiss( ncol, nz_b, u_fl_b, v_fl_b, nbv_fl_b, rho_fl_b, f_cor,  &
                 kcta, mfs_ct )
  ! OUT |  mflx_ct_XXXX ; mflx_XXXX ; mf_pos ; mf_neg
  ! OUT |  diag_spec_ctop ; diag_spec
 
  deallocate( f_cor )
  deallocate( mfs_ct )
  deallocate( u_fl_b, v_fl_b, nbv_fl_b, rho_fl_b, t_fl_b )

!=======================================================================
!  CALCULATE GRAVITY WAVE DRAG
!=======================================================================
 
  if ( l_drag_u_o .or. l_drag_v_o ) then

    allocate( p_dlev(ncol,nz_b), lnp_flev(ncol,nz_b) )
    p_dlev  (:,:) = spread(p_dlev_ref  ,1,ncol)
    lnp_flev(:,:) = spread(lnp_flev_ref,1,ncol)

    call calc_drag_lnp(ncol,nz_b,lnp_flev,p_dlev,0)
    ! IN  |  mf_pos ; mf_neg
    ! OUT |  drag_u ; drag_v
 
!    deallocate( rho_dl_b )
    deallocate( p_dlev, lnp_flev )

  end if

  deallocate( kcta )

!=======================================================================
!  DAILY, ZONAL AVERAGE
!=======================================================================

  wgt_t = 1.

  wgt_acc = wgt_acc + wgt_t

  call zonal_avg(ij_col(1:ncol,2),nx, wgt_t)
 
  hour = hour + 24/nhour

  !---------------------------------------------------------------------
  ENDDO  L_HOUR

  !---------------------------------------------------------------------
  ENDDO  L_DATE

  mon = mon + 1
  if (mon == 13) then
    year = year + 1  ;  mon = 1
  end if
  !---------------------------------------------------------------------
  ENDDO  L_MON

  nt = i_time

!  call put_vars_set
  call put_vars_zm_set

  do iv=1, nvo
    set(iv)%var_out(:,:,:,:) = set(iv)%var_out(:,:,:,:)/wgt_acc
  enddo

  deallocate( var_sc, ind_mc, var_u, var_v )

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nvo,set,'CGWP offline calculation')

! END

  call finalize

  STOP


CONTAINS


SUBROUTINE initialize

  nc = int(nc_dc(1))  ;  dc = nc_dc(2)

  nphi = 180/phi0_dphi(2)
  allocate( phi_deg(nphi) )
  do i=1, nphi
    phi_deg(i) = float(phi0_dphi(1) + (i-1)*phi0_dphi(2))
  enddo

  call set_spec_param

  call get_wm_hg2cgwp
  if (beta_wm <= 0) then
    beta_wm = beta_wm0
    write(6,*) ' beta_WM is set to the default :', beta_wm0
  end if

  call switch_defaults

!  l_spec_on = .False.

  call get_nv_output(nvo)
  allocate( set(nvo) )

  day1_0 = day1

  nmon = mm(2)  ;  nhour = hh(2)

!  problem found in the data for 00 UTC 1 each month
!  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  ndate = 30  ;  date = 2  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 3  ! for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  call get_iouter(lat,lat_rng, iy2)

  ny2 = iy2(2) - iy2(1) + 1

  x1_i = 1   ;  y1_i = iy2(1)   ;  z1_i = 1       ! for getalt
  nx_i = nx  ;  ny_i = ny2      ;  nz_i = nzu-1
  call getalt
  do k=1, nzu-1
    z_th(:,:,k) = z_th(:,:,k) - z_th(:,:,0)
  enddo
  z_th(:,:,0) = 0.

  ! define vertical coordinates
  allocate( p_dlev_ref(nz_b), zp_dlev_ref(nz_b), zp_flev_ref(nz_b),      &
            lnp_flev_ref(nz_b) )

  call switch_para_in

  allocate( p_grd(nz_b), zp_grd(nz_b) )
  call read_erai('coord')

  p_dlev_ref (:) = p_grd (:)
  zp_dlev_ref(:) = zp_grd(:)

  zp_flev_ref(1:nz_b-1) = 0.5*(zp_dlev_ref(1:nz_b-1) +                   &
                               zp_dlev_ref(2:nz_b  ))
  zp_flev_ref(nz_b) = 2.*zp_dlev_ref(nz_b) - zp_flev_ref(nz_b-1)
  lnp_flev_ref(:) = (-zp_flev_ref(:)/7.e3) + log(1.e5)

  allocate( lat_grd(ny2) )
!  lat_grd(:) = lat(iy2(1):iy2(2))
  lat_grd(:) = 0.

  allocate( u_grd(nx,nz_b), v_grd(nx,nz_b), t_grd(nx,nz_b) )

  ncolm = nx*ny2

  allocate( ij_col(ncolm, 2) )
  allocate( var_sc(nx,ny2,15), ind_mc(nx,ny2,1) )
  allocate( var_u(nx,ny2,nzu), var_v(nx,ny2,nzu) )

  call alloc_zonal_avg(ny2,nz_b,nc,nphi)

  wgt_acc = 0.

END subroutine initialize

SUBROUTINE read_erai(c_coord)

  implicit none

  character(len=*), intent(in), optional ::  c_coord

  integer, parameter                ::  nvi_b = 3
  real, dimension(nx,nz_b,nvi_b)    ::  var_grd
  real, dimension(nx)               ::  ix_grd
  real, dimension(nx_b0,nz_b,nvi_b) ::  var_ra
  real, dimension(nx_b0)            ::  ix_ra, lon_ra
  real, dimension(nx_b0,nz_b)       ::  tmp_ra
  real, dimension(0:nk_b)           ::  cc, cs, ifc

  integer ::  ii, ik

  real, parameter ::  twopi = 6.283185 !3

  day1 = -999

  ex0 = .TRUE.
  do iv_i=1, nvi_b  ! for get_ifilename
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  iv_i = 1
  call read_erai_lon_data(20,file_i(iv_i),lon_ra,zp_grd,                 &
                          var_ra(:,:,iv_i),nx_b0,nz_b)
 
  zp_grd(:) = zp_grd(:)*1.e3  ! [m]

  p_grd(:) = 1.e5*exp(-zp_grd(:)/7.e3)

  if ( present(c_coord) ) then
    day1 = day1_0  ;  RETURN
  end if

  do iv_i=2, nvi_b
    call read_erai_lon_data(20,file_i(iv_i),lon_ra,zp_grd,               &
                            var_ra(:,:,iv_i),nx_b0,nz_b)
  enddo
  if (lon_ra(1) /= 0.) then
    if (lon_ra(1) == -180.) then
      do iv_i=1, nvi_b
        tmp_ra(:,:) = var_ra(:,:,iv_i)
        var_ra(1:nx_b0/2,:,iv_i) = tmp_ra(nx_b0/2+1:,:)
        var_ra(nx_b0/2+1:,:,iv_i) = tmp_ra(1:nx_b0/2,:)
      enddo
      tmp_ra(:,1) = lon_ra(:)
      lon_ra(1:nx_b0/2) = tmp_ra(nx_b0/2+1:,1)
      lon_ra(nx_b0/2+1:) = tmp_ra(1:nx_b0/2,1)
    else
      print*, 'Check longitudes in ERA-I.'  ;  STOP
    end if
  end if

  ix_ra(:) = (/ ( float(ii-1), ii=1, nx_b0 ) /)
  ix_grd(:) = (/ ( float(ii-1), ii=1, nx ) /)
  ifc(:) = (/ ( float(ik), ik=0, nk_b ) /)

  if ( nx == nx_b0 .and. ( nk_b < 0 .or. nk_b >= nx_b0/2 ) ) then
    var_grd(:,:,:) = var_ra(:,:,:)
  else
    do iv_i=1, nvi_b
    do k=1, nz_b
      do ik=0, nk_b
        cc(ik) = sum(var_ra(:,k,iv_i)*cos(ix_ra(:)*ifc(ik)/float(nx_b0)*twopi))
        cs(ik) = sum(var_ra(:,k,iv_i)*sin(ix_ra(:)*ifc(ik)/float(nx_b0)*twopi))
      enddo
      cc(:) = cc(:)/float(nx_b0)
      cs(:) = cs(:)/float(nx_b0)
      do ii=1, nx
        var_grd(ii,k,iv_i) = cc(0) + 2.*(                                &
            sum(cc(1:)*cos(ix_grd(ii)*ifc(1:)/float(nx)*twopi)) +        &
            sum(cs(1:)*sin(ix_grd(ii)*ifc(1:)/float(nx)*twopi)) )
      enddo
    enddo
    enddo
  end if

  u_grd(:,:) = var_grd(:,:,1)
  v_grd(:,:) = var_grd(:,:,2)
  t_grd(:,:) = var_grd(:,:,3)

  day1 = day1_0

END subroutine read_erai

SUBROUTINE sc05vars_um

  integer, parameter              ::  nvi_sc = 15
  real, dimension(nx,nz_b,nvi_sc) ::  var_grd

  character(len=128) ::  file_i1
  integer            ::  ncid

  ex0 = .TRUE.
  do iv_i=1, 5  ! for get_ifilename
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  ! read var.
  print*, trim(file_i(1))

  iv_i = 1  ;  var_sc(:,:,1:1) = get_ivara3d(1,nx,iy2(1),ny2,1,1)

  iv_i = 2  ;  var_sc(:,:,2:11) = get_ivara3d(1,nx,iy2(1),ny2,1,10)

  iv_i = 3  ;  var_u(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,1,nzu)
  iv_i = 4  ;  var_v(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,1,nzu)

  iv_i = 5  ;  ind_mc(:,:,1:1) = get_ivara3d(1,nx,iy2(1),ny2,1,1)

  print*, 'time index :', it_i(1)

  l = 0
  do j=1, ny2
  do i=1, nx
    if (var_sc(i,j,1) == missv)  CYCLE
    if (var_sc(i,j,2) == var_sc(i,j,3))  CYCLE  ! happens rarely !
!    if (ind_mc(i,j,1) /= 0.)  CYCLE
    if ( ind_mc(i,j,1) /= 0. .and. var_sc(i,j,2) .gt. 4.e3 )  CYCLE
    l = l + 1
    ij_col(l,1) = i
    ij_col(l,2) = j
  enddo
  enddo
  ncol = l

  if ( allocated(u_ct) )  deallocate( u_ct, v_ct, u_cb, v_cb, t_ct,      &
                                      n_q, n_ct, rho_ct, zcta, zcba,     &
                                      cqx, cqy, heatmax )
  allocate( u_ct(ncol), v_ct(ncol), u_cb(ncol), v_cb(ncol), t_ct(ncol),  &
            n_q(ncol), n_ct(ncol), rho_ct(ncol), zcta(ncol), zcba(ncol), &
            cqx(ncol), cqy(ncol), heatmax(ncol) )

  allocate( u_sfc(ncol), v_sfc(ncol) )

  do l=1, ncol
    i = ij_col(l,1)
    j = ij_col(l,2)
    heatmax(l) = var_sc(i,j,1 )
    zcba   (l) = var_sc(i,j,2 )
    zcta   (l) = var_sc(i,j,3 )
    rho_ct (l) = var_sc(i,j,4 )
    n_q    (l) = var_sc(i,j,5 )
    n_ct   (l) = var_sc(i,j,6 )
    t_ct   (l) = var_sc(i,j,7 )
    cqx    (l) = var_sc(i,j,8 )
    cqy    (l) = var_sc(i,j,9 )
    u_ct   (l) = var_sc(i,j,10)
    v_ct   (l) = var_sc(i,j,11)
    u_sfc  (l) = var_u (i,j,1 )
    v_sfc  (l) = var_v (i,j,1 )
    k = minloc(abs(zcba(l) - z_th(i,j,1:)),1)
    u_cb(l) = 0.5*(var_u(i,j,k) + var_u(i,j,k+1))
    v_cb(l) = 0.5*(var_v(i,j,k) + var_v(i,j,k+1))
  enddo

END subroutine sc05vars_um

SUBROUTINE put_vars_set

  real, dimension(ncol) ::  coln
  real, dimension(nz_b) ::  zo

  coln(:) = (/ ( float(l), l=1, ncol ) /)

  iv = 0

  axisname1d_put = 'case'
  ndim1d_put = ncol
  if ( l_mflx_u_ctop_o ) then
    call putset(iv,'mflx_ct_east',mflx_ct_east, coln)
    call putset(iv,'mflx_ct_west',mflx_ct_west, coln)
  end if
  if ( l_mflx_v_ctop_o ) then
    call putset(iv,'mflx_ct_north',mflx_ct_north, coln)
    call putset(iv,'mflx_ct_south',mflx_ct_south, coln)
  end if

  axisname2d_put = (/'z_f ','case'/)
  ndim2d_put = (/nz_b,ncol/)
  zo(:) = zp_flev_ref(:)
  if ( l_mflx_u_o ) then
    call putset(iv,'mflx_east',transpose(mflx_east), zo,coln)
    call putset(iv,'mflx_west',transpose(mflx_west), zo,coln)
  end if
  if ( l_mflx_v_o ) then
    call putset(iv,'mflx_north',transpose(mflx_north), zo,coln)
    call putset(iv,'mflx_south',transpose(mflx_south), zo,coln)
  end if

  axisname2d_put = (/'z_d ','case'/)
  ndim2d_put = (/nz_b,ncol/)
  zo(:) = zp_dlev_ref(:)
  if ( l_drag_u_o ) then
    call putset(iv,'drag_u',transpose(drag_u), zo,coln)
  end if
  if ( l_drag_v_o ) then
    call putset(iv,'drag_v',transpose(drag_v), zo,coln)
  end if

  axisname3d_put = (/'c_ph','dir ','case'/)
  ndim3d_put = (/nc*2+1,nphi*2,ncol/)
  if ( l_spec_ctop_o ) then
    call putset(iv,'mflx_ct_spec',diag_spec_ctop, c_phase,phi_deg2,coln)
  end if
 
  axisname4d_put = (/'c_ph','dir ','z_f ','case'/)
  ndim4d_put = (/nc*2+1,nphi*2,nz_b,ncol/)
  zo(:) = zp_flev_ref(:)
  if ( l_spec_o ) then
    call putset(iv,'mflx_spec',diag_spec, c_phase,phi_deg2,zo,coln)
  end if
 
END subroutine put_vars_set

SUBROUTINE put_vars_zm_set

  real, dimension(ny2 ) ::  lato
  real, dimension(nz_b) ::  zo

  lato(:) = lat(iy2(1):iy2(2))

  iv = 0

  axisname1d_put = 'lat'
  ndim1d_put = ny2
  if ( l_mflx_u_ctop_o ) then
    call putset(iv,'mflx_ct_east',zm_mflx_ct_east, lato)
    call putset(iv,'mflx_ct_west',zm_mflx_ct_west, lato)
  end if
  if ( l_mflx_v_ctop_o ) then
    call putset(iv,'mflx_ct_north',zm_mflx_ct_north, lato)
    call putset(iv,'mflx_ct_south',zm_mflx_ct_south, lato)
  end if

  axisname2d_put = (/'lat','z_f'/)
  ndim2d_put = (/ny2,nz_b/)
  zo(:) = zp_flev_ref(:)
  if ( l_mflx_u_o ) then
    call putset(iv,'mflx_east',zm_mflx_east, lato,zo)
    call putset(iv,'mflx_west',zm_mflx_west, lato,zo)
  end if
  if ( l_mflx_v_o ) then
    call putset(iv,'mflx_north',zm_mflx_north, lato,zo)
    call putset(iv,'mflx_south',zm_mflx_south, lato,zo)
  end if

  axisname2d_put = (/'lat','z_d'/)
  ndim2d_put = (/ny2,nz_b/)
  zo(:) = zp_dlev_ref(:)
  if ( l_drag_u_o ) then
    call putset(iv,'drag_u',zm_drag_u, lato,zo)
  end if
  if ( l_drag_v_o ) then
    call putset(iv,'drag_v',zm_drag_v, lato,zo)
  end if

  axisname3d_put = (/'c_ph','dir ','lat '/)
  ndim3d_put = (/nc*2+1,nphi*2,ny2/)
  if ( l_spec_ctop_o ) then
    call putset(iv,'mflx_ct_spec',zm_diag_spec_ctop,                     &
                c_phase,phi_deg2,lato)
  end if
 
  axisname4d_put = (/'c_ph','dir ','z_f ','lat '/)
  ndim4d_put = (/nc*2+1,nphi*2,nz_b,ny2/)
  zo(:) = zp_flev_ref(:)
  if ( l_spec_o ) then
    call putset(iv,'mflx_spec',zm_diag_spec, c_phase,phi_deg2,zo,lato)
  end if
 
END subroutine put_vars_zm_set

SUBROUTINE finalize

  deallocate( p_dlev_ref, zp_dlev_ref, zp_flev_ref, lnp_flev_ref )
  deallocate( p_grd, zp_grd, u_grd, v_grd, t_grd )
  deallocate( lat_grd )
  deallocate( ij_col )
  call dealloc_zonal_avg
  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo
  deallocate( set )

END subroutine finalize


END program cgwp


PROGRAM QBOI_DYN_UM

  use hadgem
  use netio

  implicit none

  include 'c_math.inc'
  include 'c_phys.inc'

  integer, parameter ::  nv = 12
  integer, parameter ::  k_start = 4  ! > 1; max(k_t1) = 5 for h1 = 200 m
  real,    parameter ::  h1_realistic_temp = 200.  ! [m]

  integer ::  days_avrg
  real    ::  p_intp(999)

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ P_INTP, LAT_RNG, DAYS_AVRG, MISSV
  namelist /FILEIO/ DAY1, NDAY_I, FID, FILE_I_HEAD, FILE_I_FORM,         &
                    FILE_I_XXXX, VAR_I_NAME, FILE_ALT, VAR_ALT, FILE_O

  integer ::  nv_in, nv_intp, npi
  integer ::  imon, ihour, i_time, i_time_last
  integer ::  iy2(2), iz2(2), ny2, nz2, ipi1, ipi2
  integer ::  i,j,k, kpi, k2, k0
  integer ::  iv_t, iv_u, iv_v, iv_w                            ! 1--4
  integer ::  iv_vu, iv_dtdz, iv_dudz                           ! 5--20
  integer ::  iv_p0l, iv_t0l, iv_w0l                            ! 21--30
  integer ::  ivm_dvtisdz, ivm_dudz, ivm_fhat                   ! 31--40
  integer ::  iv_vs, iv_ws, iv_fy, iv_fz, iv_epd, iv_psis, iv_dudt,      &
              iv_um, iv_tm, iv_advy, iv_advz, iv_dyn            ! 51--99
  real    ::  hlnp1_min, hlnp9_max, c_ghri, c_lrgi, ntavg
  logical ::  l_w_in, l_thlev(20), l_belowgrd, l_abovetop, l_k_t1
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:),  allocatable ::  var_i, var_p, var4d
  real, dimension(:,:,:),    allocatable ::  var3d
  real, dimension(:,:,:),    allocatable ::  hlnp_r, inv_dzp_r, p_r
  real, dimension(:,:,:),    allocatable ::  hlnp_t, inv_dzp_t, p_t
  real, dimension(:,:,:),    allocatable ::  tmp3d, var_i0l
  real, dimension(:,:,:),    allocatable ::  coef_rho1, coef_rho2
  real, dimension(:,:),      allocatable ::  coslat, f, rho0, ritanlat
  real, dimension(:,:),      allocatable ::  hlnpi_l2
  real, dimension(:),        allocatable ::  hlnpi, inv_2dy, dzp_deriv
  integer, dimension(:,:,:), allocatable ::  k_low
  integer, dimension(:,:),   allocatable ::  k_t1, k_t2
  integer, dimension(:),     allocatable ::  iv_avg, iv_dif

  type(vset) ::  set(nv)

  real, dimension(2), parameter ::  c_sgn2 = (/-1.,1./)

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  year = yyyy
  mon  = mm(1)

  call initialize

  i_time_last = 0  ;  i_time = 0

  L_MON:  DO imon=1, nmon+1
  !---------------------------------------------------------------------
  if (opt_30d == 0)  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------
  if (days_avrg /= 0) then
    if ( mod(date,days_avrg) == 1 .or. days_avrg == 1 )                  &
       i_time = i_time + 1
  end if

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
  if (days_avrg == 0)  i_time = i_time + 1

  day_from_ref = get_dayfromref(year,mon,date,hour)

! READ VARIABLES

  call get_vars

  call calc_req_vars

! INTERPOLATION 1 and COVARIANCE SHEAR TERMS

  var3d(:,:,:) = 0.  ;  tmp3d(:,:,:) = 0.

  ! indices for tmp3d: ivm_*
  ivm_dvtisdz = 31
  ivm_dudz = 32  ;  ivm_fhat = 33

  ! z = zi -+ dz/2
  do k2=1, 2

    call vintp_z2p(3,(/iv_u,iv_v,iv_dtdz/),hlnp_r,inv_dzp_r,hlnpi_l2(:,k2))
    call vintp_z2p(2,(/iv_t,iv_w/),hlnp_t,inv_dzp_t,hlnpi_l2(:,k2))
 
    var3d(:,:,iv_um) = sum(var_p(:,:,:,iv_u), dim=1)/float(nx)
    var3d(:,:,iv_vs) = sum(var_p(:,:,:,iv_v), dim=1)/float(nx)
    var3d(:,:,iv_ws) = sum(var_p(:,:,:,iv_w), dim=1)/float(nx)
    var3d(:,:,iv_tm) = sum(var_p(:,:,:,iv_t), dim=1)/float(nx)
    do k=1, npi
    do j=1, ny2
      var_p(:,j,k,iv_u) = var_p(:,j,k,iv_u) - var3d(j,k,iv_um)
      var_p(:,j,k,iv_v) = var_p(:,j,k,iv_v) - var3d(j,k,iv_vs)
      var_p(:,j,k,iv_w) = var_p(:,j,k,iv_w) - var3d(j,k,iv_ws)
      var_p(:,j,k,iv_t) = var_p(:,j,k,iv_t) - var3d(j,k,iv_tm)
    enddo
    enddo

    tmp3d(:,:,ivm_dvtisdz) = tmp3d(:,:,ivm_dvtisdz) + c_sgn2(k2)*        &
        ( sum(var_p(:,:,:,iv_v)*var_p(:,:,:,iv_t), dim=1)/               &
          ( sum(var_p(:,:,:,iv_dtdz), dim=1) +                           &
            var3d(:,:,iv_tm)*((kappa/h_scale)*float(nx)) ) )
    var3d(:,:,iv_epd) = var3d(:,:,iv_epd) + c_sgn2(k2)*                  &
        ( sum(var_p(:,:,:,iv_w)*var_p(:,:,:,iv_u), dim=1)/float(nx) )

  enddo

  do k=1, npi
    ! d([v'T']/S)/dz,  used for v_res, epd
    tmp3d(:,k,ivm_dvtisdz) = tmp3d(:,k,ivm_dvtisdz)/dzp_deriv(k)
    ! iv_epd: -d[w'u']/dz,  used for epd
    var3d(:,k,iv_epd) = var3d(:,k,iv_epd)/(dzp_deriv(k)*(-1.))
  enddo

! INTERPOLATION 2 and TEM DIAGNOSTICS CALCULATION

  ! z = zi
  call vintp_z2p(4,(/iv_u,iv_v,iv_vu,iv_dtdz/),hlnp_r,inv_dzp_r,hlnpi)
  call vintp_z2p(3,(/iv_t,iv_w,iv_dudz/),hlnp_t,inv_dzp_t,hlnpi)

  var3d(:,:,iv_um) = sum(var_p(:,:,:,iv_u), dim=1)/float(nx)
  var3d(:,:,iv_vs) = sum(var_p(:,:,:,iv_v), dim=1)/float(nx)
  var3d(:,:,iv_ws) = sum(var_p(:,:,:,iv_w), dim=1)/float(nx)
  var3d(:,:,iv_tm) = sum(var_p(:,:,:,iv_t), dim=1)/float(nx)
  do k=1, npi
  do j=1, ny2
    var_p(:,j,k,iv_u) = var_p(:,j,k,iv_u) - var3d(j,k,iv_um)
    var_p(:,j,k,iv_v) = var_p(:,j,k,iv_v) - var3d(j,k,iv_vs)
    var_p(:,j,k,iv_w) = var_p(:,j,k,iv_w) - var3d(j,k,iv_ws)
    var_p(:,j,k,iv_t) = var_p(:,j,k,iv_t) - var3d(j,k,iv_tm)
  enddo
  enddo

  tmp3d(:,:,ivm_dudz) = sum(var_p(:,:,:,iv_dudz), dim=1)/float(nx)

  ! iv_psis: [v'T']/S,  used for psi_res, v_res, w_res, f_y, f_z, epd
  var3d(:,:,iv_psis) = sum(var_p(:,:,:,iv_v)*var_p(:,:,:,iv_t), dim=1)/  &
      ( sum(var_p(:,:,:,iv_dtdz), dim=1) +                               &
        var3d(:,:,iv_tm)*((kappa/h_scale)*float(nx)) )

  ! iv_fz: -[w'u'],  used for f_z, epd
  var3d(:,:,iv_fz) = sum(var_p(:,:,:,iv_w)*var_p(:,:,:,iv_u), dim=1)/    &
                     (float(nx)*(-1.))

  ! iv_fy: -[v'u'],  used for f_y, epd
  var3d(:,:,iv_fy) = ( sum(var_p(:,:,:,iv_vu), dim=1)/float(nx) -        &
                       var3d(:,:,iv_vs)*var3d(:,:,iv_um) )*(-1.)

  ! v_res,  used for adv_y, u_force_dyn
  var3d(:,:,iv_vs) = var3d(:,:,iv_vs) + ( var3d(:,:,iv_psis)/h_scale -   &
                     tmp3d(:,:,ivm_dvtisdz) )

  ! w_res,  used for adv_z
  var3d(:,:,iv_ws) = var3d(:,:,iv_ws) + ( grady(var3d(:,:,iv_psis)) -    &
                     var3d(:,:,iv_psis)*ritanlat(:,:) )

  ! iv_fy: f_y/(rho0*a*coslat),  used for f_y, epd
  var3d(:,:,iv_fy) = var3d(:,:,iv_fy) + var3d(:,:,iv_psis)*              &
                                        tmp3d(:,:,ivm_dudz)

  ! f_hat,  used for f_z, epd, adv_y, u_force_dyn
  tmp3d(:,:,ivm_fhat) = f(:,:) - ( grady(var3d(:,:,iv_um)) -             &
                        var3d(:,:,iv_um)*ritanlat(:,:) )

  ! iv_fz: f_z/(rho0*a*coslat),  used for f_z, epd
  var3d(:,:,iv_fz) = var3d(:,:,iv_fz) + var3d(:,:,iv_psis)*              &
                                        tmp3d(:,:,ivm_fhat)

  ! epd  [ = ( epd_z ) + ( epd_y ) ],  used for u_force_dyn
  var3d(:,:,iv_epd) =                                                    &
      ( var3d(:,:,iv_epd) +                                              &
        ( tmp3d(:,:,ivm_dvtisdz)*tmp3d(:,:,ivm_fhat) +                   &
          var3d(:,:,iv_psis)*( tmp3d(:,:,ivm_dudz)*ritanlat(:,:) -       &
                               grady(tmp3d(:,:,ivm_dudz)) )              &
        ) - var3d(:,:,iv_fz)/h_scale ) +                                 &
      ( grady(var3d(:,:,iv_fy)) - 2.*var3d(:,:,iv_fy)*ritanlat(:,:) )
 
  ! f_y, f_z
  var3d(:,:,iv_fy) = var3d(:,:,iv_fy)*(rho0(:,:)*r_earth*coslat(:,:))
  var3d(:,:,iv_fz) = var3d(:,:,iv_fz)*(rho0(:,:)*r_earth*coslat(:,:))

  ! psi_res
  var3d(:,:,iv_psis) = var3d(:,:,iv_psis)*(rho0(:,:)*coslat(:,:))

  ! adv_z,  used for u_force_dyn
  var3d(:,:,iv_advz) = -var3d(:,:,iv_ws)*tmp3d(:,:,ivm_dudz)

  ! adv_y (including the curvature term)
  var3d(:,:,iv_advy) = var3d(:,:,iv_vs)*(tmp3d(:,:,ivm_fhat) - f(:,:))

  ! u_force_dyn
  var3d(:,:,iv_dyn) = var3d(:,:,iv_epd) + var3d(:,:,iv_advz) +           &
                      var3d(:,:,iv_vs)*(tmp3d(:,:,ivm_fhat))

! TIME AVERAGING / DIFFERENCING

  ntavg = float(nhour*days_avrg)
  if (days_avrg == 0)  ntavg = 1.

  if (i_time == i_time_last) then
    var4d(:,:,i_time,iv_avg(:)) = var4d(:,:,i_time,iv_avg(:)) +          &
                                  var3d(:,:,iv_avg(:))/ntavg
    t(i_time) = t(i_time) + day_from_ref/ntavg
  else
    var3d(:,:,iv_dudt) = var3d(:,:,iv_um)
    if (days_avrg == 0)  var3d(:,:,iv_dif(:)) =                          &
                         var3d(:,:,iv_dif(:))*(float(nhour)/86400.)
    if (i_time /= 1) then
      var4d(:,:,i_time_last,iv_avg(:)) =                                 &
          var4d(:,:,i_time_last,iv_avg(:)) +                             &
          var3d(:,:,iv_avg(:))/(ntavg*2.)
      t(i_time_last) = t(i_time_last) + day_from_ref/(ntavg*2.)
      var4d(:,:,i_time_last,iv_dif(:)) =                                 &
          var4d(:,:,i_time_last,iv_dif(:)) +                             &
          var3d(:,:,iv_dif(:))/max(1., float(days_avrg)*86400.)
    end if
    if (imon /= nmon+1) then
      var4d(:,:,i_time,iv_avg(:)) = var4d(:,:,i_time,iv_avg(:)) +        &
                                    var3d(:,:,iv_avg(:))/(ntavg*2.)
      t(i_time) = t(i_time) + day_from_ref/(ntavg*2.)
      var4d(:,:,i_time,iv_dif(:)) = var4d(:,:,i_time,iv_dif(:)) -        &
          var3d(:,:,iv_dif(:))/max(1., float(days_avrg)*86400.)
    end if
  end if

  hour = hour + 24/nhour

  i_time_last = i_time

  if (imon == nmon+1)  EXIT L_MON
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

  nt = i_time - 1

  deallocate( hlnp_r, inv_dzp_r, p_r, hlnp_t, inv_dzp_t, p_t )
  deallocate( var_i, k_low, tmp3d, var_i0l, var_p )
  deallocate( k_t1, k_t2 )
  deallocate( coef_rho1, coef_rho2, z_th, z_rho )

  call setbdy


  nd1a = NY2
  nd2a = NPI
  nd3a = NT 
  nd4a = 1

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,1) = var4d(:,:,1:nd3a,iv+50)
  enddo


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'TEM diagnostics',missv=missv)

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  do iv=3, 99
    if ( trim(var_i_name(iv)) == '-999' )  nv_in = iv - 2  ! excluding p
  enddo
  do k=2, 999
    if ( p_intp(k) <= 0. .or. p_intp(k) == 999. .or.                     &
         p_intp(k) > 1100. ) then
      npi = k - 1  ;  EXIT
    end if
  enddo

  allocate( hlnpi(npi) )
  hlnpi(:) = h_scale*log( p_intp(1:npi)*100. )

  iv_t = -999  ;  iv_u = -999  ;  iv_v = -999  ;  iv_w = -999
  do iv=1, nv_in
    if (trim(var_i_name(iv+1)) == 'theta')  iv_t = iv
    if (trim(var_i_name(iv+1)) == 'u'    )  iv_u = iv
    if (trim(var_i_name(iv+1)) == 'v'    )  iv_v = iv
    if (trim(var_i_name(iv+1)) == 'dz_dt')  iv_w = iv
  enddo
  if ( any( (/iv_t, iv_u, iv_v/) == -999 ) ) then
    write(6,*)  'Error: no T, u, v.'  ;  STOP
  end if
  l_thlev(:) = .False.
  l_thlev(iv_t) = .True.
  if (iv_w /= -999) then
    l_w_in = .True.  ;  l_thlev(iv_w) = .True.
  else
    l_w_in = .False.  ;  iv_w = nv_in + 1
  end if
  ! iv_u, iv_v, iv_t, iv_w : from 1 to 4
  iv = max(nv_in, iv_w)
  iv_vu   = iv + 1
  iv_dtdz = iv + 2
  iv_dudz = iv + 3
  l_thlev(iv_dudz) = .True.

  nv_intp = iv_dudz   ! CAUTION ! the last value above: 7, currently

  iv_p0l = 21  ;  iv_t0l = 22  ;  iv_w0l = 23

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 1  ! p at rho  ! for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  allocate( dzp_deriv(npi), hlnpi_l2(npi,2) )
  do k=1, npi
    k0 = minval(minloc(abs( ht(:) - 7.e3*log(1.e3/p_intp(k)) )))
    ! 7 km is used here regardless of the actual h_scale value, as a reference
    dzp_deriv(k) = max(100., 0.5*(ht_th(k0) - ht_th(k0-1)))
  enddo
  hlnpi_l2(:,1) = hlnpi(:) + 0.5*dzp_deriv(:)
  hlnpi_l2(:,2) = hlnpi(:) - 0.5*dzp_deriv(:)

  call get_iouter(lat,lat_rng, iy2)
  ny2 = iy2(2) - iy2(1) + 1

  x1_i = 1   ;  y1_i = iy2(1)   ;  z1_i = 1    ! for getalt
  nx_i = nx  ;  ny_i = ny2      ;  nz_i = nz
  call getalt
  allocate( coef_rho1(nx,ny2,nz-1), coef_rho2(nx,ny2,nz-1) )
  coef_rho1(:,:,:) = ( z_rho(:,:,2:nz) - z_th (:,:,1:nz-1) )/            &
                     ( z_rho(:,:,2:nz) - z_rho(:,:,1:nz-1) )
  coef_rho2(:,:,:) = 1. - coef_rho1(:,:,:)

  allocate( hlnp_r(nx,ny2,nz), inv_dzp_r(nx,ny2,nz), p_r(nx,ny2,nz) )
  allocate( hlnp_t(nx,ny2,nz), inv_dzp_t(nx,ny2,nz), p_t(nx,ny2,nz) )
  allocate( var_i(nx,ny2,nz,nv_intp), var_i0l(nx,ny2,21:23) )
  allocate( var_p(nx,ny2,npi,nv_intp) )
  allocate( k_low(nx,ny2,npi), tmp3d(ny2,npi,31:33) )
  allocate( k_t1(nx,ny2), k_t2(nx,ny2) )

  if (days_avrg == 0) then
    allocate( var4d(ny2,npi,nmon*31*nhour,51:50+nv), t(nmon*31*nhour) )
  else
    allocate( var4d(ny2,npi,nmon*31,51:50+nv), t(nmon*31) )
  end if
  var4d(:,:,:,:) = 0.  ;  t(:) = 0.

  ! some constants
  c_ghri = g*h_scale/rd
  c_lrgi = lapse_rate_sfc*rd/g
  l_k_t1 = .False.

  allocate( coslat(ny2,npi), f(ny2,npi), rho0(ny2,npi),                  &
            ritanlat(ny2,npi) )
  allocate( inv_2dy(ny2) )
!  zp(:) = -h_scale*log(p_intp(1:npi)/1.e3)
  coslat(:,:) = spread(cos(lat(:)*deg2rad),2,nz)
  if (abs(lat(1 )) == 90.)  coslat(1 ,:) = 0.
  if (abs(lat(ny)) == 90.)  coslat(ny,:) = 0.
  ritanlat(:,:) = spread(sin(lat(:)*deg2rad)/( cos(lat(:)*deg2rad)*      &
                         r_earth ),2,nz)
  if (abs(lat(1 )) == 90.)  ritanlat(1 ,:) = 0.
  if (abs(lat(ny)) == 90.)  ritanlat(ny,:) = 0.
 
  f   (:,:) = spread(2.*omega_earth*sin(lat(:)*deg2rad),2,nz)
  rho0(:,:) = spread((p_intp(1:npi)*100.)/g/h_scale  ,1,ny)
  ! Ts is set to gH/R
  inv_2dy(2:ny2-1) = 1./((lat(3:ny2)-lat(1:ny2-2))*deg2rad*r_earth)
  inv_2dy(1  ) = 0.5/((lat(2  )-lat(1    ))*deg2rad*r_earth)
  inv_2dy(ny2) = 0.5/((lat(ny2)-lat(ny2-1))*deg2rad*r_earth)

  if (p_intp(1) > 300.) then  ! [hPa]
    k_t1(:,:) = k_start
    ! searching up
    do k=k_start, nz
    do j=1, ny2
    do i=1, nx
      if ( (z_th(i,j,k) - z_th(i,j,0)) <= h1_realistic_temp )          &
         k_t1(i,j) = k
    enddo
    enddo
    enddo
    l_k_t1 = .True.
  end if

  iv_vs = 51  ;  iv_ws = 52  ;  iv_fy = 53  ;  iv_fz = 54
  iv_epd = 55  ;  iv_psis = 56
  iv_um = 57  ;  iv_tm = 58
  iv_dudt = 59
  iv_advy = 60  ;  iv_advz = 61  ;  iv_dyn = 62
  allocate( iv_dif(1), iv_avg(nv-1) )
  iv_dif = (/iv_dudt/)
  iv_avg = (/iv_vs, iv_ws, iv_fy, iv_fz, iv_epd, iv_psis, iv_um, iv_tm,  &
             iv_advy, iv_advz, iv_dyn/)
 
  allocate( var3d(ny2,npi,51:50+nv) )
  ovarname(iv_vs  -50) = 'v_res'
  ovarname(iv_ws  -50) = 'w_res'
  ovarname(iv_fy  -50) = 'f_y'
  ovarname(iv_fz  -50) = 'f_z'
  ovarname(iv_epd -50) = 'epd'
  ovarname(iv_psis-50) = 'psi_res'
  ovarname(iv_dudt-50) = 'u_tend'
  ovarname(iv_um  -50) = 'u'
  ovarname(iv_tm  -50) = 'T'
  ovarname(iv_advy-50) = 'adv_y'
  ovarname(iv_advz-50) = 'adv_z'
  ovarname(iv_dyn -50) = 'u_force_dyn'

  END subroutine initialize

  SUBROUTINE get_vars

  ex0 = .TRUE.
  do iv_i=1, nv_in+1  ! for get_ifilename
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  ! read p at rho level
  print*, trim(file_i(1))
  iv_i = 1  ;  p_r(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,1,nz)
  hlnp_r(:,:,:) = h_scale*log(p_r(:,:,:))

  ! interpolate p and its log to theta levels
  ! (The Exner ftn is a better choice for vertical interpolation in the
  !  troposphere, but we use log p throughout this program.)
  hlnp_t(:,:,1:nz-1) = hlnp_r(:,:,1:nz-1)*coef_rho1(:,:,:) +             &
                       hlnp_r(:,:,2:nz  )*coef_rho2(:,:,:)
  hlnp_t(:,:,nz) = hlnp_r(:,:,nz)*2. - hlnp_t(:,:,nz-1)
  p_t(:,:,:) = exp(hlnp_t(:,:,:)/h_scale)
 
  iz2(1) = k_start  ;  iz2(2) = nz
  do k=k_start+1, nz
    hlnp1_min = minval(hlnp_r(:,:,k))
    if (hlnp1_min < maxval(hlnpi_l2)) then
      iz2(1) = k-2  ;  EXIT
    end if
  enddo
  do k=nz-1, k_start, -1
    hlnp9_max = maxval(hlnp_r(:,:,k))
    if (hlnp9_max > minval(hlnpi_l2)) then
      iz2(2) = k+2  ;  EXIT
    end if
  enddo
  l_belowgrd = .False.  ;  l_abovetop = .False.
  ipi1 = 1  ;  ipi2 = npi
  if ( iz2(1) < k_start .or. iz2(2) > nz ) then
    if (iz2(1) < k_start) then
      l_belowgrd = .True.
      do kpi=1, npi
        if (hlnp1_min < hlnpi_l2(kpi,1))  ipi1 = kpi+1
      enddo
    end if
    if (iz2(2) > nz) then
      l_abovetop = .True.
      do kpi=npi, 1, -1
        if (hlnp9_max > hlnpi_l2(kpi,2))  ipi2 = kpi-1
      enddo
    end if
    write(6,*) 'Warning: out of domain,', iz2(:)
    iz2(1) = max(k_start, iz2(1))
    iz2(2) = min(nz     , iz2(2))
    write(6,*) '         bounded,      ', iz2(:)
  end if
  if ( .not. l_w_in )  iz2(2) = nz   ! to utilize continuity Eq.
  nz2 = iz2(2) - iz2(1) + 1

  var_i0l(:,:,iv_p0l) = p_t(:,:,iz2(1)-1)

  ! reduce vertical levels to be used
  var_i(:,:,1:nz2,1) = hlnp_r(:,:,iz2(1):iz2(2))
  hlnp_r(:,:,1:nz2) = var_i(:,:,1:nz2,1)

  var_i(:,:,1:nz2,1) = p_r(:,:,iz2(1):iz2(2))
  p_r(:,:,1:nz2) = var_i(:,:,1:nz2,1)

  var_i(:,:,1:nz2,1) = hlnp_t(:,:,iz2(1):iz2(2))
  hlnp_t(:,:,1:nz2) = var_i(:,:,1:nz2,1)

  var_i(:,:,1:nz2,1) = p_t(:,:,iz2(1):iz2(2))
  p_t(:,:,1:nz2) = var_i(:,:,1:nz2,1)

  if (nz2 /= nz) then
    hlnp_r(:,:,nz2+1:nz) = 0.  ;  p_r(:,:,nz2+1:nz) = 0.
    hlnp_t(:,:,nz2+1:nz) = 0.  ;  p_t(:,:,nz2+1:nz) = 0.
  end if
 
  ! -1/d(Hlnp)
  inv_dzp_r(:,:,1:nz2-1) = 1./(hlnp_r(:,:,1:nz2-1) - hlnp_r(:,:,2:nz2))
  inv_dzp_r(:,:,nz2:nz) = 0.

  inv_dzp_t(:,:,2:nz2) = 1./(hlnp_t(:,:,1:nz2-1) - hlnp_t(:,:,2:nz2))
  inv_dzp_t(:,:,1) = 1./(h_scale*log(var_i0l(:,:,iv_p0l)) - hlnp_t(:,:,1))
  if (nz2 /= nz)  inv_dzp_t(:,:,nz2+1:nz) = 0.

  ! read u, v, T, w
  var_i(:,:,:,:) = 0.
  do iv=1, nv_in
    iv_i = iv + 1  ! for get_ivara3d
    print*, trim(file_i(iv_i))
    var_i(:,:,1:nz2,iv) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  enddo
  iv_i = iv_t + 1  ! for get_ivara3d
  var_i0l(:,:,iv_t0l) = reshape(get_ivara3d(1,nx,iy2(1),ny2,iz2(1)-1,1), &
                                (/nx,ny2/))
  if ( l_w_in ) then
    iv_i = iv_w + 1  ! for get_ivara3d
    var_i0l(:,:,iv_w0l) = reshape(get_ivara3d(1,nx,iy2(1),ny2,           &
                                  iz2(1)-1,1), (/nx,ny2/))
  end if
  print*, 'time index :', it_i(2:nv_in)
  if (abs(lat(1  )) == 90.)  var_i(:,1  ,:,iv_v) = 0.
  if (abs(lat(ny2)) == 90.)  var_i(:,ny2,:,iv_v) = 0.

  ! potential temperature to T
  var_i(:,:,1:nz2,iv_t) = var_i(:,:,1:nz2,iv_t)*                         &
                          ( (p_t(:,:,1:nz2)/1.e5)**kappa )
  var_i0l(:,:,iv_t0l) = var_i0l(:,:,iv_t0l)*                             &
                          ( (var_i0l(:,:,iv_p0l)/1.e5)**kappa )

  if ( l_w_in ) then
    ! dz_dt to w_omega
    ! w_omega = -omega*H/p, where omega ~= -w*g*rho = -w*g*p/Rd/Tv
    ! Here, Tv is approximated to T.
    var_i(:,:,1:nz2,iv_w) = var_i(:,:,1:nz2,iv_w)*c_ghri/                &
                            var_i(:,:,1:nz2,iv_t)
    var_i0l(:,:,iv_w0l) = var_i0l(:,:,iv_w0l)*c_ghri/var_i0l(:,:,iv_t0l)
    ! w at theta
  else
    ! using continuity Eq.
    print*,'not coded yet'  ;  STOP
  end if

  if (nz2 /= nz) then
    hlnp_r(:,:,nz2+1:nz) = 0.
    hlnp_t(:,:,nz2+1:nz) = 0.
  end if
 
  END subroutine get_vars

  SUBROUTINE calc_req_vars

  ! vu at rho
  var_i(:,:,1:nz2,iv_vu) = var_i(:,:,1:nz2,iv_v)*var_i(:,:,1:nz2,iv_u)

  ! dT/dz at rho
  var_i(:,:,2:nz2,iv_dtdz) = ( var_i(:,:,2:nz2  ,iv_t) -                 &
                               var_i(:,:,1:nz2-1,iv_t) )*                &
                             inv_dzp_t(:,:,2:nz2)
  var_i(:,:,1,iv_dtdz) = ( var_i(:,:,1,iv_t) - var_i0l(:,:,iv_t0l) )*    &
                         inv_dzp_t(:,:,1)
 
  ! du/dz at theta
  var_i(:,:,1:nz2-1,iv_dudz) = ( var_i(:,:,2:nz2  ,iv_u) -               &
                                 var_i(:,:,1:nz2-1,iv_u) )*              &
                               inv_dzp_r(:,:,1:nz2-1)
  var_i(:,:,nz2,iv_dudz) = var_i(:,:,nz2-1,iv_dudz)

  if (nz2 /= nz)  var_i(:,:,nz2+1:nz,:) = 0.
 
  END subroutine calc_req_vars

  SUBROUTINE vintp_z2p(nvs,ivs,hlnp,inv_dzp,hlnp_req)

  integer,                     intent(in) ::  nvs
  integer, dimension(:),       intent(in) ::  ivs
  real, dimension(nx,ny2,nz2), intent(in) ::  hlnp, inv_dzp
  real, dimension(:),          intent(in) ::  hlnp_req

  real, dimension(nx,ny2,nz2,nvs) ::  vars
  real, dimension(nx,ny2,npi,nvs) ::  vars_p

  integer ::  klast, kbuf, k1add, iv_ti, iv_dtdzi

  if ( size(ivs) /= nvs ) then
    write(6,*) 'Error: nvs/ivs in subroutine vintp_z2p'  ;  STOP
  end if
  if ( size(hlnp_req) /= npi ) then
    write(6,*) 'Error: npi/hlnp_req in subroutine vintp_z2p'  ;  STOP
  end if

  iv_ti = -999  ;  iv_dtdzi = -999
  do iv=1, nvs
    vars(:,:,:,iv) = var_i(:,:,:,ivs(iv))
    if (ivs(iv) == iv_t   )  iv_ti    = iv
    if (ivs(iv) == iv_dtdz)  iv_dtdzi = iv
  enddo

  k_low(:,:,:) = -1   ! if above the top

  L_PLEV1:  DO kpi=1, npi

  klast = 1
  do j=1, ny2
  do i=1, nx
    if ( hlnp(i,j,klast) >= hlnp_req(kpi) ) then
      ! searching up
      do k=klast, nz2-1
        if ( hlnp(i,j,k+1) <= hlnp_req(kpi) ) then
          k_low(i,j,kpi) = k  ;  EXIT
        end if
      enddo
    else
      ! searching down
      k_low(i,j,kpi) = 0   ! if below the bottom
      do k=klast-1, 1, -1
        if ( hlnp(i,j,k) >= hlnp_req(kpi) ) then
          k_low(i,j,kpi) = k  ;  EXIT
        end if
      enddo
    end if
    klast = max(1, k_low(i,j,kpi))
  enddo
  enddo

  ENDDO  L_PLEV1

!+ for vector platforms
!  L_PLEV1v:  DO kpi=1, npi
!
!  ! searching up
!  do k=1, nz2
!  do j=1, ny2
!  do i=1, nx
!    if ( hlnp(i,j,k) < hlnp_req(kpi) .and. k_low(i,j,kpi) == -1 )        &
!       k_low(i,j,kpi) = k-1
!  enddo
!  enddo
!  enddo
!  ! k_low = 0, if below the bottom
!
!  ENDDO  L_PLEV1v
!-

  if ( all( l_thlev(ivs(:)) ) ) then
    kbuf = 1  ;  k1add = 1
  else if ( all( .not. l_thlev(ivs(:)) ) ) then
    kbuf = 0  ;  k1add = 0
  else
    write(6,*)  'Error: levels are indicated wrongly.'
  end if

  if (ipi1 <= ipi2-kbuf) then

    L_PLEV2:  DO kpi=ipi1, ipi2-kbuf  ! add kbuf for theta level variables

    do j=1, ny2
    do i=1, nx
      vars_p(i,j,kpi,:) = inv_dzp(i,j,k_low(i,j,kpi)+k1add)*             &
          ( vars(i,j,k_low(i,j,kpi),:)*                                  &
              ( hlnp_req(kpi) - hlnp(i,j,k_low(i,j,kpi)+1) ) +           &
            vars(i,j,k_low(i,j,kpi)+1,:)*                                &
              ( hlnp(i,j,k_low(i,j,kpi)) - hlnp_req(kpi) ) )
    enddo
    enddo

    ENDDO  L_PLEV2

  end if

  if ( l_belowgrd ) then  ! based on rho level
                          ! safe also for theta level variables

    if ( .not. l_k_t1 ) then
      write(6,*) 'Error: k_t1 has not been obtained.'  ;  STOP
    end if

    k_t2(:,:) = k_t1(:,:) - (iz2(1)-1)

    L_PLEV2l1:  DO kpi=1, ipi1-1

    do j=1, ny2
    do i=1, nx
      if (k_low(i,j,kpi) == 0) then
        vars_p(i,j,kpi,:) = vars(i,j,1,:)
        ! following UM (vn6.6), use the nearest level below 200 m if the
        ! first level is below 200 m
        ! It is chosen as a realistic physical value for the atmosphere.
        if (iv_ti /= -999)  vars_p(i,j,kpi,iv_ti) =                      &
           vars(i,j,k_t2(i,j),iv_ti)*                                    &
             ( exp((hlnp_req(kpi) - hlnp(i,j,k_t2(i,j)))/h_scale)        &
             )**c_lrgi
        if (iv_dtdzi /= -999)  vars_p(i,j,kpi,iv_dtdzi) =                &
           vars(i,j,k_t2(i,j)+1,iv_dtdzi)
      else
        vars_p(i,j,kpi,:) = inv_dzp(i,j,k_low(i,j,kpi)+k1add)*           &
            ( vars(i,j,k_low(i,j,kpi),:)*                                &
                ( hlnp_req(kpi) - hlnp(i,j,k_low(i,j,kpi)+1) ) +         &
              vars(i,j,k_low(i,j,kpi)+1,:)*                              &
                ( hlnp(i,j,k_low(i,j,kpi)) - hlnp_req(kpi) ) )
      end if
    enddo
    enddo

    ENDDO  L_PLEV2l1

  end if

  if ( l_abovetop .or. kbuf > 0 ) then  ! based on rho level
                          ! add kbuf for theta level variables

    L_PLEV2l2:  DO kpi=ipi2-kbuf+1, npi

    do j=1, ny2
    do i=1, nx
      if (k_low(i,j,kpi) == -1) then
        ! assume constant
        vars_p(i,j,kpi,:) = vars(i,j,nz2,:)
      else
        vars_p(i,j,kpi,:) = inv_dzp(i,j,k_low(i,j,kpi)+k1add)*           &
            ( vars(i,j,k_low(i,j,kpi),:)*                                &
                ( hlnp_req(kpi) - hlnp(i,j,k_low(i,j,kpi)+1) ) +         &
              vars(i,j,k_low(i,j,kpi)+1,:)*                              &
                ( hlnp(i,j,k_low(i,j,kpi)) - hlnp_req(kpi) ) )
      end if
    enddo
    enddo

    ENDDO  L_PLEV2l2

  end if

  do iv=1, nvs
    var_p(:,:,:,ivs(iv)) = vars_p(:,:,:,iv)
  enddo

  END subroutine vintp_z2p

  FUNCTION grady(var0)

  real, dimension(:,:), intent(in)  ::  var0
  real, dimension(ny2,size(var0,2)) ::  grady

  integer ::  nz_gr, k_gr

  nz_gr = size(var0,2)

  do k_gr=1, nz_gr
    grady(2:ny2-1,k_gr) = (var0(3:ny2,k_gr)-var0(1:ny2-2,k_gr))*inv_2dy(2:ny2-1)
  enddo
  grady(1  ,:) = 0.
  grady(ny2,:) = 0.

  END function grady

  SUBROUTINE setbdy

  if (abs(lat(1)) == 90.) then
    var4d(1,:,:,iv_vs  ) = 0.
    var4d(1,:,:,iv_ws  ) = missv
    var4d(1,:,:,iv_fy  ) = 0.
    var4d(1,:,:,iv_fz  ) = 0.
    var4d(1,:,:,iv_epd ) = 0.
    var4d(1,:,:,iv_psis) = 0.
    var4d(1,:,:,iv_dudt) = 0.
    var4d(1,:,:,iv_um  ) = 0.
    var4d(1,:,:,iv_advy) = 0.
    var4d(1,:,:,iv_advz) = 0.
    var4d(1,:,:,iv_dyn ) = 0.
  end if
  if (abs(lat(ny2)) == 90.) then
    var4d(ny2,:,:,iv_vs  ) = 0.
    var4d(ny2,:,:,iv_ws  ) = missv
    var4d(ny2,:,:,iv_fy  ) = 0.
    var4d(ny2,:,:,iv_fz  ) = 0.
    var4d(ny2,:,:,iv_epd ) = 0.
    var4d(ny2,:,:,iv_psis) = 0.
    var4d(ny2,:,:,iv_dudt) = 0.
    var4d(ny2,:,:,iv_um  ) = 0.
    var4d(ny2,:,:,iv_advy) = 0.
    var4d(ny2,:,:,iv_advz) = 0.
    var4d(ny2,:,:,iv_dyn ) = 0.
  end if

  END subroutine setbdy

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lat  ','p','time',' '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat(iy2(1):iy2(2))
  set(iv)%axis2 = p_intp(1:nd2a)
  set(iv)%axis3 = t(1:nd3a)
  set(iv)%axis4 = -999.
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var4d, hlnpi, dzp_deriv, hlnpi_l2 )
  deallocate( lon, lat, ht, ht_th, t )
  deallocate( coslat, f, rho0, ritanlat, inv_2dy )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program QBOI_DYN_UM


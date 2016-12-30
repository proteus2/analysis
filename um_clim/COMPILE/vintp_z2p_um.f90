PROGRAM VINTP_UM_z2p

  use hadgem
  use netio

  implicit none

  include 'c_phys.inc'

  real, parameter  ::  h1_realistic_temp = 200.  ! [m]
  real, parameter  ::  h1_safecalc_hgt = 2000.  ! [m]

  integer          ::  days_avrg
  real             ::  p_intp(999)
  character(len=8) ::  var_i_vert_grid

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ P_INTP, LAT_RNG, DAYS_AVRG
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, VAR_I_NAME, FILE_ALT, VAR_I_VERT_GRID,  &
                    VAR_ALT, FILE_O

  integer ::  nv, npi
  integer ::  imon, ihour, i_time, i_time_last, ipi1, ipi2
  integer ::  iy2(2), iz2(2), ny2, nz2, klast
  integer ::  i,j,k, kpi, iv_ht, iv_t, iv_w
  real    ::  lnp1_min, lnp9_max, c_ghri, c_lrgi, ntavg
  logical ::  l_thlev, l_belowgrd, l_abovetop, l_k_z1, l_k_t1
  character(len=32), dimension(:), allocatable ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  var5d
  real, dimension(:,:,:,:),   allocatable ::  var_i, var4d
  real, dimension(:,:,:),     allocatable ::  lnp3d, inv_dlnpi3d, p3d
  real, dimension(:,:,:),     allocatable ::  p3d_th_ht
  real, dimension(:,:,:),     allocatable ::  coef_rho0, coef_rho1
  real, dimension(:),         allocatable ::  lnpi
  integer, dimension(:,:,:),  allocatable ::  k_low
  integer, dimension(:,:),    allocatable ::  k_z1, k_t1

  type(vset) ::  set(10)

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

! INTERPOLATE VARIABLES VERTICALLY

  k_low(:,:,:) = -1   ! if above the top

  L_PLEV1:  DO kpi=1, npi

  klast = 1
  do j=1, ny2
  do i=1, nx
    if ( lnp3d(i,j,klast) >= lnpi(kpi) ) then
      ! searching up
      do k=klast, nz2-1
        if ( lnp3d(i,j,k+1) <= lnpi(kpi) ) then
          k_low(i,j,kpi) = k  ;  EXIT
        end if
      enddo
    else
      ! searching down
      k_low(i,j,kpi) = 0   ! if below the bottom
      do k=klast-1, 1, -1
        if ( lnp3d(i,j,k) >= lnpi(kpi) ) then
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
!    if ( lnp3d(i,j,k) < lnpi(kpi) .and. k_low(i,j,kpi) == -1 )           &
!       k_low(i,j,kpi) = k-1
!  enddo
!  enddo
!  enddo
!  ! k_low = 0, if below the bottom
!
!  ENDDO  L_PLEV1v
!-

  inv_dlnpi3d(:,:,1:nz2-1) = 1./(lnp3d(:,:,1:nz2-1) - lnp3d(:,:,2:nz2))
  inv_dlnpi3d(:,:,nz2:nz) = 0.

  if (ipi1 <= ipi2) then

    L_PLEV2:  DO kpi=ipi1, ipi2

    do j=1, ny2
    do i=1, nx
      var4d(i,j,kpi,:) = inv_dlnpi3d(i,j,k_low(i,j,kpi))*                &
          ( var_i(i,j,k_low(i,j,kpi),:)*                                 &
              ( lnpi(kpi) - lnp3d(i,j,k_low(i,j,kpi)+1) ) +              &
            var_i(i,j,k_low(i,j,kpi)+1,:)*                               &
              ( lnp3d(i,j,k_low(i,j,kpi)) - lnpi(kpi) ) )
    enddo
    enddo

    ENDDO  L_PLEV2

  end if

  if ( l_belowgrd ) then

    if ( iv_ht /= -999 .and. ( .not. l_k_z1 ) ) then
      k_z1(:,:) = 1
      ! searching up
      do k=1, nz  ! nz, not nz2, as h1_safecalc_hgt is quite high
      do j=1, ny2
      do i=1, nx
        if ( (z_th(i,j,k) - z_th(i,j,0)) <= h1_safecalc_hgt )            &
           k_z1(i,j) = k
      enddo
      enddo
      enddo
      l_k_z1 = .True.
    end if

    if ( iv_t /= -999 .and. ( .not. l_k_t1 ) ) then
      k_t1(:,:) = 1
      ! searching up
      do k=1, nz2
      do j=1, ny2
      do i=1, nx
        if ( (z_th(i,j,k) - z_th(i,j,0)) <= h1_realistic_temp )          &
           k_t1(i,j) = k
      enddo
      enddo
      enddo
      l_k_t1 = .True.
    end if

    L_PLEV2l1:  DO kpi=1, ipi1-1

    do j=1, ny2
    do i=1, nx
      if (k_low(i,j,kpi) == 0) then
        var4d(i,j,kpi,:) = var_i(i,j,1,:)
        ! following UM (vn6.6), use the nearest level below 2 km
        if (iv_ht /= -999)  var4d(i,j,kpi,iv_ht) =                       &
           var_i(i,j,1,iv_ht) +                                          &
           ( z_th(i,j,k_z1(i,j)) - var_i(i,j,1,iv_ht) )*                 &
             ( 1. - ( exp(lnpi(kpi) - lnp3d(i,j,1)) )**c_lrgi )/         &
             ( 1. - ( p3d_th_ht(i,j,k_z1(i,j))/p3d(i,j,1) )**c_lrgi )
        ! following UM (vn6.6), use the nearest level below 200 m if the
        ! first level is below 200 m
        ! It is chosen as a realistic physical value for the atmosphere.
        if (iv_t /= -999)  var4d(i,j,kpi,iv_t) =                         &
           var_i(i,j,k_t1(i,j),iv_t)*                                    &
             ( exp(lnpi(kpi) - lnp3d(i,j,k_t1(i,j))) )**c_lrgi
      else
        var4d(i,j,kpi,:) = inv_dlnpi3d(i,j,k_low(i,j,kpi))*              &
            ( var_i(i,j,k_low(i,j,kpi),:)*                               &
                ( lnpi(kpi) - lnp3d(i,j,k_low(i,j,kpi)+1) ) +            &
              var_i(i,j,k_low(i,j,kpi)+1,:)*                             &
                ( lnp3d(i,j,k_low(i,j,kpi)) - lnpi(kpi) ) )
      end if
    enddo
    enddo

    ENDDO  L_PLEV2l1

  end if

  if ( l_abovetop ) then

    L_PLEV2l2:  DO kpi=ipi2+1, npi

    do j=1, ny2
    do i=1, nx
      if (k_low(i,j,kpi) == -1) then
        ! assume constant T, winds, etc. except for Z
        ! simple linear extrapolation for Z (assuming dZ/dlnp = const.)
        var4d(i,j,kpi,:) = var_i(i,j,nz2,:)
        if (iv_ht /= -999)  var4d(i,j,kpi,iv_ht) =                       &
           var_i(i,j,nz2,iv_ht) +                                        &
           ( var_i(i,j,nz2,iv_ht) - var_i(i,j,nz2-1,iv_ht) )*            &
             inv_dlnpi3d(i,j,nz2-1)*(lnp3d(i,j,nz2) - lnpi(kpi))
      else
        var4d(i,j,kpi,:) = inv_dlnpi3d(i,j,k_low(i,j,kpi))*              &
            ( var_i(i,j,k_low(i,j,kpi),:)*                               &
                ( lnpi(kpi) - lnp3d(i,j,k_low(i,j,kpi)+1) ) +            &
              var_i(i,j,k_low(i,j,kpi)+1,:)*                             &
                ( lnp3d(i,j,k_low(i,j,kpi)) - lnpi(kpi) ) )
      end if
    enddo
    enddo

    ENDDO  L_PLEV2l2

  end if


  if (days_avrg == 0) then
    var5d(:,:,:,i_time,:) = var4d(:,:,:,:)
    t(i_time) = day_from_ref
  else
    ntavg = float(nhour*days_avrg)
    if (i_time == i_time_last) then
      var5d(:,:,:,i_time,:) = var5d(:,:,:,i_time,:) +                    &
                              var4d(:,:,:,:)/ntavg
      t(i_time) = t(i_time) + day_from_ref/ntavg
    else
      if (i_time /= 1) then
        var5d(:,:,:,i_time_last,:) = var5d(:,:,:,i_time_last,:) +        &
                                     var4d(:,:,:,:)/(ntavg*2.)
        t(i_time_last) = t(i_time_last) + day_from_ref/(ntavg*2.)
      end if
      if (imon /= nmon+1) then
        var5d(:,:,:,i_time,:) = var5d(:,:,:,i_time,:) +                  &
                                var4d(:,:,:,:)/(ntavg*2.)
        t(i_time) = t(i_time) + day_from_ref/(ntavg*2.)
      end if
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

  deallocate( lnp3d, inv_dlnpi3d, p3d, var_i, k_low )
  if (iv_ht /= -999)  deallocate( k_z1, p3d_th_ht )
  if (iv_t  /= -999)  deallocate( k_t1 )
  if ( allocated(coef_rho0) )  deallocate( coef_rho0, coef_rho1 )
  if ( allocated(z_rho    ) )  deallocate( z_th, z_rho )


  nd1a = NX
  nd2a = NY2
  nd3a = NPI
  nd4a = NT

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = var5d(:,:,:,1:nd4a,iv)
  enddo


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'interpolated to pressure levels')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  do iv=3, 99
    if ( trim(var_i_name(iv)) == '-999' )  nv = iv - 2  ! excluding p
  enddo
  do k=2, 999
    if ( p_intp(k) <= 0. .or. p_intp(k) == 999. .or.                     &
         p_intp(k) > 1100. ) then
      npi = k - 1  ;  EXIT
    end if
  enddo

  allocate( lnpi(npi) )
  lnpi(:) = log( p_intp(1:npi)*100. )

  l_thlev = .False.
  if ( trim(var_i_vert_grid) == 'theta' .or.                             &
       trim(var_i_vert_grid) == 'zt' )  l_thlev = .True.

  iv_ht = -999  ;  iv_t = -999  ;  iv_w = -999
  do iv=1, nv
    if ( trim(var_i_name(iv+1)) == 'ht' .or.                             &
         trim(var_i_name(iv+1)) == 'ht_1' )  iv_ht = iv
    if (trim(var_i_name(iv+1)) == 'theta')  iv_t = iv
    if (trim(var_i_name(iv+1)) == 'dz_dt')  iv_w = iv
  enddo
 
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

  call get_iouter(lat,lat_rng, iy2)
  ny2 = iy2(2) - iy2(1) + 1

  if ( l_thlev .or. iv_ht /= -999 ) then
    x1_i = 1   ;  y1_i = iy2(1)   ;  z1_i = 1    ! for getalt
    nx_i = nx  ;  ny_i = ny2      ;  nz_i = nz
    call getalt
    allocate( coef_rho0(nx,ny2,nz-1), coef_rho1(nx,ny2,nz-1) )
    coef_rho0(:,:,:) = ( z_rho(:,:,2:nz) - z_th (:,:,1:nz-1) )/          &
                       ( z_rho(:,:,2:nz) - z_rho(:,:,1:nz-1) )
    coef_rho1(:,:,:) = 1. - coef_rho0(:,:,:)
  end if

  allocate( lnp3d(nx,ny2,nz), inv_dlnpi3d(nx,ny2,nz), p3d(nx,ny2,nz) )
  allocate( var_i(nx,ny2,nz,nv) )
  allocate( k_low(nx,ny2,npi) )
  if (iv_ht /= -999)  allocate( k_z1(nx,ny2), p3d_th_ht(nx,ny2,nz) )
  if (iv_t  /= -999)  allocate( k_t1(nx,ny2) )

  if (days_avrg == 0) then
    allocate( var5d(nx,ny2,npi,nmon*31*nhour,nv), t(nmon*31*nhour) )
  else
    allocate( var5d(nx,ny2,npi,nmon*31,nv), t(nmon*31) )
  end if
  var5d(:,:,:,:,:) = 0.  ;  t(:) = 0.
  allocate( var4d(nx,ny2,npi,nv) )

  allocate( ovarname(nv) )
  do iv=1, nv
    ovarname(iv) = trim(var_i_name(iv+1))
  enddo
  if (iv_ht /= -999)  ovarname(iv_ht) = 'ht'
  if (iv_t  /= -999)  ovarname(iv_t ) = 'T'
  if (iv_w  /= -999)  ovarname(iv_w ) = 'w_omega'

  ! some constants
  c_ghri = g*h_scale/rd
  c_lrgi = lapse_rate_sfc*rd/g
  l_k_z1 = .False.  ;  l_k_t1 = .False.

  END subroutine initialize

  SUBROUTINE get_vars

  ex0 = .TRUE.
  do iv_i=1, nv+1  ! for get_ifilename
    if (iv_i == iv_ht+1)  CYCLE
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  ! read p at rho level
  print*, trim(file_i(1))
  iv_i = 1  ;  p3d(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,1,nz)
  lnp3d(:,:,:) = log(p3d(:,:,:))

  if ( l_thlev ) then
    ! interpolate p and its log to theta levels
    var_i(:,:,:,1) = lnp3d(:,:,:)
    lnp3d(:,:,1:nz-1) = var_i(:,:,1:nz-1,1)*coef_rho0(:,:,:) +           &
                        var_i(:,:,2:nz  ,1)*coef_rho1(:,:,:)
    lnp3d(:,:,nz) = var_i(:,:,nz,1)*2. - lnp3d(:,:,nz-1)
    p3d(:,:,:) = exp(lnp3d(:,:,:))
  else if (iv_ht /= -999) then
    ! interpolate p to theta levels
    ! Exner ftn is a better choice than log p in the troposphere
    var_i(:,:,:,1) = (p3d(:,:,:)/1.e5)**kappa
    p3d_th_ht(:,:,1:nz-1) = var_i(:,:,1:nz-1,1)*coef_rho0(:,:,:) +       &
                            var_i(:,:,2:nz  ,1)*coef_rho1(:,:,:)
    p3d_th_ht(:,:,nz) = var_i(:,:,nz,1)*2. - p3d_th_ht(:,:,nz-1)
    p3d_th_ht(:,:,:) = 1.e5*(p3d_th_ht(:,:,:)**(1./kappa))
  end if
 
  iz2(1) = 1  ;  iz2(2) = nz
  do k=1, nz
    lnp1_min = minval(lnp3d(:,:,k))
    if (lnp1_min < maxval(lnpi)) then
      iz2(1) = k-1  ;  EXIT
    end if
  enddo
  do k=nz, 1, -1
    lnp9_max = maxval(lnp3d(:,:,k))
    if (lnp9_max > minval(lnpi)) then
      iz2(2) = k+1  ;  EXIT
    end if
  enddo
  l_belowgrd = .False.  ;  l_abovetop = .False.
  ipi1 = 1  ;  ipi2 = npi
  if ( iz2(1) < 1 .or. iz2(2) > nz ) then
    if (iz2(1) < 1) then
      l_belowgrd = .True.
      do kpi=1, npi
        if (lnp1_min < lnpi(kpi))  ipi1 = kpi+1
      enddo
    end if
    if (iz2(2) > nz) then
      l_abovetop = .True.
      do kpi=npi, 1, -1
        if (lnp9_max > lnpi(kpi))  ipi2 = kpi-1
      enddo
    end if
    write(6,*) 'Warning: out of domain,', iz2(:)
    iz2(1) = max(1 , iz2(1))
    iz2(2) = min(nz, iz2(2))
    write(6,*) '         bounded,      ', iz2(:)
  end if
  nz2 = iz2(2) - iz2(1) + 1

  ! reduce vertical levels to be used
  var_i(:,:,1:nz2,1) = lnp3d(:,:,iz2(1):iz2(2))
  lnp3d(:,:,1:nz2) = var_i(:,:,1:nz2,1)

  var_i(:,:,1:nz2,1) = p3d(:,:,iz2(1):iz2(2))
  p3d(:,:,1:nz2) = var_i(:,:,1:nz2,1)

  ! read other var.s
  var_i(:,:,:,:) = 0.
  it_i(2:nv+1) = -999
  do iv=1, nv
    if (iv == iv_ht) then
      var_i(:,:,1:nz2,iv_ht) = z_rho(:,:,iz2(1):iz2(2))
      CYCLE
    end if
    iv_i = iv + 1  ! for get_ivara3d
    print*, trim(file_i(iv_i))
    var_i(:,:,1:nz2,iv) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  enddo
  print*, 'time index :', it_i(2:nv+1)

  ! potential temperature to T
  if (iv_t /= -999)  var_i(:,:,1:nz2,iv_t) =                             &
     var_i(:,:,1:nz2,iv_t)*( (p3d(:,:,1:nz2)/1.e5)**kappa )

  ! dz_dt to w_omega
  ! w_omega = -omega*H/p, where omega ~= -w*g*rho = -w*g*p/Rd/Tv
  ! Here, Tv is approximated to T.
  if (iv_w /= -999) then
    if (iv_t /= -999) then
      var_i(:,:,1:nz2,iv_w) = var_i(:,:,1:nz2,iv_w)*c_ghri/              &
                              var_i(:,:,1:nz2,iv_t)
    else
      write(6,*) ' Error: omega needs T.'  ;  STOP
    end if
  end if
 
  if (nz2 /= nz) then
    lnp3d(:,:,nz2+1:nz) = 0.  ;  var_i(:,:,nz2+1:nz,:) = 0.
  end if

  ! keep z_th(:,:,1:nz) to obtain k_z1 correctly

  END subroutine get_vars

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lon  ','lat','p','time'/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lon
  set(iv)%axis2 = lat(iy2(1):iy2(2))
  set(iv)%axis3 = p_intp(1:nd3a)
  set(iv)%axis4 = t(1:nd4a)
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var5d, var4d, lnpi )
  deallocate( ovarname )
  deallocate( lon, lat, ht, ht_th, t )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program VINTP_UM_z2p


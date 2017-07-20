PROGRAM TEM3d_REANALYSIS
! attributes (scale facter, add_offset, "_FillValue" or "missing_value") must be considered.

  use tem3d
  use util,  only: lowpass_k, filter121_y
  use reanal
  use netio
  use const_glob,  only: g, rd

  implicit none

  integer, parameter ::  nv = 26, nv2 = 0
  integer, parameter ::  nwgth = 28  ! 28*2+1 days
  real   , parameter ::  h_s = 7.e3  ! 7 km scale height
  real   , parameter ::  lapse_rate_sfc = 6.5e-3  ! for missing data

  integer ::  k_max

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH
  namelist /PARAM/ LAT_RNG, P_RNG, K_MAX
  namelist /FILEIO/ NT_F4, MISSV, FILE_I_HEAD, FILE_I_FORM, FILE_I_XXXX, &
                    VAR_I, VAR_I_NAME, FILE_O

  integer ::  imon, ihour, date_s, date_e, i_time, i_save, i_target
  integer ::  iy2(2), iz2(2), ny2, nz2, iy3(2), iz3(2)
  integer ::  iy2b(2), iz2b(2), iy2o(2), iz2o(2), nbuf_y2(2), nbuf_z2(2)
  integer ::  i,j,k,n, kk
  real    ::  ntavg, dlon
  real, dimension(-nwgth:nwgth) ::  lanczos, lanczos_shift, ltmp
  character(len=32), dimension(nv+nv2) ::  ovarname

  real, dimension(:,:,:,:), allocatable ::  var4d, var4ds
!nv2  real, dimension(:,:,:),   allocatable ::  var3d, var3ds
  real, dimension(:,:,:,:), allocatable ::  u, v, te, gp, w
  real, dimension(:,:,:),   allocatable ::  u1, v1, te1, gp1, w1
  real, dimension(:,:,:),   allocatable ::  us, vs, tes, gps, ws
  real, dimension(:,:,:),   allocatable ::  te_prt, gp_prt
  real, dimension(:),       allocatable ::  lat0, p0a

  type(vset), dimension(nv+nv2) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  call initialize

  var4ds(:,:,:,:) = 0.

  i_time = 0
  i_save = -nwgth
  i_target = 0

  L_MON:  DO imon=0, nmon+1
  !---------------------------------------------------------------------

  ndate = get_ndate()
  date_s = 1  ;  date_e = ndate

  if (imon == 0) then
    date_s = ndate+1 - nwgth
  else if (imon == nmon+1) then
    date_e = nwgth
  end if

  L_DATE:  DO date=date_s, date_e
  !---------------------------------------------------------------------
  print*, year, mon, date

  u(:,:,:,i_save) = 0.  ;  v(:,:,:,i_save) = 0.
  te(:,:,:,i_save) = 0.  ;  gp(:,:,:,i_save) = 0.
  w(:,:,:,i_save) = 0.

print*,'read'

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------

  ! get variable
  call get_4var

  u (:,:,:,i_save) = u (:,:,:,i_save) + u1 (:,:,:)
  v (:,:,:,i_save) = v (:,:,:,i_save) + v1 (:,:,:)
  te(:,:,:,i_save) = te(:,:,:,i_save) + te1(:,:,:)
  gp(:,:,:,i_save) = gp(:,:,:,i_save) + gp1(:,:,:)
  w (:,:,:,i_save) = w (:,:,:,i_save) + w1 (:,:,:)

  hour = hour + 24/nhour

  !---------------------------------------------------------------------
  ENDDO  L_HOUR

  u (:,:,:,i_save) = u (:,:,:,i_save)/float(nhour)
  v (:,:,:,i_save) = v (:,:,:,i_save)/float(nhour)
  te(:,:,:,i_save) = te(:,:,:,i_save)/float(nhour)
  gp(:,:,:,i_save) = gp(:,:,:,i_save)/float(nhour)
  w (:,:,:,i_save) = w (:,:,:,i_save)/float(nhour)

print*, 'filter'
  ! filter out small-scale waves
  call lowpass_k(u (:,:,:,i_save),k_max)
  call lowpass_k(v (:,:,:,i_save),k_max)
  call lowpass_k(te(:,:,:,i_save),k_max)
  call lowpass_k(gp(:,:,:,i_save),k_max)
  call lowpass_k(w (:,:,:,i_save),k_max)
  call filter121_y(u (:,:,:,i_save),'U')
  call filter121_y(v (:,:,:,i_save),'V')
  call filter121_y(te(:,:,:,i_save),'S')
  call filter121_y(gp(:,:,:,i_save),'S')
  call filter121_y(w (:,:,:,i_save),'S')

  if ( i_time == 0 .and. i_save < nwgth ) then
    i_save = i_save + 1
    CYCLE  L_DATE
  end if

  i_time = i_time + 1

print*,'mean and pert.'
  us(:,:,:) = 0.  ;  vs(:,:,:) = 0.
  tes(:,:,:) = 0.  ;  gps(:,:,:) = 0.
  ws(:,:,:) = 0.
  do n=-nwgth, nwgth
    us (:,:,:) = us (:,:,:) + lanczos_shift(n)*u (:,:,:,n)
    vs (:,:,:) = vs (:,:,:) + lanczos_shift(n)*v (:,:,:,n)
    tes(:,:,:) = tes(:,:,:) + lanczos_shift(n)*te(:,:,:,n)
    gps(:,:,:) = gps(:,:,:) + lanczos_shift(n)*gp(:,:,:,n)
    ws (:,:,:) = ws (:,:,:) + lanczos_shift(n)*w (:,:,:,n)
  enddo
  te_prt(:,:,:) = te(:,:,:,i_target) - tes(:,:,:)
  gp_prt(:,:,:) = gp(:,:,:,i_target) - gps(:,:,:)

print*,'TEM'
  ! calculate zonal mean
  call waf3d_nons_qg( nx,ny2,nz2,lat0,p0a*100.,                          &
      gp_prt,te_prt,gps,tes,us,vs,ws,dlon,h_s,1.e32,                     &
      var4d(:,:,:,1),var4d(:,:,:,2),var4d(:,:,:,3),var4d(:,:,:,4),       &
      var4d(:,:,:,5),var4d(:,:,:,6),var4d(:,:,:,7),var4d(:,:,:,8),       &
      var4d(:,:,:,9),var4d(:,:,:,10),var4d(:,:,:,11),var4d(:,:,:,12),    &
      var4d(:,:,:,13),var4d(:,:,:,14),var4d(:,:,:,15),var4d(:,:,:,16),   &
      var4d(:,:,:,17),var4d(:,:,:,18),var4d(:,:,:,19),var4d(:,:,:,20),   &
      var4d(:,:,:,21),var4d(:,:,:,22) )

  var4ds(:,:,:,:22) = var4ds(:,:,:,:22) + var4d(:,:,:,:22)
  var4ds(:,:,:,23) = var4ds(:,:,:,23) + us (:,:,:)
  var4ds(:,:,:,24) = var4ds(:,:,:,24) + vs (:,:,:)
  var4ds(:,:,:,25) = var4ds(:,:,:,25) + tes(:,:,:)
  if (i_time == 1)  var4ds(:,:,:,26) = var4d(:,:,:,15)

  i_save = i_save + 1
  if (i_save == nwgth+1)  i_save = -nwgth
  i_target = i_target + 1
  if (i_target == nwgth+1)  i_target = -nwgth

  ltmp(-nwgth) = lanczos_shift(nwgth)
  ltmp(-nwgth+1:nwgth) = lanczos_shift(-nwgth:nwgth-1)
  lanczos_shift(:) = ltmp(:)

  !---------------------------------------------------------------------
  ENDDO  L_DATE

  mon = mon + 1
  if (mon == 13) then
    year = year + 1
    mon = 1
  end if

  !---------------------------------------------------------------------
  ENDDO  L_MON

  nt = i_time

  deallocate( u, v, te, gp, u1, v1, te1, gp1, te_prt, gp_prt )
  deallocate( us, vs, tes, gps )
  deallocate( w, w1, ws )

  var4ds(:,:,:,:25) = var4ds(:,:,:,:25)/float(nt)
  var4ds(:,:,:,26) = (var4d(:,:,:,15) - var4ds(:,:,:,26))/float(nt-1)
print*,nt

  ! missing

  nd1a = nx
  nd2a = ny2 - (nbuf_y2(1) + nbuf_y2(2))
  nd3a = nz2 - (nbuf_z2(1) + nbuf_z2(2))
  nd4a = 1

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,1) = var4ds(:,iy3(1):iy3(2),iz3(1):iz3(2),iv)
  enddo
!nv2  do iv=nv+1, nv+nv2
!nv2    call setdim2
!nv2    allocate( set(iv)%var_out(nd2a,nd3a,nd4a,1) )
!nv2    set(iv)%var_out(:,:,:,:) = 1.e32
!nv2
!nv2    set(iv)%var_out(:,:,1,1) = var3ds(iy3(1):iy3(2),iz3(1):iz3(2),iv-nv)
!nv2  enddo


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv+nv2,set,'transformed Eulerian mean eqn.')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  if ( nt_f4(4) /= 1 .or. ( nt_f4(3) /= 1 .and. nt_f4(3) /= 12 ) ) then
    print*, 'NT_I(3:4) should be (/1,1/) or (/12,1/).'  ;  STOP
  end if

  year = yyyy
  mon  = mm(1) - 1
  if (mon == 0) then
    year = yyyy - 1
    mon  = 12
  end if

  nmon = mm(2)  ;  nhour = hh(2)

  ex0 = .TRUE.
  do iv_i=1, 5
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP
  print*, trim(file_i(3))

  iv_i = 3  ! temperature
  call getdim(file_i(iv_i),var_i_name(iv_i))
  dlon = lon(2) - lon(1)

  ovarname(1:22) = varname_waf3d_nons_qg(1:22)
  ovarname(23) = 'u_s'
  ovarname(24) = 'v_s'
  ovarname(25) = 'T_s'
  ovarname(26) = 'tend_A'

  call wgt_lanczos30d

  call get_iouter(lat,lat_rng, iy2o)
  iy2b(1) = max(1 ,iy2o(1)-3)
  iy2b(2) = min(ny,iy2o(2)+3)
  nbuf_y2(:) = abs(iy2b(:) - iy2o(:))
  iy2(:) = iy2b(:)
  if ( l_rev(2) ) then
    iy2(1) = ny + 1 - iy2b(2)
    iy2(2) = ny + 1 - iy2b(1)
  end if
  call get_iouter(p*(-1.),p_rng*(-1.), iz2o)
  iz2b(1) = max(1 ,iz2o(1)-3)
  iz2b(2) = min(nz,iz2o(2)+3)
  nbuf_z2(:) = abs(iz2b(:) - iz2o(:))
  iz2(:) = iz2b(:)
  if ( l_rev(3) ) then
    iz2(1) = nz + 1 - iz2b(2)
    iz2(2) = nz + 1 - iz2b(1)
  end if

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  iy3(1) = 1 + nbuf_y2(1)  ;  iy3(2) = ny2 - nbuf_y2(2)
  iz3(1) = 1 + nbuf_z2(1)  ;  iz3(2) = nz2 - nbuf_z2(2)

  allocate( lat0(ny2), p0a(nz2) )
  lat0(:) = lat(iy2b(1):iy2b(2))  ;  p0a(:) = p(iz2b(1):iz2b(2))

  allocate( var4ds(nx,ny2,nz2,nv) )
!nv2  allocate( var3ds(ny2,nz2,nv2) )
  var4ds(:,:,:,:) = 0.
!nv2  var3ds(:,:,:) = 0.
  allocate( var4d(nx,ny2,nz2,nv) )
!nv2  allocate( var3d(ny2,nz2,nv2) )

  allocate( u1(nx,ny2,nz2), v1(nx,ny2,nz2), te1(nx,ny2,nz2),             &
            gp1(nx,ny2,nz2) )
  allocate( w1(nx,ny2,nz2) )
  allocate( u(nx,ny2,nz2,-nwgth:nwgth), v(nx,ny2,nz2,-nwgth:nwgth),      &
            te(nx,ny2,nz2,-nwgth:nwgth), gp(nx,ny2,nz2,-nwgth:nwgth) )
  allocate( w(nx,ny2,nz2,-nwgth:nwgth) )
  allocate( us(nx,ny2,nz2), vs(nx,ny2,nz2), tes(nx,ny2,nz2),             &
            gps(nx,ny2,nz2) )
  allocate( ws(nx,ny2,nz2) )
  allocate( te_prt(nx,ny2,nz2), gp_prt(nx,ny2,nz2) )
 
  END subroutine initialize

  SUBROUTINE get_4var

  real ::  ztmp1, ztmp2, tmp

  ex0 = .TRUE.
  do iv_i=1, 5
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  print*, trim(file_i(3))
  if (.not. ex0)  STOP

  ! read 4 var.s
  iv_i = 1  ;  u1 (:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 2  ;  v1 (:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 3  ;  te1(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 4  ;  gp1(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 5  ;  w1 (:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)

  ! check whether the input gp1 is geopotential or GPH
  if ( trim(unit_h) == '' ) then
    ztmp1 = h_s*log(1.e3/p0a(nz2))
    Z_IDENT:  do k=nz2, 1, -1
    do j=1, ny2
    do i=1, nx
      if (gp1(i,j,k) /= missv) then
        ztmp2 = gp1(i,j,k)  ;  EXIT Z_IDENT
      end if
    enddo
    enddo
    enddo  Z_IDENT
    if ( abs(ztmp2 - ztmp1) < abs(ztmp2/g - ztmp1) ) then
      unit_h = 'm'
    else
      unit_h = 'm**2 s**-2'
    end if
  end if

  ! GPH to GP  (The missing value also changes, if it exists.)
  if ( trim(unit_h) == 'm' )  gp1(:,:,:) = gp1(:,:,:)*g

  ! omega to w(zp)
  do k=1, nz2
    w1(:,:,k) = -w1(:,:,k)/(p0a(k)*100.)*h_s
  enddo

  if (missv /= 1.0) then
    tmp = lapse_rate_sfc*rd/g
    do j=1, ny2
    do i=1, nx
    do k=nz2, 1, -1
      if (te1(i,j,k) == missv) then
        u1(i,j,1:k) = u1(i,j,k+1)
        v1(i,j,1:k) = v1(i,j,k+1)
        w1(i,j,1:k) = 0.
        do kk=k, 1, -1
          te1(i,j,kk) = te1(i,j,kk+1)*(p0a(kk)/p0a(kk+1))**tmp
!          gp1(i,j,kk) = gp1(i,j,kk+1)*
        enddo
        EXIT
      endif
    enddo
    enddo
    enddo
  end if

  END subroutine get_4var

  SUBROUTINE setdim

  if ( .not. allocated(t) ) then
    allocate( t(nmon) )
    t(:) = (/ (n, n=1,nmon) /)
  end if

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lon  ','lat  ','p ',' '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lon
  set(iv)%axis2 = lat0(iy3(1):iy3(2))
  set(iv)%axis3 = p0a(iz3(1):iz3(2))
  set(iv)%axis4 = t
    
  END subroutine setdim

  SUBROUTINE setdim2

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lat  ','p ',' ',' '/) 
  set(iv)%nd(:) = (/nd2a,nd3a,nd4a,1/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat0(iy3(1):iy3(2))
  set(iv)%axis2 = p0a(iz3(1):iz3(2))
  set(iv)%axis3 = t
  set(iv)%axis4 = -999.
    
  END subroutine setdim2

  SUBROUTINE finalize

  deallocate( var4ds, var4d )
!nv2  deallocate( var3ds, var3d )
  deallocate( lon, lat, p, t, t2pt, lat0, p0a )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize

  SUBROUTINE wgt_lanczos30d

    lanczos = (/ -0.0001626153, -0.0005024954, -0.00101763, -0.001686738, -0.002467858, -0.003298386, -0.004096661, -0.00476516, -0.005195176, -0.005272799, -0.00488586, -0.003931446, -0.002323499, -2.108048e-09, 0.003070771, 0.006885407, 0.01140231, 0.01654075, 0.02218224, 0.02817415, 0.03433555, 0.04046498, 0.04634967, 0.05177583, 0.0565393, 0.06045602, 0.06337157, 0.06516935, 0.06577682, 0.06516935, 0.06337157, 0.06045602, 0.0565393, 0.05177583, 0.04634967, 0.04046498, 0.03433555, 0.02817415, 0.02218224, 0.01654075, 0.01140231, 0.006885407, 0.003070771, -2.108048e-09, -0.002323499, -0.003931446, -0.00488586, -0.005272799, -0.005195176, -0.00476516, -0.004096661, -0.003298386, -0.002467858, -0.001686738, -0.00101763, -0.0005024954, -0.0001626153 /)

    lanczos_shift(:) = lanczos(:)

  END subroutine wgt_lanczos30d

END program TEM3d_REANALYSIS


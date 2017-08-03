PROGRAM TEM3d_REANALYSIS
! attributes (scale facter, add_offset, "_FillValue" or "missing_value") must be considered.

  use tem3d
  use util,  only: lowpass_k
  use reanal
  use netio
  use const_glob,  only: g, rd

  implicit none

  integer, parameter ::  nv = 9, nv2 = 2
  real   , parameter ::  h_s = 7.e3  ! 7 km scale height
  real   , parameter ::  lapse_rate_sfc = 6.5e-3  ! for missing data

  integer ::  k_max

  namelist /ANALCASE/ EXPNAME, YYYY, MM
  namelist /PARAM/ LAT_RNG, P_RNG, K_MAX
  namelist /FILEIO/ NT_F4, MISSV, FILE_I_HEAD, FILE_I_FORM, FILE_I_XXXX, &
                    VAR_I, VAR_I_NAME, FILE_O

  integer ::  imon
  integer ::  iy2(2), iz2(2), ny2, nz2, iy3(2), iz3(2)
  integer ::  iy2b(2), iz2b(2), iy2o(2), iz2o(2), nbuf_y2(2), nbuf_z2(2)
  integer ::  i,j,k,n, kk
  real    ::  ntavg, dlon
  character(len=32), dimension(nv+nv2) ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  var5d
  real, dimension(:,:,:,:), allocatable ::  var4d, var4d2
  real, dimension(:,:,:),   allocatable ::  var3d2
  real, dimension(:,:,:),   allocatable ::  u, v, te, gp, w
  real, dimension(:,:),     allocatable ::  wm
  real, dimension(:),       allocatable ::  lat0, p0a

  type(vset), dimension(nv+nv2) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  call initialize

  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------

  ! get variable
  allocate( w(nx,ny2,nz2) )
  call get_4var

  wm(:,:) = sum(w, dim=1)/float(nx)
  wm(:,:) = -wm(:,:)*spread(h_s/(p0a(:)*100.),1,ny2)

  deallocate( w )

  ! filter out small-scale waves
  call lowpass_k(gp,k_max)
  call lowpass_k(te,k_max)

  ! calculate zonal mean
  call waf3d_s_qg_gp(                                                    &
      nx,ny2,nz2,lat0,p0a*100.,gp,te,dlon,h_s,1.e32,                     &
      var4d(:,:,:,1),var4d(:,:,:,2),var4d(:,:,:,3),var4d(:,:,:,4),       &
      var4d(:,:,:,5),var4d(:,:,:,6),var4d(:,:,:,7),                      &
      var3d2(:,:,1),var3d2(:,:,2) )

  var5d(:,:,:,imon,:7) = var4d(:,:,:,:7)
  var5d(:,:,:,imon,8) = u(:,:,:)
  var5d(:,:,:,imon,9) = v(:,:,:)
  var4d2(:,:,imon,:2) = var3d2(:,:,:2)

  mon = mon + 1
  if (mon == 13) then
    year = year + 1
    mon = 1
  end if

  !---------------------------------------------------------------------
  ENDDO  L_MON

  deallocate( u, v, te, gp, wm )

  ! missing

  nd1a = nx
  nd2a = ny2 - (nbuf_y2(1) + nbuf_y2(2))
  nd3a = nz2 - (nbuf_z2(1) + nbuf_z2(2))
  nd4a = nmon

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = var5d(:,iy3(1):iy3(2),iz3(1):iz3(2),1:nd4a,iv)
  enddo
  do iv=nv+1, nv+nv2
    call setdim2
    allocate( set(iv)%var_out(nd2a,nd3a,nd4a,1) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,1) = var4d2(iy3(1):iy3(2),iz3(1):iz3(2),1:nd4a,iv-nv)
  enddo


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
  mon  = mm(1)

  nmon = mm(2)

  ex0 = .TRUE.
  do iv_i=1, 4
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

  ovarname(1:7) = varname_waf3d_s_qg(1:7)
  ovarname(8) = 'u'
  ovarname(9) = 'v'
  ovarname(10:11) = varname_waf3d_s_qg(18:19)

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

  allocate( var5d(nx,ny2,nz2,nmon,nv) )
  allocate( var4d2(ny2,nz2,nmon,nv2) )
  var5d(:,:,:,:,:) = 0.
  var4d2(:,:,:,:) = 0.
  allocate( var4d(nx,ny2,nz2,nv) )
  allocate( var3d2(ny2,nz2,nv2) )

  allocate( u(nx,ny2,nz2), v(nx,ny2,nz2), te(nx,ny2,nz2), gp(nx,ny2,nz2) )
  allocate( wm(ny2,nz2) )
 
  END subroutine initialize

  SUBROUTINE get_4var

  real ::  ztmp1, ztmp2, tmp

  ex0 = .TRUE.
  do iv_i=1, 4
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  print*, trim(file_i(3))
  if (.not. ex0)  STOP

  ! read 4 var.s
  iv_i = 1  ;  u (:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 2  ;  v (:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 3  ;  te(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 4  ;  gp(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
!  iv_i = 5  ;  w (:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)

  ! check whether the input gp is geopotential or GPH
  if ( trim(unit_h) == '' ) then
    ztmp1 = h_s*log(1.e3/p0a(nz2))
    Z_IDENT:  do k=nz2, 1, -1
    do j=1, ny2
    do i=1, nx
      if (gp(i,j,k) /= missv) then
        ztmp2 = gp(i,j,k)  ;  EXIT Z_IDENT
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
  if ( trim(unit_h) == 'm' )  gp(:,:,:) = gp(:,:,:)*g

  if (missv /= 1.0) then
    tmp = lapse_rate_sfc*rd/g
    do j=1, ny2
    do i=1, nx
    do k=nz2, 1, -1
      if (u(i,j,k) == missv) then  ! do not use gp whose missing value is changed
        u(i,j,1:k) = u(i,j,k+1)
        v(i,j,1:k) = v(i,j,k+1)
        w(i,j,1:k) = w(i,j,k+1)
        do kk=k, 1, -1
          te(i,j,kk) = te(i,j,kk+1)*(p0a(kk)/p0a(kk+1))**tmp
!          gp(i,j,kk) = gp(i,j,kk+1) ...
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
  set(iv)%axis = (/'lon  ','lat  ','p ','time'/) 
  if (nmon == 1)  set(iv)%axis(4) = ' '
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
  set(iv)%axis = (/'lat  ','p ','time',' '/) 
  if (nmon == 1)  set(iv)%axis(3) = ' '
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

  deallocate( var5d, var4d, var4d2, var3d2 )
  deallocate( lon, lat, p, t, t2pt, lat0, p0a )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program TEM3d_REANALYSIS


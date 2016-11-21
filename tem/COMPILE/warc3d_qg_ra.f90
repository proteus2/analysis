PROGRAM TEM3d_REANALYSIS
! attributes (scale facter, add_offset, "_FillValue" or "missing_value") must be considered.

  use warc3d
  use reanal
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 9, nv2 = 2

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH
  namelist /PARAM/ OPT_AVRG, LAT_RNG, P_RNG
  namelist /FILEIO/ NT_F4, MISSV, FILE_I_HEAD, FILE_I_FORM, FILE_I_XXXX, &
                    VAR_I, VAR_I_NAME, UNIT_H, FILE_O

  integer ::  imon, ihour, i_time, i_time_last
  integer ::  iy2(2), iz2(2), ny2, nz2, iy3(2), iz3(2)
  integer ::  iy2b(2), iz2b(2), iy2o(2), iz2o(2), nbuf_y2(2), nbuf_z2(2)
  integer ::  i,j,k,n, kk
  integer ::  tag_exit
  real    ::  ntavg, dlon
  character(len=32), dimension(nv+nv2) ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  var5d
  real, dimension(:,:,:,:), allocatable ::  var4d, var4d2
  real, dimension(:,:,:),   allocatable ::  var3d2
  real, dimension(:,:,:),   allocatable ::  u, v, pt, gp, w
  real, dimension(:,:),     allocatable ::  wm
  real, dimension(:),       allocatable ::  lat0, p0a, t2pt0

  type(vset), dimension(nv+nv2) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  if ( nt_f4(4) /= 1 .or. ( nt_f4(3) /= 1 .and. nt_f4(3) /= 12 ) ) then
    print*, 'NT_I(3:4) should be (/1,1/) or (/12,1/).'  ;  STOP
  end if

  call initialize

  tag_exit = 0

  i_time_last = 0
  i_time      = 0

  L_MON:  DO imon=1, nmon+1
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT
if (date == 3)  exit  !yh

  if (opt_avrg /= 0)  i_time = i_time + 1

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

  if (opt_avrg == 0)  i_time = i_time + 1

  ! get variable
  allocate( w(nx,ny2,nz2) )
  call get_4var

  if (missv /= 1.0) then
do j=1, ny2
do i=1, nx
do k=nz2, 1, -1
  if (u(i,j,k) == missv) then  ! CAUTION, if pt or gp, of which
                               ! the missing value is changed
    u(i,j,1:k) = u(i,j,k+1)
    v(i,j,1:k) = v(i,j,k+1)
    w(i,j,1:k) = w(i,j,k+1)
    do kk=k, 1, -1
      pt(i,j,kk) = pt(i,j,kk+1)*exp( -log(pt(i,j,kk+2)/pt(i,j,kk+1))*  &
                   log(p0a(kk)/p0a(kk+1))/log(p0a(kk+1)/p0a(kk+2)) )
!      gp(i,j,kk) = gp(i,j,kk+1) ...
    enddo
    exit
  endif
enddo
enddo
enddo
  end if

  wm(:,:) = sum(w, dim=1)/float(nx)
  wm(:,:) = -wm(:,:)*spread(h_scale/(p0a(:)*100.),1,ny2)

  deallocate( w )

  ! calculate zonal mean
  call warc_s_qg(                                                        &
      nx,ny2,nz2,lat0,p0a*100.,u,v,pt,gp,dlon,h_scale,1.e32,             &
      var4d(:,:,:,1),var4d(:,:,:,2),var4d(:,:,:,3),var4d(:,:,:,4),       &
      var4d(:,:,:,5),var4d(:,:,:,6),var4d(:,:,:,7),var3d2(:,:,1) )

! if tadv_z exists among output from the above subroutine
!  var4d(:,:,:,XX) = var4d(:,:,:,XX)*spread(spread((p0a(:)*100./p0)**kappa,1,nx),2,ny2)

  ntavg = float(nhour)
  if (opt_avrg == 0)  ntavg = 1.

  if (i_time == i_time_last) then
    var5d(:,:,:,i_time,:7) = var5d(:,:,:,i_time,:7) +                    &
                             var4d(:,:,:,:7)/ntavg
    var4d2(:,:,i_time,:1) = var4d2(:,:,i_time,:1) +                      &
                            var3d2(:,:,:1)/ntavg
  else
    var4d(:,:,:,8) = u(:,:,:)
    var4d(:,:,:,9) = pt(:,:,:)*spread(spread((p0a(:)*100./p0)**kappa,    &
                     1,nx),2,ny2)
    var3d2(:,:,2) = var3d2(:,:,1)
    if (opt_avrg == 0)  var4d(:,:,:,8:9) = var4d(:,:,:,8:9)*float(nhour)
    if (opt_avrg == 0)  var3d2(:,:,2) = var3d2(:,:,2)*float(nhour)

    if (i_time /= 1) then
      var5d(:,:,:,i_time_last,:7) = var5d(:,:,:,i_time_last,:7) +        &
                                    var4d(:,:,:,:7)/(ntavg*2.)
      var5d(:,:,:,i_time_last,8:9) = var5d(:,:,:,i_time_last,8:9) +      &
                                     var4d(:,:,:,8:9)
      var4d2(:,:,i_time_last,:1) = var4d2(:,:,i_time_last,:1) +          &
                                   var3d2(:,:,:1)/(ntavg*2.)
      var4d2(:,:,i_time_last,2) = var4d2(:,:,i_time_last,2) +            &
                                  var3d2(:,:,2)
    end if
    if (imon /= nmon+1) then
      var5d(:,:,:,i_time,:7) = var5d(:,:,:,i_time,:7) +                  &
                               var4d(:,:,:,:7)/(ntavg*2.)
      var5d(:,:,:,i_time,8:9) = var5d(:,:,:,i_time,8:9) -                &
                                var4d(:,:,:,8:9)
      var4d2(:,:,i_time,:1) = var4d2(:,:,i_time,:1) +                    &
                              var3d2(:,:,:1)/(ntavg*2.)
      var4d2(:,:,i_time,2) = var4d2(:,:,i_time,2) - var3d2(:,:,2)
    end if
  end if

  hour = hour + 24/nhour

  i_time_last = i_time

  if (imon == nmon+1)  tag_exit = 1
  !---------------------------------------------------------------------
  ENDDO  L_HOUR

  !---------------------------------------------------------------------
  ENDDO  L_DATE

  mon = mon + 1
  if (mon == 13) then
    year = year + 1
    mon = 1
  end if

  !---------------------------------------------------------------------
  ENDDO  L_MON

  nt = i_time - 1

  deallocate( u, v, pt, gp, wm )

  ! pt2t
  ! missing

  nd1a = nx
  nd2a = ny2 - (nbuf_y2(1) + nbuf_y2(2))
  nd3a = nz2 - (nbuf_z2(1) + nbuf_z2(2))
  nd4a = nt

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

  year = yyyy
  mon  = mm(1)

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = get_ndate()  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename

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

  ovarname(1:7) = varname_warc_qg(1:7)
  ovarname(8) = 'u_tend'
  ovarname(9) = 't_tend'
  ovarname(10) = varname_warc_qg(8)
  ovarname(11) = 'U0_tend'

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

  allocate( lat0(ny2), p0a(nz2), t2pt0(nz2) )
  lat0(:) = lat(iy2b(1):iy2b(2))  ;  p0a(:) = p(iz2b(1):iz2b(2))
  t2pt0(:) = t2pt(iz2b(1):iz2b(2))

  if (opt_avrg == 0) then
    allocate( var5d(nx,ny2,nz2,nmon*31*nhour,nv) )
    allocate( var4d2(ny2,nz2,nmon*31*nhour,nv2) )
  else
    allocate( var5d(nx,ny2,nz2,nmon*31,nv) )
    allocate( var4d2(ny2,nz2,nmon*31,nv2) )
  end if
  var5d(:,:,:,:,:) = 0.
  var4d2(:,:,:,:) = 0.
  allocate( var4d(nx,ny2,nz2,nv) )
  allocate( var3d2(ny2,nz2,nv2) )

  allocate( u(nx,ny2,nz2), v(nx,ny2,nz2), pt(nx,ny2,nz2), gp(nx,ny2,nz2) )
  allocate( wm(ny2,nz2) )
 
  END subroutine initialize

  SUBROUTINE get_4var

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
  iv_i = 3  ;  pt(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 4  ;  gp(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  iv_i = 5  ;  w (:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)

  ! t to pt  (The missing value changes, if it exists.)
  do k=1, nz2
    pt(:,:,k) = pt(:,:,k)*t2pt0(k)
  enddo

  ! gph to gp  (The missing value changes, if it exists.)
  if ( unit_h == "m" .or. unit_h == "M" )  gp(:,:,:) = gp(:,:,:)*g

  END subroutine get_4var

  SUBROUTINE setdim

  if ( .not. allocated(t) ) then
    allocate( t(nt) )
    do n=1, nt
      t(n) = float(n-1)+0.5
    enddo
    if (opt_avrg == 0)  t(:) = t(:)/nhour
  end if

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lon  ','lat  ','p ','time'/) 
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
  deallocate( lon, lat, p, t, lat0, p0a, t2pt, t2pt0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program TEM3d_REANALYSIS


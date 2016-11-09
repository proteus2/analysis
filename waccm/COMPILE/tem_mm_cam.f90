PROGRAM TEM_CAM

  use tem
  use cam
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 12
  real,    parameter ::  zb_out = 650.

  integer ::  days_avrg

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_FEBLEAP
  namelist /PARAM/ DAYS_AVRG
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, VAR_I_NAME, FILE_ALT, VAR_ALT, FILE_O

  integer ::  iz, imon, ihour, i_time, i_time_last
  integer ::  tag_exit
  real    ::  ntavg
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:), allocatable ::  var4d
  real, dimension(:,:,:),   allocatable ::  var3d
  real, dimension(:,:,:),   allocatable ::  u, v, w, pt, pt_th
  real, dimension(:,:),     allocatable ::  ptm, ptm_th

  type(vset), dimension(nv) ::  set

  ovarname(:11) = varname_tem_p(:11)
  ovarname(12) = 'u_tend'

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  year = yyyy
  mon  = mm(1)

  call initialize

  tag_exit = 0
  i_time_last = 0  ;  i_time = 0

  L_MON:  DO imon=1, nmon+1
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

!mm  ndate = get_ndate()
  ndate = 1  !mm

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

  if (days_avrg /= 0) then
    if (mod(date,days_avrg) == 1)  i_time = i_time + 1
  end if

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

  if (days_avrg == 0)  i_time = i_time + 1

!mm  day_from_ref = get_dayfromref(year,mon,date,hour)
  day_from_ref = dim4(1)  !mm
  if (imon == nmon+1)  EXIT  !mm

  ex0 = .TRUE.
  do iv_i=1, 4
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  ! get variable
  allocate( u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), pt(nx,ny,nz) )
  allocate( dptmdz_from_thlev(ny,nz) )

  allocate( pt_th(nx,ny,0:nz) )

  call get_4var
  call v_interpol2z_5var

  allocate( ptm_th(ny,0:nz), ptm(ny,nz) )
  ptm_th(:,:) = sum(pt_th, dim=1)/float(nx)
  ptm   (:,:) = sum(pt   , dim=1)/float(nx)
  do k=1, nz
    dptmdz_from_thlev(:,k) = ptm(:,k)*                                   &
        log(ptm_th(:,k)/ptm_th(:,k-1))/(h_scale*log(ht_th(k-1)/ht_th(k)))
  enddo
  deallocate( ptm, ptm_th )

  deallocate( pt_th )

  ! calculate zonal mean
  call tem_hydro_p(                                                      &
      nx,ny,nz,lat,ht*1.e2,u,v,w,pt,h_scale,1.e32,                       &
      var3d(:,:,1 ),var3d(:,:,2 ),var3d(:,:,3 ),var3d(:,:,4 ),           &
      var3d(:,:,5 ),var3d(:,:,6 ),var3d(:,:,7 ),var3d(:,:,8 ),           &
      var3d(:,:,9 ),var3d(:,:,10),var3d(:,:,11) )

  deallocate( dptmdz_from_thlev )

  ntavg = float(nhour*days_avrg)
  if (days_avrg == 0)  ntavg = 1.
  ntavg = 0.5  !mm

  if (i_time == i_time_last) then
    var4d(:,:,i_time,:11) = var4d(:,:,i_time,:11) +                      &
                            var3d(:,:,:11)/ntavg
    t(i_time) = t(i_time) + day_from_ref/ntavg
  else
    var3d(:,:,12) = sum(u, dim=1)/float(nx)
    if (days_avrg == 0)  var3d(:,:,12) = var3d(:,:,12)*float(nhour)

    if (i_time /= 1) then
      var4d(:,:,i_time_last,:11) = var4d(:,:,i_time_last,:11) +          &
                                   var3d(:,:,:11)/(ntavg*2.)
      t(i_time_last) = t(i_time_last) + day_from_ref/(ntavg*2.)
      var4d(:,:,i_time_last,12) = var4d(:,:,i_time_last,12) +            &
                                     var3d(:,:,12)/float(days_avrg)
    end if
    if (imon /= nmon+1) then
      var4d(:,:,i_time,:11) = var4d(:,:,i_time,:11) +                    &
                              var3d(:,:,:11)/(ntavg*2.)
      t(i_time) = t(i_time) + day_from_ref/(ntavg*2.)
      var4d(:,:,i_time,12) = var4d(:,:,i_time,12) -                      &
                                var3d(:,:,12)/float(days_avrg)
    end if
  end if

  deallocate( pt )
  deallocate( u, v, w )

  hour = hour + 24/nhour

  i_time_last = i_time

  if (imon == nmon+1)  tag_exit = 1
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

  nd1a = NY
  nd2a = NZ
  nd3a = NT
  nd4a = 1

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,1) = var4d(:,:,1:nd3a,iv)

    set(iv)%var_out(:,:2   ,:,:) = 1.e32
    set(iv)%var_out(:,nz-1:,:,:) = 1.e32
    set(iv)%var_out(1:4,nz-3:nz,:,:) = 1.e32
  enddo

  if (zb_out /= 0.) then
    iz = nd2a + 1
    do k=2, nd2a
      if (ht(k) > zb_out) then  ;  iz = k - 1  ;  EXIT  ;  end if
    enddo
    do iv=1, nv
      set(iv)%var_out(:,iz+1:,:,1) = 1.e32
    enddo
  end if


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'transformed Eulerian mean eqn.')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nmon = mm(2)  ;  nhour = hh(2)

  date = 1  ;  hour = hh(1)  ! for get_ifilename
  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 1  ! u
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz   ! for getalt
!  call getalt

  if (days_avrg == 0) then
    allocate( var4d(ny,nz,nmon*31*nhour,nv), t(nmon*31*nhour) )
  else
    allocate( var4d(ny,nz,nmon*31,nv), t(nmon*31) )
  end if
  var4d(:,:,:,:) = 0.  ;  t(:) = 0.
  allocate( var3d(ny,nz,nv) )

  END subroutine initialize

  SUBROUTINE get_4var

  real ::  tmp(3)

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz   ! for get_ivar3d

  print*, trim(file_i(1))
  iv_i = 1  ;  u    (:,:,:) = get_ivar3d()
  print*, trim(file_i(2))
  iv_i = 2  ;  v    (:,:,:) = get_ivar3d()
  print*, trim(file_i(3))
  iv_i = 3  ;  w    (:,:,:) = get_ivar3d()

  w(:,:,:) = -h_scale/spread(spread(ht(:)*1.e2,1,nx),2,ny)*w(:,:,:)

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz+1   ! for get_ivar3d

  print*, trim(file_i(4))
  iv_i = 4  ;  pt_th(:,:,:) = get_ivar3d()

  do k=2, nz
    if ( minval(pt_th(:,:,k)) /= 0. ) then
      tmp(1) = (ht(k-1)/1.e3)**0.286
      tmp(2) = (ht(k  )/1.e3)**0.286
      tmp(3) = 6.5e-3*h_scale*log(ht(k)/ht(k-1))
      where ( pt_th(:,:,k) == 0. )
        pt_th(:,:,k) = ( pt_th(:,:,k-1)*tmp(1) + tmp(3) )/tmp(2)
      end where
    end if
  enddo

  print*, 'time index :', it_i(1:4)

  END subroutine get_4var

  SUBROUTINE v_interpol2z_5var

  real, dimension(:,:,:), allocatable ::  lnvar

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz   ! for v_linintp_th2z
  k_const_rho = 0  !!
  allocate( lnvar(nx_i,ny_i,0:nz_i) )
  lnvar(:,:,:) = log(pt_th(:,:,:))
  call v_linintp_0th2z(lnvar)
  pt(:,:,:) = exp(lnvar(:,:,1:))
  deallocate( lnvar )

!  nx_i = nx  ;  ny_i = ny  ;  nz_i = k_const_rho+1   ! for v_cubintp_*
!  call v_cubintp_rho2z(u(:,:,:nz_i))
!  call v_cubintp_rho2z(v(:,:,:nz_i))
!  call v_cubintp_rho2z(w(:,:,:nz_i))
!
!  allocate( lnvar(nx_i,ny_i,0:nz_i) )
!
!  lnvar(:,:,:) = log(pt_th(:,:,:nz_i))
!  call v_cubintp_th2z(lnvar(:,:,:),log(pt(:,:,nz_i+1)))
!  pt(:,:,:nz_i) = exp(lnvar(:,:,:))
!
!  lnvar(:,:,:) = log(pt_th(:,:,:nz_i))
!  call v_cubintp_th2zth(lnvar)
!  pt_th(:,:,:nz_i) = exp(lnvar(:,:,:))
!
!  deallocate( lnvar )

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz

  END subroutine v_interpol2z_5var

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lat  ','p ','time',' '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat
  set(iv)%axis2 = ht
  set(iv)%axis3 = t
  set(iv)%axis4 = -999.
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var4d, var3d )
  deallocate( lon, lat, ht, t )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program TEM_CAM


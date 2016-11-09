PROGRAM ZME2

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  nv = 4 !5
  integer, parameter ::  nvo = nv*4+1
  real,    parameter ::  missv = 1.e32, v_small = 1.e-5
  character(len=64)  ::  ivarname(nv), ovarname(nvo)
  character(len=128) ::  f_namelist

  integer ::  yyyymm, date(2), utc(12), xy_itv(2)
  real    ::  fct(2), fct_itv
  character(len=10)  ::  expname, exception(99)
  character(len=128) ::  file_i, file_o

  namelist /ANALCASE/ EXPNAME, YYYYMM
  namelist /PARAM/ DATE, UTC, FCT, FCT_ITV, EXCEPTION, XY_ITV
  namelist /FILEIO/ FILE_I, FILE_O

  ! input variables
  data  ivarname(1) /'u'/
  data  ivarname(2) /'v'/
  data  ivarname(3) /'temp'/
  data  ivarname(4) /'ht'/
!  data  ivarname(5) /'p_1'/
  data  f_namelist  /'./namelist/nl.input'/
  ! output variables
  data  ovarname(1 :4 ) /'ZME_U','ZMSE_U','cZME_U','cZMSE_U'/
  data  ovarname(5 :8 ) /'ZME_V','ZMSE_V','cZME_V','cZMSE_V'/
  data  ovarname(9 :12) /'ZME_T','ZMSE_T','cZME_T','cZMSE_T'/
  data  ovarname(13:16) /'ZME_Z','ZMSE_Z','cZME_Z','cZMSE_Z'/
!  data  ovarname(17:20) /'ZME_P','ZMSE_P','cZME_P','cZMSE_P'/
  data  ovarname(17)    /'coslat'/

  integer ::  nutc, nt, nfct, nx, ny, nz, nxu, nyv, nt0
  integer ::  nxr, nyr, nzr, ntr, nxn, nyn, nd1a, nd2a, nd3a, nd4a
  integer ::  i,j,k,n, iv, idat, iutc, ifct, i_time, it, ncid, ncid2, nn, iv2
  integer ::  year, month, tyear, tmonth, tdate, thour, fct_h
  logical ::  l_getdim, ex1, ex2
  character(len=10)  ::  timec, ttimec
  character(len=32)  ::  c_axis(3,2)
  character(len=256) ::  fname, ftrue
  character(len=3)   ::  hhh

  real, dimension(:), allocatable ::  t, t_fct, lon, lat, p, lonu, latv
  real, dimension(:), allocatable ::  lato, cosphi, time0
  real, dimension(:,:,:,:), allocatable ::  var, vart, var2, var3

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvo) ::  set

! READ NAMELISTS

  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! DEFINE AXIS

  do n=1, 12
    if (utc(n) < 0) then  ;  nutc = n - 1  ;  EXIT  ;  end if
  enddo

  call gettaxis

  l_getdim = .TRUE.

! LOOP

  L_VAR:  DO iv=1, nv


  i_time = 0


  L_DAT:  DO idat=date(1), date(2)
  L_UTC:  DO iutc=1, nutc


  year  = int(yyyymm/100)
  month = yyyymm - year*100

  i_time = i_time + 1


  L_FCT:  DO ifct=1, nfct


  fct_h = int(t_fct(ifct)*24.)

  write(hhh,'(i3.3)') fct_h
  call fct_target(year,month,idat,utc(iutc),fct_h, tyear,tmonth,tdate,thour)

  write(timec ,'(i4.4,3i2.2)') year, month, idat, utc(iutc)
  write(ttimec,'(i4.4,3i2.2)') tyear, tmonth, tdate, thour

  write(6,*)
  write(6,'(a,i3.3,a)') ' '//trim(ivarname(iv))//' : '//timec//' + ',fct_h,' h'

  call getfname
  call check_ex(fname, ex1)  ;  call check_ex(ftrue, ex2)
  if (.not. (ex1 .and. ex2))  CYCLE

  call opennc(fname,ncid )
  call opennc(ftrue,ncid2)

  if (l_getdim)  call getdim
  if (iv >= 5)  nz = 1

  ! allocate output
  nxr = NX   ;   nd1a = NYN
  nyr = NY   ;   nd2a = NZ
  nzr = NZ   ;   nd3a = NFCT
  ntr = NT   ;   nd4a = NT

  do nn=1, 4
    iv2 = (iv-1)*4 + nn
    if (.not. allocated(set(iv2)%var_out)) then
      allocate( set(iv2)%var_out(nd1a,nd2a,nd3a,nd4a) )
      set(iv2)%var_out(:,:,:,:) = missv
    end if
  enddo

  ! get var
  allocate( var(nxr,nyr,nzr,1), vart(nxr,nyr,nzr,1) )
  select case ( iv )
  case ( 1 )
    allocate( var2(nxr,nyr,nzr,1), var3(nxr,nyr,nzr,1) )
    call geta4d(ncid, trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,1,1, var2)
    call geta4d(ncid2,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,1,1, var3)
    do k=1, nzr
    do j=1, nyr
      do i=2, nxr
        var (i,j,k,1) = 0.5*(var2(i-1,j,k,1)+var2(i,j,k,1))
        vart(i,j,k,1) = 0.5*(var3(i-1,j,k,1)+var3(i,j,k,1))
      enddo
      var (1,j,k,1) = 0.5*(var2(nxr,j,k,1)+var2(1,j,k,1))
      vart(1,j,k,1) = 0.5*(var3(nxr,j,k,1)+var3(1,j,k,1))
    enddo
    enddo
    deallocate( var2, var3 )
  case ( 2 )
    allocate( var2(nxr,nyv,nzr,1), var3(nxr,nyv,nzr,1) )
    call geta4d(ncid, trim(ivarname(iv)),1,nxr,1,nyv,1,nzr,1,1, var2)
    call geta4d(ncid2,trim(ivarname(iv)),1,nxr,1,nyv,1,nzr,1,1, var3)
    do k=1, nzr
    do j=2, nyr-1
      var (:,j,k,1) = 0.5*(var2(:,j-1,k,1)+var2(:,j,k,1))
      vart(:,j,k,1) = 0.5*(var3(:,j-1,k,1)+var3(:,j,k,1))
    enddo
    enddo
    var (:,1,:,1) = 0.  ;  var (:,nyr,:,1) = 0.
    vart(:,1,:,1) = 0.  ;  vart(:,nyr,:,1) = 0.
    deallocate( var2, var3 )
  case default
    call geta4d(ncid ,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,1,1, var )
    call geta4d(ncid2,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,1,1, vart)
  end select

  call closenc(ncid )
  call closenc(ncid2)

  ! calculate ZME, ZMSE
  allocate( var2(nyn,nzr,2,1) )
  do k=1, nzr
    call zm_err2(nxr,nyr,1,xy_itv(1),xy_itv(2),1,var(:,:,k,1),vart(:,:,k,1), &
                  var2(:,k,1,1),var2(:,k,2,1))
  enddo

  deallocate( var, vart )

  ! save output
  iv2 = (iv-1)*4
  set(iv2+1)%var_out(:,:,ifct,i_time) = var2(:,:,1,1)
  set(iv2+2)%var_out(:,:,ifct,i_time) = var2(:,:,2,1)
  do k=1, nd2a
    set(iv2+3)%var_out(:,k,ifct,i_time) =  &
    set(iv2+1)%var_out(:,k,ifct,i_time)*cosphi(:)
    set(iv2+4)%var_out(:,k,ifct,i_time) =  &
    set(iv2+2)%var_out(:,k,ifct,i_time)*cosphi(:)
  enddo

  deallocate ( var2 )


  ENDDO  L_FCT


  ENDDO  L_UTC
  ENDDO  L_DAT


  do nn=1, 4
    iv2 = (iv-1)*4 + nn

    set(iv2)%vname = ovarname(iv2)
    set(iv2)%axis  = (/'lat  ','p','fcst','t'/)
    set(iv2)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set(iv2)%axis1(set(iv2)%nd(1)) )
    allocate( set(iv2)%axis2(set(iv2)%nd(2)) )
    allocate( set(iv2)%axis3(set(iv2)%nd(3)) )
    allocate( set(iv2)%axis4(set(iv2)%nd(4)) )

    set(iv2)%axis1 = lato
!    set(iv2)%axis2 = p
    set(iv2)%axis3 = t_fct
    set(iv2)%axis4 = t
    if (iv < 5) then
      set(iv2)%axis2 = p
    else
      set(iv2)%axis2 = -999.
      set(iv2)%axis(2) = '    '
    end if
  enddo


  ENDDO  L_VAR


  allocate( set(nvo)%var_out(nd1a,1,1,1) )
  set(nvo)%var_out(:,1,1,1) = cosphi(:)
  set(nvo)%vname = ovarname(nvo)
  set(nvo)%axis  = (/'lat  ',' ',' ',' '/)
  set(nvo)%nd(:) = (/nd1a,1,1,1/)
  allocate( set(nvo)%axis1(set(nvo)%nd(1)) )
  allocate( set(nvo)%axis2(set(nvo)%nd(2)) )
  allocate( set(nvo)%axis3(set(nvo)%nd(3)) )
  allocate( set(nvo)%axis4(set(nvo)%nd(4)) )
  set(nvo)%axis1 = lato
  set(nvo)%axis2 = -999.
  set(nvo)%axis3 = -999.
  set(nvo)%axis4 = -999.

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nvo,set,'ZME2')

! END

  deallocate( t_fct, t )
  deallocate( lon, lat, p, lonu, latv )
  deallocate( lato, cosphi )
  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3, set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  STOP


  CONTAINS

  SUBROUTINE gettaxis

    ! t
    nt = (date(2)-date(1)+1)*nutc
    allocate( t(nt) )
    n = 0
    do idat=date(1), date(2)
    do iutc=1, nutc
      n = n + 1
      t(n) = float(idat) + float(utc(iutc))/24.
    enddo
    enddo

    ! fct
    nfct = int((fct(2)-fct(1))/fct_itv)+1
    allocate( t_fct(nfct) )
    do ifct=1, nfct
      t_fct(ifct) = fct(1) + fct_itv*(ifct-1)
    enddo

  END subroutine gettaxis

  SUBROUTINE getfname

    if (iv < 5) then
      fname = trim(file_i)//'/'//timec//'/'// &
              trim(expname)//'.std_p.'//timec//'+'//hhh//'.nc'

      ftrue = trim(file_i)//'/'//ttimec//'/'// &
              trim(expname)//'.std_p.'//ttimec//'+000.nc'
    else
!      fname = trim(file_i)//'/'//timec//'/'//timec//'.sfc.nc'
!      ftrue = trim(file_i)//'/'//ttimec//'/'//ttimec//'.sfc.nc'
    end if

  END subroutine getfname

  SUBROUTINE getdim

    call diminfop(ncid,.TRUE., nx,ny,nz,nxu,nyv,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( p  (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)

    nxn = nx/xy_itv(1)
    nyn = (ny-1)/xy_itv(2)+1
    allocate( lato(nyn), cosphi(nyn) )
    do j=1, nyn
      lato(j) = lat((j-1)*xy_itv(2)+1)
    enddo
    cosphi(:)   = cos(lato(:)*deg2rad)
    cosphi(1)   = 0.
    cosphi(nyn) = 0.

    l_getdim = .FALSE.

  END subroutine getdim

END program ZME2


SUBROUTINE fct_target(year0,month0,date0,hour0,fct, year,month,date,hour)

  implicit none

  integer, intent(in)  ::  year0, month0, date0, hour0, fct
  integer, intent(out) ::  year, month, date, hour

  integer, dimension(12) ::  mondate = (/31,28,31,30,31,30,31,31,30,31,30,31/)


  year  = year0
  month = month0
  date  = date0
  hour  = hour0 + fct

  do while (hour < 0)
    hour = hour + 24
    date = date - 1
  enddo
  do while (hour >= 24)
    hour = hour - 24
    date = date + 1
  enddo
  do while (date < 1)
    if (month == 1) then
      date = date + mondate(12)
    else if ( month == 2 .and. mod(year,4) == 0 ) then
      date = date + 29
    else
      date = date + mondate(month-1)
    end if
    month = month - 1
    if (month == 0) then
      month = 12
      year = year - 1
    end if
  enddo
  do while (date > mondate(month))
    if ( month == 2 .and. mod(year,4) == 0 ) then
      date = date - 29
    else
      date = date - mondate(month)
    end if
    month = month + 1
    if (month == 13) then
      month = 1
      year = year + 1
    end if
  enddo

  RETURN

END subroutine fct_target


SUBROUTINE check_ex(fname,existence)

  implicit none

  character(len=*), intent(in) ::  fname

  logical ::  existence


  inquire(file=trim(fname), exist=existence)
  if (.not. existence)  print*, '    ',trim(fname),' not found. - passed'

  RETURN

END subroutine check_ex


SUBROUTINE zm_err2(nx,ny,nt,itvx,itvy,itvt,var,var0, zme,zmse)

  implicit none

  integer,                   intent(in) ::  nx, ny, nt, itvx, itvy, itvt
  real, dimension(nx,ny,nt), intent(in) ::  var, var0

  real, dimension((ny-1)/itvy+1), intent(out) ::  zme, zmse

  integer                  ::  i,j,n, ii,jj,nn
  integer                  ::  nxo, nyo, nto
  real, dimension(nx/itvx) ::  err, err2


  nyo = (ny-1)/itvy+1
  nxo = nx/itvx  ;  nto = nt/itvt

  zme (:) = 0.
  zmse(:) = 0.
  do n=1, nto
    nn = n*itvt
    do j=1, nyo
      jj = (j-1)*itvy+1
      do i=1, nxo
        ii = (i-1)*itvx+1
        err(i) = var(ii,jj,nn) - var0(ii,jj,nn)
      enddo
      err2(:) = err(:)*err(:)
      zme (j) = zme (j) + sum(err (:))/nxo
      zmse(j) = zmse(j) + sum(err2(:))/nxo
    enddo
  enddo
  zme (:) = zme (:)/nto
  zmse(:) = zmse(:)/nto

  RETURN

END subroutine zm_err2


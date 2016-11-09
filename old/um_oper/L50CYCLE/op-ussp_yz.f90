! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2010, month =  1
  integer, parameter ::  dstart =  6, dend = 31
  integer, parameter ::  nutc  =  1, ustart = 0
  real,    parameter ::  fstart = 1.0, fend = 5.0, fitv = 1.0   ! [day]
  integer, parameter ::  nvar = 6
  real,    parameter ::  missv = 1.e32
  real,    parameter ::  zh_int = 12.0
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070812'/
  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'gwdc'/
  data                   vartype /'ussp'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/ 
  data                   outname /'ussp_yz'/
  ! input variables
  data                   ivarname(1) /'field424'/
  data                   ivarname(2) /'field425'/
  data                   ivarname(3) /'field420'/
  data                   ivarname(4) /'field422'/
  data                   ivarname(5) /'field423'/
  data                   ivarname(6) /'field421'/
  ! output variables
  data                   ovarname(1) /'drag_x'/
  data                   ovarname(2) /'drag_y'/
  data                   ovarname(3) /'mf_e'/
  data                   ovarname(4) /'mf_w'/
  data                   ovarname(5) /'mf_n'/
  data                   ovarname(6) /'mf_s'/


  real, dimension(:,:,:,:),   allocatable ::  var, var_za
  real, dimension(:),         allocatable ::  lon, lat, zh, lonu, latv, zht
  real, dimension(:),         allocatable ::  time0

  integer ::  nfct, i_time, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid, ncid2, tag
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=6)   ::  cmonth0
  character(len=10)  ::  cutc0
  character(len=128) ::  fname

  real, parameter ::  v_small = 1.e-5

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  nfct = int((fend-fstart)/fitv)+1

  nt = (dend-dstart+1)*nutc*nfct

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month


  N_VAR:   DO iv=1, nvar

  i_time = 0


  N_DAT:   DO idt=1, dend-dstart+1
  N_UTC:   DO iut=1, nutc


  date = dstart + (idt-1)
  utc  = ustart + (iut-1)*24/nutc


  N_FCT:   DO ifc=1, nfct


  i_time = i_time + 1

  fct = int((fstart + (ifc-1)*fitv)*24.)
  call fct0_time(year,month,date,utc,fct, year0,month0,date0,utc0)

  write(cmonth0,'(i4.4,i2.2)') year0, month0
  write(cutc0,  '(a6,2i2.2)' ) cmonth0, date0, utc0

  write(6,*)
  write(6,'(a10,a3,i4.4)') cutc0, ' + ', fct
  write(6,*)

  tag = 0
  do n=1, nexcept
    if (cutc0 == excdate(n))  tag = 1
  enddo
  if ( tag ) then
    write(6,*) ' Exception - passed.'
    i_time = i_time - 1
    CYCLE
  end if

  fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          trim(vartype)//'.'//cutc0//'.a24h.nc'
  if (fitv == 0.5) then
    fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
            trim(vartype)//'.'//cutc0//'.a12h.nc'
  end if

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (iv*idt*iut*ifc == 1) then
    call diminfo(ncid,.TRUE.,zh_int, nx,ny,nz,nxu,nyv,nzt,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( zh (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    if (trim(c_axis(3,2)) /= empty)  allocate( zht (nzt) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,nzt,c_axis, lon,lat,zh,lonu,latv,zht)
  end if

  call timeinfo(ncid, nt0)
  allocate( time0(nt0) )
  call get1d(ncid,'t',nt0, time0)
  do n=1, nt0
    if ( abs(time0(n) - (real(fct+utc0)/24.-0.5*fitv)) < v_small ) then
      it = n
    end if
  enddo
  deallocate( time0 )


  ! get var
  nxr = NX   ;   nxa =  1
  nyr = NY   ;   nya = NY
  nzr = NZ   ;   nza = NZ  ! ;   iz = 28
  ntr = NT   ;   nta =  1

  if ( iv == 2 .or. iv == 5 .or. iv == 6 ) then
    nyr = NYV   ;   nya = NYV
  else
    nyr = NY   ;   nya = NY
  end if

  if (idt*iut*ifc == 1) then
    allocate( set(iv)%var_out(nxa,nya,nza,nta) )
    allocate( var_za(1,nyr,nzr,ntr) )
  end if

  allocate( var(nxr,nyr,nzr,1) )
  call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var)

  call zonal_avg(nxr,nyr,nzr,1,missv,var, var_za(:,:,:,i_time))
  deallocate ( var )

  if (idt*iut*ifc == 1) then
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis4 = (year-2000)*100.+month
  end if

  call closenc(ncid)


  ENDDO  N_FCT


  ENDDO  N_UTC
  ENDDO  N_DAT


  call tempo_avg(1,nyr,nzr,i_time,missv,var_za(:,:,:,:i_time), set(iv)%var_out)
  deallocate( var_za )


  ENDDO  N_VAR


  set(1)%axis = (/'     ','lat ','zh ','t'/)
  set(2)%axis = (/'     ','latv','zh ','t'/)
  set(3)%axis = (/'     ','lat ','zht','t'/)
  set(4)%axis = (/'     ','lat ','zht','t'/)
  set(5)%axis = (/'     ','latv','zht','t'/)
  set(6)%axis = (/'     ','latv','zht','t'/)

  set(1)%axis1 = -999.
  set(1)%axis2 = lat
  set(1)%axis3 = zh
  set(2)%axis1 = -999.
  set(2)%axis2 = latv
  set(2)%axis3 = zh
  do iv=3, 4
    set(iv)%axis1 = -999.
    set(iv)%axis2 = lat
    set(iv)%axis3 = zht
  enddo
  do iv=5, 6
    set(iv)%axis1 = -999.
    set(iv)%axis2 = latv
    set(iv)%axis3 = zht
  enddo

  set(1)%var_out = set(1)%var_out * 86400.   ! [m/s/day]
  set(2)%var_out = set(2)%var_out * 86400.


  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'

  write(6,*)
  write(6,*) trim(fname)
  write(6,*)

  call outnc(trim(fname),nvar,set,'')

  do iv=1, nvar
    deallocate( set(iv)%var_out )
  enddo


  STOP

END program


SUBROUTINE fct0_time(year,month,date,utc,fct, year0,month0,date0,utc0)

  implicit none

  integer, intent(in)  ::  year, month, date, utc, fct
  integer, intent(out) ::  year0, month0, date0, utc0

  integer, dimension(12) ::  mondate = (/31,28,31,30,31,30,31,31,30,31,30,31/)


  year0  = year
  month0 = month
  date0  = date
  utc0   = utc - fct

  do while (utc0 < 0)
    utc0 = utc0 + 24
    date0 = date0 - 1
  enddo
  do while (date0 < 1)
    if (month0 == 1) then
      date0 = date0 + mondate(12)
    else if ( month0 == 2 .and. mod(year0,4) == 0 ) then
      date0 = date0 + 29
    else
      date0 = date0 + mondate(month0-1)
    end if
    month0 = month0 - 1
    if (month0 == 0) then
      month0 = 12
      year0 = year0 - 1
    end if
  enddo


  RETURN

END subroutine fct0_time


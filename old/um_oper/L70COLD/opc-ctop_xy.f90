! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2010, month = 1
  integer, parameter ::  dstart = 1, dend = 26
  integer, parameter ::  nutc  =  1, ustart = 0
  real,    parameter ::  fstart = 5.125, fend = 10.0, fitv = 0.125   ! [day]
  integer, parameter ::  nvar = 5
  real,    parameter ::  missv = 1.e32
  real,    parameter ::  zh_int = 12.0
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070812'/
  ! files
  data                   ifdir   /'/data4/kyh/UM_cold/UM_OUT/'/
  data                   expname /'gwdc'/
  data                   vartype /'gwdc'/
  data                   ofdir   /'/data14/kyh/UM_cold/UM_OUT/'/ 
  data                   outname /'ctop_xy'/
  ! input variables
  data                   ivarname( 1) /'unspecified_2'/
  data                   ivarname( 2) /'unspecified_3'/
  data                   ivarname( 3) /'unspecified_4'/
  data                   ivarname( 4) /'unspecified_5'/
  data                   ivarname( 5) /'unspecified_6'/
  ! output variables
  data                   ovarname( 1) /'mf_ct_e '/
  data                   ovarname( 2) /'mf_ct_w '/
  data                   ovarname( 3) /'mf_ct_n '/
  data                   ovarname( 4) /'mf_ct_s '/
  data                   ovarname( 5) /'heatmax'/


  real, dimension(:,:,:,:),   allocatable ::  var
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
  character(len=3)   ::  hhh

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
  year0  = year
  month0 = month
  date0  = date
  utc0   = utc

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

  write(hhh,'(i3.3)') fct
  fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          trim(expname)//'.'//trim(vartype)//'.'//cutc0//'+'//hhh//'.nc'

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
  it = 1


  ! get var
  nxr = NX   ;   nxa = NX
  nyr = NY   ;   nya = NY
  nzr =  1   ;   nza =  1   ;   iz = 1 !29
  ntr = NT   ;   nta =  1

  if (idt*iut*ifc == 1) then
    allocate( set(iv)%var_out(nxa,nya,nza,nta) )
    allocate( var(nxr,nyr,nzr,ntr) )
  end if

!  if ( iv >= 9 .and. iv <= 12 ) then
!    call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,iz,nzr,it,1, var(:,:,:,i_time))
!  else
    call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var(:,:,:,i_time))
!  end if

  if (idt*iut*ifc == 1) then
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%vname = ovarname(iv)
    set(iv)%axis  = (/'lon  ','lat ','   ','t'/)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = lon
    set(iv)%axis2 = lat
    set(iv)%axis3 = -999.
    set(iv)%axis4 = (year-2000)*100.+month
  end if

  call closenc(ncid)


  ENDDO  N_FCT


  ENDDO  N_UTC
  ENDDO  N_DAT


  call tempo_avg(nxr,nyr,nzr,i_time,missv,var(:,:,:,:i_time), set(iv)%var_out)
  deallocate( var )


  ENDDO  N_VAR


  set(5)%var_out = set(5)%var_out * 86400.   ! [K/day]


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


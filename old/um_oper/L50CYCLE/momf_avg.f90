! temporal average of monthly results
! use the results of "momf".

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year1  = 2008, year2  = 2010
  integer, parameter ::  month1 =    7, month2 =    6
  integer, parameter ::  nvar = 9
  real,    parameter ::  missv = 2.e20
  real,    parameter ::  zh_int = 12.0
  character(len=128) ::  ifdir, expname, vartype, inpname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data5/kyh/umres/gwdcx_f200/anal/'/
  data                   expname /'gwdcx'/
  data                   inpname /'f200.momf'/
  ! input variables
  data                   ivarname(1) /'mf_ct_e '/
  data                   ivarname(2) /'mf_ct_w '/
  data                   ivarname(3) /'mf_ct_n '/
  data                   ivarname(4) /'mf_ct_s '/
  data                   ivarname(5) /'mf_z16_e'/
  data                   ivarname(6) /'mf_z16_w'/
  data                   ivarname(7) /'mf_z16_n'/
  data                   ivarname(8) /'mf_z16_s'/
  data                   ivarname(9) /'heatmax'/
  ! output variables
  data                   ovarname(1) /'mf_ct_e '/
  data                   ovarname(2) /'mf_ct_w '/
  data                   ovarname(3) /'mf_ct_n '/
  data                   ovarname(4) /'mf_ct_s '/
  data                   ovarname(5) /'mf_z16_e'/
  data                   ovarname(6) /'mf_z16_w'/
  data                   ovarname(7) /'mf_z16_n'/
  data                   ovarname(8) /'mf_z16_s'/
  data                   ovarname(9) /'heatmax'/


  real, dimension(:,:,:,:,:), allocatable ::  var
  real, dimension(:,:,:,:,:), allocatable ::  var_avg
  real, dimension(:),         allocatable ::  lon, lat, zh, lonu, latv, zht

  integer ::  year, month, nmonth, imonth
  integer ::  nx, ny, nz, nxu, nyv, nzt
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz
  integer ::  i,j,k,n,iv, ncid
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear, cyear1, cyear2
  character(len=2)   ::  cmonth, cmonth1, cmonth2
  character(len=128) ::  fname


  write(cyear1, '(i4.4)') year1
  write(cmonth1,'(i2.2)') month1
  write(cyear2, '(i4.4)') year2
  write(cmonth2,'(i2.2)') month2

  nmonth = (year2-year1)*12 + month2-month1+1

  year  = year1
  month = month1 - 1
  imonth = 0


  MONTHLY:  DO WHILE ( year /= year2 .or. month /= month2 )


  imonth = imonth + 1

  month = month + 1
  if (month == 13) then
    month = 1
    year = year + 1
  end if

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month

  fname = trim(ifdir)//trim(expname)//'.'//trim(inpname)//'.'//cyear//'.'//cmonth//'.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (imonth == 1) then
    call diminfo(ncid,.FALSE.,zh_int, nx,ny,nz,nxu,nyv,nzt,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( zh (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    if (trim(c_axis(3,2)) /= empty)  allocate( zht (nzt) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,nzt,c_axis, lon,lat,zh,lonu,latv,zht)
  end if


  ! get var
  nxr = NX       ;   nxa = NX
  nyr = NY       ;   nya = NY
  nzr =  1       ;   nza =  1
  ntr = NMONTH   ;   nta =  1

  if (imonth == 1)  allocate( var(nxr,nyr,nzr,ntr,nvar) )

  do iv=1, nvar
    call get4d(ncid,trim(ivarname(iv)),nxr,nyr,nzr,  1, var(:,:,:,imonth,iv))
  enddo

  call closenc(ncid)


  ENDDO  MONTHLY


  if (imonth /= nmonth)  print*, ' ERROR !!!'


  allocate( var_avg(nxa,nya,nza,nta,nvar) )
  do iv=1, nvar
    call tempo_avg(nxr,nyr,nzr,ntr,missv,var(:,:,:,:,iv), var_avg(:,:,:,:,iv))
  enddo

  ! dump
  fname = 'res/'//trim(expname)//'.'//trim(inpname)//'.'// &
                cyear1//cmonth1//'-'//cyear2//cmonth2//'.nc'
  call out3d(trim(fname),nvar,ovarname,var_avg,                        &
             'lon',nxa,lon,'lat',nya,lat,                              &
             't',1,(/(year1+year2)*0.5/),'')



  STOP

END program


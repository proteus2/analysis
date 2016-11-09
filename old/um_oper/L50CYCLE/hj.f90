PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  nmonth =  1, mstart =  7
  integer, parameter ::  nyear  =  1, ystart = 2008
  integer, parameter ::  nvar = 8
  real,    parameter ::  missv = 2.e20
  real,    parameter ::  zh_int = 12.0
  character(len=128) ::  ifdir, expname, vartype, outname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! files
  data                   ifdir   /'/export31/chj/isgwdclm2-veres/monthly_new/sgwdclm2_veres.waccm.mon.'/
  data                   expname /'xf20'/
  data                   vartype /'gwdc'/
  data                   outname /'gwdc_yz'/
  ! input variables
  data                   ivarname(1) /'TAUCTXP '/
  data                   ivarname(2) /'TAUCTXN '/
  data                   ivarname(3) /'TAUCTYP '/
  data                   ivarname(4) /'TAUCTYN '/
  data                   ivarname(5) /'TAUCTXPF'/
  data                   ivarname(6) /'TAUCTXNF'/
  data                   ivarname(7) /'TAUCTYPF'/
  data                   ivarname(8) /'TAUCTYNF'/
  ! output variables
  data                   ovarname(1) /'TAUCTXP '/
  data                   ovarname(2) /'TAUCTXN '/
  data                   ovarname(3) /'TAUCTYP '/
  data                   ovarname(4) /'TAUCTYN '/
  data                   ovarname(5) /'TAUCTXPF'/
  data                   ovarname(6) /'TAUCTXNF'/
  data                   ovarname(7) /'TAUCTYPF'/
  data                   ovarname(8) /'TAUCTYNF'/


  real, dimension(:,:,:,:),   allocatable ::  var
  real, dimension(:),         allocatable ::  lon, lat, zh, lonu, latv, zht

  integer ::  year, month
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, imn, iyr
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz
  integer ::  i,j,k,n,iv, ncid
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=128) ::  fname

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  N_MON:   DO imn=1, nmonth
  N_YEAR:  DO iyr=1, nyear


  year  = ystart + (iyr-1)
  month = mstart + (imn-1)
  if (month > 12)  month = month - 12
  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month

  fname = trim(ifdir)//cyear//'-'//cmonth//'.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis

  if (imn == 1 .and. iyr == 1) then
    call diminfo(ncid,.TRUE.,zh_int, nx,ny,nz,nxu,nyv,nzt,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( zh (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    if (trim(c_axis(3,2)) /= empty)  allocate( zht (nzt) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,nzt,c_axis, lon,lat,zh,lonu,latv,zht)
  end if
  call timeinfo(ncid, nt)


  ! get var
  nxr = NX   ;   nxa =  1
  nyr = NY   ;   nya = NY
  nzr = NZT  ;   nza = NZT ! ;   iz = 28
  ntr = NT   ;   nta =  1

  do iv=1, nvar

    if (iv == 2) then
      nyr = NYV   ;   nya = NYV
    else
      nyr = NY   ;   nya = NY
    end if

    allocate( var(nxr,nyr,nzr,ntr) )
    call get4d(ncid,trim(ivarname(iv)),nxr,nyr,nzr,ntr, var)

    allocate( set(iv)%var_out(nxa,nya,nza,nta) )
    call zonal_tempo_avg(nxr,nyr,nzr,ntr,missv,var, set(iv)%var_out)

    deallocate ( var )

    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta

  enddo

  set(1)%axis = (/'     ','lat ','zh ','t'/)
  set(2)%axis = (/'     ','latv','zh ','t'/)
  set(3)%axis = (/'     ','lat ','zht','t'/)
  set(4)%axis = (/'     ','lat ','zht','t'/)
  set(5)%axis = (/'     ','lat ','zht','t'/)
  set(6)%axis = (/'     ','lat ','zht','t'/)

  do iv=1, nvar
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    if (iv == 1) then
      set(iv)%axis1 = -999.
      set(iv)%axis2 = lat
      set(iv)%axis3 = zh
    else if (iv == 2) then
      set(iv)%axis1 = -999.
      set(iv)%axis2 = latv
      set(iv)%axis3 = zh
    else
      set(iv)%axis1 = -999.
      set(iv)%axis2 = lat
      set(iv)%axis3 = zht
    end if
    set(1)%axis4 = (year-2000)*100.+month
  enddo
  call closenc(ncid)

  set(1)%var_out = set(1)%var_out * 86400.   ! [m/s/day]
  set(2)%var_out = set(2)%var_out * 86400.

  ! dump
  fname = 'res/'//trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
  call outnc(trim(fname),nvar,set,'drag [m/s/day], flux [N/m2]')


  ENDDO  N_YEAR
  ENDDO  N_MON


  STOP

END program


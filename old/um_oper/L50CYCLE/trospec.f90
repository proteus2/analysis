! temporal average of monthly results
! use the results of "uvt_yz".

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis
  use pwrspd

  implicit none

  integer, parameter ::  year1  = 2010, year2  = 2015
  integer, parameter ::  month1 =    1, month2 =   12
  real,    parameter ::  lat1   = -12.5, lat2   =  12.5, dt = 1.
  integer, parameter ::  nvar = 1
  real,    parameter ::  missv = 2.e20
  real,    parameter ::  zh_int = 12.0
  character(len=128) ::  ifdir, expname, vartype, inpname, outname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data5/kyh/umres/'/
  data                   expname /'ctl'/
  data                   inpname /'uvt_yz'/
  data                   outname /'uspec_oz'/
  ! input variables
  data                   ivarname(1) /'u'/
  ! output variables
  data                   ovarname(1) /'psd_u'/


  real, dimension(:,:,:,:,:), allocatable ::  var
  real, dimension(:,:,:,:,:), allocatable ::  var_avg, psd
  real, dimension(:),         allocatable ::  lon, lat, zh, lonu, latv, zht
  real, dimension(:),         allocatable ::  time, oo

  integer ::  year, month, nmonth, imonth
  integer ::  nx, ny, nz, nxu, nyv, nzt
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  js, cnt, iy
  integer ::  i,j,k,n,iv,noo, ncid
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear, cyear1, cyear2
  character(len=2)   ::  cmonth, cmonth1, cmonth2
  character(len=128) ::  fname


  write(cyear1, '(i4.4)') year1
  write(cmonth1,'(i2.2)') month1
  write(cyear2, '(i4.4)') year2
  write(cmonth2,'(i2.2)') month2

  nmonth = (year2-year1)*12 + month2-month1+1
  allocate( time(nmonth) )

  year  = year1
  month = month1 - 1
  imonth = 0


  MONTHLY:  DO WHILE ( year /= year2 .or. month /= month2 )


  imonth = imonth + 1
  time(imonth) = float(imonth)*dt
  month = month + 1
  if (month == 13) then
    month = 1
    year = year + 1
  end if

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month

  fname = trim(ifdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(inpname)//'.'//cyear//'.'//cmonth//'.nc'

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
  cnt = 0
  do j=1, ny
    if ( lat(j) >= lat1 .and. lat(j) <= lat2 ) then
      cnt = cnt + 1
      if (cnt == 1)  js = j
    end if
  enddo

  nxr =  1       ;   nxa =  1
  nyr = CNT      ;   nya =  1   ;   iy = JS
  nzr = NZ       ;   nza = NZ
  ntr = NMONTH   ;   nta = NMONTH

  if (imonth == 1)  allocate( var(nxr,nyr,nzr,ntr,nvar) )

  do iv=1, nvar
    call geta3d(ncid,trim(ivarname(iv)),iy,nyr,1,nzr,1,1, var(:,:,:,imonth,iv))
  enddo

  call closenc(ncid)


  ENDDO  MONTHLY


  if (imonth /= nmonth)  print*, ' ERROR !!!'


  allocate( var_avg(nxa,nya,nza,nta,nvar) )
  do iv=1, nvar
    call merid_avg(nxr,nyr,nzr,ntr,missv,lat(iy:iy+nyr-1), &
                   var(:,:,:,:,iv), var_avg(:,:,:,:,iv))
  enddo

  deallocate( var )


  noo = nta/2+1
  allocate( oo(noo) )
  allocate( psd(noo,nxa,nya,nza,nvar) )

  do iv=1, nvar
  do k=1, nza
  do j=1, nya
  do i=1, nxa
    call psd1d(nta,3,time,dt,var_avg(i,j,k,:,iv),oo,psd(:,i,j,k,iv),1,0)
  enddo
  enddo
  enddo
  enddo


  ! dump
  fname = 'res/'//trim(expname)//'.'//trim(outname)//'.'// &
                cyear1//cmonth1//'-'//cyear2//cmonth2//'.nc'
  print*, fname
  call out2d(trim(fname),nvar,ovarname,psd, &
             'freq',noo,oo,'zh',nza,zh,'starting from '//cyear1//'.'//cmonth1)



  STOP

END program


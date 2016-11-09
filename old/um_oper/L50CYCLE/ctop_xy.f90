PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  nmonth =  1, mstart =  7
  integer, parameter ::  nyear  =  1, ystart = 2008
  integer, parameter ::  nvar = 13
  real,    parameter ::  missv = 2.e20
  real,    parameter ::  zh_int = 12.0
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data5/kyh/umres/'/
  data                   expname /'test'/
  data                   vartype /'gwdc'/
  data                   ofdir   /'/data5/kyh/umres/'/
  data                   outname /'ctop_xy'/
  ! input variables
  data                   ivarname( 1) /'unspecified_6'/
  data                   ivarname( 2) /'unspecified_7'/
  data                   ivarname( 3) /'unspecified_8'/
  data                   ivarname( 4) /'unspecified_9'/
  data                   ivarname( 5) /'unspecified_2'/
  data                   ivarname( 6) /'unspecified_3'/
  data                   ivarname( 7) /'unspecified_4'/
  data                   ivarname( 8) /'unspecified_5'/
  data                   ivarname( 9) /'unspecified_2'/
  data                   ivarname(10) /'unspecified_3'/
  data                   ivarname(11) /'unspecified_4'/
  data                   ivarname(12) /'unspecified_5'/
  data                   ivarname(13) /'unspecified_10'/
  ! output variables
  data                   ovarname( 1) /'mf_ct_e '/
  data                   ovarname( 2) /'mf_ct_w '/
  data                   ovarname( 3) /'mf_ct_n '/
  data                   ovarname( 4) /'mf_ct_s '/
  data                   ovarname( 5) /'mf_ct0_e'/
  data                   ovarname( 6) /'mf_ct0_w'/
  data                   ovarname( 7) /'mf_ct0_n'/
  data                   ovarname( 8) /'mf_ct0_s'/
  data                   ovarname( 9) /'mf_z17_e'/
  data                   ovarname(10) /'mf_z17_w'/
  data                   ovarname(11) /'mf_z17_n'/
  data                   ovarname(12) /'mf_z17_s'/
  data                   ovarname(13) /'heatmax'/


  real, dimension(:,:,:,:),   allocatable ::  var
  real, dimension(:,:,:,:,:), allocatable ::  var_avg
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


  N_MON:   DO imn=1, nmonth
  N_YEAR:  DO iyr=1, nyear


  year  = ystart + (iyr-1)
  month = mstart + (imn-1)
  if (month > 12)  month = month - 12
  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month

  fname = trim(ifdir)//trim(expname)//'/'//cyear//'/'//                     &
          trim(expname)//'.'//cyear//'.'//cmonth//'.'//trim(vartype)//'.nc'

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
  nxr = NX   ;   nxa = NX
  nyr = NY   ;   nya = NY
  nzr =  1   ;   nza =  1   ;   iz = 29
  ntr = NT   ;   nta =  1

  allocate( var_avg(nxa,nya,nza,nta,nvar) )

  do iv=1, nvar

    allocate( var(nxr,nyr,nzr,ntr) )

    if ( iv >= 5 .and. iv <= 8 ) then
      call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,1,ntr, var)
    else if ( iv >= 9 .and. iv <= 12 ) then
      call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,iz,nzr,1,ntr, var)
    else
      call get4d(ncid,trim(ivarname(iv)),nxr,nyr,nzr,ntr, var)
    end if

    call tempo_avg(nxr,nyr,nzr,ntr,missv,var, var_avg(:,:,:,:,iv))

    deallocate ( var )

  enddo

  call closenc(ncid)


  var_avg(:,:,:,:,13) = var_avg(:,:,:,:,13) * 86400.   ! [K/day]

  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
  print*, fname
  call out3d(trim(fname),nvar,ovarname,var_avg,                        &
             'lon',nxa,lon,'lat',nya,lat,                              &
             't',1,(/(year-2000)*100.+month/),'heatmax [K/day]')


  ENDDO  N_YEAR
  ENDDO  N_MON


  STOP

END program


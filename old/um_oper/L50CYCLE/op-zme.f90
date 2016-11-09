! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2010, month = 1
  integer, parameter ::  dstart = 1, dend = 26
  integer, parameter ::  nutc  =  2, ustart = 0
  real,    parameter ::  fstart = 0.5, fend = 5.0, fitv = 0.5   ! [day]
  integer, parameter ::  xitv = 4, yitv = 6
  integer, parameter ::  nivar = 5
  integer, parameter ::  nvar = nivar*4
  real,    parameter ::  missv = 1.e32
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  character(len=128) ::  ifdir, expname, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070000'/
  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'gwdc'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/ 
  data                   outname /'zme'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'temp'/
  data                   ivarname(4) /'ht'/
  data                   ivarname(5) /'p'/
  ! output variables
  data         ovarname(1 :4 ) /'ZME_U','ZMSE_U','cZME_U','cZMSE_U'/
  data         ovarname(5 :8 ) /'ZME_V','ZMSE_V','cZME_V','cZMSE_V'/
  data         ovarname(9 :12) /'ZME_T','ZMSE_T','cZME_T','cZMSE_T'/
  data         ovarname(13:16) /'ZME_Z','ZMSE_Z','cZME_Z','cZMSE_Z'/
  data         ovarname(17:20) /'ZME_P','ZMSE_P','cZME_P','cZMSE_P'/


  real, dimension(:,:,:,:),   allocatable ::  var, var2, var_za, vara, var3
  real, dimension(:),         allocatable ::  var_x, var_x2
  real, dimension(:),         allocatable ::  lon, lat, p, lonu, latv, pt, lato, cosphi
  real, dimension(:),         allocatable ::  time0, t_fct, t

  integer ::  nfct, i_time, year1, month1, date1, utc1, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nd1a, nd2a, nd3a, nd4a, nxn, nyn
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid, ncid2, tag, ii, jj, nn, iv2
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=6)   ::  cmonth0, cmonth1
  character(len=10)  ::  cutc0, cutc1
  character(len=128) ::  fname, fanal

  real, parameter ::  v_small = 1.e-5

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  nfct = int((fend-fstart)/fitv)+1
  allocate( t_fct(nfct) )
  do ifc=1, nfct
    t_fct(ifc) = fstart + fitv*(ifc-1)
  enddo

  nt = (dend-dstart+1)*nutc
  allocate( t(nt) )
  i_time = 0
  do idt=dstart, dend
  do iut=1, nutc
    i_time = i_time + 1
    t(i_time) = float(idt) + float(ustart+(iut-1)*24/nutc)/24.
  enddo
  enddo

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month


  N_VAR:   DO iv=1, nivar

  i_time = 0


  N_DAT:   DO idt=1, dend-dstart+1
  N_UTC:   DO iut=1, nutc


  date = dstart + (idt-1)
  utc  = ustart + (iut-1)*24/nutc

  i_time = i_time + 1


  N_FCT:   DO ifc=1, nfct


  fct = int((fstart + (ifc-1)*fitv)*24.)

!  call fct0_time(year,month,date,utc,fct, year0,month0,date0,utc0)
  year0 = year  ;  month0 = month  ;  date0 = date  ;  utc0 = utc
  call fct0_time(year0,month0,date0,utc0,-fct, year1,month1,date1,utc1)

  write(cmonth0,'(i4.4,i2.2)') year0, month0
  write(cutc0,  '(a6,2i2.2)' ) cmonth0, date0, utc0
  write(cmonth1,'(i4.4,i2.2)') year1, month1
  write(cutc1,  '(a6,2i2.2)' ) cmonth1, date1, utc1

  write(6,*)
  write(6,'(a,a10,a3,i4.4)') trim(ivarname(iv))//' : ',cutc0,' + ',fct
  write(6,*)

  tag = 0
  do n=1, nexcept
    if (cutc0 == excdate(n))  tag = 1
  enddo
  if ( tag ) then
    write(6,*) ' Exception - passed.'
    CYCLE
  end if

  if (iv < 5) then
    fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
            cutc0//'.press1.nc'
    if (fct > 24)  &
       fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
               cutc0//'.press2.nc'
    fanal = trim(ifdir)//trim(expname)//'/'//cmonth1//'/'//cutc1//'/'// &
            cutc1//'.press1.nc'
  else
    fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
            cutc0//'.sfc.nc'
    fanal = trim(ifdir)//trim(expname)//'/'//cmonth1//'/'//cutc1//'/'// &
            cutc1//'.sfc.nc'
  end if


  call opennc(trim(fname),ncid)
  call opennc(trim(fanal),ncid2)

  ! get dim. sizes and axis
  if (iv*idt*iut*ifc == 1) then
    call diminfop(ncid,.TRUE., nx,ny,nz,nxu,nyv,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( p  (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)
  end if

  call timeinfo(ncid, nt0)
  allocate( time0(nt0) )
  call get1d(ncid,'t',nt0, time0)
  do n=1, nt0
    if ( abs(time0(n) - real(fct)/24.) < v_small )  it = n
  enddo
  deallocate( time0 )

  if (iv*idt*iut*ifc == 1) then
    nxn = nx/xitv
    nyn = (ny-1)/yitv+1
    allocate( lato(nyn), cosphi(nyn) )
    do j=1, nyn
      lato(j) = lat((j-1)*yitv+1)
    enddo
    cosphi(:)   = cos(lato(:)*deg2rad)
    cosphi(1)   = 0.
    cosphi(nyn) = 0.
  end if


  ! get var
  nxr = NX   ;   nd1a = NYN
  nyr = NY   ;   nd2a = NZ
  nzr = NZ   ;   nd3a = NFCT  ! ;   iz = 28
  ntr = NT   ;   nd4a = NT
  if (iv >= 5) then
    nzr = 1   ;   nd2a = 1
  end if

  if (idt*iut*ifc == 1) then
    do nn=1, 4
      iv2 = (iv-1)*4 + nn
      allocate( set(iv2)%var_out(nd1a,nd2a,nd3a,nd4a) )
    enddo
  end if

  allocate( var(nxr,nyr,nzr,1), vara(nxr,nyr,nzr,1) )
  if (iv == 2) then
    allocate( var2(nxr,nyv,nzr,1), var3(nxr,nyv,nzr,1) )
    call geta4d(ncid, trim(ivarname(iv)),1,nxr,1,nyv,1,nzr,it,1, var2)
    call geta4d(ncid2,trim(ivarname(iv)),1,nxr,1,nyv,1,nzr, 1,1, var3)
    do k=1, nzr
    do j=2, nyr-1
    do i=1, nxr
      var (i,j,k,1) = 0.5*(var2(i,j-1,k,1)+var2(i,j,k,1))
      vara(i,j,k,1) = 0.5*(var3(i,j-1,k,1)+var3(i,j,k,1))
    enddo
    enddo
    enddo
    var (:,  1,:,:) = 0.
    var (:,nyr,:,:) = 0.
    vara(:,  1,:,:) = 0.
    vara(:,nyr,:,:) = 0.
    deallocate( var2, var3 )
  else
    call geta4d(ncid ,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var )
    call geta4d(ncid2,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr, 1,1, vara)
  end if

  call closenc(ncid )
  call closenc(ncid2)

  if (idt*iut*ifc == 1) then
    allocate( var_x(nxn), var_x2(nxn) )
    allocate( var_za(1,nyn,nzr,2) )
  end if

  do k=1, nzr
  do j=1, nyn
    jj = (j-1)*yitv+1
    do i=1, nxn
      ii = (i-1)*xitv+1
      var_x(i) = var(ii,jj,k,1) - vara(ii,jj,k,1)
    enddo
    var_x2(:) = var_x(:)*var_x(:)
    call zonal_avg(nxn,1,1,1,missv,var_x , var_za(1,j,k,1))
    call zonal_avg(nxn,1,1,1,missv,var_x2, var_za(1,j,k,2))
  enddo
  enddo

  iv2 = (iv-1)*4
  set(iv2+1)%var_out(:,:,ifc,i_time) = var_za(1,:,:,1)
  set(iv2+2)%var_out(:,:,ifc,i_time) = var_za(1,:,:,2)
  do k=1, nd2a
    set(iv2+3)%var_out(:,k,ifc,i_time) =  &
    set(iv2+1)%var_out(:,k,ifc,i_time)*cosphi(:)
    set(iv2+4)%var_out(:,k,ifc,i_time) =  &
    set(iv2+2)%var_out(:,k,ifc,i_time)*cosphi(:)
  enddo

  deallocate ( var, vara )

  if (idt*iut*ifc == 1) then
    do nn=1, 4
      iv2 = (iv-1)*4 + nn
      set(iv2)%nd(1) = nd1a
      set(iv2)%nd(2) = nd2a
      set(iv2)%nd(3) = nd3a
      set(iv2)%nd(4) = nd4a
      set(iv2)%axis  = (/'lat  ','p   ','fcst','t'/)
      set(iv2)%vname = ovarname(iv2)
      allocate( set(iv2)%axis1(set(iv2)%nd(1)) )
      allocate( set(iv2)%axis2(set(iv2)%nd(2)) )
      allocate( set(iv2)%axis3(set(iv2)%nd(3)) )
      allocate( set(iv2)%axis4(set(iv2)%nd(4)) )
      set(iv2)%axis1 = lato
      set(iv2)%axis2 = p
      set(iv2)%axis3 = t_fct
      set(iv2)%axis4 = t
      if (iv < 5) then
        set(iv2)%axis2 = p
      else
        set(iv2)%axis2 = -999.
        set(iv2)%axis(2) = '    '
      end if
    enddo
  end if


  ENDDO  N_FCT


  ENDDO  N_UTC
  ENDDO  N_DAT


  deallocate( var_za, var_x, var_x2 )


  ENDDO  N_VAR


  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/pp/'// &
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
  do while (utc0 >= 24)
    utc0 = utc0 - 24
    date0 = date0 + 1
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
  do while (date0 > mondate(month0))
    if ( month0 == 2 .and. mod(year0,4) == 0 ) then
      date0 = date0 - 29
    else
      date0 = date0 - mondate(month0)
    end if
    month0 = month0 + 1
    if (month0 == 13) then
      month0 = 1
      year0 = year0 + 1
    end if
  enddo


  RETURN

END subroutine fct0_time


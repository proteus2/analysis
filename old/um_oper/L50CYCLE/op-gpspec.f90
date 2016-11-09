! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis
  use fft

  implicit none

  integer, parameter ::  year = 2010, month = 1
  integer, parameter ::  dstart = 12, dend =  12
  integer, parameter ::  nutc  =  1, ustart = 0
  real,    parameter ::  fstart = 5., fend = 5., fitv = 0.5   ! [day]
  integer, parameter ::  wn_tail = 15
  integer, parameter ::  nivar = 1, nvar = 1
  real,    parameter ::  missv = 1.e32
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'ctl'/
  data                   vartype /'press1'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/ 
  data                   outname /'gphspec'/
  ! input variables
  data                   ivarname(1) /'ht'/
  ! output variables
  data                   ovarname(1) /'amp_z'/


  real, dimension(:,:,:,:),     allocatable ::  var, var2, amp
  real, dimension(:),           allocatable ::  lon, lat, p, lonu, latv, pt
  real, dimension(:),           allocatable ::  time0, wn
  double complex, dimension(:), allocatable ::  coef

  integer ::  nfct, i_time, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz, it
  integer ::  i,j,k,n,iv,iv2, iw, ncid, ncid2, tag
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

  allocate( wn(wn_tail) )
  do iw=1, wn_tail
    wn(iw) = iw*1.
  enddo


  N_VAR:   DO iv=1, nivar

  i_time = 0


  N_DAT:   DO idt=1, dend-dstart+1
  N_UTC:   DO iut=1, nutc


  date = dstart + (idt-1)
  utc  = ustart + (iut-1)*24/nutc


  N_FCT:   DO ifc=1, nfct


  i_time = i_time + 1

  fct = int((fstart + (ifc-1)*fitv)*24.)
  call fct0_time(year,month,date,utc,fct, year0,month0,date0,utc0)

  write(cmonth0,'(i4.4,i2.2)' ) year0, month0
  write(cutc0,  '(a6,2i2.2)') cmonth0, date0, utc0

  write(6,*)
  write(6,'(a10,a3,i4.4)') cutc0, ' + ', fct
  write(6,*)

  fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          cutc0//'.'//trim(vartype)//'.nc'

  if ( trim(vartype) == 'press1' .and. fct > 24 )  &
     fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
             cutc0//'.press2.nc'

  call opennc(trim(fname),ncid)

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


  ! get var
  nxr = NX   ;   nxa = WN_TAIL
  nyr = NY   ;   nya = NY
  nzr = NZ   ;   nza = NZ  ! ;   iz = 28
  ntr = NT   ;   nta =  1

  if (idt*iut*ifc == 1) then
    do iv2=1, nvar
      allocate( set(iv2)%var_out(nxa,nya,nza,nta) )
    enddo
  end if

  allocate( var(nxr,nyr,nzr,1) )
  call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var)

  if (idt*iut*ifc == 1)  allocate( amp(nxa,nyr,nzr,ntr) )
  allocate( coef(nxr) )
  do k=1, nzr
  do j=1, nyr
    tag = 0
    do i=1, nxr
      if (var(i,j,k,1) == missv)  tag = 1
    enddo
    if (tag == 0) then
      call fft1df(nxr,var(:,j,k,1),coef)
      do iw=1, nxa
        amp(iw,j,k,i_time) = 2.*cdabs(coef(iw+1))/nxr
      enddo
    else
      amp(:,j,k,i_time) = missv
    end if
  enddo
  enddo
  deallocate ( var )

  if (idt*iut*ifc == 1) then
    do iv2=1, nvar
      set(iv2)%nd(1) = nxa
      set(iv2)%nd(2) = nya
      set(iv2)%nd(3) = nza
      set(iv2)%nd(4) = nta
      set(iv2)%axis = (/'wn   ','lat ','p ','t'/)
      set(iv2)%vname = ovarname(iv2)
      allocate( set(iv2)%axis1(set(iv2)%nd(1)) )
      allocate( set(iv2)%axis2(set(iv2)%nd(2)) )
      allocate( set(iv2)%axis3(set(iv2)%nd(3)) )
      allocate( set(iv2)%axis4(set(iv2)%nd(4)) )
      set(iv2)%axis1 = wn
      set(iv2)%axis2 = lat
      set(iv2)%axis3 = p
      set(iv2)%axis4 = (year-2000)*100.+month
    enddo
  end if

  call closenc(ncid)


  ENDDO  N_FCT


  ENDDO  N_UTC
  ENDDO  N_DAT


  do iv2=1, nvar
    call tempo_avg(nxa,nyr,nzr,ntr,missv,amp, set(iv2)%var_out)
  enddo


  ENDDO  N_VAR


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


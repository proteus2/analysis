! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2010, month =  1
  integer, parameter ::  dstart = 2, dend = 31
  integer, parameter ::  nutc = 2, ustart = 0
  real,    parameter ::  fstart = 1.0, fend = 5.0, fitv = 1.0   ! [day]
  integer, parameter ::  xitv = 4, yitv = 6
  integer, parameter ::  nvar = 4
  real,    parameter ::  missv = 1.e32
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070000'/
  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'ctl'/
  data                   vartype /'press1'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/ 
  data                   outname /'rmse_xa'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'temp'/
  data                   ivarname(4) /'ht'/
  ! output variables
  data                   ovarname(1) /'u'/
  data                   ovarname(2) /'v'/
  data                   ovarname(3) /'temp'/
  data                   ovarname(4) /'z'/


  real, dimension(:,:,:,:),   allocatable ::  var, var2, vara, var3
  real, dimension(:),         allocatable ::  lon, lat, p, lonu, latv, pt, lato, cosphi
  real, dimension(:),         allocatable ::  time0, fcst, time

  integer ::  nfct, i_time, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta, nxn, nyn
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid, ncid2, tag, ii, jj
  real    ::  sumcosphi
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=6)   ::  cmonth0
  character(len=10)  ::  cutc0, cutc
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
  allocate( fcst(nfct) )
  do ifc=1, nfct
    fcst(ifc) = fstart + (ifc-1)*fitv
  enddo

  nt = (dend-dstart+1)*nutc
  allocate( time(nt) )
  do n=1, nt
    time(n) = dstart + (n-1)*(1./nutc) + ustart/24.
  enddo
 

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month


  N_VAR:   DO iv=1, nvar
  N_FCT:   DO ifc=1, nfct


  fct = int((fstart + (ifc-1)*fitv)*24.)

  i_time = 0


  N_DAT:   DO idt=1, dend-dstart+1
  N_UTC:   DO iut=1, nutc


  date = dstart + (idt-1)
  utc  = ustart + (iut-1)*24/nutc


  i_time = i_time + 1

  call fct0_time(year,month,date,utc,fct, year0,month0,date0,utc0)

  write(cmonth0,'(i4.4,i2.2)' ) year0, month0
  write(cutc0,  '(a6,2i2.2)'  ) cmonth0, date0, utc0
  write(cutc,   '(a4,a2,2i2.2)') cyear, cmonth, date, utc

  write(6,*)
  write(6,'(a10,a3,i4.4)') cutc0, ' + ', fct
  write(6,*)

  tag = 0
  do n=1, nexcept
    if (cutc0 == excdate(n))  tag = 1
  enddo
  if (month0 /= month)  tag = 2
  if ( tag /= 0 ) then
    if (tag == 1)  write(6,*) ' Exception - passed.'
    if (tag == 2)  write(6,*) ' No data.'
!    i_time = i_time - 1
    CYCLE
  end if

  fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          cutc0//'.'//trim(vartype)//'.nc'

  if ( trim(vartype) == 'press1' .and. fct > 24 )  &
     fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
             cutc0//'.press2.nc'

  fanal = trim(ifdir)//trim(expname)//'/'//cyear//cmonth//'/'//cutc//'/'// &
          cutc//'.'//trim(vartype)//'.nc'


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
    nyn = ny/yitv+1
    allocate( lato(nyn), cosphi(nyn) )
    do j=1, nyn
      lato(j) = lat((j-1)*yitv+1)
    enddo
    cosphi(:)   = cos(lato(:)*deg2rad)
    cosphi(1)   = 0.
    cosphi(nyn) = 0.
  end if


  ! get var
  nxr = NX   ;   nxa = NFCT
  nyr = NY   ;   nya =  1
  nzr = NZ   ;   nza = NZ  ! ;   iz = 28
  ntr = NT   ;   nta = NT

  if (idt*iut*ifc == 1) then
    allocate( set(iv)%var_out(nxa,nya,nza,nta) )
    set(iv)%var_out(:,:,:,:) = 1.e32
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


  allocate( var2(1,nyn,nzr,1), var3(1,nyn,nzr,1) )

  do j=1, nyn
    jj = (j-1)*yitv+1
    call zonal_avg(nxr,1,nzr,1,missv,var (:,jj,:,1), var2(1,j,:,1))
    call zonal_avg(nxr,1,nzr,1,missv,vara(:,jj,:,1), var3(1,j,:,1))
  enddo

  deallocate( var, vara )

  set(iv)%var_out(ifc,1,:,i_time) = 0.
  sumcosphi = sum(cosphi(:))
  do k=1, nzr
  do j=1, nyn
    set(iv)%var_out(ifc,1,k,i_time) = set(iv)%var_out(ifc,1,k,i_time) + &
                             cosphi(j)*(var2(1,j,k,1)-var3(1,j,k,1))**2
  enddo
  enddo
  set(iv)%var_out(ifc,1,:,i_time) = sqrt(set(iv)%var_out(ifc,1,:,i_time)/sumcosphi)

  deallocate( var2,var3 )


  if (idt*iut*ifc == 1) then
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%axis = (/'FCST ','    ','p ','time'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = fcst
    set(iv)%axis2 = -999.
    set(iv)%axis3 = p
    set(iv)%axis4 = time
  end if


  ENDDO  N_UTC
  ENDDO  N_DAT


  ENDDO  N_FCT


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


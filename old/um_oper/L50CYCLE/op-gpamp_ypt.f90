! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis
  use fft

  implicit none

  integer, parameter ::  year = 2010, month = 1
  integer, parameter ::  dstart = 1, dend = 31
  integer, parameter ::  nutc  =  2, ustart = 0
  real,    parameter ::  fstart = 3.0, fend = 5.0, fitv = 0.5   ! [day]
  integer, parameter ::  wn_tail = 64
  integer, parameter ::  nvar = 1
  real,    parameter ::  missv = 1.e32
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nvar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'0009070812'/
  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'ctl'/
  data                   vartype /'press1'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/ 
  data                   outname /'gpamp_ypt'/
  ! input variables
  data                   ivarname(1) /'ht'/
  ! output variables
  data                   ovarname(1) /'amp_z'/


  real, dimension(:,:,:,:),   allocatable ::  var, amp
  real, dimension(:),         allocatable ::  lon, lat, p, lonu, latv, pt
  real, dimension(:),         allocatable ::  time0, t_fct, t
  double complex, dimension(:), allocatable ::  coef

  integer ::  nfct, i_time, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nd1a, nd2a, nd3a, nd4a
  integer ::  iz, it, iw
  integer ::  i,j,k,n,iv, ncid, ncid2, tag
  real    ::  avgcoef
  real, dimension(wn_tail) ::  wn
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

  do iw=1, wn_tail
    wn(iw) = iw*1.
  enddo

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month


  N_VAR:   DO iv=1, nvar

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

  write(cmonth0,'(i4.4,i2.2)' ) year0, month0
  write(cutc0,  '(a6,2i2.2)') cmonth0, date0, utc0

  write(6,*)
  write(6,'(a10,a3,i4.4)') cutc0, ' + ', fct
  write(6,*)

  avgcoef = 1./float(nfct-1)
  if ( ifc == 1 .or. ifc == nfct )  avgcoef = avgcoef*0.5

  tag = 0
  do n=1, nexcept
    if (cutc0 == excdate(n))  tag = 1
  enddo
  if ( tag ) then
    write(6,*) ' Exception - passed.'
    CYCLE
  end if

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
  nxr = NX   ;   nd1a = NY
  nyr = NY   ;   nd2a = NZ
  nzr = NZ   ;   nd3a = WN_TAIL
  ntr = NT   ;   nd4a = NT

  if (idt*iut*ifc == 1) then
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 0.
  end if

  allocate( var(nxr,nyr,nzr,1) )
  call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var)

  if (idt*iut*ifc == 1)  allocate( amp(nd3a,nyr,nzr,1), coef(nxr) )
  do k=1, nzr
  do j=1, nyr
    tag = 0
    do i=1, nxr
      if ( var(i,j,k,1) == missv .or. amp(1,j,k,1) == missv )  tag = 1
    enddo
    if (tag == 0) then
      call fft1df(nxr,var(:,j,k,1),coef)
      do iw=1, nd3a
        amp(iw,j,k,1) = 2.*cdabs(coef(iw+1))/nxr
      enddo
    else
      amp(:,j,k,1) = missv
    end if
  enddo
  enddo
  deallocate ( var )

  do iw=1, nd3a
    set(iv)%var_out(:,:,iw,i_time) = set(iv)%var_out(:,:,iw,i_time) + &
       avgcoef*amp(iw,:,:,1)
  enddo
  amp(:,:,:,:) = 0.


  if (idt*iut*ifc == 1) then
    set(iv)%nd(1) = nd1a
    set(iv)%nd(2) = nd2a
    set(iv)%nd(3) = nd3a
    set(iv)%nd(4) = nd4a
    set(iv)%axis = (/'lat  ','p   ','wn','t'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = lat
    set(iv)%axis2 = p
    set(iv)%axis3 = wn
    set(iv)%axis4 = t
  end if

  call closenc(ncid)


  ENDDO  N_FCT


  ENDDO  N_UTC
  ENDDO  N_DAT
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


! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2009, month = 7
  integer, parameter ::  dstart = 1, dend = 31   ! 6, 31
  integer, parameter ::  nutc  = 2, ustart = 0
  real,    parameter ::  fstart = 1.0, fend = 5.0, fitv = 1.0   ! [day]
  integer, parameter ::  nvar = 4
  real,    parameter ::  h_scale = 7.e3
  real,    parameter ::  missv = 1.e32
  real,    parameter ::  zh_int = 12.0
  integer, parameter ::  nexcept = 1
  integer, parameter ::  np = 21
  real               ::  p(np)
  character(len=10)  ::  excdate(nexcept)
  character(len=128) ::  ifdir, expname, ofdir, outname
  character(len=64)  ::  vartype(nvar), ivarname(nvar), ovarname(nvar)

  ! pressure levels
  data  p / 1000., 925., 850., 700., 600., 500., 400., 300., 250., 200., &
            150., 100., 70., 50., 30., 20., 10., 5., 3., 2., 1. /
  ! exceptions
  data                   excdate(1) /'2009070812'/
  ! files
  data                   ifdir   /'/data1/kyh/portal/UM_OPER/L50CYCLE/'/
  data                   expname /'GWDC'/
  data                   ofdir   /'/data1/kyh/portal/UM_OPER/anal/L50CYCLE/'/
  data                   outname /'zgwd_ypft'/
  ! input file types
  data                   vartype(1) /'ussp'/
  data                   vartype(2) /'gwdo'/
  data                   vartype(3) /'bldo'/
  data                   vartype(4) /'gwdc'/
  ! input variables
  data                   ivarname(1) /'field424'/
  data                   ivarname(2) /'field68'/
  data                   ivarname(3) /'field1596'/
  data                   ivarname(4) /'unspecified'/
  ! output variables
  data                   ovarname(1) /'fu_ussp'/
  data                   ovarname(2) /'fu_gwdo'/
  data                   ovarname(3) /'fu_bldo'/
  data                   ovarname(4) /'fu_gwdc'/


  real, dimension(:,:,:,:),   allocatable ::  var, var_i, exner
  real, dimension(:),         allocatable ::  lon, lat, zh, lonu, latv, zht
  real, dimension(:),         allocatable ::  fcst, time, time0

  integer ::  nfct, t_cnt, fct
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid, ncid2, tag
  real    ::  temp
  real, dimension(np) ::  zp, exner_o
  character(len=32)   ::  c_axis(3,2)
  character(len=4)    ::  cyear
  character(len=2)    ::  cmonth
  character(len=6)    ::  cmonth0
  character(len=10)   ::  cutc0
  character(len=128)  ::  fname, fnamep

  real, parameter ::  v_small = 1.e-5

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  zp(:) = -h_scale*(log(p(:)/1.e3))
  exner_o(:) = (p(:)/1.e3)**kappa

  nfct = int((fend-fstart)/fitv)+1
  allocate( fcst(nfct) )
  do n=1, nfct
    fcst(n) = fstart + fitv*(real(n-1)-0.5)
  enddo

  nt = (dend-dstart+1)*nutc
  allocate( time(nt) )
  do n=1, nt
    time(n) = dstart + (n-1)*(1./nutc) + ustart/24.
  enddo

  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month


  N_VAR:   DO iv=1, nvar


  t_cnt = 0


  N_DAT:   DO idt=dstart, dend
  N_UTC:   DO iut=ustart, ustart+12*(nutc-1), 12


  t_cnt = t_cnt + 1

  write(cmonth0,'(i4.4,i2.2)' ) year, month
  write(cutc0,  '(a6,2i2.2)') cmonth0, idt, iut

  tag = 0
  do n=1, nexcept
    if (cutc0 == excdate(n))  tag = 1
  enddo
  if ( tag ) then
    write(6,*) ' Exception - passed.'
    set(iv)%var_out(:,:,:,t_cnt) = missv
    CYCLE
  end if

  fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          trim(expname)//'.'//trim(vartype(iv))//'.'//cutc0//'+120_av24h.nc'

  fnamep = trim(ifdir)//trim(expname)//'/'//cmonth0//'/diag/'//cutc0//'/'// &
           'std2.'//cutc0//'.a24h.nc'

  call opennc(trim(fname) ,ncid )
  call opennc(trim(fnamep),ncid2)

  ! get dim. sizes and axis
  if (iv*t_cnt == 1) then
    call diminfo(ncid,.TRUE.,zh_int, nx,ny,nz,nxu,nyv,nzt,c_axis)
    allocate( lat(ny) )
    allocate( lon(nx), zh(nz), lonu(nxu), latv(nyv), zht(nzt) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,nzt,c_axis, lon,lat,zh,lonu,latv,zht)
  end if

  call timeinfo(ncid, nt0)
  allocate( time0(nt0) )
  call get1d(ncid,'t',nt0, time0)


  nxr = NXU  ;   nxa = NY
  nyr = NY   ;   nya = NP
  nzr = NZ   ;   nza = NFCT  ! ;   iz = 28
  ntr = NT   ;   nta = NT

  if (t_cnt == 1)  allocate( set(iv)%var_out(nxa,nya,nza,nta) )


  N_FCT:   DO ifc=1, nfct


  fct = int((fstart + (ifc-1)*fitv)*24.)

  write(6,*)
  write(6,'(a10,a3,i4.4)') cutc0, ' + ', fct
  write(6,*)

  do n=1, nt0
    if ( abs(time0(n) - (real(fct)/24.-0.5*fitv)) < v_small )  it = n
  enddo


  ! get var
  allocate( var(nxr,nyr,nzr,1), var_i(nxr,nyr,np,1), exner(nxr,nyr,nzr,1) )
  call geta4d(ncid ,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var  )
  call geta4d(ncid2,'p_1'             ,1,nx ,1,nyr,1,nzr,it,1, exner)
  exner(:,:,:,:) = (exner(:,:,:,:)/1.e5)**kappa

  do k=1, nz
  do j=1, ny
    temp = exner(1,j,k,1)
    do i=1, nx-1
      exner(i,j,k,1) = 0.5*(exner(i,j,k,1) + exner(i+1,j,k,1))
    enddo
    exner(nx,j,k,1) = 0.5*(exner(nx,j,k,1) + temp)
  enddo
  enddo

  call vintp_p(nxr,ny,nz,exner,var,np,exner_o,1., var_i)

  call zonal_avg(nxr,nyr,np,1,missv,var_i, set(iv)%var_out(:,:,ifc,t_cnt))

  deallocate ( var, var_i, exner )

  if (t_cnt*ifc == 1) then
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%axis = (/'lat ','zp ','fcst','t'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = lat
    set(iv)%axis2 = zp
    set(iv)%axis3 = fcst
    set(iv)%axis4 = time
  end if


  ENDDO  N_FCT


  deallocate( time0 )
  call closenc(ncid )
  call closenc(ncid2)


  ENDDO  N_UTC
  ENDDO  N_DAT


  do n=1, nta
    if (set(iv)%var_out(3,3,1,n) /= missv)  &
       set(iv)%var_out(:,:,:,n) = set(iv)%var_out(:,:,:,n) * 86400.
  enddo


  ENDDO  N_VAR


  ! dump
  fname = trim(ofdir)//trim(expname)//'/tmp/'// &
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


! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2010, month = 1
  integer, parameter ::  dstart =  1, dend = 31
  integer, parameter ::  nutc  =  2, ustart = 0
  real,    parameter ::  fstart = 0.0, fend = 5.0, fitv = 0.5   ! [day]
  integer, parameter ::  xitv = 4, yitv = 6
  integer, parameter ::  nxe = 240, nye = 121, nze = 37
  integer, parameter ::  nivar = 5
  integer, parameter ::  nvar = nivar*4
  real,    parameter ::  missv = 1.e32
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  character(len=128) ::  ifdir, expname, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar), ivare(nivar)

  ! exceptions
  data                   excdate(1) /'2009070000'/
  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'ctl'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/ 
  data                   outname /'zme-era'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'temp'/
  data                   ivarname(4) /'ht'/
  data                   ivarname(5) /'p'/
  data                   ivare   (1) /'u'/
  data                   ivare   (2) /'v'/
  data                   ivare   (3) /'t'/
  data                   ivare   (4) /'z'/
  data                   ivare   (5) /'msl'/
  ! output variables
  data         ovarname(1 :4 ) /'ZME_U','ZMSE_U','cZME_U','cZMSE_U'/
  data         ovarname(5 :8 ) /'ZME_V','ZMSE_V','cZME_V','cZMSE_V'/
  data         ovarname(9 :12) /'ZME_T','ZMSE_T','cZME_T','cZMSE_T'/
  data         ovarname(13:16) /'ZME_Z','ZMSE_Z','cZME_Z','cZMSE_Z'/
  data         ovarname(17:20) /'ZME_P','ZMSE_P','cZME_P','cZMSE_P'/


  real, dimension(:,:,:,:),   allocatable ::  var, var2, var_za
  real, dimension(:,:),       allocatable ::  vare
  real, dimension(:),         allocatable ::  var_x, var_x2
  real, dimension(:),         allocatable ::  lon, lat, p, lonu, latv, pt, lato, cosphi
  real, dimension(:),         allocatable ::  time0, t_fct, t
  integer, dimension(:),      allocatable ::  ple, iz

  integer ::  nfct, i_time, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nd1a, nd2a, nd3a, nd4a, nxn, nyn
  integer ::  it, itera
  integer ::  i,j,k,n,iv, ncid, ncid2, tag, ii, jj, nn, iv2, ke
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=6)   ::  cmonth0
  character(len=10)  ::  cutc0
  character(len=128) ::  fname, fera

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

  write(cmonth0,'(i4.4,i2.2)') year0, month0
  write(cutc0,  '(a6,2i2.2)' ) cmonth0, date0, utc0

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
    fera = '/data10/kyh/ERAinter/ERA.'//cmonth0//'.press.nc'
  else
    fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
            cutc0//'.sfc.nc'
    fera = '/data10/kyh/ERAinter/ERA.'//cmonth0//'.pmsl.nc'
  end if


  call opennc(trim(fname),ncid)
  call opennc(trim(fera),ncid2)

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

  ! for era
  if (iv*idt*iut*ifc == 1) then
    allocate( iz(nz), ple(nze) )
    call iget1d(ncid2,'levelist',nze, ple)
    do k=1, nz
      iz(k) = 0
      do ke=1, nze
        if (abs(p(k)-float(ple(ke))) < min(p(k),float(ple(ke)))*1.e-3)  &
           iz(k) = nze-ke+1   ! ke-index will be reversed in get_era4d
      enddo
    enddo
    deallocate( ple )
  end if

  itera = date0*2-1 + (utc0+fct)/12
  if (itera > 62)  EXIT

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
  nzr = NZ   ;   nd3a = NFCT
  ntr = NT   ;   nd4a = NT
  if (iv >= 5) then
    nzr = 1   ;   nd2a = 1
    iz(:) = 0   ;   iz(1) = 1
  end if

  if (idt*iut*ifc == 1) then
    do nn=1, 4
      iv2 = (iv-1)*4 + nn
      allocate( set(iv2)%var_out(nd1a,nd2a,nd3a,nd4a) )
      set(iv2)%var_out(:,:,:,:) = 1.e32
    enddo
  end if

  allocate( var(nxr,nyr,nzr,1) )
  if (iv == 2) then
    allocate( var2(nxr,nyv,nzr,1) )
    call geta4d(ncid, trim(ivarname(iv)),1,nxr,1,nyv,1,nzr,it,1, var2)
    do k=1, nzr
    do j=2, nyr-1
    do i=1, nxr
      var(i,j,k,1) = 0.5*(var2(i,j-1,k,1)+var2(i,j,k,1))
    enddo
    enddo
    enddo
    var(:,  1,:,:) = 0.
    var(:,nyr,:,:) = 0.
    deallocate( var2 )
  else
    call geta4d(ncid ,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var)
  end if

  ! for era
  allocate( var2(nxe,nye,nze,1) )
  call get_era4d(ncid2,ivare(iv),1,nxe,1,nye,1,nze,itera,1, var2)
  if (iv == 4)  var2 = var2 / g

  allocate( vare(nxn,nyn) )

  if (idt*iut*ifc == 1) then
    allocate( var_x(nxn), var_x2(nxn) )
    allocate( var_za(1,nyn,2,1) )
  end if


  Z_UM:  DO k=1, nzr


  IF (iz(k) /= 0) then

  ! assume that the point (lon1,lat1) in ERA is same as it in UM.
  call interp(nxe,nye,var2(:,:,iz(k),1),nxn,nyn, vare)

  do j=1, nyn
    jj = (j-1)*yitv+1
    do i=1, nxn
      ii = (i-1)*xitv+1
      var_x(i) = var(ii,jj,k,1) - vare(i,j) 
    enddo
    var_x2(:) = var_x(:)*var_x(:)
    call zonal_avg(nxn,1,1,1,missv,var_x , var_za(1,j,1,1))
    call zonal_avg(nxn,1,1,1,missv,var_x2, var_za(1,j,2,1))
  enddo

  iv2 = (iv-1)*4
  set(iv2+1)%var_out(:,k,ifc,i_time) = var_za(1,:,1,1)
  set(iv2+2)%var_out(:,k,ifc,i_time) = var_za(1,:,2,1)
  set(iv2+3)%var_out(:,k,ifc,i_time) =  &
  set(iv2+1)%var_out(:,k,ifc,i_time)*cosphi(:)
  set(iv2+4)%var_out(:,k,ifc,i_time) =  &
  set(iv2+2)%var_out(:,k,ifc,i_time)*cosphi(:)


  END if  ! iz(k) /= 0


  ENDDO  Z_UM


  deallocate( var2 )
  deallocate( vare )
  deallocate( var )

  call closenc(ncid )
  call closenc(ncid2)

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


SUBROUTINE get_era4d(ncid,varname,ix,nx,iy,ny,iz,nz,it,nt, dera)

  use netio

  implicit none

  integer         , intent(in) ::  ncid, ix,nx,iy,ny,iz,nz,it,nt
  character(len=*), intent(in) ::  varname

  real, dimension(nx,ny,nz,nt), intent(out) ::  dera

  integer ::  st, varid, j,k,n
  real*8  ::  sf, ao

  integer*2, dimension(nx,ny,nz,nt) ::  dera4, derar

  
  st = nf_inq_varid(ncid,trim(varname), varid)

  st = nf_get_att_double(ncid,varid,'scale_factor', sf)
  st = nf_get_att_double(ncid,varid,'add_offset', ao)
  if (st /= 0) then
    print*, 'ERROR in get_era4d'
    STOP
  end if

  call sgeta4d(ncid,trim(varname),ix,nx,iy,ny,iz,nz,it,nt, derar)

  do n=1, nt
  do k=1, nz
  do j=1, ny
    dera4(:,j,k,n) = derar(:,ny-j+1,nz-k+1,n)
  enddo
  enddo
  enddo

  dera = real(float(dera4)*sf+ao)

  RETURN

END subroutine


SUBROUTINE interp(nxi,nyi,vari,nx,ny, var)

  implicit none

  integer,                          intent(in)  ::  nxi, nyi, nx, ny
  real, dimension(0:nxi-1,0:nyi-1), intent(in)  ::  vari
  real, dimension(0:nx-1 ,0:ny-1 ), intent(out) ::  var

  integer                        ::  i,j


if ( float(nxi)/nx /= 1.5 .or. float(nyi-1)/(ny-1) /= 1.5 ) then

  print*, 'SUBROUTINE interp is not complete.'
  STOP

else

  do j=0, ny-1, 2
  do i=0, nx-2, 2
    var(i,j) = vari(i/2*3,j/2*3)
  enddo
  enddo

  do j=0, ny-1, 2
  do i=0, nx-2, 2
    var(i+1,j) = 0.5*( vari(i/2*3+1,j/2*3) + vari(i/2*3+2,j/2*3) )
  enddo
  enddo

  do j=0, ny-3, 2
  do i=0, nx-2, 2
    var(i,j+1) = 0.5*( vari(i/2*3,j/2*3+1) + vari(i/2*3,j/2*3+2) )
  enddo
  enddo

  do j=0, ny-3, 2
  do i=0, nx-2, 2
    var(i+1,j+1) = 0.25*( vari(i/2*3+1,j/2*3+1) + vari(i/2*3+2,j/2*3+1) &
                        + vari(i/2*3+1,j/2*3+2) + vari(i/2*3+2,j/2*3+2) )
  enddo
  enddo

  RETURN

end if


END subroutine


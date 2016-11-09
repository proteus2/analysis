! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  dcase = 10010712
  real,    parameter ::  fstart = 0.0, fend = 5.0, fitv = 2./24. ! [day]
  integer, parameter ::  nvar = 4
  real,    parameter ::  missv = 1.e32
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivfile(nvar), ivarname(nvar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data11/kyh/UM_case/'/
  data                   expname /'gwdc'/
  data                   vartype /'std'/
  data                   ofdir   /'/data11/kyh/UM_case/'/ 
  data                   outname /'uvtz_ypf'/
  ! var. names in input files
  data                   ivfile(1) /'u_p'/
  data                   ivfile(2) /'v_p'/
  data                   ivfile(3) /'t_p'/
  data                   ivfile(4) /'gph_p'/
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


  real, dimension(:,:,:,:),   allocatable ::  var, var2, var_za
  real, dimension(:),         allocatable ::  lon, lat, p, lonu, latv
  real, dimension(:),         allocatable ::  time0, t_fct, t

  integer ::  nfct, fct
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, ifc
  integer ::  nxr, nyr, nzr, ntr, nd1a, nd2a, nd3a, nd4a
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid
  character(len=32)  ::  c_axis(3,2)
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


  N_VAR:   DO iv=1, nvar
  N_FCT:   DO ifc=1, nfct


  fct = int((fstart + (ifc-1)*fitv)*24.)

  write(6,'(a2,i8.8,a3,i4.4)') '20', dcase, ' + ', fct

  write(fname,'(a,i8.8,a,i8.8,a,i3.3,a)')                  &
        trim(ifdir)//trim(expname)//'/20',dcase,'/'//      &
        trim(vartype)//'.'//trim(ivfile(iv))//'.20',dcase, &
        '+',24*(int(t_fct(ifc)-0.001)+1),'.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (ifc == 1) then
    call diminfop(ncid,.TRUE., nx,ny,nz,nxu,nyv,c_axis)
    if (nx == 1)  nx = nxu
    if (ny == 1)  ny = nyv + 1
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
  nxr = NX   ;   nd1a = 1 
  nyr = NY   ;   nd2a = NY
  nzr = NZ   ;   nd3a = NZ    ! ;   iz = 28
  ntr = 1    ;   nd4a = NFCT

  if (ifc == 1)  allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )

  allocate( var(nxr,nyr,nzr,1) )
  if (iv == 2) then
    allocate( var2(nxr,nyv,nzr,1) )
    call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyv,1,nzr,it,1, var2)
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
    call geta4d(ncid,trim(ivarname(iv)),1,nxr,1,nyr,1,nzr,it,1, var)
  end if

  call zonal_avg(nxr,nyr,nzr,1,missv,var, set(iv)%var_out(1,:,:,ifc))
  deallocate ( var )

  if (ifc == 1) then
    set(iv)%nd(1) = nd1a
    set(iv)%nd(2) = nd2a
    set(iv)%nd(3) = nd3a
    set(iv)%nd(4) = nd4a
    set(iv)%axis = (/'    ','lat  ','p   ','fcst'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = 0.
    set(iv)%axis2 = lat
    set(iv)%axis3 = p
    set(iv)%axis4 = t_fct
  end if

  call closenc(ncid)


  ENDDO  N_FCT


  if ( allocated(lon ) )  deallocate( lon  )
  if ( allocated(lat ) )  deallocate( lat  )
  if ( allocated(p   ) )  deallocate( p    )
  if ( allocated(lonu) )  deallocate( lonu )
  if ( allocated(latv) )  deallocate( latv )


  ENDDO  N_VAR


  ! dump
  write(fname,'(a,i8.8,a)') trim(ofdir)//trim(expname)//'/anal/'// &
        trim(expname)//'.'//trim(outname)//'.20',dcase,'.nc'

  write(6,*)
  write(6,*) trim(fname)
  write(6,*)

  call outnc(trim(fname),nvar,set,'')

  do iv=1, nvar
    deallocate( set(iv)%var_out )
  enddo


  STOP

END program


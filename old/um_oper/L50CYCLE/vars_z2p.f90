PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

!  integer, parameter  ::  np = 30
!  real, dimension(np) ::  p = (/1000.,925.,850.,700.,500.,400.,300.,250.,200.,150.,100., &
!                                70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,0.7,0.5,0.3,0.2,0.1,  &
!                                0.07,0.05,0.04,0.03/)
  integer, parameter  ::  np = 50
  real, dimension(np) ::  p = (/1000.,925.,850.,750.,600.,500.,420.,350., &
                                300.,250.,200.,160.,125.,100., &
                                80.,65.,55.,45.,37.,30.,25.,20.,16.,12.5,10., &
                                8.0,6.5,5.5,4.5,3.7,3.0,2.5,2.0,1.6,1.25,1.0, &
                                0.8,0.65,0.55,0.45,0.37,0.3,0.25,0.2,0.16,0.125,0.1, &
                                0.08,0.065,0.055/)
  real,    parameter  ::  h_scale = 7.0e3
  integer, parameter  ::  nmonth =  1, mstart =  7
  integer, parameter  ::  nyear  =  1, ystart = 2008
  integer, parameter  ::  nvar = 5
  real,    parameter  ::  missv = 1.e32
  real,    parameter  ::  zh_int = 12.0
  character(len=128)  ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)   ::  ivarname(nvar), ovarname(0:nvar)
  character(len=64)   ::  c_prho, c_pth, c_psfc, c_zrho, c_zth, c_zsfc

  ! files
  data                   ifdir   /'/data5/kyh/umres/'/
  data                   expname /'test'/
  data                   vartype /'std'/
  data                   ofdir   /'/data5/kyh/umres/'/ 
  data                   outname /'var_xypt'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'dz_dt'/
  data                   ivarname(4) /'temp'/
  data                   ivarname(5) /'field83'/  ! PV on theta lev.
  data                   c_prho/'p_2'/
  data                   c_pth /'p_3'/
  data                   c_psfc/'p_4'/
  data                   c_zrho/'ht_1'/
  data                   c_zth /'ht'/
  data                   c_zsfc/'ht_2'/
  ! output variables
  data                   ovarname(0) /'z'/
  data                   ovarname(1) /'u'/
  data                   ovarname(2) /'v'/
  data                   ovarname(3) /'wp'/
  data                   ovarname(4) /'temp'/
  data                   ovarname(5) /'PV'/


  real, dimension(:,:,:,:),   allocatable ::  var, var2
  real, dimension(:,:,:,:),   allocatable ::  p_rho, p_th, p_all, z_all
  real, dimension(:,:,:),     allocatable ::  z_rho, z_th, u1, v1
  real, dimension(:,:),       allocatable ::  dhdx, dhdy, h_side

  real, dimension(:),         allocatable ::  lon, lat, zh, lonu, latv, zht, t

  integer ::  year, month
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, imn, iyr
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz, ind0, ind1, ind2
  integer ::  i,j,k,n,iv, ncid, ncid2
  real    ::  temp
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
  type(vset), dimension(0:nvar) ::  set


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
  allocate( t(nt) )
  call get1d(ncid,'t',nt, t)

! ==============================================
  nxr = NX   ;   nxa = NX
  nyr = NY   ;   nya = NY
  nzr = NZ   ;   nza = NP  ! ;   iz = 28
  ntr = NT   ;   nta = NT
!===============================================
  do iv=0, nvar
    allocate( set(iv)%var_out(nxa,nya,nza,nta) )
  enddo

  ! get pressure data
  allocate( p_rho(nxr,nyr,0:nzr,ntr) )
  allocate( p_th (nxr,nyr,0:nzr,ntr) )
  allocate( p_all(nxr,nyr,0:nzr*2,ntr) )

  call get4d(ncid,trim(c_prho),nxr,nyr,nzr,ntr, p_rho(:,:,1:nzr,:))
  call get4d(ncid,trim(c_pth ),nxr,nyr,nzr,ntr, p_th (:,:,1:nzr,:))
  call get4d(ncid,trim(c_psfc),nxr,nyr,  1,ntr, p_all(:,:,0    ,:))
  p_rho(:,:,0,:) = p_all(:,:,0,:)
  p_th (:,:,0,:) = p_all(:,:,0,:)
  do n=1, ntr
  do k=1, nzr
    p_all(:,:,k*2-1,n) = p_rho(:,:,k,n)
    p_all(:,:,k*2  ,n) = p_th (:,:,k,n)
  enddo
  enddo
  p_rho(:,:,:,:) = p_rho(:,:,:,:)/100.
  p_th (:,:,:,:) = p_th (:,:,:,:)/100.
  p_all(:,:,:,:) = p_all(:,:,:,:)/100.

  ! get (geometric) height data (from initial field)
  allocate( z_all(nxr,nyr,0:nzr*2,ntr) )
  allocate( z_rho(nxr,nyr,nzr) )
  allocate( z_th (nxr,nyr,nzr) )

  fname = trim(ifdir)//'ctl/ctl.grid.nc'
  call opennc(trim(fname),ncid2)
  call get4d(ncid2,trim(c_zrho),nxr,nyr,nzr,1, z_rho)
  call get4d(ncid2,trim(c_zth ),nxr,nyr,nzr,1, z_th )
  call get4d(ncid2,trim(c_zsfc),nxr,nyr,  1,1, z_all(:,:,0,1))
  call closenc(ncid2)
  do k=1, nzr
    z_all(:,:,k*2-1,1) = z_rho(:,:,k)
    z_all(:,:,k*2  ,1) = z_th (:,:,k)
  enddo
  deallocate( z_rho, z_th )
  do n=2, ntr
    z_all(:,:,:,n) = z_all(:,:,:,1)
  enddo

  ! save terrain-height gradient to calculate w_sfc
  allocate( dhdx(nxr,nyr), dhdy(nxr,nyr), h_side(nxr,2) )
  allocate( u1(nxr,nyr,ntr), v1(nxr,nyr,ntr) )
  call gradx_2nd(nxr,nyr,1,1,z_all(:,:,0,1),lon,lat,missv, dhdx)
  call grady_2nd(nxr,nyr,1,1,z_all(:,:,0,1),lat,missv, dhdy)
  h_side(:,1) = z_all(:,2    ,0,1)
  h_side(:,2) = z_all(:,nyr-1,0,1)

  ! vertical interpolation for height
  call int_z2p(nxr,nyr,nzr*2,ntr,0,z_all,p_all,nza,p,missv, set(0)%var_out)
  deallocate( z_all )


  allocate( var(nxr,nyr,0:nzr,ntr) )
 
  do iv=1, nvar

    ! horizontal interpolation for u, v
    if (iv == 1) then
      allocate( var2(NXU,nyr,nzr,ntr) )
      call get4d(ncid,trim(ivarname(iv)),NXU,nyr,nzr,ntr, var2)
      do n=1, ntr
      do k=1, nzr
      do j=1, nyr
        var(1,j,k,n) = 0.5*(var2(1,j,k,n)+var2(nxu,j,k,n))
        do i=2, nxr
          var(i,j,k,n) = 0.5*(var2(i-1,j,k,n)+var2(i,j,k,n))
        enddo
      enddo
      enddo
      enddo
      deallocate( var2 )
      u1(:,:,:) = var(:,:,1,:)
    else if (iv == 2) then
      allocate( var2(nxr,NYV,nzr,ntr) )
      call get4d(ncid,trim(ivarname(iv)),nxr,NYV,nzr,ntr, var2)
      do n=1, ntr
      do k=1, nzr
        do j=2, nyr-1
        do i=1, nxr
          var(i,j,k,n) = 0.5*(var2(i,j-1,k,n)+var2(i,j,k,n))
        enddo
        enddo
        var(:,  1,k,n) = 0.
        var(:,nyr,k,n) = 0.
      enddo
      enddo
      deallocate( var2 )
      v1(:,:,:) = var(:,:,1,:)
    else if (iv == 3) then
      call get4d(ncid,trim(ivarname(iv)),nxr,nyr,nzr,ntr, var(:,:,1:nzr,:))
      ! calculate w_sfc
      do n=1, ntr
        do j=2, nyr-1
          var(:,j,0,n) = u1(:,j,n)*dhdx(:,j) + v1(:,j,n)*dhdy(:,j)
        enddo
        ind1 = maxloc(u1(:,1,n),1)
        ind0 = ind1 - nx/4
        ind2 = ind1 + nx/4
        if (ind0 < 1)   ind0 = ind0 + nx
        if (ind2 > nx)  ind2 = ind2 - nx
        var(:,1,0,n) = u1(ind1,1,n)*(h_side(ind2,1)-h_side(ind0,1)) &
                       /((lat(3)-lat(1))*deg2rad*r_earth)
        ind1 = maxloc(u1(:,nyr,n),1)
        ind0 = ind1 - nx/4
        ind2 = ind1 + nx/4
        if (ind0 < 1)   ind0 = ind0 + nx
        if (ind2 > nx)  ind2 = ind2 - nx
        var(:,nyr,0,n) = u1(ind1,nyr,n)*(h_side(ind2,2)-h_side(ind0,2)) &
                         /((lat(3)-lat(1))*deg2rad*r_earth)
      enddo
      deallocate( u1, v1, dhdx, dhdy, h_side )
    else if (iv == 4) then
      call get4d(ncid,trim(ivarname(iv)),nxr,nyr,nzr,ntr, var(:,:,1:nzr,:))
      call get4d(ncid,trim(ivarname(iv))//'_1',nxr,nyr,1,ntr, var(:,:,0,:))
    else
      call get4d(ncid,trim(ivarname(iv)),nxr,nyr,nzr,ntr, var(:,:,1:nzr,:))
    end if

    ! calculate d(zp)/dt (vertical velocity in pressure-height coord.)
    if (iv == 3) then
      allocate( var2(nxr,nyr,0:nzr,ntr) )
      call get4d(ncid,'temp',nxr,nyr,nzr,ntr, var2(:,:,1:nzr,:))
      call get4d(ncid,'temp_1',nxr,nyr,1,ntr, var2(:,:,0,:))
      temp = g*h_scale/rd
      var(:,:,:,:) = var(:,:,:,:)/var2(:,:,:,:)*temp
      deallocate( var2 )
    end if

    ! vertical interpolation into p-coordinate
    if ( iv == 1 .or. iv == 2 ) then
      call int_z2p(nxr,nyr,nzr,ntr,1,var(:,:,1:nzr,:),p_rho,nza,p,missv, &
                   set(iv)%var_out)
    else if ( iv == 3 .or. iv == 4 ) then
      call int_z2p(nxr,nyr,nzr,ntr,0,var(:,:,0:nzr,:),p_th,nza,p,missv, &
                   set(iv)%var_out)
    else
      call int_z2p(nxr,nyr,nzr,ntr,1,var(:,:,1:nzr,:),p_th,nza,p,missv, &
                   set(iv)%var_out)
    end if

  enddo

  call closenc(ncid)

  deallocate( var )
  deallocate( p_rho, p_th, p_all )


  do iv=0, nvar
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%axis = (/'lon','lat','p','t'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = lon
    set(iv)%axis2 = lat
    set(iv)%axis3 = p
    set(iv)%axis4 = t
  enddo

  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
  call outnc(trim(fname),nvar+1,set,'')

  do iv=0, nvar
    deallocate( set(iv)%var_out )
  enddo
  deallocate( t )


  ENDDO  N_YEAR
  ENDDO  N_MON


  STOP

END program


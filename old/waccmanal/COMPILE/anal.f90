! for the operational output

PROGRAM analysis_UM

  ! Hydrostatic TEM eq. in p-coord.
  ! (Andrews et al., 1987)

  use netio
  use um_anal
  use um_axis
  use fft

  implicit none

  integer, parameter ::  year1 = 1980, nyear = 10
  integer, parameter ::  month1 =  1, nmon = 12
  integer, parameter ::  wn1 = 1, wn2 = 10                      ! wn2 < nx/2
  real,    parameter ::  h_scale = 7.e3
  real,    parameter ::  missv = 1.e32, zbottom = 0. !5.e3
  integer, parameter ::  nivar = 4, nvar = 3
  character(len=128) ::  ifdir, ifname, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! files
!  data                   ifdir   /'/export6/sis/control1_result/daily'/
!  data                   ifname  /'control1.waccm.day'/
  data                   ifdir   /'/data10/BACK/sis/ex7-8/lsgwdc_result/clim/daily1'/
  data                   ifname  /'lsgwdclm.waccm.day1'/
  data                   ofdir   /'./'/
  data                   outname /'epf_yp'/
  ! input variables
  data                   ivarname(1) /'U'/
  data                   ivarname(2) /'V'/
  data                   ivarname(3) /'OMEGA'/
  data                   ivarname(4) /'T'/
  ! output variables
  data                   ovarname(1) /'epd'/
  data                   ovarname(2) /'f_y'/
  data                   ovarname(3) /'f_z'/


  real, dimension(:,:,:),   allocatable ::  u, v, w, pt, var_in
  real, dimension(:,:),     allocatable ::  um, vm, wm, ptm, rho0
  real, dimension(:,:),     allocatable ::  epd, fy, fz
  real, dimension(:),       allocatable ::  lon, lat, zp, p, p_in
  real, dimension(:),       allocatable ::  mon, wn

  real, dimension(:,:),     allocatable ::  phi, vpt, wpt, vu, wu
  real, dimension(:,:),     allocatable ::  grady, gradz, temp
  real, dimension(:),       allocatable ::  cosphi, f

  double complex, dimension(:,:,:), allocatable ::  coefu, coefv, coefw, coefpt

  integer, dimension(:,:), allocatable ::  nepd, nfy, nfz

  integer ::  nx, ny, nz, iyr
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta, nwn, nd
  integer ::  iz, it
  integer ::  i,j,k,n,iv,iwn, ncid, tag, iwn1, iwn2, imon, iday
  logical ::  l_err
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=128) ::  fname

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  allocate( mon(nmon) )

  nwn = wn2 - wn1 + 1
  allocate( wn(nwn+1) )
  do iwn=1, nwn
    wn(iwn) = iwn*1.
  enddo
  wn(nwn+1) = 999.


  it = 0


  N_MON:   DO imon=month1, month1+nmon-1


  it = it + 1
  mon(it) = imon


  N_YR:    DO iyr=year1, year1+nyear-1


  write(cyear, '(i4.4)') iyr
  write(cmonth,'(i2.2)') imon

  fname = trim(ifdir)//'/'//trim(ifname)//'.'//cyear//'-'//cmonth//'.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if ( imon == month1 .and. iyr == year1 ) then
    call dilen(ncid,'lon' , nx,l_err)
    call dilen(ncid,'lat' , ny,l_err)
    call dilen(ncid,'lev' , nz,l_err)

    allocate( lon(nx), lat(ny), p(nz), p_in(nz) )

    call get1d(ncid,'lon',nx, lon )
    call get1d(ncid,'lat',ny, lat )
    call get1d(ncid,'lev',nz, p_in)
    do k=1, nz
      p(k) = p_in(nz+1-k)*100.
    enddo

    allocate( zp(nz) )
    zp(:) = -h_scale*(log(p(:)/1.e5))
  end if

  call dilen(ncid,'time', nd,l_err)

  ! allocate var.s
  if ( imon == month1 .and. iyr == year1 ) then
    allocate( cosphi(ny), f(ny) )
    f(:) = 2.*ome_earth*sin(lat(:)*deg2rad)
    cosphi(:) = cos(lat(:)*deg2rad)
  end if

  nxr = NX      ;   nxa = NY
  nyr = NY      ;   nya = NZ
  nzr = NZ      ;   nza = NWN+1 ! ;   iz = 28
  ntr = ND      ;   nta = NMON

  if ( imon == month1 .and. iyr == year1 ) then
    do iv=1, nvar
      allocate( set(iv)%var_out(nxa,nya,nza,nta) )
      set(iv)%var_out(:,:,:,:) = 0.
    enddo
  end if


  N_DAY:   DO iday=1, nd

 
  allocate( um(ny,nz), vm(ny,nz), wm(ny,nz) )
  allocate( ptm(ny,nz) )
  allocate( rho0(ny,nz) )

  allocate( epd(ny,nz) )
  allocate( fy(ny,nz), fz(ny,nz) )

  allocate( phi(ny,nz) )
  allocate( vpt(ny,nz), wpt(ny,nz) )
  allocate( vu(ny,nz), wu(ny,nz) )

  allocate( grady(ny,nz), gradz(ny,nz) )
  allocate( temp(ny,nz) )

  allocate( u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz) )
  allocate( pt(nx,ny,nz) )
  allocate( var_in(nx,ny,nz) )

  ! get var.s
  call geta4d(ncid,trim(ivarname(1)),1,nx,1,ny,1,nz,iday,1, var_in)
  do k=1, nz
    u(:,:,k) = var_in(:,:,nz+1-k)
  enddo
  call geta4d(ncid,trim(ivarname(2)),1,nx,1,ny,1,nz,iday,1, var_in)
  do k=1, nz
    v(:,:,k) = var_in(:,:,nz+1-k)
  enddo
  call geta4d(ncid,trim(ivarname(3)),1,nx,1,ny,1,nz,iday,1, var_in)
  do k=1, nz
    w(:,:,k) = -var_in(:,:,nz+1-k)*h_scale/p(k)
  enddo
  call geta4d(ncid,trim(ivarname(4)),1,nx,1,ny,1,nz,iday,1, var_in)
  do k=1, nz
    pt(:,:,k) = var_in(:,:,nz+1-k)*(1.e5/p(k))**kappa
  enddo

  deallocate( var_in )

  ! mean and perturbation
  call zonal_avg(nx,ny,nz,1,missv,u ,um )
  call zonal_avg(nx,ny,nz,1,missv,v ,vm )
  call zonal_avg(nx,ny,nz,1,missv,w ,wm )
  call zonal_avg(nx,ny,nz,1,missv,pt,ptm)

  do k=1, nz
    rho0(:,k) = p(k)/rd/(ptm(:,k)*(p(k)/1.e5)**kappa)
  enddo

  do k=1, nz
  do j=1, ny
    u (:,j,k) = u (:,j,k) - um (j,k)
    v (:,j,k) = v (:,j,k) - vm (j,k)
    w (:,j,k) = w (:,j,k) - wm (j,k)
    pt(:,j,k) = pt(:,j,k) - ptm(j,k)
  enddo
  enddo

  deallocate( vm, wm )


  ! mean flux of waves
  iwn1 = wn1 + 1
  iwn2 = wn2 + 1

  allocate( coefu(nx,ny,nz), coefv(nx,ny,nz), coefw(nx,ny,nz), coefpt(nx,ny,nz) )

  do k=1, nz
  do j=1, ny
    call fft1df(nx,u (:,j,k), coefu (:,j,k))
    call fft1df(nx,v (:,j,k), coefv (:,j,k))
    call fft1df(nx,w (:,j,k), coefw (:,j,k))
    call fft1df(nx,pt(:,j,k), coefpt(:,j,k))
  enddo
  enddo

  deallocate( pt )
  deallocate( u, v, w )


  N_WN:  DO iwn=iwn1, iwn2+1


  if (iwn /= iwn2+1) then
    do k=1, nz
    do j=1, ny
      vpt(j,k) = 2.*(real(coefv(iwn,j,k))*real(coefpt(iwn,j,k)) + &
                 aimag(coefv(iwn,j,k))*aimag(coefpt(iwn,j,k)))/(nx*nx)
      wpt(j,k) = 2.*(real(coefw(iwn,j,k))*real(coefpt(iwn,j,k)) + &
                 aimag(coefw(iwn,j,k))*aimag(coefpt(iwn,j,k)))/(nx*nx)
      vu(j,k) = 2.*(real(coefv(iwn,j,k))*real(coefu(iwn,j,k)) + &
                aimag(coefv(iwn,j,k))*aimag(coefu(iwn,j,k)))/(nx*nx)
      wu(j,k) = 2.*(real(coefw(iwn,j,k))*real(coefu(iwn,j,k)) + &
                aimag(coefw(iwn,j,k))*aimag(coefu(iwn,j,k)))/(nx*nx)
    enddo
    enddo
  else
    do k=1, nz
    do j=1, ny
      vpt(j,k) = 2.*sum(real(coefv(iwn1:iwn2,j,k))*real(coefpt(iwn1:iwn2,j,k)) + &
                 aimag(coefv(iwn1:iwn2,j,k))*aimag(coefpt(iwn1:iwn2,j,k)))/(nx*nx)
      wpt(j,k) = 2.*sum(real(coefw(iwn1:iwn2,j,k))*real(coefpt(iwn1:iwn2,j,k)) + &
                 aimag(coefw(iwn1:iwn2,j,k))*aimag(coefpt(iwn1:iwn2,j,k)))/(nx*nx)
      vu(j,k) = 2.*sum(real(coefv(iwn1:iwn2,j,k))*real(coefu(iwn1:iwn2,j,k)) + &
                aimag(coefv(iwn1:iwn2,j,k))*aimag(coefu(iwn1:iwn2,j,k)))/(nx*nx)
      wu(j,k) = 2.*sum(real(coefw(iwn1:iwn2,j,k))*real(coefu(iwn1:iwn2,j,k)) + &
                aimag(coefw(iwn1:iwn2,j,k))*aimag(coefu(iwn1:iwn2,j,k)))/(nx*nx)
    enddo
    enddo
  end if


  ! phi
!  call grady_2nd(1,ny,nz,1,ptm,lat,missv, grady)
  grady = 0.   ! hydrostatic TEM in p-coord. (Andrews et al.)

  temp(:,:) = log(ptm(:,:))
  call gradz_2nd(1,ny,nz,1,temp,zp,0., gradz)
  gradz(:,:) = ptm(:,:)*gradz(:,:)

  phi(:,:) = (vpt(:,:)*gradz(:,:) - wpt(:,:)*grady(:,:))  &
            / (grady(:,:)**2 + gradz(:,:)**2)

  do k=1, nz
    temp(:,k) = um(:,k)*cosphi(:)
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,0., grady)
  do k=1, nz
  do j=2, ny-1
    grady(j,k) = grady(j,k)/cosphi(j)
  enddo
  enddo
  grady( 1,:) = missv
  grady(ny,:) = missv

  call gradz_2nd(1,ny,nz,1,um,zp,0., gradz)

  ! epf, epd
  fy(:,:) = missv
  fz(:,:) = missv

  do k=1, nz
    fy(:,k) = (vu(:,k) - phi(:,k)*gradz(:,k))  &
              *(-rho0(:,k)*r_earth*cosphi(:))
  enddo
  do k=1, nz
  do j=2, ny-1
    fz(j,k) = (wu(j,k) + phi(j,k)*(grady(j,k)-f(j)))  &
              *(-rho0(j,k)*r_earth*cosphi(j))
  enddo
  enddo
  fz( 1,:) = missv
  fz(ny,:) = missv

  call div_yz(ny,nz,1,fy,fz,lat,zp,missv, epd)
  do k=1, nz
  do j=2, ny-1
    epd(j,k) = epd(j,k)/(rho0(j,k)*r_earth*cosphi(j)) * 86400.
  enddo
  enddo
  epd( 1,:) = missv
  epd(ny,:) = missv


  do k=1, nz
  do j=2, ny-1
    set(1)%var_out(j,k,iwn-iwn1+1,it) = set(1)%var_out(j,k,iwn-iwn1+1,it) + epd(j,nz+1-k)/nd
    set(2)%var_out(j,k,iwn-iwn1+1,it) = set(2)%var_out(j,k,iwn-iwn1+1,it) + fy (j,nz+1-k)/nd
    set(3)%var_out(j,k,iwn-iwn1+1,it) = set(3)%var_out(j,k,iwn-iwn1+1,it) + fz (j,nz+1-k)/nd
  enddo
  enddo


  ENDDO  N_WN


  deallocate( coefu, coefv, coefw, coefpt )
  deallocate( grady, gradz )
  deallocate( epd )
  deallocate( fy, fz )
  deallocate( vpt, wpt )
  deallocate( vu, wu )
  deallocate( phi )
  deallocate( ptm, um )
  deallocate( temp )
  deallocate( rho0 )


  ENDDO  N_DAY


  call closenc(ncid)


  ENDDO  N_YR


  do iv=1, nvar
    do k=1, nz
    do j=2, ny-1
      set(iv)%var_out(j,k,:,it) = set(iv)%var_out(j,k,:,it) / nyear
    enddo
    enddo
    set(iv)%var_out( 1,:,:,it) = missv
    set(iv)%var_out(ny,:,:,it) = missv
  enddo


  ENDDO  N_MON


  do iv=1, nvar
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%axis = (/'lat  ','lev  ','wn   ','month'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = lat
    set(iv)%axis2 = p_in
    set(iv)%axis3 = wn
    set(iv)%axis4 = mon 
  enddo

  ! dump
  fname = trim(ofdir)//trim(outname)//'.nc'
  call outnc(trim(fname),nvar,set,'')

  do iv=1, nvar
    deallocate( set(iv)%var_out )
  enddo


  STOP

END program


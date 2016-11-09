! for the operational output

PROGRAM analysis_UM

  ! Hydrostatic TEM eq. in p-coord.
  ! (Andrews et al., 1987)

  use netio
  use um_anal
  use um_axis
  use fft

  implicit none

  integer, parameter ::  year = 2009, month =  7
  integer, parameter ::  dstart = 6, dend = 23   ! 6, 31
  integer, parameter ::  nutc  =  1, ustart = 0
  real,    parameter ::  fstart = 5.0, fend = 5.0, fitv = 1.0   ! [day]
  real,    parameter ::  avg_itv = 0.5
  integer, parameter ::  wn1 = 4, wn2 = 8                       ! wn2 < nx/2
  real,    parameter ::  h_scale = 7.e3
  real,    parameter ::  missv = 1.e32, zbottom = 0. !5.e3
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  integer, parameter ::  nivar = 4, nvar = 3
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070200'/
  ! files
  data                   ifdir   /'/data6/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'gwdc'/
  data                   vartype /'press1'/
  data                   ofdir   /'/data6/kyh/UM_OPER/UM_OUT/'/
  data                   outname /'epf_yp'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'dz_dt'/
  data                   ivarname(4) /'temp'/
  ! output variables
  data                   ovarname(1) /'epd'/
  data                   ovarname(2) /'f_y'/
  data                   ovarname(3) /'f_z'/


  real, dimension(:,:,:,:), allocatable ::  tem_in, tem_in2
  real, dimension(:,:,:),   allocatable ::  u, v, w, pt
  real, dimension(:,:),     allocatable ::  um, vm, wm, ptm, rho0
  real, dimension(:,:),     allocatable ::  epd, fy, fz
  real, dimension(:),       allocatable ::  lon, lat, zp, p, lonu, latv, zpt
  real, dimension(:),       allocatable ::  time0

  real, dimension(:,:),     allocatable ::  phi, vpt, wpt, vu, wu
  real, dimension(:,:),     allocatable ::  grady, gradz, temp
  real, dimension(:),       allocatable ::  cosphi, f

  double complex, dimension(:), allocatable ::  coefu, coefv, coefw, coefpt

  integer, dimension(:,:), allocatable ::  nepd, nfy, nfz
  integer, dimension(:),   allocatable ::  it

  integer ::  nfct, i_time, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta, ntavg
  integer ::  iz
  integer ::  i,j,k,n,iv, ncid, tag, na, iwn1, iwn2
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

  tag = 0
  do n=1, nexcept
    if (cutc0 == excdate(n))  tag = 1
  enddo
  if ( tag ) then
    write(6,*) ' Exception - passed.'
    i_time = i_time - 1
    CYCLE
  end if

  fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          cutc0//'.'//trim(vartype)//'.nc'

  if ( trim(vartype) == 'press1' .and. fct > 24 )  &
     fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
             cutc0//'.press2.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (idt*iut*ifc == 1) then
    call diminfop(ncid,.TRUE., nx,ny,nz,nxu,nyv,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( p  (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)
    p(:) = p(:)*100.
    allocate( zp(nz) )
    zp(:) = -h_scale*(log(p(:)/1.e5))
  end if

  ntavg = nint(fitv/avg_itv)
  allocate( it(ntavg) )

  call timeinfo(ncid, nt0)
  allocate( time0(nt0) )
  call get1d(ncid,'t',nt0, time0)
  do na=1, ntavg
  do n=1, nt0
    if ( abs(time0(n) - (real(fct)/24.-avg_itv*(ntavg-na))) < v_small ) then
      it(na) = n
    end if
  enddo
  enddo
  deallocate( time0 )


  ! allocate var.s
  if (idt*iut*ifc == 1) then
    allocate( cosphi(ny), f(ny) )
    f(:) = 2.*ome_earth*sin(lat(:)*deg2rad)
    cosphi(:) = cos(lat(:)*deg2rad)
  end if

  nxr = NX      ;   nxa =  1
  nyr = NY      ;   nya = NY
  nzr = NZ      ;   nza = NZ  ! ;   iz = 28
  ntr = NT      ;   nta =  1
  if (idt*iut*ifc == 1) then
    do iv=1, nvar
      allocate( set(iv)%var_out(nxa,nya,nza,nta) )
      set(iv)%var_out(:,:,:,:) = 0.
    enddo
    allocate( nepd(ny,nz) )
    allocate( nfy(ny,nz), nfz(ny,nz) )
    nepd(:,:) = 0
    nfy (:,:) = 0
    nfz (:,:) = 0
  end if

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

  ! get var.s
  allocate( tem_in(nx,ny,nz,ntavg) )
  allocate( tem_in2(nx,ny,nz,ntavg) )

  do na=1, ntavg
    call geta4d(ncid,trim(ivarname(1)),1,nx,1,ny,1,nz,it(na),1, tem_in(:,:,:,na))
  enddo
  call tempo_avg(nx,ny,nz,ntavg,missv,tem_in, u)

  do na=1, ntavg
    call geta4d(ncid,trim(ivarname(2)),1,nx,1,nyv,1,nz,it(na),1, tem_in(:,2:ny,:,na))
  enddo
  call tempo_avg(nx,nyv,nz,ntavg,missv,tem_in(:,2:ny,:,:), v(:,2:ny,:))
  do k=1, nzr
  do j=2, nyv
    v(:,j,k) = 0.5*(v(:,j,k)+v(:,j+1,k))
  enddo
  enddo
  v(:,  1,:) = 0.
  v(:,nyr,:) = 0.
  do na=1, ntavg
    call geta4d(ncid,trim(ivarname(3)),1,nx,1,ny,1,nz,it(na),1, tem_in(:,:,:,na))
  enddo
  do na=1, ntavg
    call geta4d(ncid,trim(ivarname(4)),1,nx,1,ny,1,nz,it(na),1, tem_in2(:,:,:,na))
  enddo
  do na=1, ntavg
  do k=1, nzr
    tem_in(:,:,k,na) = tem_in(:,:,k,na)*g*h_scale/rd/(tem_in2(:,:,k,na)*(p(k)/1.e5)**kappa)
  enddo
  enddo
  call tempo_avg(nx,ny,nz,ntavg,missv,tem_in, w)
  call tempo_avg(nx,ny,nz,ntavg,missv,tem_in2, pt)

  deallocate( tem_in2 )
  deallocate( tem_in )
  deallocate( it )

  ! mean and perturbation
  call zonal_avg(nx,ny,nz,1,missv,u ,um )
  call zonal_avg(nx,ny,nz,1,missv,v ,vm )
  call zonal_avg(nx,ny,nz,1,missv,w ,wm )
  call zonal_avg(nx,ny,nz,1,missv,pt,ptm)


  rho0(:,:) = missv
  do k=1, nz
  do j=1, ny
    if (ptm(j,k) /= missv)  rho0(j,k) = p(k)/rd/ptm(j,k)
  enddo
  enddo
  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (pt(i,j,k) /= missv)  pt(i,j,k) = pt(i,j,k)*(1.e5/p(k))**kappa
  enddo
  enddo
  enddo
  do k=1, nz
  do j=1, ny
    if (ptm(j,k) /= missv)  ptm(j,k) = ptm(j,k)*(1.e5/p(k))**kappa
  enddo
  enddo

  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (u (i,j,k) /= missv)  u (i,j,k) = u (i,j,k) - um (j,k)
    if (v (i,j,k) /= missv)  v (i,j,k) = v (i,j,k) - vm (j,k)
    if (w (i,j,k) /= missv)  w (i,j,k) = w (i,j,k) - wm (j,k)
    if (pt(i,j,k) /= missv)  pt(i,j,k) = pt(i,j,k) - ptm(j,k)
  enddo
  enddo
  enddo

  deallocate( vm, wm )


  ! mean flux of waves
  iwn1 = wn1 + 1
  iwn2 = wn2 + 1

  allocate( coefu(nx), coefv(nx), coefw(nx), coefpt(nx) )

  vpt(:,:) = missv
  wpt(:,:) = missv
  vu (:,:) = missv
  wu (:,:) = missv

  do k=1, nz
  do j=1, ny

    coefu (:) = missv
    coefv (:) = missv
    coefw (:) = missv
    coefpt(:) = missv

    tag = 1
    do i=1, nx
      if (u(i,j,k) == missv)  tag = 0
    enddo
    if ( tag )  call fft1df(nx,u(:,j,k), coefu)
    tag = 1
    do i=1, nx
      if (v(i,j,k) == missv)  tag = 0
    enddo
    if ( tag )  call fft1df(nx,v(:,j,k), coefv)
    tag = 1
    do i=1, nx
      if (w(i,j,k) == missv)  tag = 0
    enddo
    if ( tag )  call fft1df(nx,w(:,j,k), coefw)
    tag = 1
    do i=1, nx
      if (pt(i,j,k) == missv)  tag = 0
    enddo
    if ( tag )  call fft1df(nx,pt(:,j,k), coefpt)

    if (coefpt(1) /= missv) then
      if (coefv(1) /= missv)  &
         vpt(j,k) = 2.*sum( real(coefv(iwn1:iwn2))*real(coefpt(iwn1:iwn2)) + &
                    aimag(coefv(iwn1:iwn2))*aimag(coefpt(iwn1:iwn2)) )/(nx*nx)
      if (coefw(1) /= missv)  &
         wpt(j,k) = 2.*sum( real(coefw(iwn1:iwn2))*real(coefpt(iwn1:iwn2)) + &
                    aimag(coefw(iwn1:iwn2))*aimag(coefpt(iwn1:iwn2)) )/(nx*nx)
    end if

    if (coefu(1) /= missv) then
      if (coefv(1) /= missv)  &
         vu(j,k) = 2.*sum( real(coefv(iwn1:iwn2))*real(coefu(iwn1:iwn2)) + &
                   aimag(coefv(iwn1:iwn2))*aimag(coefu(iwn1:iwn2)) )/(nx*nx)
      if (coefw(1) /= missv)  &
         wu(j,k) = 2.*sum( real(coefw(iwn1:iwn2))*real(coefu(iwn1:iwn2)) + &
                   aimag(coefw(iwn1:iwn2))*aimag(coefu(iwn1:iwn2)) )/(nx*nx)
    end if

  enddo
  enddo

  deallocate( coefu, coefv, coefw, coefpt )
  deallocate( pt )
  deallocate( u, v, w )

  ! phi
!  call grady_2nd(1,ny,nz,1,ptm,lat,missv, grady)
  grady = 0.   ! hydrostatic TEM in p-coord. (Andrews et al.)

  temp(:,:) = missv
  do k=1, nz
  do j=1, ny
    if (ptm(j,k) /= missv)  temp(j,k) = log(ptm(j,k))
  enddo
  enddo
  call gradz_2nd(1,ny,nz,1,temp,zp,missv, gradz)
  do k=1, nz
  do j=1, ny
    if (gradz(j,k) /= missv)  gradz(j,k) = ptm(j,k)*gradz(j,k)
  enddo
  enddo

  deallocate( ptm )

  phi(:,:) = missv
  do k=1, nz
  do j=1, ny
    if ( grady(j,k) /= missv .and. gradz(j,k) /= missv .and.   &
         vpt(j,k) /= missv .and. wpt(j,k) /= missv )           &
       phi(j,k) = (vpt(j,k)*gradz(j,k) - wpt(j,k)*grady(j,k))  &
                 / (grady(j,k)**2 + gradz(j,k)**2)
  enddo
  enddo

  temp(:,:) = missv
  do k=1, nz
  do j=1, ny
    if (um(j,k) /= missv)  temp(j,k) = um(j,k)*cosphi(j)
  enddo
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)

  deallocate( temp )

  do k=1, nz
  do j=2, ny-1
    if (grady(j,k) /= missv)  grady(j,k) = grady(j,k)/cosphi(j)
  enddo
  enddo
  grady( 1,:) = missv
  grady(ny,:) = missv

  call gradz_2nd(1,ny,nz,1,um,zp,missv, gradz)

  deallocate( um )

  ! epf, epd
  fy(:,:) = missv
  fz(:,:) = missv
  do k=1, nz
  do j=1, ny
    if ( phi(j,k) /= missv ) then
      if ( gradz(j,k) /= missv )                      &
         fy(j,k) = (vu(j,k) - phi(j,k)*gradz(j,k))  &
                   *(-rho0(j,k)*r_earth*cosphi(j))
      if ( grady(j,k) /= missv )                             &
         fz(j,k) = (wu(j,k) + phi(j,k)*(grady(j,k)-f(j)))  &
                   *(-rho0(j,k)*r_earth*cosphi(j))
    end if
  enddo
  enddo

  deallocate( grady, gradz )
  deallocate( vpt, wpt )
  deallocate( vu, wu )
  deallocate( phi )

  call div_yz(ny,nz,1,fy,fz,lat,zp,missv, epd)
  do k=1, nz
  do j=1, ny
    if (epd(j,k) /= missv)  &
       epd(j,k) = epd(j,k)/(rho0(j,k)*r_earth*cosphi(j)) * 86400.
  enddo
  enddo

  deallocate( rho0 )


  do k=1, nza
  do j=1, nya
  do i=1, nxa
    if (epd(j,k) /= missv) then
      nepd(j,k) = nepd(j,k) + 1
      set(1)%var_out(i,j,k,1) = epd(j,k) !set(1)%var_out(i,j,k,1) + epd(j,k)
    end if
    if (fy(j,k) /= missv) then
      nfy(j,k) = nfy(j,k) + 1
      set(2)%var_out(i,j,k,1) = fy(j,k) !set(2)%var_out(i,j,k,1) + fy(j,k)
    end if
    if (fz(j,k) /= missv) then
      nfz(j,k) = nfz(j,k) + 1
      set(3)%var_out(i,j,k,1) = fz(j,k) !set(3)%var_out(i,j,k,1) + fz(j,k)
    end if
  enddo
  enddo
  enddo

  deallocate( epd )
  deallocate( fy, fz )

  call closenc(ncid)

  if (zbottom /= 0.) then
    iz = nza + 1
    do k=2, nza
      if (zp(k) > zbottom) then
        iz = k - 1
        EXIT
      end if
    enddo
    do iv=1, nvar
      set(iv)%var_out(:,:,:iz-1,:) = missv
    enddo
  end if

  do iv=1, nvar
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%axis = (/'     ','lat ','zp ','t'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = -999.
    set(iv)%axis2 = lat
    set(iv)%axis3 = zp
    set(iv)%axis4 = (year-2000)*100.+month
  enddo

  ! dump
  write(fname,'(a,i3.3,a,i3.3,a,i1,a)') trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname),wn1,'-',wn2,'.'//cutc0//'+',int(fct/24.),'.nc'
  call outnc(trim(fname),nvar,set,'')


  ENDDO  N_FCT


  ENDDO  N_UTC
  ENDDO  N_DAT


!  do k=1, nza
!  do j=1, nya
!    if (nepd(j,k) /= 0) then
!      set(1)%var_out(:,j,k,1) = set(1)%var_out(:,j,k,1)/nepd(j,k)
!    else
!      set(1)%var_out(:,j,k,1) = missv
!    end if
!    if (nfy(j,k) /= 0) then
!      set(2)%var_out(:,j,k,1) = set(2)%var_out(:,j,k,1)/nfy(j,k)
!    else
!      set(2)%var_out(:,j,k,1) = missv
!    end if
!    if (nfz(j,k) /= 0) then
!      set(3)%var_out(:,j,k,1) = set(3)%var_out(:,j,k,1)/nfz(j,k)
!    else
!      set(3)%var_out(:,j,k,1) = missv
!    end if
!  enddo
!  enddo
!
!  if (zbottom /= 0.) then
!    iz = nza + 1
!    do k=2, nza
!      if (zp(k) > zbottom) then
!        iz = k - 1
!        EXIT
!      end if
!    enddo
!    do iv=1, nvar
!      set(iv)%var_out(:,:,:iz-1,:) = missv
!    enddo
!  end if
!
!  do iv=1, nvar
!    set(iv)%nd(1) = nxa
!    set(iv)%nd(2) = nya
!    set(iv)%nd(3) = nza
!    set(iv)%nd(4) = nta
!    set(iv)%axis = (/'     ','lat ','zp ','t'/)
!    set(iv)%vname = ovarname(iv)
!    allocate( set(iv)%axis1(set(iv)%nd(1)) )
!    allocate( set(iv)%axis2(set(iv)%nd(2)) )
!    allocate( set(iv)%axis3(set(iv)%nd(3)) )
!    allocate( set(iv)%axis4(set(iv)%nd(4)) )
!    set(iv)%axis1 = -999.
!    set(iv)%axis2 = lat
!    set(iv)%axis3 = zp
!    set(iv)%axis4 = (year-2000)*100.+month
!  enddo
!
!  ! dump
!  fname = trim(ofdir)//trim(expname)//'/anal/'// &
!          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
!  call outnc(trim(fname),nvar,set,'')

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


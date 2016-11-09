! for the operational output

PROGRAM analysis_UM

  ! Non-hydrostatic TEM eq. in z-coord. (zonally symmetric density)
  ! (Andrews and McIntyre, 1978)

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2009, month =  7
  integer, parameter ::  dstart = 6, dend = 31   ! 6, 31
  integer, parameter ::  nutc  =  1, ustart = 0
  real,    parameter ::  fstart = 1.0, fend = 1.0, fitv = 1.0   ! [day]
  real,    parameter ::  missv = 1.e32, zbottom = 5.e3
  real,    parameter ::  zh_int = 12.0
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  integer, parameter ::  nivar = 5, nvar = 9
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070000'/
  ! files
  data                   ifdir   /'/data6/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'ctl'/
  data                   vartype /'std1'/
  data                   ofdir   /'/data6/kyh/UM_OPER/UM_OUT/'/
  data                   outname /'temu_yz'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'dz_dt'/
  data                   ivarname(4) /'theta'/
  data                   ivarname(5) /'temp'/
  ! output variables
  data                   ovarname(1) /'v_res'/
  data                   ovarname(2) /'w_res'/
  data                   ovarname(3) /'dudt'/
  data                   ovarname(4) /'cor'/
  data                   ovarname(5) /'adv_y'/
  data                   ovarname(6) /'adv_z'/
  data                   ovarname(7) /'epd'/
  data                   ovarname(8) /'f_y'/
  data                   ovarname(9) /'f_z'/


  real, dimension(:,:,:),   allocatable ::  u_in, v_in, w_in, pt_in, tem_in
  real, dimension(:,:,:),   allocatable ::  u, v, w, pt, rho
  real, dimension(:,:),     allocatable ::  um, vres, wres, ptm, rho0
  real, dimension(:,:),     allocatable ::  cor, advy, advz, epd, fy, fz
  real, dimension(:),       allocatable ::  lon, lat, zh, lonu, latv, zht
  real, dimension(:),       allocatable ::  time0

  real, dimension(:,:),     allocatable ::  phi, vflx, wflx
  real, dimension(:,:),     allocatable ::  grady, gradz, temp
  real, dimension(:),       allocatable ::  cosphi, f

  integer, dimension(:,:), allocatable ::  nvres, nwres, ncor, nadvy, nadvz
  integer, dimension(:,:), allocatable ::  nepd, nfy, nfz

  integer ::  nfct, i_time, date, utc, fct
  integer ::  year0, month0, date0, utc0
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, iyr, idt, iut, ifc
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid, tag
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=6)   ::  cmonth0
  character(len=10)  ::  cutc0
  character(len=128) ::  fname
  real, parameter ::  kappa_i = 1./kappa

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
          trim(vartype)//'.'//cutc0//'.a24h.nc'

  if ( trim(vartype) == 'std1' .and. fitv == 0.5 )  &
     fname = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
             'std0.'//cutc0//'.a12h.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (idt*iut*ifc == 1) then
    call diminfo(ncid,.TRUE.,zh_int, nx,ny,nz,nxu,nyv,nzt,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( zh (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    if (trim(c_axis(3,2)) /= empty)  allocate( zht (nzt) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,nzt,c_axis, lon,lat,zh,lonu,latv,zht)
  end if

  call timeinfo(ncid, nt0)
  allocate( time0(nt0) )
  call get1d(ncid,'t',nt0, time0)
  do n=1, nt0
    if ( abs(time0(n) - (real(fct+utc0)/24.-0.5*fitv)) < v_small ) then
      it = n
    end if
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
    allocate( nvres(ny,nz), nwres(ny,nz) )
    allocate( ncor(ny,nz) )
    allocate( nadvy(ny,nz), nadvz(ny,nz) )
    allocate( nepd(ny,nz) )
    allocate( nfy(ny,nz), nfz(ny,nz) )
    nvres(:,:) = 0
    nwres(:,:) = 0
    ncor (:,:) = 0
    nadvy(:,:) = 0
    nadvz(:,:) = 0
    nepd (:,:) = 0
    nfy  (:,:) = 0
    nfz  (:,:) = 0
  end if

  allocate( um(ny,nz), vres(ny,nz), wres(ny,nz) )
  allocate( ptm(ny,nz) )
  allocate( rho0(ny,nz) )

  allocate( cor(ny,nz), advy(ny,nz), advz(ny,nz) )
  allocate( epd(ny,nz) )
  allocate( fy(ny,nz), fz(ny,nz) )

  allocate( phi(ny,nz) )
  allocate( vflx(ny,nz), wflx(ny,nz) )

  allocate( grady(ny,nz), gradz(ny,nz) )
  allocate( temp(ny,nz) )

  allocate( u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz) )
  allocate( pt(nx,ny,nz) )
  allocate( rho(nx,ny,nz) )

  ! get var.s
  allocate( u_in(nxu,ny,nz) )
  call geta4d(ncid,trim(ivarname(1)),1,nxu,1,ny,1,nz,it,1, u_in)
  call u2rho(nx,ny,nz,1,u_in,missv, u)
  deallocate( u_in )
  allocate( v_in(nx,nyv,nz) )
  call geta4d(ncid,trim(ivarname(2)),1,nx,1,nyv,1,nz,it,1, v_in)
  call v2rho(nx,ny,nz,1,v_in,missv, v)
  deallocate( v_in )
  allocate( w_in(nx,ny,nzt) )
  call geta4d(ncid,trim(ivarname(3)),1,nx,1,ny,1,nzt,it,1, w_in)
  call t2rho(nx,ny,nz,1,w_in,missv, w)
  deallocate( w_in )

  allocate( pt_in(nx,ny,nzt) )
  allocate( tem_in(nx,ny,nzt) )
  call geta4d(ncid,trim(ivarname(4)),1,nx,1,ny,1,nzt,it,1, pt_in)
  call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nzt,it,1, tem_in)
  tem_in(:,:,:) = (tem_in(:,:,:)/pt_in(:,:,:))**kappa_i/tem_in(:,:,:)*1.e5/rd
  pt_in(:,:,:) = log(pt_in(:,:,:))
  call t2rho(nx,ny,nz,1,pt_in,missv, pt)
  pt(:,:,:) = exp(pt(:,:,:))
  tem_in(:,:,:) = log(tem_in(:,:,:))
  call t2rho(nx,ny,nz,1,tem_in,missv, rho)
  rho(:,:,:) = exp(rho(:,:,:))
  deallocate( tem_in )
  deallocate( pt_in )

  ! mean and perturbation
  call zonal_avg(nx,ny,nz,1,missv,u  ,um  )
  call zonal_avg(nx,ny,nz,1,missv,v  ,vres)
  call zonal_avg(nx,ny,nz,1,missv,w  ,wres)
  call zonal_avg(nx,ny,nz,1,missv,pt ,ptm )
  call zonal_avg(nx,ny,nz,1,missv,rho,rho0)
  deallocate( rho )

  do k=1, nz
  do j=1, ny
    u (:,j,k) = u (:,j,k) - um  (j,k)
    v (:,j,k) = v (:,j,k) - vres(j,k)
  enddo
  enddo
  do k=2, nz
  do j=1, ny
    w (:,j,k) = w (:,j,k) - wres(j,k)
    pt(:,j,k) = pt(:,j,k) - ptm (j,k)
  enddo
  enddo
  w (:,:,1) = missv
  pt(:,:,1) = missv

  ! phi
  call corr_zonal_avg(nx,ny,nz,1,missv,v,pt, vflx)
  call corr_zonal_avg(nx,ny,nz,1,missv,w,pt, wflx)
  deallocate( pt )

  call grady_2nd(1,ny,nz,1,ptm,lat,missv, grady)

  temp(:,1) = missv
  temp(:,2:) = log(ptm(:,2:))
  call gradz_2nd(1,ny,nz,1,temp,zh,missv, gradz)
  gradz(:,2:) = ptm(:,2:)*gradz(:,2:)

  phi(:,:) = missv
  do k=2, nz
  do j=2, ny-1
    phi(j,k) = (vflx(j,k)*gradz(j,k) - wflx(j,k)*grady(j,k)) &
                / (grady(j,k)**2 + gradz(j,k)**2)
  enddo
  enddo

  ! residual mean meridional circulations
  temp(:,:) = missv
  do k=2, nz
  do j=2, ny-1
    temp(j,k) = rho0(j,k)*phi(j,k)
  enddo
  enddo
  call gradz_2nd(1,ny,nz,1,temp,zh,missv, gradz)
  do k=2, nz
  do j=2, ny-1
    vres(j,k) = vres(j,k) - gradz(j,k)/rho0(j,k)
  enddo
  enddo
  vres( 1,:) = missv
  vres(ny,:) = missv
  vres( :,1) = missv

  do k=2, nz
  do j=2, ny-1
    temp(j,k) = temp(j,k)*cosphi(j)
  enddo
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)
  do k=2, nz
  do j=2, ny-1
    wres(j,k) = wres(j,k) + grady(j,k)/rho0(j,k)/cosphi(j)
  enddo
  enddo
  wres( 1,:) = missv
  wres(ny,:) = missv
  wres( :,1) = missv

  ! coriolis term
  cor(:,:) = missv
  do k=2, nz
  do j=2, ny-1
    cor(j,k) = vres(j,k)*f(j) * 86400.
  enddo
  enddo
  cor( 1,:) = missv
  cor(ny,:) = missv
  cor( :,1) = missv

  ! advy
  do k=1, nz
  do j=1, ny
    temp(j,k) = um(j,k)*cosphi(j)
  enddo
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)
  do k=1, nz
  do j=2, ny-1
    grady(j,k) = grady(j,k)/cosphi(j)
  enddo
  enddo
  grady( 1,:) = missv
  grady(ny,:) = missv

  do k=2, nz
  do j=2, ny-1
    advy(j,k) = -vres(j,k)*grady(j,k) * 86400.
  enddo
  enddo
  advy( 1,:) = missv
  advy(ny,:) = missv
  advy( :,1) = missv

  ! advz
  call gradz_2nd(1,ny,nz,1,um,zh,missv, gradz)
  do k=2, nz
  do j=2, ny-1
    advz(j,k) = -wres(j,k)*gradz(j,k) * 86400.
  enddo
  enddo
  advz( 1,:) = missv
  advz(ny,:) = missv
  advz( :,1) = missv

  ! epf, epd
  call corr_zonal_avg(nx,ny,nz,1,missv,v,u, vflx)
  call corr_zonal_avg(nx,ny,nz,1,missv,w,u, wflx)
  deallocate( u, v, w )

  do k=2, nz
  do j=2, ny-1
    fy(j,k) = (vflx(j,k) - phi(j,k)*gradz(j,k)) &
               *(-rho0(j,k)*r_earth*cosphi(j))
  enddo
  enddo
  fy( 1,:) = missv
  fy(ny,:) = missv
  fy( :,1) = missv

  do k=2, nz
  do j=2, ny-1
    fz(j,k) = (wflx(j,k) + phi(j,k)*(grady(j,k)-f(j))) &
               *(-rho0(j,k)*r_earth*cosphi(j))
  enddo
  enddo
  fz( 1,:) = missv
  fz(ny,:) = missv
  fz( :,1) = missv

  deallocate( grady, gradz )

  call div_yz(ny,nz,1,fy,fz,lat,zh,missv, epd)
  do k=2, nz
  do j=2, ny-1
    epd(j,k) = epd(j,k)/(rho0(j,k)*r_earth*cosphi(j)) * 86400.
  enddo
  enddo
  epd(:,1) = missv


  deallocate( temp )
  deallocate( vflx, wflx )
  deallocate( phi )
  deallocate( rho0 )
  deallocate( ptm )
  deallocate( um )


  do k=1, nza
  do j=1, nya
  do i=1, nxa
    if (vres(j,k) /= missv) then
      nvres(j,k) = nvres(j,k) + 1
      set(1)%var_out(i,j,k,1) = set(1)%var_out(i,j,k,1) + vres(j,k)
    end if
    if (wres(j,k) /= missv) then
      nwres(j,k) = nwres(j,k) + 1
      set(2)%var_out(i,j,k,1) = set(2)%var_out(i,j,k,1) + wres(j,k)
    end if
    if (cor(j,k) /= missv) then
      ncor(j,k) = ncor(j,k) + 1
      set(4)%var_out(i,j,k,1) = set(4)%var_out(i,j,k,1) + cor(j,k)
    end if
    if (advy(j,k) /= missv) then
      nadvy(j,k) = nadvy(j,k) + 1
      set(5)%var_out(i,j,k,1) = set(5)%var_out(i,j,k,1) + advy(j,k)
    end if
    if (advz(j,k) /= missv) then
      nadvz(j,k) = nadvz(j,k) + 1
      set(6)%var_out(i,j,k,1) = set(6)%var_out(i,j,k,1) + advz(j,k)
    end if
    if (epd(j,k) /= missv) then
      nepd(j,k) = nepd(j,k) + 1
      set(7)%var_out(i,j,k,1) = set(7)%var_out(i,j,k,1) + epd(j,k)
    end if
    if (fy(j,k) /= missv) then
      nfy(j,k) = nfy(j,k) + 1
      set(8)%var_out(i,j,k,1) = set(8)%var_out(i,j,k,1) + fy(j,k)
    end if
    if (fz(j,k) /= missv) then
      nfz(j,k) = nfz(j,k) + 1
      set(9)%var_out(i,j,k,1) = set(9)%var_out(i,j,k,1) + fz(j,k)
    end if
  enddo
  enddo
  enddo

  deallocate( vres, wres )
  deallocate( cor )
  deallocate( advy, advz )
  deallocate( epd )
  deallocate( fy, fz )

  call closenc(ncid)


  ENDDO  N_FCT


  ENDDO  N_UTC
  ENDDO  N_DAT


  do k=1, nza
  do j=1, nya
    if (nvres(j,k) /= 0) then
      set(1)%var_out(:,j,k,1) = set(1)%var_out(:,j,k,1)/nvres(j,k)
    else
      set(1)%var_out(:,j,k,1) = missv
    end if
    if (nwres(j,k) /= 0) then
      set(2)%var_out(:,j,k,1) = set(2)%var_out(:,j,k,1)/nwres(j,k)
    else
      set(2)%var_out(:,j,k,1) = missv
    end if
    if (ncor(j,k) /= 0) then
      set(4)%var_out(:,j,k,1) = set(4)%var_out(:,j,k,1)/ncor(j,k)
    else
      set(4)%var_out(:,j,k,1) = missv
    end if
    if (nadvy(j,k) /= 0) then
      set(5)%var_out(:,j,k,1) = set(5)%var_out(:,j,k,1)/nadvy(j,k)
    else
      set(5)%var_out(:,j,k,1) = missv
    end if
    if (nadvz(j,k) /= 0) then
      set(6)%var_out(:,j,k,1) = set(6)%var_out(:,j,k,1)/nadvz(j,k)
    else
      set(6)%var_out(:,j,k,1) = missv
    end if
    if (nepd(j,k) /= 0) then
      set(7)%var_out(:,j,k,1) = set(7)%var_out(:,j,k,1)/nepd(j,k)
    else
      set(7)%var_out(:,j,k,1) = missv
    end if
    if (nfy(j,k) /= 0) then
      set(8)%var_out(:,j,k,1) = set(8)%var_out(:,j,k,1)/nfy(j,k)
    else
      set(8)%var_out(:,j,k,1) = missv
    end if
    if (nfz(j,k) /= 0) then
      set(9)%var_out(:,j,k,1) = set(9)%var_out(:,j,k,1)/nfz(j,k)
    else
      set(9)%var_out(:,j,k,1) = missv
    end if
  enddo
  enddo

  do k=1, nza
  do j=1, nya
  do i=1, nxa
    if ( set(4)%var_out(i,j,k,1) /= missv .and. &
         set(5)%var_out(i,j,k,1) /= missv .and. &
         set(6)%var_out(i,j,k,1) /= missv .and. &
         set(7)%var_out(i,j,k,1) /= missv ) then
      set(3)%var_out(i,j,k,1) = set(4)%var_out(i,j,k,1) + &
                                set(5)%var_out(i,j,k,1) + &
                                set(6)%var_out(i,j,k,1) + &
                                set(7)%var_out(i,j,k,1)
    else
      set(3)%var_out(i,j,k,1) = missv
    end if
  enddo
  enddo
  enddo

  if (zbottom /= 0.) then
    iz = nza + 1
    do k=2, nza
      if (zh(k) > zbottom) then
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
    set(iv)%axis = (/'     ','lat ','zh ','t'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = -999.
    set(iv)%axis2 = lat
    set(iv)%axis3 = zh
    set(iv)%axis4 = (year-2000)*100.+month
  enddo

  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
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


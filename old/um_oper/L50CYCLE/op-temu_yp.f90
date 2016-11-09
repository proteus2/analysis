! for the operational output

PROGRAM analysis_UM

  ! Hydrostatic TEM eq. in p-coord.
  ! (Andrews et al., 1987)

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2010, month =  1
  integer, parameter ::  dstart = 6, dend = 31
  integer, parameter ::  nutc  =  1, ustart = 0
  real,    parameter ::  fstart = 1.0, fend = 5.0, fitv = 1.0   ! [day]
  real,    parameter ::  h_scale = 7.e3
  real,    parameter ::  missv = 1.e32, zbottom = 5.e3
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  integer, parameter ::  nivar = 4, nvar = 9
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070000'/
  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'gwdc'/
  data                   vartype /'press1'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   outname /'temu_yp'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'dz_dt'/
  data                   ivarname(4) /'temp'/
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


  real, dimension(:,:,:,:), allocatable ::  tem_in, tem_in2
  real, dimension(:,:,:),   allocatable ::  u, v, w, pt
  real, dimension(:,:),     allocatable ::  um, vres, wres, ptm, rho0
  real, dimension(:,:),     allocatable ::  cor, advy, advz, epd, fy, fz
  real, dimension(:),       allocatable ::  lon, lat, zp, p, lonu, latv, zpt
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

  call timeinfo(ncid, nt0)
  allocate( time0(nt0) )
  call get1d(ncid,'t',nt0, time0)
  do n=1, nt0
    if ( abs(time0(n) - real(fct)/24.) < v_small ) then
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

  ! get var.s
  allocate( tem_in(nx,ny,nz,8) )
  allocate( tem_in2(nx,ny,nz,8) )

  call geta4d(ncid,trim(ivarname(1)),1,nx,1,ny,1,nz,2,8, tem_in)
  call tempo_avg(nx,ny,nz,8,missv,tem_in, u)

  call geta4d(ncid,trim(ivarname(2)),1,nx,1,nyv,1,nz,2,8, tem_in(:,2:ny,:,:) )
  call tempo_avg(nx,nyv,nz,8,missv,tem_in(:,2:ny,:,:), v(:,2:ny,:))
  do k=1, nzr
  do j=2, nyv
    v(:,j,k) = 0.5*(v(:,j,k)+v(:,j+1,k))
  enddo
  enddo
  v(:,  1,:) = 0.
  v(:,nyr,:) = 0.
  call geta4d(ncid,trim(ivarname(3)),1,nx,1,ny,1,nz,2,8, tem_in)
  call geta4d(ncid,trim(ivarname(4)),1,nx,1,ny,1,nz,2,8, tem_in2)
  do n=1, 8
  do k=1, nzr
    tem_in(:,:,k,n) = tem_in(:,:,k,n)*g*h_scale/rd/tem_in2(:,:,k,n)
  enddo
  enddo
  call tempo_avg(nx,ny,nz,8,missv,tem_in, w)
  call tempo_avg(nx,ny,nz,8,missv,tem_in2, pt)

  deallocate( tem_in2 )
  deallocate( tem_in )

  ! mean and perturbation
  call zonal_avg(nx,ny,nz,1,missv,u ,um  )
  call zonal_avg(nx,ny,nz,1,missv,v ,vres)
  call zonal_avg(nx,ny,nz,1,missv,w ,wres)
  call zonal_avg(nx,ny,nz,1,missv,pt,ptm )


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
    if (u (i,j,k) /= missv)  u (i,j,k) = u (i,j,k) - um  (j,k)
    if (v (i,j,k) /= missv)  v (i,j,k) = v (i,j,k) - vres(j,k)
    if (w (i,j,k) /= missv)  w (i,j,k) = w (i,j,k) - wres(j,k)
    if (pt(i,j,k) /= missv)  pt(i,j,k) = pt(i,j,k) - ptm (j,k)
  enddo
  enddo
  enddo


  ! phi
  call corr_zonal_avg(nx,ny,nz,1,missv,v,pt, vflx)
  call corr_zonal_avg(nx,ny,nz,1,missv,w,pt, wflx)
  deallocate( pt )

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


  phi(:,:) = missv
  do k=1, nz
  do j=1, ny
    if ( grady(j,k) /= missv .and. gradz(j,k) /= missv .and.  &
         vflx(j,k) /= missv .and. wflx(j,k) /= missv )        &
       phi(j,k) = (vflx(j,k)*gradz(j,k) - wflx(j,k)*grady(j,k))  &
                 / (grady(j,k)**2 + gradz(j,k)**2)
  enddo
  enddo

  ! residual mean meridional circulations
  temp(:,:) = missv
  do k=1, nz
  do j=1, ny
    if (phi(j,k) /= missv)  temp(j,k) = rho0(j,k)*phi(j,k)
  enddo
  enddo
  call gradz_2nd(1,ny,nz,1,temp,zp,missv, gradz)
  do k=1, nz
  do j=1, ny
    if ( vres(j,k) /= missv .and. gradz(j,k) /= missv ) then
      vres(j,k) = vres(j,k) - gradz(j,k)/rho0(j,k)
    else
      vres(j,k) = missv
    end if
  enddo
  enddo

  do k=1, nz
  do j=1, ny
    if (temp(j,k) /= missv)  temp(j,k) = temp(j,k)*cosphi(j)
  enddo
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)
  do k=1, nz
  do j=2, ny-1
    if ( wres(j,k) /= missv .and. grady(j,k) /= missv ) then
      wres(j,k) = wres(j,k) + grady(j,k)/rho0(j,k)/cosphi(j)
    else
      wres(j,k) = missv
    end if
  enddo
  enddo
  wres( 1,:) = missv
  wres(ny,:) = missv

  ! coriolis term
  cor(:,:) = missv
  do k=1, nz
  do j=1, ny
    if (vres(j,k) /= missv)  cor(j,k) = vres(j,k)*f(j) * 86400.
  enddo
  enddo

  ! advy
  temp(:,:) = missv
  do k=1, nz
  do j=1, ny
    if (um(j,k) /= missv)  temp(j,k) = um(j,k)*cosphi(j)
  enddo
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)
  do k=1, nz
  do j=2, ny-1
    if (grady(j,k) /= missv)  grady(j,k) = grady(j,k)/cosphi(j)
  enddo
  enddo
  grady( 1,:) = missv
  grady(ny,:) = missv

  advy(:,:) = missv
  do k=1, nz
  do j=2, ny-1
    if ( vres(j,k) /= missv .and. grady(j,k) /= missv )  &
       advy(j,k) = -vres(j,k)*grady(j,k) * 86400.
  enddo
  enddo

  ! advz
  call gradz_2nd(1,ny,nz,1,um,zp,missv, gradz)

  advz(:,:) = missv
  do k=1, nz
  do j=2, ny-1
    if ( wres(j,k) /= missv .and. gradz(j,k) /= missv )  &
       advz(j,k) = -wres(j,k)*gradz(j,k) * 86400.
  enddo
  enddo

  ! epf, epd
  call corr_zonal_avg(nx,ny,nz,1,missv,v,u, vflx)
  call corr_zonal_avg(nx,ny,nz,1,missv,w,u, wflx)
  deallocate( u, v, w )

  fy(:,:) = missv
  do k=1, nz
  do j=1, ny
    if ( phi(j,k) /= missv .and. gradz(j,k) /= missv )  &
       fy(j,k) = (vflx(j,k) - phi(j,k)*gradz(j,k))      &
                 *(-rho0(j,k)*r_earth*cosphi(j))
  enddo
  enddo

  fz(:,:) = missv
  do k=1, nz
  do j=1, ny
    if ( phi(j,k) /= missv .and. grady(j,k) /= missv )     &
       fz(j,k) = (wflx(j,k) + phi(j,k)*(grady(j,k)-f(j)))  &
                 *(-rho0(j,k)*r_earth*cosphi(j))
  enddo
  enddo

  deallocate( grady, gradz )

  call div_yz(ny,nz,1,fy,fz,lat,zp,missv, epd)
  do k=1, nz
  do j=1, ny
    if (epd(j,k) /= missv)  &
       epd(j,k) = epd(j,k)/(rho0(j,k)*r_earth*cosphi(j)) * 86400.
  enddo
  enddo

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


! for the operational output

PROGRAM analysis_UM

  ! Hydrostatic TEM eq. in p-coord.
  ! (Andrews et al., 1987)

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  year = 2009, month =  7
  integer, parameter ::  dstart =  1, dend = 31  ! 1, 26
  integer, parameter ::  nutc  =  2, ustart = 0
  real,    parameter ::  fstart = 1.0, fend = 5.0, fitv = 1.0    ! [day]
  integer            ::  nfavg(5)              ! (fend-fstart)/fitv+1
  data                   nfavg/8,4,4,2,2/
  real,    parameter ::  h_scale = 7.e3, ts = 240.
  real,    parameter ::  missv = 1.e32, zb_out = 0.e3
  integer, parameter ::  nexcept = 1
  character(len=10)  ::  excdate(nexcept)
  integer, parameter ::  nvar = 15
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ovarname(nvar)

  ! exceptions
  data                   excdate(1) /'2009070000'/
  ! files
  data                   ifdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   expname /'ctl'/
  data                   vartype /'press1'/
  data                   ofdir   /'/data10/kyh/UM_OPER/UM_OUT/'/
  data                   outname /'tem_ypt'/
  ! output variables
  data                   ovarname(1 ) /'v_res'/
  data                   ovarname(2 ) /'w_res'/
  data                   ovarname(3 ) /'u_force'/
  data                   ovarname(4 ) /'t_force'/
  data                   ovarname(5 ) /'cor'/
  data                   ovarname(6 ) /'uadv_y'/
  data                   ovarname(7 ) /'uadv_z'/
  data                   ovarname(8 ) /'epd'/
  data                   ovarname(9 ) /'f_y'/
  data                   ovarname(10) /'f_z'/
  data                   ovarname(11) /'tadv_y'/
  data                   ovarname(12) /'tadv_z'/
  data                   ovarname(13) /'t_eddy'/
  data                   ovarname(14) /'u_tend'/
  data                   ovarname(15) /'t_tend'/

  real, dimension(:,:,:,:),   allocatable ::  u, v, w, pt
  real, dimension(:,:),       allocatable ::  utend, ttend
  real, dimension(:),         allocatable ::  lon, lat, zp, p, lonu, latv, t2pt
  real, dimension(:),         allocatable ::  fday, fcst, time

  integer ::  nfct, date, utc, fct, t_cnt
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, idt, iut, ifc
  integer ::  nf1, nf2, nf
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid1, ncid2, tag
  character(len=32)  ::  c_axis(3,2)
  character(len=4)   ::  cyear
  character(len=2)   ::  cmonth
  character(len=6)   ::  cmonth0
  character(len=10)  ::  cutc0
  character(len=128) ::  fname1, fname2, fname

  real, parameter ::  v_small = 1.e-5

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  nfavg(:) = nfavg(:) + 1

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
    CYCLE
  end if

  fname1 = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          cutc0//'.'//trim(vartype)//'.nc'

  fname2 = trim(ifdir)//trim(expname)//'/'//cmonth0//'/'//cutc0//'/'// &
          cutc0//'.press2.nc'

  ! get dim. sizes and axis
  if ( t_cnt == 1 ) then

    call opennc(trim(fname1),ncid1)
    call opennc(trim(fname2),ncid2)

    call diminfop(ncid1,.TRUE., nx,ny,nz,nxu,nyv,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( p  (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    call axisinfo(ncid1,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)
    p(:) = p(:)*100.
    allocate( zp(nz) )
    zp(:) = -h_scale*(log(p(:)/1.e5))
    allocate( t2pt(nz) )
    t2pt(:) = (1.e5/p(:))**kappa

    call timeinfo(ncid1, nf1)
    call timeinfo(ncid2, nf2)
    nf = nf1 + nf2
    allocate( fday(nf) )
    call get1d(ncid1,'t',nf1, fday(1:nf1))
    call get1d(ncid2,'t',nf2, fday(nf1+1:nf))

    call closenc(ncid1)
    call closenc(ncid2)

    nxr = NX      ;   nxa = NY
    nyr = NY      ;   nya = NZ
    nzr = NZ      ;   nza = NFCT  ! ;   iz = 28
    ntr =  1      ;   nta = NT

    do iv=1, nvar
      allocate( set(iv)%var_out(nxa,nya,nza,nta) )
      set(iv)%var_out(:,:,:,:) = missv
    enddo

    do iv=1, nvar
      set(iv)%nd(1) = nxa
      set(iv)%nd(2) = nya
      set(iv)%nd(3) = nza
      set(iv)%nd(4) = nta
      set(iv)%axis = (/'lat  ','zp ','fcst','t'/)
      set(iv)%vname = ovarname(iv)
      allocate( set(iv)%axis1(set(iv)%nd(1)) )
      allocate( set(iv)%axis2(set(iv)%nd(2)) )
      allocate( set(iv)%axis3(set(iv)%nd(3)) )
      allocate( set(iv)%axis4(set(iv)%nd(4)) )
      set(iv)%axis1 = lat
      set(iv)%axis2 = zp
      set(iv)%axis3 = fcst
      set(iv)%axis4 = time
    enddo

  end if


  N_FCT:   DO ifc=1, nfct


  fct = int((fstart + (ifc-1)*fitv)*24.)

  write(6,*)
  write(6,'(a10,a3,i4.4)') cutc0, ' + ', fct
  write(6,*)


  allocate( u(nx,ny,nz,nfavg(ifc)), v(nx,ny,nz,nfavg(ifc)), w(nx,ny,nz,nfavg(ifc)) )
  allocate( pt(nx,ny,nz,nfavg(ifc)) )
  allocate( utend(ny,nz), ttend(ny,nz) )

  ! get var.s
  call getvars(fname1,fname2,nx,ny,nz,nf,nfavg(ifc),fday,fct,fitv,h_scale,t2pt, &
                u,v,w,pt,utend,ttend)

  set(14)%var_out(:,:,ifc,t_cnt) = utend(:,:)
  set(15)%var_out(:,:,ifc,t_cnt) = ttend(:,:)

  deallocate( utend, ttend )

  ! calculate forcing
  call tem_force(nx,ny,nz,nfavg(ifc),lat,p,u,v,w,pt,ts,h_scale,missv,           &
                 set(1 )%var_out(:,:,ifc,t_cnt),set(2 )%var_out(:,:,ifc,t_cnt), &
                 set(5 )%var_out(:,:,ifc,t_cnt),set(6 )%var_out(:,:,ifc,t_cnt), &
                 set(7 )%var_out(:,:,ifc,t_cnt),set(8 )%var_out(:,:,ifc,t_cnt), &
                 set(9 )%var_out(:,:,ifc,t_cnt),set(10)%var_out(:,:,ifc,t_cnt), &
                 set(11)%var_out(:,:,ifc,t_cnt),set(12)%var_out(:,:,ifc,t_cnt), &
                 set(13)%var_out(:,:,ifc,t_cnt))

  deallocate( u, v, w, pt )


  do j=1, nya
  do i=2, nxa-1
    set(3)%var_out(i,j,ifc,t_cnt) = set(5)%var_out(i,j,ifc,t_cnt) + &
                                    set(6)%var_out(i,j,ifc,t_cnt) + &
                                    set(7)%var_out(i,j,ifc,t_cnt) + &
                                    set(8)%var_out(i,j,ifc,t_cnt)
  enddo
  enddo

  do j=1, nya
  do i=2, nxa-1
    set(4)%var_out(i,j,ifc,t_cnt) = set(11)%var_out(i,j,ifc,t_cnt) + &
                                    set(12)%var_out(i,j,ifc,t_cnt) + &
                                    set(13)%var_out(i,j,ifc,t_cnt) 
  enddo
  enddo

  if (zb_out /= 0.) then
    iz = nza + 1
    do k=2, nza
      if (zp(k) > zb_out) then
        iz = k - 1
        EXIT
      end if
    enddo
    do iv=1, nvar
      set(iv)%var_out(:,:iz-1,ifc,t_cnt) = missv
    enddo
  end if


  ENDDO  N_FCT
  ENDDO  N_UTC
  ENDDO  N_DAT


  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
  call outnc(trim(fname),nvar,set,'')


  deallocate( fday )
  deallocate( lon, lat, p, lonu, latv, zp, t2pt )
  deallocate( fcst )
  do iv=1, nvar
    deallocate( set(iv)%var_out )
  enddo


  STOP

END program


SUBROUTINE getvars(fname1,fname2,nx,ny,nz,nf,nfavg,fday,fct,fitv,h_scale,t2pt, &
                    u,v,w,pt,u_tend,t_tend)

  use netio
  use um_anal

  implicit none

  integer,             intent(in) ::  nx, ny, nz, nf, nfavg, fct
  real,                intent(in) ::  fitv, h_scale
  real, dimension(nz), intent(in) ::  t2pt
  real, dimension(nf), intent(in) ::  fday
  character(len=128),  intent(in) ::  fname1, fname2

  real, dimension(nx,ny,nz,nfavg), intent(out) ::  u, v, w, pt
  real, dimension(ny,nz),          intent(out) ::  u_tend, t_tend

  integer, parameter          ::  nt_p1 = 9  ! nt for .press1 file
  integer                     ::  i,j,k,n,na
  integer                     ::  it, ncid, ncid2, order_dt, nyh
  real                        ::  dayavg, c_dt, dt, dt1, dt2, dy
  real, dimension(nx,ny,nz,3) ::  gph
  real, dimension(nx,ny,nz)   ::  d_gph, temp
  real, dimension(3)          ::  t, dtcoef
  real, dimension(ny)         ::  dx

  character(len=128) ::  fname
  character(len=32 ) ::  ivarname(5)

  real, parameter ::  v_small = 1.e-5

  data  ivarname(1) /'u'/
  data  ivarname(2) /'v'/
  data  ivarname(3) /'dz_dt'/
  data  ivarname(4) /'temp'/
  data  ivarname(5) /'ht'/


  c_dt = fitv  ! [day]


  NA:   DO na=1, nfavg

  dayavg = real(fct)/24. + fitv*real(na-nfavg)/(nfavg-1)

  do n=1, nf
    if ( abs(fday(n) - dayavg) < v_small )  it = n
  enddo

  fname = fname1

  if (dayavg > 1.0) then
    fname = fname2
    it = it - nt_p1
  end if

  call opennc(trim(fname),ncid)

  call geta4d(ncid,trim(ivarname(1)),1,nx,1,ny  ,1,nz,it,1, u  (:, :  ,:,na))
  call geta4d(ncid,trim(ivarname(2)),1,nx,1,ny-1,1,nz,it,1, v  (:,2:ny,:,na))
  call geta4d(ncid,trim(ivarname(3)),1,nx,1,ny  ,1,nz,it,1, w  (:, :  ,:,na))
  call geta4d(ncid,trim(ivarname(4)),1,nx,1,ny  ,1,nz,it,1, pt (:, :  ,:,na))
  call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny  ,1,nz,it,1, gph(:, :  ,:,2 ))
  call geta1d(ncid,'t',it,1,t(2))

  order_dt = 2
  if (it /= 1) then
    call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,it-1,1, gph(:,:,:,1))
    call geta1d(ncid,'t',it-1,1,t(1))
  else if (dayavg <= 1.0) then
    order_dt = 1
    gph(:,:,:,1) = gph(:,:,:,2)
    t(1) = t(2)
  else
    call opennc(trim(fname1),ncid2)
    call geta4d(ncid2,trim(ivarname(5)),1,nx,1,ny,1,nz,nt_p1,1, gph(:,:,:,1))
    call geta1d(ncid2,'t',nt_p1,1,t(1))
    call closenc(ncid2)
  end if
  if ( (dayavg <= 1.0 .and. it < nt_p1) .or. &
       (dayavg > 1.0 .and. it < nf-nt_p1) ) then
    call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,it+1,1, gph(:,:,:,3))
    call geta1d(ncid,'t',it+1,1,t(3))
  else if (dayavg > 1.0) then
    order_dt = 1
    gph(:,:,:,3) = gph(:,:,:,2)
    t(3) = t(2)
  else
    call opennc(trim(fname2),ncid2)
    call geta4d(ncid2,trim(ivarname(5)),1,nx,1,ny,1,nz,1,1, gph(:,:,:,3))
    call geta1d(ncid2,'t',1,1,t(3))
    call closenc(ncid2)
  end if

  call closenc(ncid)

  ! dzdt
  if (order_dt == 2) then
    dt1 = (t(2) - t(1))*86400.
    dt2 = (t(3) - t(2))*86400.
    if ( abs(t(2)-0.5*(t(1)+t(3))) < v_small) then
      dt1 = (t(3) - t(1))*86400.
      dt2 = dt1
    end if
    dt = dt1 + dt2

    dtcoef(1) = -dt2/(dt1*dt)
    dtcoef(2) = (dt2-dt1)/(dt1*dt2)
    dtcoef(3) = dt1/(dt2*dt)

    d_gph(:,:,:) = 0.
    do n=1, 3
      d_gph(:,:,:) = d_gph(:,:,:) + dtcoef(n)*gph(:,:,:,n)
    enddo
  else
    dt = (t(3) - t(1))*86400.
    d_gph(:,:,:) = (gph(:,:,:,3) - gph(:,:,:,1))/dt
  end if
print*, sum(abs(d_gph))

  ! udzdx at u-grid
  nyh = (ny-1)/2
  dx(nyh+1) = r_earth*2.*pi/nx
  do j=1, nyh
    dx(nyh+1+j) = dx(nyh+1)*cos(0.5*pi/float(nyh)*j)
  enddo
  do j=1, nyh
    dx(nyh+1-j) = dx(nyh+1+j)
  enddo

  do k=1, nz
    do j=2, ny-1
      do i=1, nx-1
        temp(i,j,k) = u(i,j,k,na)*(gph(i+1,j,k,2) - gph(i,j,k,2))/dx(j)
      enddo
      temp(nx,j,k) = u(nx,j,k,na)*(gph(1,j,k,2) - gph(nx,j,k,2))/dx(j)
    enddo
    temp(1,1 ,k) = sum(temp(:,2   ,k))
    temp(1,ny,k) = sum(temp(:,ny-1,k))
    temp(:,1 ,k) = temp(1,1 ,k)
    temp(:,ny,k) = temp(1,ny,k)
  enddo

  d_gph(1,:,:) = d_gph(1,:,:) + 0.5*(temp(nx,:,:)+temp(1,:,:))
  do k=1, nz
  do j=1, ny
  do i=2, nx
    d_gph(i,j,k) = d_gph(i,j,k) + 0.5*(temp(i-1,j,k)+temp(i,j,k))
  enddo
  enddo
  enddo
print*, sum(abs(temp))

  ! vdzdy at v-grid :  temp(:,2:ny,:)
  temp(:,:,:) = 0.
  dy = r_earth*pi/(ny-1)
  do k=1, nz
  do j=2, ny
    temp(:,j,k) = v(:,j,k,na)*(gph(:,j,k,2) - gph(:,j-1,k,2))/dy
  enddo
  enddo

  do k=1, nz
  do j=2, ny-1
    d_gph(:,j,k) = d_gph(:,j,k) + 0.5*(temp(:,j,k)+temp(:,j+1,k))
  enddo
  enddo
print*, sum(abs(temp))
print*, sum(abs(w(:,:,:,na)))

  w(:,:,:,na) = g*h_scale/rd/pt(:,:,:,na) * ( w(:,:,:,na) - d_gph(:,:,:) )

  do k=1, nz
    pt(:,:,k,na) = pt(:,:,k,na) * t2pt(k)
  enddo

  temp(:,:,:) = u(:,:,:,na)
  u(1,:,:,na) = 0.5*(temp(nx,:,:)+temp(1,:,:))
  do k=1, nz
  do j=1, ny
  do i=2, nx
    u(i,j,k,na) = 0.5*(temp(i-1,j,k)+temp(i,j,k))
  enddo
  enddo
  enddo
  do k=1, nz
  do j=2, ny-1
    v(:,j,k,na) = 0.5*(v(:,j,k,na)+v(:,j+1,k,na))
  enddo
  enddo
  v(:, 1,:,na) = 0.
  v(:,ny,:,na) = 0.


  ENDDO  NA


  do k=1, nz
  do j=1, ny
    u_tend(j,k) = (sum(u(:,j,k,nfavg))-sum(u(:,j,k,1)))/nx / c_dt
  enddo
  enddo
  do k=1, nz
  do j=1, ny
    t_tend(j,k) = (sum(pt(:,j,k,nfavg))-sum(pt(:,j,k,1)))/nx / c_dt
  enddo
  enddo


  RETURN

END subroutine getvars


SUBROUTINE tem_force(nx,ny,nz,nt,lat,p,u,v,w,pt,ts,h_scale,missv, &
                      vres1,wres1,cor1,uadvy1,uadvz1,epd1,fy1,fz1,tadvy1,tadvz1,teddy1)

  use um_anal

  implicit none

  integer,                      intent(in) ::  nx, ny, nz, nt
  real,                         intent(in) ::  ts, h_scale, missv
  real, dimension(ny),          intent(in) ::  lat
  real, dimension(nz),          intent(in) ::  p
  real, dimension(nx,ny,nz,nt), intent(in) ::  u, v, w, pt

  real, dimension(ny,nz), intent(out) ::  vres1, wres1, fy1, fz1
  real, dimension(ny,nz), intent(out) ::  cor1, uadvy1, uadvz1, epd1
  real, dimension(ny,nz), intent(out) ::  tadvy1, tadvz1, teddy1

  real, dimension(ny,nz) ::  vres, wres, fy, fz
  real, dimension(ny,nz) ::  cor, uadvy, uadvz, epd
  real, dimension(ny,nz) ::  tadvy, tadvz, teddy

  real                   ::  coef
  real, dimension(ny)    ::  cosphi, f
  real, dimension(nz)    ::  zp, t2pt
  real, dimension(ny,nz) ::  um, ptm, rho0
  real, dimension(ny,nz) ::  phi, rvpt, rwpt, rvu, rwu
  real, dimension(ny,nz) ::  dptmdy, dptmdz, dptmd2z, divy_um, dumdz
  real, dimension(ny,nz) ::  grady, gradz, temp

  real, dimension(:,:,:), allocatable ::  prt, v_prt, w_prt

  integer ::  j,k,n


  cosphi(:) = cos(lat(:)*deg2rad)
  f(:) = 2.*ome_earth*sin(lat(:)*deg2rad)
  zp(:) = -h_scale*(log(p(:)/1.e5))
  do k=1, nz
    rho0(:,k) = p(k)/rd/ts
  enddo
  t2pt(:) = (1.e5/p(:))**kappa

  vres1  = 0.
  wres1  = 0.
  fy1    = 0.
  fz1    = 0.
  cor1   = 0.
  uadvy1 = 0.
  uadvz1 = 0.
  epd1   = 0.
  tadvy1 = 0.
  tadvz1 = 0.
  teddy1 = 0.


  DO n=1, nt


  coef = 1./real(nt-1)
  if ( n == 1 .or. n == nt )  coef = 0.5*coef

  ! mean and perturbation
  call zonal_avg(nx,ny,nz,1,0.,u (:,:,:,n),um  )
  call zonal_avg(nx,ny,nz,1,0.,v (:,:,:,n),vres)
  call zonal_avg(nx,ny,nz,1,0.,w (:,:,:,n),wres)
  call zonal_avg(nx,ny,nz,1,0.,pt(:,:,:,n),ptm )

  allocate( prt(nx,ny,nz), v_prt(nx,ny,nz), w_prt(nx,ny,nz) )

  do k=1, nz
  do j=1, ny
    v_prt(:,j,k) = v(:,j,k,n) - vres(j,k)
    w_prt(:,j,k) = w(:,j,k,n) - wres(j,k)
  enddo
  enddo

  do k=1, nz
  do j=1, ny
    prt(:,j,k) = pt(:,j,k,n) - ptm(j,k)
  enddo
  enddo
  call corr_zonal_avg(nx,ny,nz,1,0.,v_prt,prt, rvpt)
  call corr_zonal_avg(nx,ny,nz,1,0.,w_prt,prt, rwpt)

  do k=1, nz
  do j=1, ny
    prt(:,j,k) = u(:,j,k,n) - um(j,k)
  enddo
  enddo
  call corr_zonal_avg(nx,ny,nz,1,0.,v_prt,prt, rvu )
  call corr_zonal_avg(nx,ny,nz,1,0.,w_prt,prt, rwu )

  deallocate( prt, v_prt, w_prt )

  rvpt(:,:) = rho0(:,:)*rvpt(:,:)
  rwpt(:,:) = rho0(:,:)*rwpt(:,:)
  rvu (:,:) = rho0(:,:)*rvu (:,:)
  rwu (:,:) = rho0(:,:)*rwu (:,:)

  ! grad, pt
  call grady_2nd(1,ny,nz,1,ptm,lat,0., dptmdy)
  dptmdy(1 ,:) = missv
  dptmdy(ny,:) = missv

  temp(:,:) = log(ptm(:,:))
  call gradz_2nd_irr(1,ny,nz,1,temp,zp, dptmdz)
!a  call grad2z_2nd_irr(1,ny,nz,1,temp,zp,missv, dptmd2z)
  call gradz_2nd_irr(1,ny,nz,1,dptmdz,zp, dptmd2z)
  dptmd2z(:,:) = dptmd2z(:,:) + dptmdz(:,:)*dptmdz(:,:)
  dptmdz (:,:) = ptm(:,:)*dptmdz (:,:)
  dptmd2z(:,:) = ptm(:,:)*dptmd2z(:,:)

!  dptmdz (:,1 ) = missv
!  dptmdz (:,nz) = missv
!  dptmd2z(:,1 ) = missv
!  dptmd2z(:,nz) = missv
!a  dptmd2z(:,1 ) = dptmdz(:,1 )*dptmdz(:,1 )/ptm(:,1 ) ! assume dN/dz = 0
!a  dptmd2z(:,nz) = dptmdz(:,nz)*dptmdz(:,nz)/ptm(:,nz)

  ! grad, u
  do k=1, nz
    temp(:,k) = um(:,k)*cosphi(:)
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,0., divy_um)
  do k=1, nz
  do j=2, ny-1
    divy_um(j,k) = divy_um(j,k)/cosphi(j)
  enddo
  enddo
  divy_um(1 ,:) = missv
  divy_um(ny,:) = missv

  call gradz_2nd_irr(1,ny,nz,1,um,zp, dumdz)
!  dumdz(:,1 ) = missv
!  dumdz(:,nz) = missv

  ! phi
  do k=1, nz
    phi(:,k) = rvpt(:,k)/dptmdz(:,k)
  enddo
!  phi(:,1 ) = missv
!  phi(:,nz) = missv

  ! residual mean meridional circulations
  call gradz_2nd_irr(1,ny,nz,1,rvpt,zp, gradz)
  do k=1, nz
    vres(:,k) = vres(:,k) - ( gradz(:,k)-rvpt(:,k)*dptmd2z(:,k)/dptmdz(:,k) ) &
                /dptmdz(:,k)/rho0(:,k)
  enddo
!  vres(:,1 ) = missv
!  vres(:,nz) = missv

  do k=1, nz
    temp(:,k) = phi(:,k)*cosphi(:)
  enddo
!  call grady_2nd(1,ny,nz-2,1,temp(:,2:nz-1),lat,0., grady(:,2:nz-1))
  call grady_2nd(1,ny,nz,1,temp,lat,0., grady)
  do k=1, nz
  do j=2, ny-1
    wres(j,k) = wres(j,k) + grady(j,k)/rho0(j,k)/cosphi(j)
  enddo
  enddo
  wres(1 ,: ) = missv
  wres(ny,: ) = missv
!  wres(: ,1 ) = missv
!  wres(: ,nz) = missv

  ! coriolis term
  do k=1, nz
    cor(:,k) = vres(:,k)*f(:) * 86400.
  enddo
!  cor(:,1 ) = missv
!  cor(:,nz) = missv

  ! advy
  do k=1, nz
  do j=2, ny-1
    uadvy(j,k) = -vres(j,k)*divy_um(j,k) * 86400.
  enddo
  enddo
  uadvy(1 ,: ) = missv
  uadvy(ny,: ) = missv
!  uadvy(: ,1 ) = missv
!  uadvy(: ,nz) = missv

  do k=1, nz
  do j=2, ny-1
    tadvy(j,k) = -vres(j,k)*dptmdy(j,k)/t2pt(k) * 86400.
  enddo
  enddo
  tadvy(1 ,: ) = missv
  tadvy(ny,: ) = missv
!  tadvy(: ,1 ) = missv
!  tadvy(: ,nz) = missv

  ! advz
  do k=1, nz
  do j=2, ny-1
    uadvz(j,k) = -wres(j,k)*dumdz(j,k) * 86400.
  enddo
  enddo
  uadvz(1 ,: ) = missv
  uadvz(ny,: ) = missv
!  uadvz(: ,1 ) = missv
!  uadvz(: ,nz) = missv

  do k=1, nz
  do j=2, ny-1
    tadvz(j,k) = -wres(j,k)*dptmdz(j,k)/t2pt(k) * 86400.
  enddo
  enddo
  tadvz(1 ,: ) = missv
  tadvz(ny,: ) = missv
!  tadvz(: ,1 ) = missv
!  tadvz(: ,nz) = missv

  ! epf, epd
  do k=1, nz
    fy(:,k) = -r_earth*cosphi(:)*(rvu(:,k) - phi(:,k)*dumdz(:,k))
  enddo
!  fy(:,1 ) = missv
!  fy(:,nz) = missv

  do k=1, nz
  do j=2, ny-1
    fz(j,k) = -r_earth*cosphi(j)*(rwu(j,k) + phi(j,k)*(divy_um(j,k)-f(j)))
  enddo
  enddo
  fz(1 ,: ) = missv
  fz(ny,: ) = missv
!  fz(: ,1 ) = missv
!  fz(: ,nz) = missv

  do k=1, nz
    temp(:,k) = fy(:,k)*cosphi(:)
  enddo
!  call grady_2nd(1,ny,nz-2,1,temp(:,2:nz-1),lat,0., grady(:,2:nz-1))
  call grady_2nd(1,ny,nz,1,temp,lat,0., grady)
  do k=1, nz
  do j=2, ny-1
    epd(j,k) = grady(j,k)/r_earth/(cosphi(j)*cosphi(j))
  enddo
  enddo

  do k=1, nz
  do j=2, ny-1
    temp(j,k) = (divy_um(j,k)-f(j))*rvpt(j,k)
  enddo
  enddo
  temp(1 ,:) = 0.
  temp(ny,:) = 0.
  call gradz_2nd_irr(1,ny,nz,1,temp,zp, gradz)
  do k=1, nz
    temp(:,k) = (gradz(:,k) - temp(:,k)*dptmd2z(:,k)/dptmdz(:,k))/dptmdz(:,k)
  enddo
  call gradz_2nd_irr(1,ny,nz,1,rwu,zp, gradz)
  do k=1, nz
  do j=2, ny-1
    epd(j,k) = ( epd(j,k) - (gradz(j,k)+temp(j,k)) ) / rho0(j,k) * 86400.
  enddo
  enddo
  epd(1 ,: ) = missv
  epd(ny,: ) = missv
!  epd(: ,1 ) = missv
!  epd(: ,nz) = missv

  ! t_eddy
  do k=1, nz
  do j=2, ny-1
    temp(j,k) = dptmdy(j,k)*rvpt(j,k)
  enddo
  enddo
  temp(1 ,:) = 0.
  temp(ny,:) = 0.
  call gradz_2nd_irr(1,ny,nz,1,temp,zp, gradz)
  do k=1, nz
    temp(:,k) = (gradz(:,k) - temp(:,k)*dptmd2z(:,k)/dptmdz(:,k))/dptmdz(:,k)
  enddo
  call gradz_2nd_irr(1,ny,nz,1,rwpt,zp, gradz)
  do k=1, nz
    teddy(:,k) = -(gradz(:,k)+temp(:,k)) / rho0(:,k) / t2pt(k) * 86400.
  enddo
!  teddy(:,1 ) = missv
!  teddy(:,nz) = missv

  ! sum
  vres1  = vres1  + coef*vres
  wres1  = wres1  + coef*wres
  fy1    = fy1    + coef*fy
  fz1    = fz1    + coef*fz
  cor1   = cor1   + coef*cor
  uadvy1 = uadvy1 + coef*uadvy
  uadvz1 = uadvz1 + coef*uadvz
  epd1   = epd1   + coef*epd
  tadvy1 = tadvy1 + coef*tadvy
  tadvz1 = tadvz1 + coef*tadvz
  teddy1 = teddy1 + coef*teddy


  ENDDO  ! n


  RETURN

END subroutine tem_force


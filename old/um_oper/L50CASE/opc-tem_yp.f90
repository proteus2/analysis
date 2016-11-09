! for the operational output
! forcing with interval of tavg is averaged.

PROGRAM analysis_UM

  ! Hydrostatic TEM eq. in p-coord.
  ! (Andrews et al., 1987)

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  dcase = 10010712
  real,    parameter ::  fstart = 0.0, fend = 5.0, fitv = 2./24. ! [day]
  real,    parameter ::  h_scale = 7.e3, ts = 240.
  real,    parameter ::  missv = 1.e32, zb_out = 0.e3
  integer, parameter ::  nvar = 13
  character(len=128) ::  ifdir, expname, ofdir, outname
  character(len=64)  ::  ovarname(nvar)

  ! files
  data                   ifdir   /'/data11/kyh/UM_case/'/
  data                   expname /'ctl'/
  data                   ofdir   /'/data11/kyh/UM_case/'/
  data                   outname /'tem_ypf'/
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

  real, dimension(:,:,:), allocatable ::  u, v, w, pt
  real, dimension(:),     allocatable ::  lon, lat, zp, p, lonu, latv, t2pt
  real, dimension(:),     allocatable ::  t_fct

  integer            ::  nfct, nx, ny, nz, nxu, nyv
  integer            ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer            ::  ifc, iz
  integer            ::  i,j,k,n,iv, ncid
  character(len=10)  ::  datec
  character(len=32)  ::  c_axis(3,2)
  character(len=256) ::  fname, fdir

  real, parameter ::  v_small = 1.e-5

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvar) ::  set


  write(datec,'(i2.2,i8.8)') 20, dcase

  nfct = int((fend-fstart)/fitv)+1
  allocate( t_fct(nfct) )
  do ifc=1, nfct
    t_fct(ifc) = fstart + fitv*(ifc-1)
  enddo

  fdir = trim(ifdir)//trim(expname)//'/'//trim(datec)//'/'

  ! get axis
  write(fname,'(a,i3.3,a)') trim(fdir)//                            &
        'std.t_p.'//trim(datec)//'+',24*(int(fstart-v_small)+1),'.nc'

  call opennc(trim(fname),ncid)
  call diminfop(ncid,.TRUE., nx,ny,nz,nxu,nyv,c_axis)
  allocate( lon(nx), lat(ny), p(nz) )
  allocate( lonu(nxu), latv(nyv) )   ! not used
  call axisinfo(ncid,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)
  p(:) = p(:)*100.
  allocate( zp(nz) )
  zp(:) = -h_scale*(log(p(:)/1.e5))
  allocate( t2pt(nz) )
  t2pt(:) = (1.e5/p(:))**kappa
  call closenc(ncid)

  ! set output variables
  nxr = NX      ;   nxa =  1
  nyr = NY      ;   nya = NY
  nzr = NZ      ;   nza = NZ  ! ;   iz = 28
  ntr =  1      ;   nta = NFCT

  do iv=1, nvar
    allocate( set(iv)%var_out(nxa,nya,nza,nta) )
  enddo

  do iv=1, nvar
    set(iv)%nd(1) = nxa
    set(iv)%nd(2) = nya
    set(iv)%nd(3) = nza
    set(iv)%nd(4) = nta
    set(iv)%axis = (/'     ','lat ','zp ','fcst'/)
    set(iv)%vname = ovarname(iv)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = -999.
    set(iv)%axis2 = lat
    set(iv)%axis3 = zp
    set(iv)%axis4 = t_fct
  enddo

  ! allocate for input variables
  allocate( u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz) )
  allocate( pt(nx,ny,nz) )


  N_FCT:   DO ifc=1, nfct


  write(6,'(a,i3.3)') datec//' + ', int(t_fct(ifc)*24.)

  call getvars(datec,t_fct(ifc),nx,ny,nz,h_scale,t2pt,fdir, u,v,w,pt)

  ! calculate forcing
  call tem_force(nx,ny,nz,lat,p,u,v,w,pt,ts,h_scale,missv,              &
                 set(1 )%var_out(1,:,:,ifc),set(2 )%var_out(1,:,:,ifc), &
                 set(5 )%var_out(1,:,:,ifc),set(6 )%var_out(1,:,:,ifc), &
                 set(7 )%var_out(1,:,:,ifc),set(8 )%var_out(1,:,:,ifc), &
                 set(9 )%var_out(1,:,:,ifc),set(10)%var_out(1,:,:,ifc), &
                 set(11)%var_out(1,:,:,ifc),set(12)%var_out(1,:,:,ifc), &
                 set(13)%var_out(1,:,:,ifc))


  ENDDO  N_FCT


  deallocate( u, v, w, pt )


  set(3)%var_out(:,:,:,:) = missv
  do ifc=1, nfct
  do k=2, nza-1
  do j=2, nya-1
    set(3)%var_out(:,j,k,ifc) = set(5)%var_out(:,j,k,ifc) + &
                                set(6)%var_out(:,j,k,ifc) + &
                                set(7)%var_out(:,j,k,ifc) + &
                                set(8)%var_out(:,j,k,ifc)
  enddo
  enddo
  enddo

  set(4)%var_out(:,:,:,:) = missv
  do ifc=1, nfct
  do k=2, nza-1
  do j=2, nya-1
    set(4)%var_out(:,j,k,ifc) = set(11)%var_out(:,j,k,ifc) + &
                                set(12)%var_out(:,j,k,ifc) + &
                                set(13)%var_out(:,j,k,ifc) 
  enddo
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
      set(iv)%var_out(:,:,:iz-1,:) = missv
    enddo
  end if


  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//trim(datec)//'.nc'
  print*, fname
  call outnc(trim(fname),nvar,set,'')


  deallocate( lon, lat, p, lonu, latv, zp, t2pt )
  deallocate( t_fct )
  do iv=1, nvar
    deallocate( set(iv)%var_out )
  enddo


  STOP

END program


SUBROUTINE getvars(datec,t_fct,nx,ny,nz,h_scale,t2pt,fdir,u,v,w,pt)

  use netio
  use um_axis
  use um_anal

  implicit none

  integer,             intent(in) ::  nx, ny, nz
  real,                intent(in) ::  t_fct, h_scale
  real, dimension(nz), intent(in) ::  t2pt
  character(len=10),   intent(in) ::  datec
  character(len=256),  intent(in) ::  fdir

  real, dimension(nx,ny,nz), intent(out) ::  u, v, w, pt

  integer            ::  i,j,k,n
  integer            ::  it, nt, ncid, st, order_dt, it2, nt2, nyh
  real               ::  dt1, dt2, dt, dy
  character(len=256) ::  fname(5)
  character(len=32 ) ::  vartype, ivfile(5), ivarname(5)

  real, dimension(nx,ny,nz,3) ::  gph
  real, dimension(nx,ny,nz)   ::  d_gph, temp
  real, dimension(3)          ::  t, dtcoef
  real, dimension(ny)         ::  dx

  real, dimension(:), allocatable ::  fcst, fcst2

  real, parameter    ::  v_small = 1.e-5

  data  vartype     /'std'/
  ! var. names in input files
  data  ivfile(1)   /'u_p'/
  data  ivfile(2)   /'v_p'/
  data  ivfile(3)   /'w_p'/
  data  ivfile(4)   /'t_p'/
  data  ivfile(5)   /'gph_p'/
  ! input variables
  data  ivarname(1) /'u'/
  data  ivarname(2) /'v'/
  data  ivarname(3) /'dz_dt'/
  data  ivarname(4) /'temp'/
  data  ivarname(5) /'ht'/


  do i=1, 5
    write(fname(i),'(a,i3.3,a)') trim(fdir)//                         &
          trim(vartype)//'.'//trim(ivfile(i))//'.'//trim(datec)//'+', &
          24*(int(t_fct-v_small)+1),'.nc'
  enddo

  call opennc(fname(1),ncid)
  call timeinfo(ncid, nt)
  allocate( fcst(nt) )
  call get1d(ncid,'t',nt, fcst)
  it = 0 
  do n=1, nt
    if ( abs(fcst(n)-t_fct) < v_small )  it = n
  enddo
  if (it == 0)  print*, t_fct, 'CHECK THE FCST TIME TO ANALYSIS !!', fcst
  deallocate( fcst )
  call closenc(ncid)

  call opennc(trim(fname(1)),ncid)
  call geta4d(ncid,trim(ivarname(1)),1,nx,1,ny,1,nz,it,1, u)
  call closenc(ncid)

  call opennc(trim(fname(2)),ncid)
  call geta4d(ncid,trim(ivarname(2)),1,nx,1,ny-1,1,nz,it,1, v(:,2:ny,:) )
  call closenc(ncid)

  call opennc(trim(fname(3)),ncid)
  call geta4d(ncid,trim(ivarname(3)),1,nx,1,ny,1,nz,it,1, w)
  call closenc(ncid)
  call opennc(trim(fname(4)),ncid)
  call geta4d(ncid,trim(ivarname(4)),1,nx,1,ny,1,nz,it,1, pt)
  call closenc(ncid)


  order_dt = 2
  call opennc(trim(fname(5)),ncid)
  if ( it > 1 .and. it < nt ) then
    call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,it-1,3, gph)
    call geta1d(ncid,'t',it-1,3, t)
    call closenc(ncid)
  else if (it == 1) then
    call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,it,2, gph(:,:,:,2:3))
    call geta1d(ncid,'t',it,2, t(2:3))
    call closenc(ncid)
    write(fname(5),'(a,i3.3,a)') trim(fdir)//                         &
          trim(vartype)//'.'//trim(ivfile(5))//'.'//trim(datec)//'+', &
          24*(int(t_fct-v_small)),'.nc'
    st = nf_open(trim(fname(5)),nf_nowrite,ncid)
    if (st == 0) then
      call timeinfo(ncid, nt2)
      allocate( fcst2(nt2) )
      call get1d(ncid,'t',nt2, fcst2)
      it2 = nt2
      if ( abs(fcst2(nt2)-t_fct) < v_small )  it2 = nt2 - 1
      t(1) = fcst2(it2)
      deallocate( fcst2 )
      call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,it2,1, gph(:,:,:,1))
      call closenc(ncid)
    else
      order_dt = 1
      t(1) = t(2)
      gph(:,:,:,1) = gph(:,:,:,2)
    end if
  else if (it == nt) then
    call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,it-1,2, gph(:,:,:,1:2))
    call geta1d(ncid,'t',it-1,2, t(1:2))
    call closenc(ncid)
    write(fname(5),'(a,i3.3,a)') trim(fdir)//                         &
          trim(vartype)//'.'//trim(ivfile(5))//'.'//trim(datec)//'+', &
          24*(int(t_fct-v_small)+2),'.nc'
    st = nf_open(trim(fname(5)),nf_nowrite,ncid)
    if (st == 0) then
      call timeinfo(ncid, nt2)
      allocate( fcst2(nt2) )
      call get1d(ncid,'t',nt2, fcst2)
      it2 = 1
      if ( abs(fcst2(1)-t_fct) < v_small )  it2 = 2
      t(3) = fcst2(it2)
      deallocate( fcst2 )
      call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,it2,1, gph(:,:,:,3))
      call closenc(ncid)
    else
      order_dt = 1
      t(3) = t(2)
      gph(:,:,:,3) = gph(:,:,:,2)
    end if
  end if

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
        temp(i,j,k) = u(i,j,k)*(gph(i+1,j,k,2) - gph(i,j,k,2))/dx(j)
      enddo
      temp(nx,j,k) = u(nx,j,k)*(gph(1,j,k,2) - gph(nx,j,k,2))/dx(j)
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
    temp(:,j,k) = v(:,j,k)*(gph(:,j,k,2) - gph(:,j-1,k,2))/dy
  enddo
  enddo

  do k=1, nz
  do j=2, ny-1
    d_gph(:,j,k) = d_gph(:,j,k) + 0.5*(temp(:,j,k)+temp(:,j+1,k))
  enddo
  enddo
print*, sum(abs(temp))
print*, sum(abs(w))

  w(:,:,:) = g*h_scale/rd/pt(:,:,:) * ( w(:,:,:) - d_gph(:,:,:) )

  do k=1, nz 
    pt(:,:,k) = pt(:,:,k) * t2pt(k)
  enddo

  temp(:,:,:) = u(:,:,:)
  u(1,:,:) = 0.5*(temp(nx,:,:)+temp(1,:,:))
  do k=1, nz
  do j=1, ny
  do i=2, nx
    u(i,j,k) = 0.5*(temp(i-1,j,k)+temp(i,j,k))
  enddo
  enddo
  enddo

  do k=1, nz
  do j=2, ny-1
    v(:,j,k) = 0.5*(v(:,j,k)+v(:,j+1,k))
  enddo
  enddo
  v(:, 1,:) = 0.
  v(:,ny,:) = 0.


  RETURN

END subroutine getvars


SUBROUTINE tem_force(nx,ny,nz,lat,p,u,v,w,pt,ts,h_scale,missv, &
                     vres,wres,cor,uadvy,uadvz,epd,fy,fz,tadvy,tadvz,teddy)

  use um_anal

  implicit none

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  ts, h_scale, missv
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(nx,ny,nz), intent(in) ::  u, v, w, pt

  real, dimension(ny,nz), intent(out) ::  vres, wres, fy, fz
  real, dimension(ny,nz), intent(out) ::  cor, uadvy, uadvz, epd
  real, dimension(ny,nz), intent(out) ::  tadvy, tadvz, teddy

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


  ! mean and perturbation
  call zonal_avg(nx,ny,nz,1,0.,u ,um  )
  call zonal_avg(nx,ny,nz,1,0.,v ,vres)
  call zonal_avg(nx,ny,nz,1,0.,w ,wres)
  call zonal_avg(nx,ny,nz,1,0.,pt,ptm )

  allocate( prt(nx,ny,nz), v_prt(nx,ny,nz), w_prt(nx,ny,nz) )

  do k=1, nz
  do j=1, ny
    v_prt(:,j,k) = v(:,j,k) - vres(j,k)
    w_prt(:,j,k) = w(:,j,k) - wres(j,k)
  enddo
  enddo

  do k=1, nz
  do j=1, ny
    prt(:,j,k) = pt(:,j,k) - ptm(j,k)
  enddo
  enddo
  call corr_zonal_avg(nx,ny,nz,1,0.,v_prt,prt, rvpt)
  call corr_zonal_avg(nx,ny,nz,1,0.,w_prt,prt, rwpt)

  do k=1, nz
  do j=1, ny
    prt(:,j,k) = u(:,j,k) - um(j,k)
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
  call grad2z_2nd_irr(1,ny,nz,1,temp,zp,missv, dptmd2z)
  dptmd2z(:,:) = dptmd2z(:,:) + dptmdz(:,:)*dptmdz(:,:)
  dptmdz (:,:) = ptm(:,:)*dptmdz (:,:)
  dptmd2z(:,:) = ptm(:,:)*dptmd2z(:,:)

  dptmdz (:,1 ) = missv
  dptmdz (:,nz) = missv
  dptmd2z(:,1 ) = missv
  dptmd2z(:,nz) = missv

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
  dumdz(:,1 ) = missv
  dumdz(:,nz) = missv

  ! phi
  do k=2, nz-1
    phi(:,k) = rvpt(:,k)/dptmdz(:,k)
  enddo
  phi(:,1 ) = missv
  phi(:,nz) = missv

  ! residual mean meridional circulations
  call gradz_2nd_irr(1,ny,nz,1,rvpt,zp, gradz)
  do k=2, nz-1
    vres(:,k) = vres(:,k) - ( gradz(:,k)-rvpt(:,k)*dptmd2z(:,k)/dptmdz(:,k) ) &
                /dptmdz(:,k)/rho0(:,k)
  enddo
  vres(:,1 ) = missv
  vres(:,nz) = missv


  do k=2, nz-1
    temp(:,k) = phi(:,k)*cosphi(:)
  enddo
  call grady_2nd(1,ny,nz-2,1,temp(:,2:nz-1),lat,0., grady(:,2:nz-1))
  do k=2, nz-1
  do j=2, ny-1
    wres(j,k) = wres(j,k) + grady(j,k)/rho0(j,k)/cosphi(j)
  enddo
  enddo
  wres(1 ,: ) = missv
  wres(ny,: ) = missv
  wres(: ,1 ) = missv
  wres(: ,nz) = missv

  ! coriolis term
  do k=2, nz-1
    cor(:,k) = vres(:,k)*f(:) * 86400.
  enddo
  cor(:,1 ) = missv
  cor(:,nz) = missv

  ! advy
  do k=2, nz-1
  do j=2, ny-1
    uadvy(j,k) = -vres(j,k)*divy_um(j,k) * 86400.
  enddo
  enddo
  uadvy(1 ,: ) = missv
  uadvy(ny,: ) = missv
  uadvy(: ,1 ) = missv
  uadvy(: ,nz) = missv

  do k=2, nz-1
  do j=2, ny-1
    tadvy(j,k) = -vres(j,k)*dptmdy(j,k)/t2pt(k) * 86400.
  enddo
  enddo
  tadvy(1 ,: ) = missv
  tadvy(ny,: ) = missv
  tadvy(: ,1 ) = missv
  tadvy(: ,nz) = missv

  ! advz
  do k=2, nz-1
  do j=2, ny-1
    uadvz(j,k) = -wres(j,k)*dumdz(j,k) * 86400.
  enddo
  enddo
  uadvz(1 ,: ) = missv
  uadvz(ny,: ) = missv
  uadvz(: ,1 ) = missv
  uadvz(: ,nz) = missv

  do k=2, nz-1
  do j=2, ny-1
    tadvz(j,k) = -wres(j,k)*dptmdz(j,k)/t2pt(k) * 86400.
  enddo
  enddo
  tadvz(1 ,: ) = missv
  tadvz(ny,: ) = missv
  tadvz(: ,1 ) = missv
  tadvz(: ,nz) = missv

  ! epf, epd
  do k=2, nz-1
    fy(:,k) = -r_earth*cosphi(:)*(rvu(:,k) - phi(:,k)*dumdz(:,k))
  enddo
  fy(:,1 ) = missv
  fy(:,nz) = missv

  do k=2, nz-1
  do j=2, ny-1
    fz(j,k) = -r_earth*cosphi(j)*(rwu(j,k) + phi(j,k)*(divy_um(j,k)-f(j)))
  enddo
  enddo
  fz(1 ,: ) = missv
  fz(ny,: ) = missv
  fz(: ,1 ) = missv
  fz(: ,nz) = missv

  do k=2, nz-1
    temp(:,k) = fy(:,k)*cosphi(:)
  enddo
  call grady_2nd(1,ny,nz-2,1,temp(:,2:nz-1),lat,0., grady(:,2:nz-1))
  do k=2, nz-1
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
  do k=2, nz-1
    temp(:,k) = (gradz(:,k) - temp(:,k)*dptmd2z(:,k)/dptmdz(:,k))/dptmdz(:,k)
  enddo
  call gradz_2nd_irr(1,ny,nz,1,rwu,zp, gradz)
  do k=2, nz-1
  do j=2, ny-1
    epd(j,k) = ( epd(j,k) - (gradz(j,k)+temp(j,k)) ) / rho0(j,k) * 86400.
  enddo
  enddo
  epd(1 ,: ) = missv
  epd(ny,: ) = missv
  epd(: ,1 ) = missv
  epd(: ,nz) = missv

  ! t_eddy
  do k=1, nz
  do j=2, ny-1
    temp(j,k) = dptmdy(j,k)*rvpt(j,k)
  enddo
  enddo
  temp(1 ,:) = 0.
  temp(ny,:) = 0.
  call gradz_2nd_irr(1,ny,nz,1,temp,zp, gradz)
  do k=2, nz-1
    temp(:,k) = (gradz(:,k) - temp(:,k)*dptmd2z(:,k)/dptmdz(:,k))/dptmdz(:,k)
  enddo
  call gradz_2nd_irr(1,ny,nz,1,rwpt,zp, gradz)
  do k=2, nz-1
    teddy(:,k) = -(gradz(:,k)+temp(:,k)) / rho0(:,k) / t2pt(k) * 86400.
  enddo
  teddy(:,1 ) = missv
  teddy(:,nz) = missv


  RETURN

END subroutine tem_force


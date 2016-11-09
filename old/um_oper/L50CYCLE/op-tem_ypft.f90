PROGRAM TEM

  ! Hydrostatic TEM eq. in p-coord.
  ! (Andrews et al., 1987)

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  nvo = 15
  real,    parameter ::  h_scale = 7.e3, ts = 240.
  real,    parameter ::  missv = 1.e32
  real,    parameter ::  zb_out = 0.e3
  character(len=64)  ::  ovarname(nvo)
  character(len=128) ::  f_namelist

  integer ::  yyyymm, date(2), utc(12), nf_intg(999)
  real    ::  fct(2), fct_intg
  character(len=10)  ::  expname, exception(99)
  character(len=128) ::  file_i, file_o

  namelist /ANALCASE/ EXPNAME, YYYYMM
  namelist /PARAM/ DATE, UTC, FCT, FCT_INTG, EXCEPTION, NF_INTG
  namelist /FILEIO/ FILE_I, FILE_O

  ! input variables
  data  f_namelist  /'./namelist/nl.input'/
  ! output variables
  data  ovarname(1 ) /'v_res'/
  data  ovarname(2 ) /'w_res'/
  data  ovarname(3 ) /'u_force'/
  data  ovarname(4 ) /'t_force'/
  data  ovarname(5 ) /'cor'/
  data  ovarname(6 ) /'uadv_y'/
  data  ovarname(7 ) /'uadv_z'/
  data  ovarname(8 ) /'epd'/
  data  ovarname(9 ) /'f_y'/
  data  ovarname(10) /'f_z'/
  data  ovarname(11) /'tadv_y'/
  data  ovarname(12) /'tadv_z'/
  data  ovarname(13) /'t_eddy'/
  data  ovarname(14) /'u_tend'/
  data  ovarname(15) /'t_tend'/

  integer ::  nutc, nt, nfct, nx, ny, nz, nxu, nyv, nf, nf1, nf2
  integer ::  nxr, nyr, nzr, ntr, nd1a, nd2a, nd3a, nd4a
  integer ::  i,j,k,n, iv, idat, iutc, ifct, i_time, it, ncid, ncid2, iz, ni
  integer ::  year, month
  logical ::  l_getdim, ex1, ex2
  character(len=10)  ::  timec
  character(len=32)  ::  c_axis(3,2)
  character(len=256) ::  fname, fname1, fname2

  real, dimension(:,:,:,:), allocatable ::  u, v, w, pt
  real, dimension(:,:),     allocatable ::  utend, ttend
  real, dimension(:),       allocatable ::  lon, lat, p, lonu, latv
  real, dimension(:),       allocatable ::  zp, t2pt, t, t_fct, tmed_fct, fday

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvo) ::  set

! READ NAMELISTS

  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! DEFINE AXIS

  do n=1, 12
    if (utc(n) < 0) then  ;  nutc = n - 1  ;  EXIT  ;  end if
  enddo

  call gettaxis

  l_getdim = .TRUE.

! LOOP

  i_time = 0


  L_DAT:  DO idat=date(1), date(2)
  L_UTC:  DO iutc=1, nutc


  year  = int(yyyymm/100)                    
  month = yyyymm - year*100

  i_time = i_time + 1


  L_FCT:  DO ifct=1, nfct


  write(timec ,'(i4.4,3i2.2)') year, month, idat, utc(iutc)

  write(6,*)
  write(6,'(a,i3.3,a)') ' '//timec//' + ',int(t_fct(ifct)*24.),' h'

  fname1 = trim(file_i)//'/'//timec//'/'//trim(expname)//'.std_p.'//timec//'+024.nc'
  fname2 = trim(file_i)//'/'//timec//'/'//trim(expname)//'.std_p.'//timec//'+120.nc'
  call check_ex(fname1, ex1)  ;  call check_ex(fname2, ex2)
  if (.not. (ex1 .and. ex2))  CYCLE

  if (l_getdim)  call getdim

  ! allocate output
  nxr = NX   ;   nd1a = NY
  nyr = NY   ;   nd2a = NZ
  nzr = NZ   ;   nd3a = NFCT
  ntr = NF   ;   nd4a = NT

  if (.not. allocated(set(1)%var_out)) then
    do iv=1, nvo
      allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
      set(iv)%var_out(:,:,:,:) = missv
    enddo
  end if

  allocate( u(nxr,nyr,nzr,nf_intg(ifct)), pt(nxr,nyr,nzr,nf_intg(ifct)) )
  allocate( v(nxr,nyr,nzr,nf_intg(ifct)), w (nxr,nyr,nzr,nf_intg(ifct)) )
  allocate( utend(nyr,nzr), ttend(nyr,nzr) )

  ! get var.s
  call getvars

  set(14)%var_out(:,:,ifct,i_time) = utend(:,:)
  set(15)%var_out(:,:,ifct,i_time) = ttend(:,:)

  deallocate( utend, ttend )

  ! calculate forcing
  call tem_force(nx,ny,nz,nf_intg(ifct),lat,p,u,v,w,pt,ts,h_scale,missv,           &
                 set(1 )%var_out(:,:,ifct,i_time),set(2 )%var_out(:,:,ifct,i_time), &
                 set(5 )%var_out(:,:,ifct,i_time),set(6 )%var_out(:,:,ifct,i_time), &
                 set(7 )%var_out(:,:,ifct,i_time),set(8 )%var_out(:,:,ifct,i_time), &
                 set(9 )%var_out(:,:,ifct,i_time),set(10)%var_out(:,:,ifct,i_time), &
                 set(11)%var_out(:,:,ifct,i_time),set(12)%var_out(:,:,ifct,i_time), &
                 set(13)%var_out(:,:,ifct,i_time))

  deallocate( u, pt, v, w )


  do j=1, nd2a
  do i=2, nd1a-1
    set(3)%var_out(i,j,ifct,i_time) = set(5)%var_out(i,j,ifct,i_time) + &
                                      set(6)%var_out(i,j,ifct,i_time) + &
                                      set(7)%var_out(i,j,ifct,i_time) + &
                                      set(8)%var_out(i,j,ifct,i_time)
  enddo
  enddo

  do j=1, nd2a
  do i=2, nd1a-1
    set(4)%var_out(i,j,ifct,i_time) = set(11)%var_out(i,j,ifct,i_time) + &
                                      set(12)%var_out(i,j,ifct,i_time) + &
                                      set(13)%var_out(i,j,ifct,i_time) 
  enddo
  enddo

  if (zb_out /= 0.) then
    iz = nd2a + 1
    do k=2, nd2a
      if (zp(k) > zb_out) then  ;  iz = k - 1  ;  EXIT  ;  end if
    enddo
    do iv=1, nvo
      set(iv)%var_out(:,:iz-1,ifct,i_time) = missv
    enddo
  end if


  ENDDO  L_FCT


  ENDDO  L_UTC
  ENDDO  L_DAT


  do iv=1, nvo
    set(iv)%vname = ovarname(iv)
    set(iv)%axis = (/'lat  ','zp ','fcst','t'/)
    set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = lat
    set(iv)%axis2 = zp
    set(iv)%axis3 = tmed_fct
    set(iv)%axis4 = t
  enddo

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nvo,set,'TEM')

! END

  deallocate( t_fct, tmed_fct, t, fday )
  deallocate( lon, lat, p, lonu, latv, zp, t2pt )
  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3, set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  STOP


  CONTAINS

  SUBROUTINE gettaxis

    ! t
    nt = (date(2)-date(1)+1)*nutc
    allocate( t(nt) )
    n = 0
    do idat=date(1), date(2)
    do iutc=1, nutc
      n = n + 1
      t(n) = float(idat) + float(utc(iutc))/24.
    enddo
    enddo

    ! fct
    nfct = int((fct(2)-fct(1))/fct_intg)+1
    allocate( t_fct(nfct), tmed_fct(nfct) )
    do ifct=1, nfct
      t_fct(ifct) = fct(1) + fct_intg*(ifct-1)
    enddo
    tmed_fct(:) = t_fct(:) - 0.5*fct_intg

  END subroutine gettaxis

  SUBROUTINE getdim

    call opennc(fname1,ncid )
    call opennc(fname2,ncid2)

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
    allocate( t2pt(nz) )
    t2pt(:) = (1.e5/p(:))**kappa

    call timeinfo(ncid , nf1)
    call timeinfo(ncid2, nf2)
    nf = nf1 + nf2
    allocate( fday(nf) )
    call get1d(ncid ,'t',nf1, fday(1:nf1))
    call get1d(ncid2,'t',nf2, fday(nf1+1:nf))

    call closenc(ncid )
    call closenc(ncid2)

    l_getdim = .FALSE.

  END subroutine getdim

  SUBROUTINE getvars

    implicit none

    character(len=64)              ::  ivarname(5)

    ! input variables
    data  ivarname(1) /'u'/
    data  ivarname(2) /'v'/
    data  ivarname(3) /'dz_dt'/
    data  ivarname(4) /'temp'/
    data  ivarname(5) /'ht'/

    integer                        ::  order_dt, nyh
    real                           ::  dayintg, dt, dt1, dt2, dy
    real, dimension(nxr,nyr,nzr,3) ::  gph
    real, dimension(nxr,nyr,nzr)   ::  d_gph, temp
    real, dimension(2,nzr)         ::  adv_pole
    real, dimension(3)             ::  t3, dtcoef
    real, dimension(nyr)           ::  dx
    real, parameter                ::  v_small = 1.e-5


    NI:   DO ni=1, nf_intg(ifct)


    dayintg = t_fct(ifct) + fct_intg*real(ni-nf_intg(ifct))/(nf_intg(ifct)-1)

    do n=1, ntr
      if ( abs(fday(n) - dayintg) < v_small )  it = n
    enddo

    fname = fname1

    if (dayintg > 1.0) then
      fname = fname2
      it = it - nf1
    end if

    call opennc(trim(fname),ncid)

    call geta4d(ncid,trim(ivarname(1)),1,nxr,1,nyr  ,1,nzr,it,1, u  (:, :   ,:,ni))
    call geta4d(ncid,trim(ivarname(2)),1,nxr,1,nyr-1,1,nzr,it,1, v  (:,2:nyr,:,ni))
    call geta4d(ncid,trim(ivarname(3)),1,nxr,1,nyr  ,1,nzr,it,1, w  (:, :   ,:,ni))
    call geta4d(ncid,trim(ivarname(4)),1,nxr,1,nyr  ,1,nzr,it,1, pt (:, :   ,:,ni))
    call geta4d(ncid,trim(ivarname(5)),1,nxr,1,nyr  ,1,nzr,it,1, gph(:, :   ,:,2 ))
    call geta1d(ncid,'t',it,1,t3(2))

    order_dt = 2
    if (it /= 1) then
      call geta4d(ncid,trim(ivarname(5)),1,nxr,1,nyr,1,nzr,it-1,1, gph(:,:,:,1))
      call geta1d(ncid,'t',it-1,1,t3(1))
    else if (dayintg <= 1.0) then
      order_dt = 1
      gph(:,:,:,1) = gph(:,:,:,2)
      t3(1) = t3(2)
    else
      call opennc(trim(fname1),ncid2)
      call geta4d(ncid2,trim(ivarname(5)),1,nxr,1,nyr,1,nzr,nf1,1, gph(:,:,:,1))
      call geta1d(ncid2,'t',nf1,1,t3(1))
      call closenc(ncid2)
    end if
      if ( (dayintg <= 1.0 .and. it < nf1) .or. &
         (dayintg > 1.0 .and. it < ntr-nf1) ) then
      call geta4d(ncid,trim(ivarname(5)),1,nxr,1,nyr,1,nzr,it+1,1, gph(:,:,:,3))
      call geta1d(ncid,'t',it+1,1,t3(3))
    else if (dayintg > 1.0) then
      order_dt = 1
      gph(:,:,:,3) = gph(:,:,:,2)
      t3(3) = t3(2)
    else
      call opennc(trim(fname2),ncid2)
      call geta4d(ncid2,trim(ivarname(5)),1,nxr,1,nyr,1,nzr,1,1, gph(:,:,:,3))
      call geta1d(ncid2,'t',1,1,t3(3))
      call closenc(ncid2)
    end if

    call closenc(ncid)

    ! dzdt
    dt = (t3(3) - t3(1))*86400.
    if (order_dt == 2) then
      dt1 = (t3(2) - t3(1))*86400.
      dt2 = (t3(3) - t3(2))*86400.

      dtcoef(1) = -dt2/(dt1*dt)
      dtcoef(2) = (dt2-dt1)/(dt1*dt2)
      dtcoef(3) = dt1/(dt2*dt)

      d_gph(:,:,:) = 0.
      do n=1, 3
        d_gph(:,:,:) = d_gph(:,:,:) + dtcoef(n)*gph(:,:,:,n)
      enddo
    else
      d_gph(:,:,:) = (gph(:,:,:,3) - gph(:,:,:,1))/dt
    end if
print*, sum(abs(d_gph))

    ! udzdx at u-grid
    nyh = (nyr-1)/2
    dx(nyh+1) = r_earth*2.*pi/nxr
    do j=1, nyh
      dx(nyh+1+j) = dx(nyh+1)*cos(0.5*pi/float(nyh)*j)
    enddo
    do j=1, nyh
      dx(nyh+1-j) = dx(nyh+1+j)
    enddo

    do k=1, nzr
      do j=2, nyr-1
        do i=1, nxr-1
          temp(i,j,k) = u(i,j,k,ni)*(gph(i+1,j,k,2) - gph(i,j,k,2))/dx(j)
        enddo
        temp(nxr,j,k) = u(nxr,j,k,ni)*(gph(1,j,k,2) - gph(nxr,j,k,2))/dx(j)
      enddo
      adv_pole(1,k) = sum(temp(:,2    ,k)) / nxr
      adv_pole(2,k) = sum(temp(:,nyr-1,k)) / nxr
    enddo

    do k=1, nzr
    do j=2, nyr-1
      d_gph(1,j,k) = d_gph(1,j,k) + 0.5*(temp(nxr,j,k)+temp(1,j,k))
      do i=2, nxr
        d_gph(i,j,k) = d_gph(i,j,k) + 0.5*(temp(i-1,j,k)+temp(i,j,k))
      enddo
    enddo
    enddo
print*, sum(abs(temp))

    ! vdzdy at v-grid :  temp(:,2:nyr,:)
    temp(:,:,:) = 0.
    dy = r_earth*pi/(nyr-1)
    do k=1, nzr
      do j=2, nyr
        temp(:,j,k) = v(:,j,k,ni)*(gph(:,j,k,2) - gph(:,j-1,k,2))/dy
      enddo
      adv_pole(1,k) = adv_pole(1,k) + 0.5*sum(temp(:,2  ,k)+temp(:,3    ,k)) / nxr
      adv_pole(2,k) = adv_pole(2,k) + 0.5*sum(temp(:,nyr,k)+temp(:,nyr-1,k)) / nxr
    enddo

    do k=1, nzr
    do j=2, nyr-1
      d_gph(:,j,k) = d_gph(:,j,k) + 0.5*(temp(:,j,k)+temp(:,j+1,k))
    enddo
    enddo

    ! add the term of horizontal advection at the poles
    do k=1, nzr
      d_gph(:,1  ,k) = d_gph(:,1  ,k) + adv_pole(1,k)
      d_gph(:,nyr,k) = d_gph(:,nyr,k) + adv_pole(2,k)
    enddo

print*, sum(abs(temp))
print*, sum(abs(w(:,:,:,ni)))

    w(:,:,:,ni) = g*h_scale/rd/pt(:,:,:,ni) * ( w(:,:,:,ni) - d_gph(:,:,:) )

    do k=1, nzr
      pt(:,:,k,ni) = pt(:,:,k,ni) * t2pt(k)
    enddo

    temp(:,:,:) = u(:,:,:,ni)
    u(1,:,:,ni) = 0.5*(temp(nxr,:,:)+temp(1,:,:))
    do k=1, nzr
    do j=1, nyr
    do i=2, nxr
      u(i,j,k,ni) = 0.5*(temp(i-1,j,k)+temp(i,j,k))
    enddo
    enddo
    enddo
    do k=1, nzr
    do j=2, nyr-1
      v(:,j,k,ni) = 0.5*(v(:,j,k,ni)+v(:,j+1,k,ni))
    enddo
    enddo
    v(:,1  ,:,ni) = 0.
    v(:,nyr,:,ni) = 0.


    ENDDO  NI


    do k=1, nzr
    do j=1, nyr
      utend(j,k) = (sum(u(:,j,k,nf_intg(ifct)))-sum(u(:,j,k,1)))/nxr / fct_intg
    enddo
    enddo
    do k=1, nzr
    do j=1, nyr
      ttend(j,k) = (sum(pt(:,j,k,nf_intg(ifct)))-sum(pt(:,j,k,1)))/nxr / fct_intg
    enddo
    enddo

  END subroutine getvars

END program TEM


SUBROUTINE check_ex(fname,existence)
  
  implicit none
  
  character(len=*), intent(in) ::  fname

  logical ::  existence
  
  
  inquire(file=trim(fname), exist=existence)
  if (.not. existence)  print*, '    ',trim(fname),' not found. - passed'
  
  RETURN
  
END subroutine check_ex


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


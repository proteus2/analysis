PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  nmonth =  1, mstart =  1
  integer, parameter ::  nyear  =  1, ystart = 2009
  integer, parameter ::  t_interval = 1
  real,    parameter ::  h_scale = 7.0e3
  real,    parameter ::  missv = 1.e32
  integer, parameter ::  nivar = 5, nvar = 6
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data5/kyh/umres/'/
  data                   expname /'ctl'/
  data                   vartype /'var_xypt'/
  data                   ofdir   /'/data5/kyh/umres/'/ 
  data                   outname /'tem_yp'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'wp'/
  data                   ivarname(4) /'temp'/
  data                   ivarname(5) /'z'/
  ! output variables
  data                   ovarname(1) /'v_res'/
  data                   ovarname(2) /'w_res'/
  data                   ovarname(3) /'dudt'/
  data                   ovarname(4) /'cor'/
  data                   ovarname(5) /'adv_y'/
  data                   ovarname(6) /'adv_z'/

!  data                   ovarname(4) /'temp'/
!  data                   ovarname(5) /'z'/


  real, dimension(:,:,:,:),   allocatable ::  u_in, v_in, wp_in, te_in
  real, dimension(:,:,:,:),   allocatable ::  u, v, wp, te
  real, dimension(:,:,:),     allocatable ::  um, vres, wres, tm
  real, dimension(:,:,:),     allocatable ::  cor, advy, advz, epd
  real, dimension(:),         allocatable ::  lon, lat, p, zp, lonu, latv, pt

  real, dimension(:,:,:),     allocatable ::  vtm, vum, vwm
  real, dimension(:,:,:),     allocatable ::  term1, stabz, t1, temp
  real, dimension(:),         allocatable ::  t0, cosphi, f

  integer ::  year, month
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, imn, iyr, nt_in
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz
  integer ::  i,j,k,n,iv, ncid, n_in1, n_in2
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
  type(vset), dimension(nvar) ::  set


  N_MON:   DO imn=1, nmonth
  N_YEAR:  DO iyr=1, nyear


  year  = ystart + (iyr-1)
  month = mstart + (imn-1)
  if (month > 12)  month = month - 12
  write(cyear, '(i4.4)') year
  write(cmonth,'(i2.2)') month

  fname = trim(ifdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(vartype)//'.'//cyear//'.'//cmonth//'.nc'

  call opennc(trim(fname),ncid)

  ! get dim. sizes and axis
  if (imn == 1 .and. iyr == 1) then
    call diminfop(ncid,.FALSE., nx,ny,nz,nxu,nyv,c_axis)
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(2,1)) /= empty)  allocate( lat(ny) )
    if (trim(c_axis(3,1)) /= empty)  allocate( p  (nz) )
!    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
!    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)
    allocate( zp(nz) )
    zp(:) = -h_scale*log(p(:)/1000.)
  end if

  call timeinfo(ncid, nt_in)
  nt = nt_in / t_interval

  ! allocate var.s
  if (imn == 1 .and. iyr == 1) then
    allocate( cosphi(ny), f(ny) )
    f(:) = 2.*ome_earth*sin(lat(:)*deg2rad)
    cosphi(:) = cos(lat(:)*deg2rad)
  end if

  allocate( u(nx,ny,nz,nt), v(nx,ny,nz,nt), wp(nx,ny,nz,nt) )
  allocate( te(nx,ny,nz,nt) )

  allocate( um(ny,nz,nt), vres(ny,nz,nt), wres(ny,nz,nt) )
  allocate( tm(ny,nz,nt) )

  allocate( cor(ny,nz,nt), advy(ny,nz,nt), advz(ny,nz,nt) )
  allocate( epd(ny,nz,nt) )

  allocate( vtm(ny,nz,nt), vum(ny,nz,nt), vwm(ny,nz,nt) )

  allocate( term1(ny,nz,nt), t1(ny,nz,nt), stabz(ny,nz,nt) )
  allocate( temp(ny,nz,nt) )
  allocate( t0(nz) )


  ! get var.s
  allocate( u_in(nx,ny,nz,nt_in), v_in(nx,ny,nz,nt_in), wp_in(nx,ny,nz,nt_in) )
  allocate( te_in(nx,ny,nz,nt_in) )
  call get4d(ncid,trim(ivarname(1)),nx,ny,nz,nt_in, u_in )
  call get4d(ncid,trim(ivarname(2)),nx,ny,nz,nt_in, v_in )
  call get4d(ncid,trim(ivarname(3)),nx,ny,nz,nt_in, wp_in)
  call get4d(ncid,trim(ivarname(4)),nx,ny,nz,nt_in, te_in)
  do n=1, nt
    n_in2 = n*t_interval
    n_in1 = n_in2 - t_interval + 1
    call tempo_avg(nx,ny,nz,t_interval,missv,u_in (:,:,:,n_in1:n_in2), u (:,:,:,n))
    call tempo_avg(nx,ny,nz,t_interval,missv,v_in (:,:,:,n_in1:n_in2), v (:,:,:,n))
    call tempo_avg(nx,ny,nz,t_interval,missv,wp_in(:,:,:,n_in1:n_in2), wp(:,:,:,n))
    call tempo_avg(nx,ny,nz,t_interval,missv,te_in(:,:,:,n_in1:n_in2), te(:,:,:,n))
  enddo
  deallocate( u_in, v_in, wp_in, te_in )

  ! mean and perturbation
  call zonal_avg(nx,ny,nz,nt,missv,u , um  )
  call zonal_avg(nx,ny,nz,nt,missv,v , vres)
  call zonal_avg(nx,ny,nz,nt,missv,wp, wres)
  call zonal_avg(nx,ny,nz,nt,missv,te, tm  )
  do n=1, nt
  do k=1, nz
  do j=1, ny
    u (:,j,k,n) = u (:,j,k,n) - um  (j,k,n)
    v (:,j,k,n) = v (:,j,k,n) - vres(j,k,n)
    wp(:,j,k,n) = wp(:,j,k,n) - wres(j,k,n)
    te(:,j,k,n) = te(:,j,k,n) - tm  (j,k,n)
  enddo
  enddo
  enddo

  ! t0
  call tempo_avg(1,ny,nz,nt,missv,tm, t1(:,:,1))
  call merid_avg(1,ny,nz,1,missv,lat,t1(:,:,1), t0)
  do n=1, nt
  do k=1, nz
    t1(:,k,n) = tm(:,k,n) - t0(k)
  enddo
  enddo

  ! flux terms
  call corr_zonal_avg(nx,ny,nz,nt,missv,v,te, vtm)
  call corr_zonal_avg(nx,ny,nz,nt,missv,v,u , vum)
  call corr_zonal_avg(nx,ny,nz,nt,missv,v,wp, vwm)

  ! residual mean meridional circulations
  call gradz_2nd(1,ny,nz,nt,tm,zp,missv, stabz)
  do n=1, nt
  do k=1, nz
  do j=1, ny
    if ( stabz(j,k,n) /= missv .and. t0(k) /= missv ) then
      stabz(j,k,n) = stabz(j,k,n) + kappa*t0(k)/h_scale
    else
      stabz(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  do n=1, nt
  do k=1, nz
  do j=1, ny
    if ( vtm(j,k,n) /= missv .and. stabz(j,k,n) /= missv ) then
      term1(j,k,n) = vtm(j,k,n)/stabz(j,k,n)
    else
      term1(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  call gradz_2nd(1,ny,nz,nt,term1,zp,missv, temp)
  do n=1, nt
  do k=1, nz
  do j=1, ny
    if ( vres(j,k,n) /= missv .and. temp(j,k,n) /= missv .and. &
         term1(j,k,n) /= missv ) then
      vres(j,k,n) = vres(j,k,n) - (temp(j,k,n)-term1(j,k,n)/h_scale)
    else
      vres(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  do n=1, nt
  do k=1, nz
  do j=1, ny
    if (term1(j,k,n) /= missv) then
      term1(j,k,n) = term1(j,k,n)*cosphi(j)
    end if
  enddo
  enddo
  enddo

  call grady_2nd(1,ny,nz,nt,term1,lat,missv, temp)
  do n=1, nt
  do k=1, nz
  do j=1, ny
    if ( wres(j,k,n) /= missv .and. temp(j,k,n) /= missv .and. &
         cosphi(j) /= 0. ) then
      wres(j,k,n) = wres(j,k,n) + temp(j,k,n)/cosphi(j)
    else
      wres(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  ! coriolis term
  do n=1, nt
  do k=1, nz
  do j=1, ny
    if (vres(j,k,n) /= missv) then
      cor(j,k,n) = vres(j,k,n)*f(j) * 86400.
    else
      cor(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  ! advy
  do n=1, nt
  do k=1, nz
  do j=1, ny
    if (um(j,k,n) /= missv) then
      temp(j,k,n) = um(j,k,n)*cosphi(j)
    else
      temp(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  call grady_2nd(1,ny,nz,nt,temp,lat,missv, advy)
  do n=1, nt
  do k=1, nz
  do j=1, ny
    if ( vres(j,k,n) /= missv .and. advy(j,k,n) /= missv .and. &
         cosphi(j) /= 0. ) then
      advy(j,k,n) = -vres(j,k,n)*advy(j,k,n)/cosphi(j) * 86400.
    else
      advy(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  ! advz
  call gradz_2nd(1,ny,nz,nt,um,zp,missv, temp)
  do n=1, nt
  do k=1, nz
  do j=1, ny
    if ( wres(j,k,n) /= missv .and. temp(j,k,n) /= missv ) then
      advz(j,k,n) = -wres(j,k,n)*temp(j,k,n) * 86400.
    else
      advz(j,k,n) = missv
    end if
  enddo
  enddo
  enddo

  ! epd
  epd(:,:,:) = 0.



  deallocate( um, tm )
  deallocate( u, v, wp, te )
  deallocate( vtm, vum, vwm )
  deallocate( term1, t1, stabz, temp )
  deallocate( t0 )


  nxr =  1   ;   nxa =  1
  nyr = ny   ;   nya = ny
  nzr = nz   ;   nza = nz  ! ;   iz = 28
  ntr = nt   ;   nta =  1

  allocate( set(1)%var_out(nxa,nya,nza,nta) )
  call tempo_avg(nxr,nyr,nzr,ntr,missv,vres, set(1)%var_out)
  deallocate( vres )

  allocate( set(2)%var_out(nxa,nya,nza,nta) )
  call tempo_avg(nxr,nyr,nzr,ntr,missv,wres, set(2)%var_out)
  deallocate( wres )

  allocate( set(4)%var_out(nxa,nya,nza,nta) )
  call tempo_avg(nxr,nyr,nzr,ntr,missv,cor, set(4)%var_out)
  deallocate( cor )

  allocate( set(5)%var_out(nxa,nya,nza,nta) )
  call tempo_avg(nxr,nyr,nzr,ntr,missv,advy, set(5)%var_out)
  deallocate( advy )

  allocate( set(6)%var_out(nxa,nya,nza,nta) )
  call tempo_avg(nxr,nyr,nzr,ntr,missv,advz, set(6)%var_out)
  deallocate( advz )

  allocate( set(7)%var_out(nxa,nya,nza,nta) )
  call tempo_avg(nxr,nyr,nzr,ntr,missv,epd, set(7)%var_out)
  deallocate( epd )

  allocate( set(3)%var_out(nxa,nya,nza,nta) )
  do n=1, nta
  do k=1, nza
  do j=1, nya
  do i=1, nxa
    if ( set(4)%var_out(i,j,k,n) /= missv .and. &
         set(5)%var_out(i,j,k,n) /= missv .and. &
         set(6)%var_out(i,j,k,n) /= missv .and. &
         set(7)%var_out(i,j,k,n) /= missv ) then
      set(3)%var_out(i,j,k,n) = set(4)%var_out(i,j,k,n) + &
                                set(5)%var_out(i,j,k,n) + &
                                set(6)%var_out(i,j,k,n) + &
                                set(7)%var_out(i,j,k,n)
    else
      set(3)%var_out(i,j,k,n) = missv
    end if
  enddo
  enddo
  enddo
  enddo


  if (imn == 1 .and. iyr == 1) then
    do iv=1, nvar
      set(iv)%nd(1) = nxa
      set(iv)%nd(2) = nya
      set(iv)%nd(3) = nza
      set(iv)%axis = (/'     ','lat ','p ','t'/)
      set(iv)%vname = ovarname(iv)
      allocate( set(iv)%axis1(set(iv)%nd(1)) )
      allocate( set(iv)%axis2(set(iv)%nd(2)) )
      allocate( set(iv)%axis3(set(iv)%nd(3)) )
      set(iv)%axis1 = -999.
      set(iv)%axis2 = lat
      set(iv)%axis3 = p
    enddo
  end if
  do iv=1, nvar
    set(iv)%nd(4) = nta
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis4 = (year-2000)*100.+month
  enddo
  call closenc(ncid)

  ! dump
  fname = trim(ofdir)//trim(expname)//'/anal/'// &
          trim(expname)//'.'//trim(outname)//'.'//cyear//'.'//cmonth//'.nc'
  call outnc(trim(fname),nvar,set,'')

  do iv=1, nvar
    deallocate( set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo


  ENDDO  N_YEAR
  ENDDO  N_MON


  STOP

END program


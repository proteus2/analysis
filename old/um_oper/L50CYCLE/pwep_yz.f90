PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis
  use fft

  implicit none

  integer, parameter ::  nmonth =  1, mstart =  1
  integer, parameter ::  nyear  =  7, ystart = 2009
  integer, parameter ::  t_interval = 4
  real,    parameter ::  missv = 1.e32, zbottom = 0.e3
  real,    parameter ::  zh_int = 12.0
  integer, parameter ::  nwave = 3
  integer, parameter ::  nivar = 5, nvar = nwave*3
  character(len=128) ::  ifdir, expname, vartype, ofdir, outname
  character(len=64)  ::  ivarname(nivar), ovarname(nvar)

  ! files
  data                   ifdir   /'/data5/kyh/umres/'/
  data                   expname /'ctl'/
  data                   vartype /'std'/
  data                   ofdir   /'/data5/kyh/umres/'/ 
  data                   outname /'pwep_yz'/
  ! input variables
  data                   ivarname(1) /'u'/
  data                   ivarname(2) /'v'/
  data                   ivarname(3) /'dz_dt'/
  data                   ivarname(4) /'theta_1'/
  data                   ivarname(5) /'unspecified'/
  ! output variables
  data                   ovarname(1) /'epd1'/
  data                   ovarname(2) /'f1_y'/
  data                   ovarname(3) /'f1_z'/
  data                   ovarname(4) /'epd2'/
  data                   ovarname(5) /'f2_y'/
  data                   ovarname(6) /'f2_z'/
  data                   ovarname(7) /'epd3'/
  data                   ovarname(8) /'f3_y'/
  data                   ovarname(9) /'f3_z'/


  real, dimension(:,:,:,:), allocatable ::  tem_in
  real, dimension(:,:,:),   allocatable ::  u_in, v_in, w_in, pt_in
  real, dimension(:,:,:),   allocatable ::  u, v, w, pt, rho
  real, dimension(:,:,:),   allocatable ::  uprm, vprm, wprm, ptprm
  real, dimension(:,:),     allocatable ::  um, vm, wm, ptm, rho0
  real, dimension(:,:),     allocatable ::  epd, fy, fz
  real, dimension(:),       allocatable ::  lon, lat, zh, lonu, latv, zht

  real, dimension(:,:),     allocatable ::  phi, vflx, wflx
  real, dimension(:,:),     allocatable ::  grady, gradz, temp
  real, dimension(:),       allocatable ::  cosphi, f

  integer, dimension(:,:), allocatable ::  nvm, nwm
  integer, dimension(:,:), allocatable ::  nepd, nfy, nfz

  double complex, dimension(:), allocatable ::  coef, coef2

  integer ::  year, month
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, imn, iyr, nt_in
  integer ::  nxr, nyr, nzr, ntr, nxa, nya, nza, nta
  integer ::  iz
  integer ::  i,j,k,n,nn,iv, ncid, n_in
  integer ::  iwave, nnn
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

  fname = trim(ifdir)//trim(expname)//'/'//cyear//'/'// &
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
  call timeinfo(ncid, nt_in)


  ! allocate var.s
  if (imn == 1 .and. iyr == 1) then
    allocate( cosphi(ny), f(ny) )
    f(:) = 2.*ome_earth*sin(lat(:)*deg2rad)
    cosphi(:) = cos(lat(:)*deg2rad)
  end if

  nxr =  1      ;   nxa =  1
  nyr = ny      ;   nya = ny
  nzr = nz      ;   nza = nz  ! ;   iz = 28
  ntr = nt_in   ;   nta =  1
  do iv=1, nvar
    allocate( set(iv)%var_out(nxa,nya,nza,nta) )
    set(iv)%var_out(:,:,:,:) = 0.
  enddo

  nt = ntr / t_interval / nta

  allocate( nepd(ny,nz) )
  allocate( nfy(ny,nz), nfz(ny,nz) )


  RTIME:  DO nn=1, nta


  nepd (:,:) = 0
  nfy  (:,:) = 0
  nfz  (:,:) = 0


  ATIME:  DO n=1, nt


  allocate( um(ny,nz), vm(ny,nz), wm(ny,nz) )
  allocate( ptm(ny,nz) )
  allocate( rho0(ny,nz) )

  allocate( epd(ny,nz) )
  allocate( fy(ny,nz), fz(ny,nz) )

  allocate( phi(ny,nz) )
  allocate( vflx(ny,nz), wflx(ny,nz) )

  allocate( grady(ny,nz), gradz(ny,nz) )
  allocate( temp(ny,nz) )

  allocate( u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz) )
  allocate( pt(nx,ny,nz) )
  allocate( rho(nx,ny,nz) )

  allocate( uprm(nx,ny,nz), vprm(nx,ny,nz), wprm(nx,ny,nz) )
  allocate( ptprm(nx,ny,nz) )
  allocate( coef(nx), coef2(nx) )

  ! get var.s
  n_in = (nn-1)*nt + (n-1)*t_interval + 1

  allocate( u_in(nxu,ny,nz), v_in(nx,nyv,nz), w_in(nx,ny,nzt) )
  allocate( pt_in(nx,ny,nzt) )

  allocate( tem_in(nxu,ny,nz,t_interval) )
  call geta4d(ncid,trim(ivarname(1)),1,nxu,1,ny,1,nz,n_in,t_interval, tem_in )
  call tempo_avg(nxu,ny,nz,t_interval,missv,tem_in, u_in)
  call u2rho(nx,ny,nz,1,u_in,missv, u)
  deallocate( tem_in )

  allocate( tem_in(nx,nyv,nz,t_interval) )
  call geta4d(ncid,trim(ivarname(2)),1,nx,1,nyv,1,nz,n_in,t_interval, tem_in )
  call tempo_avg(nx,nyv,nz,t_interval,missv,tem_in, v_in)
  call v2rho(nx,ny,nz,1,v_in,missv, v)
  deallocate( tem_in )

  allocate( tem_in(nx,ny,nzt,t_interval) )
  call geta4d(ncid,trim(ivarname(3)),1,nx,1,ny,1,nzt,n_in,t_interval, tem_in )
  call tempo_avg(nx,ny,nzt,t_interval,missv,tem_in, w_in)
  call t2rho(nx,ny,nz,1,w_in,missv, w)
  deallocate( tem_in )

  allocate( tem_in(nx,ny,nzt,t_interval) )
  call geta4d(ncid,trim(ivarname(4)),1,nx,1,ny,1,nzt,n_in,t_interval, tem_in )
  call tempo_avg(nx,ny,nzt,t_interval,missv,tem_in, pt_in)
  pt_in(:,:,:) = log(pt_in(:,:,:))
  call t2rho(nx,ny,nz,1,pt_in,missv, pt)
  pt(:,:,:) = exp(pt(:,:,:))
  deallocate( tem_in )

  allocate( tem_in(nx,ny,nz,t_interval) )
  call geta4d(ncid,trim(ivarname(5)),1,nx,1,ny,1,nz,n_in,t_interval, tem_in )
  call tempo_avg(nx,ny,nz,t_interval,missv,tem_in, rho)
  deallocate( tem_in )

  deallocate( u_in, v_in, w_in, pt_in )

  ! mean and perturbation
  call zonal_avg(nx,ny,nz,1,missv,u  ,um  )
  call zonal_avg(nx,ny,nz,1,missv,v  ,vm  )
  call zonal_avg(nx,ny,nz,1,missv,w  ,wm  )
  call zonal_avg(nx,ny,nz,1,missv,pt ,ptm )
  call zonal_avg(nx,ny,nz,1,missv,rho,rho0)
  deallocate( rho )

  do k=1, nz
  do j=1, ny
    u(:,j,k) = u(:,j,k) - um(j,k)
    v(:,j,k) = v(:,j,k) - vm(j,k)
  enddo
  enddo
  do k=2, nz
  do j=1, ny
    w (:,j,k) = w (:,j,k) - wm (j,k)
    pt(:,j,k) = pt(:,j,k) - ptm(j,k)
  enddo
  enddo
  w (:,:,1) = missv
  pt(:,:,1) = missv


  WAVENO:  DO iwave=1, nwave


  coef2 = (0.0,0.0)
  do k=1, nz
  do j=1, ny
    call fft1df(nx,u(:,j,k),coef)
    coef2(iwave+1   ) = coef(iwave+1   )
    coef2(nx-iwave+1) = coef(nx-iwave+1)
    call fft1db(nx,coef2,uprm(:,j,k))
  enddo
  enddo
  do k=1, nz
  do j=1, ny
    call fft1df(nx,v(:,j,k),coef)
    coef2(iwave+1   ) = coef(iwave+1   )
    coef2(nx-iwave+1) = coef(nx-iwave+1)
    call fft1db(nx,coef2,vprm(:,j,k))
  enddo
  enddo
  do k=2, nz
  do j=1, ny
    call fft1df(nx,w(:,j,k),coef)
    coef2(iwave+1   ) = coef(iwave+1   )
    coef2(nx-iwave+1) = coef(nx-iwave+1)
    call fft1db(nx,coef2,wprm(:,j,k))
  enddo
  enddo
  wprm(:,:,1) = missv
  do k=2, nz
  do j=1, ny
    call fft1df(nx,pt(:,j,k),coef)
    coef2(iwave+1   ) = coef(iwave+1   )
    coef2(nx-iwave+1) = coef(nx-iwave+1)
    call fft1db(nx,coef2,ptprm(:,j,k))
  enddo
  enddo
  ptprm(:,:,1) = missv

  if (iwave == nwave)  deallocate( coef, coef2 )


  ! phi
  call corr_zonal_avg(nx,ny,nz,1,missv,vprm,ptprm, vflx)
  call corr_zonal_avg(nx,ny,nz,1,missv,wprm,ptprm, wflx)
  if (iwave == nwave)  deallocate( ptprm )
  if (iwave == nwave)  deallocate( pt )

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
    temp(j,k) = temp(j,k)*cosphi(j)
  enddo
  enddo
  call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)

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

  call gradz_2nd(1,ny,nz,1,um,zh,missv, gradz)

  ! epf, epd
  call corr_zonal_avg(nx,ny,nz,1,missv,vprm,uprm, vflx)
  call corr_zonal_avg(nx,ny,nz,1,missv,wprm,uprm, wflx)
  if (iwave == nwave)  deallocate( uprm, vprm, wprm )
  if (iwave == nwave)  deallocate( u, v, w )

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

  if (iwave == nwave)  deallocate( grady, gradz )

  call div_yz(ny,nz,1,fy,fz,lat,zh,missv, epd)
  do k=2, nz
  do j=2, ny-1
    epd(j,k) = epd(j,k)/(rho0(j,k)*r_earth*cosphi(j)) * 86400.
  enddo
  enddo
  epd(:,1) = missv


  if (iwave == nwave) then
    deallocate( temp )
    deallocate( vflx, wflx )
    deallocate( phi )
    deallocate( rho0 )
    deallocate( ptm )
    deallocate( um )
  end if


  nnn = (iwave-1)*3
  do k=1, nza
  do j=1, nya
  do i=1, nxa
    if (epd(j,k) /= missv) then
      nepd(j,k) = nepd(j,k) + 1
      set(nnn+1)%var_out(i,j,k,nn) = set(nnn+1)%var_out(i,j,k,nn) + epd(j,k)
    end if
    if (fy(j,k) /= missv) then
      nfy(j,k) = nfy(j,k) + 1
      set(nnn+2)%var_out(i,j,k,nn) = set(nnn+2)%var_out(i,j,k,nn) + fy(j,k)
    end if
    if (fz(j,k) /= missv) then
      nfz(j,k) = nfz(j,k) + 1
      set(nnn+3)%var_out(i,j,k,nn) = set(nnn+3)%var_out(i,j,k,nn) + fz(j,k)
    end if
  enddo
  enddo
  enddo

  if (iwave == nwave) then
    deallocate( vm, wm )
    deallocate( epd )
    deallocate( fy, fz )
  end if


  ENDDO  WAVENO


  ENDDO  ATIME


  do k=1, nza
  do j=1, nya
  do i=1, nxa
  do iwave=1, nwave
    nnn = (iwave-1)*3
    if (nepd(j,k) /= 0) then
      set(nnn+1)%var_out(i,j,k,nn) = set(nnn+1)%var_out(i,j,k,nn)/nepd(j,k)
    else
      set(nnn+1)%var_out(i,j,k,nn) = missv
    end if
    if (nfy(j,k) /= 0) then
      set(nnn+2)%var_out(i,j,k,nn) = set(nnn+2)%var_out(i,j,k,nn)/nfy(j,k)
    else
      set(nnn+2)%var_out(i,j,k,nn) = missv
    end if
    if (nfz(j,k) /= 0) then
      set(nnn+3)%var_out(i,j,k,nn) = set(nnn+3)%var_out(i,j,k,nn)/nfz(j,k)
    else
      set(nnn+3)%var_out(i,j,k,nn) = missv
    end if
  enddo
  enddo
  enddo
  enddo


  ENDDO  RTIME


  if (zbottom /= 0.) then
    iz = nza + 1
    do k=1, nza
      if (zh(k) >= zbottom) then
        iz = k
        EXIT
      end if
    enddo
    do iv=1, nvar
      set(iv)%var_out(:,:,:iz-1,:) = missv
    enddo
  end if


  if (imn == 1 .and. iyr == 1) then
    do iv=1, nvar
      set(iv)%nd(1) = nxa
      set(iv)%nd(2) = nya
      set(iv)%nd(3) = nza
      set(iv)%axis = (/'     ','lat ','zh ','t'/)
      set(iv)%vname = ovarname(iv)
      allocate( set(iv)%axis1(set(iv)%nd(1)) )
      allocate( set(iv)%axis2(set(iv)%nd(2)) )
      allocate( set(iv)%axis3(set(iv)%nd(3)) )
      set(iv)%axis1 = -999.
      set(iv)%axis2 = lat
      set(iv)%axis3 = zh
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


MODULE raytrace

  use netcdfio

  public :: ray

  contains


subroutine ray(fdir,stid,no,subno,contopt,time,intopt,fini,kini,lini,lonini,latini,zini,delz,tlimit)

!=======================================================================
!
! termination condition : title of output file (NETCDF)
!
!   0, initially m2<0
!   1, reflection level
!   2, critical level
!   3, time limitation
!   9, out of domain
!
! There are no validation of (1) WKB approx. 
!   or (2) critical level interaction.
!
!   (1) When m^2 is not so small, calculation is reliable.
!   When m^2 is getting smaller enough to break the WKB assumption, 
!   the ray moves very small distance (in any direction),
!   and finally m^2 < 0 "in most cases".
!
!   (2) Must set some critical range in your results.
!   ( default : Lz < 100 m )
!
!=======================================================================

  use interpolate

  implicit none

!-----------------------------------------------------------------------
  integer,          intent(in) ::  stid
  integer,          intent(in) ::  no, subno, contopt, intopt
  real,             intent(in) ::  kini, lini, fini    ! [rad/m], [rad/s]
  real,             intent(in) ::  lonini, latini      ! [deg]
  real,             intent(in) ::  zini, delz, tlimit  ! [m]
  character(len=6), intent(in) ::  time
!-----------------------------------------------------------------------
  real,    parameter ::  g0=9.80665, Rd=287., cp=1004., cv=707.
  integer            ::  i,j,k,n, ncid
  integer            ::  nlon, nlat, nz
  integer            ::  nt, ntdn, ntup
  logical            ::  iex
  character(len=128) ::  term, termdn, termup
  character(len=255) ::  ifname, wfname, fdir
  real, dimension(:),     allocatable ::  lon, lat, z
  real, dimension(:,:),   allocatable ::  vars, varsdn, varsup
  real, dimension(:,:,:), allocatable ::  t0,u0,v0,n2,al2

!-----------------------------------------------------------------------
!  Read netcdf-type data (must be sorted in order, by regular interval)
!-----------------------------------------------------------------------

! initial temperature, wind data (background)
  write(ifname,'(a)') trim(fdir)//trim(time)//'.nc' 

  call opennc(ifname,ncid)
  call dilen(ncid,'lon',nlon)
  call dilen(ncid,'lat',nlat)
  call dilen(ncid,'Z',nz)

  allocate(lon(1:nlon))             ;  lon(:)=0.0
  allocate(lat(1:nlat))             ;  lat(:)=0.0
  allocate(z(1:nz))                 ;  z(:)=0.0
  allocate(t0(1:nlon,1:nlat,1:nz))  ;  t0(:,:,:)=0.0

  call get1d(ncid,'lon',nlon,lon)
  call get1d(ncid,'lat',nlat,lat)
  call get1d(ncid,'Z',nz, z)
  call get3d(ncid,'T',nlon,nlat,nz,t0)

!  z(:) = z(:) * 1.e3

  allocate(u0(1:nlon,1:nlat,1:nz))  ;  u0(:,:,:)=0.0
  allocate(v0(1:nlon,1:nlat,1:nz))  ;  v0(:,:,:)=0.0

  call get3d(ncid,'U',nlon,nlat,nz,u0)
  call get3d(ncid,'V',nlon,nlat,nz,v0)
  call closenc(ncid)

! obtain Brunt-Vaisala frequency, al2=(1/2H)**2
  allocate(n2(1:nlon,1:nlat,1:nz))   ;   n2(:,:,:)=0.0
  do j=1, nlat
  do i=1, nlon
    do k=2, nz-1
      n2(i,j,k) = g0/t0(i,j,k)*((t0(i,j,k+1)-t0(i,j,k-1))/(z(k+1)-z(k-1))+g0/cp)
    enddo
    n2(i,j,1) = 2.*g0/t0(i,j,1)*((t0(i,j,2)-t0(i,j,1))/(z(2)-z(1))+g0/cp)-n2(i,j,2)
    n2(i,j,nz) = 2.*g0/t0(i,j,nz)*((t0(i,j,nz)-t0(i,j,nz-1))/(z(nz)-z(nz-1))+g0/cp) &
                 - n2(i,j,nz-1)
  enddo
  enddo 

  allocate(al2(1:nlon,1:nlat,1:nz))   ;   al2(:,:,:)=0.0
  do k=1, nz
  do j=1, nlat
  do i=1, nlon
    al2(i,j,k) = ((g0/((cp/cv)*Rd*t0(i,j,k)) + n2(i,j,k)/g0)/2.)**2
  enddo
  enddo
  enddo

!  call out3d('brnal.nc',2,(/'N2','A2'/),(/n2,al2/), &
!                       'lon',nlon,lon,'lat',nlat,lat,'z',nz,z,'')

!-----------------------------------------------------------------------
! Output file and var.s
!-----------------------------------------------------------------------

  write(wfname,'(i5.5,a,i3.3,a,i1,a)') stid,'/', no,'.',subno,'.'//trim(time)//'.nc'

  nt = 10000
  allocate(vars(1:nt,1:18))   ;  vars(:,:) = 1.e+32

!-----------------------------------------------------------------------
! Ray-tracing calculation
!-----------------------------------------------------------------------

!  inquire(file='output/'//trim(wfname), exist=iex)

  print*, 'initial properties -------------------------',no,subno
  print*, 'lonr', '   latr', '   zr'
  print*, lonini, latini, zini
  print*, 'input freq.', '  k   ', '   l'
  print*, fini, kini, lini
  print*, '--------------------------------------------'

  ! downward direction
  vars(1,12) = 0.        ! time
  vars(1,1)  = kini      ! k
  vars(1,2)  = lini      ! l
  vars(1,5)  = fini      ! freq. or intrinsic freq.
  vars(1,9)  = lonini    ! lon.
  vars(1,10) = latini    ! lat.
  vars(1,11) = zini      ! z
  vars(nt,3)   = 0.      ! mwn_last

  ntdn = nt
  call ray_integ(-1,delz,tlimit,intopt,                  &
                 nlon,nlat,nz,lon,lat,z,t0,u0,v0,n2,al2, &
                 ntdn,vars,termdn)

  if (ntdn+1 .eq. nt)  print*, 'ERROR!! Set nt a large number.'

  allocate(varsdn(1:ntdn,1:18))
  do i=1, 18
  do n=1, ntdn
    varsdn(n,i) = vars(n,i)
  enddo
  enddo

  ! upward direction
  vars(1,12) = 0.        ! time
  vars(1,1)  = kini      ! k
  vars(1,2)  = lini      ! l
  vars(1,5)  = fini      ! freq. or intrinsic freq.
  vars(1,9)  = lonini    ! lon.
  vars(1,10) = latini    ! lat.
  vars(1,11) = zini      ! z
  vars(nt,3)   = varsdn(2,3)    ! mwn_last

  ntup = nt
  call ray_integ(1,delz,tlimit,intopt,                   &
                 nlon,nlat,nz,lon,lat,z,t0,u0,v0,n2,al2, &
                 ntup,vars,termup)

  if (ntup+1 .eq. nt)  print*, 'ERROR!! Set nt a large number.'

  allocate(varsup(1:ntup,1:18))
  do i=1, 18
  do n=1, ntup
    varsup(n,i) = vars(n,i)
  enddo
  enddo

  ! combine the results
  deallocate(vars)

  nt = ntup + ntdn - 1
  allocate(vars(1:nt,1:18))

  do i=1, 18
    do n=1, ntdn
      vars(n,i) = varsdn(ntdn-n+1,i)
    enddo
    do n=1, ntup
      vars(ntdn+n-1,i) = varsup(n,i)
    enddo
  enddo

  write(term,'(a,f5.1,a)') 'z0: ',zini/1.e3,'km // dn: '//trim(termdn)

!-----------------------------------------------------------------------
!  Dump
!-----------------------------------------------------------------------

  ! for some machines...
  do i=1, 18
  do n=1, nt
    if (abs(vars(n,i)) .lt. 1.e-30)  vars(n,i) = 0.
  enddo
  enddo

  call out1d(trim(wfname),18,            &
     (/'k    ',    'l',    'm',   'mi','ifreq',  'cgx',  'cgy',  'cgz',     &
         'lon',  'lat',  'hgt', 'time',   'u0',   'v0',   'N2', 'Hrho',     &
          'Ri','R_WKB'/),vars(1:nt,1:18),                                   &
       'z',nt,vars(1:nt,11),trim(term))

  deallocate(vars)


  return

end subroutine ray


subroutine ray_integ(updn,delz,tlimit,intopt,            &
                 nlon,nlat,nz,lon,lat,z,t0,u0,v0,n2,al2, &
                 nt,vars,term)

  use interpolate

  implicit none

  integer, intent(in) ::  updn, intopt
  integer, intent(in) ::  nlon, nlat, nz
  real,    intent(in) ::  delz, tlimit
  real,    intent(in) ::  lon(nlon), lat(nlat), z(nz)
  real,    intent(in), dimension(nlon,nlat,nz) ::  t0, u0, v0, n2, al2
  integer, intent(inout) ::  nt
  real,    intent(inout) ::  vars(nt,18)
  character(len=128), intent(out) ::  term

  real, parameter ::  dtlimit = 3600.
  real, parameter ::  pi=3.141593
  real, parameter ::  g0=9.80665, r_ear=6.37e6, ome_ear=7.292e-5
  integer            ::  i,j,k, n_iter, count, endcond
  integer            ::  m2sign
  real               ::  dlon, dlat, dz, delzs

  integer ::  indx,indy,indz,indxp1,indyp1,indzp1
  real*8  ::  lonr, latr, zr
  real*8  ::  dt, tsum
  real*8  ::  kwn, lwn, mwn, mwni, ifreq, ifreq_last, f0, beta, rwkb
  real*8  ::  sig2, freq, mwn_dis, mwn_geq
  real*8  ::  dx8, dy8, dz8, deltax, deltay, deltaz
  real*8  ::  cgx, cgy, cgz, delk, dell, delm
  real*8  ::  uu, vv, nn, aa, hrho, ri
  real*8  ::  dudx, dudy, dudz, dvdx, dvdy, dvdz
  real*8  ::  dn2dx, dn2dy, dn2dz, da2dx, da2dy, da2dz
  real*8  ::  mwn_last, mwn_last2, delm1, delm2


  delzs = sign(delz,real(updn))

  dlon = lon(2) - lon(1)
  dlat = lat(2) - lat(1)
  dz   = z(2) - z(1)

  dy8 = r_ear*dlat*pi/180.d0
  dz8 = dble(dz)


  print *, 'RAY-INTEGRATION STARTS.'

  tsum = dble(vars(1,12))
  lonr = dble(vars(1,9))
  latr = dble(vars(1,10))
  zr   = dble(vars(1,11))
  freq = dble(vars(1,5))
  kwn  = dble(vars(1,1))
  lwn  = dble(vars(1,2))

  mwn_last   = dble(vars(nt,3))
  mwn_last2  = 0.d0

  vars(:,:) = 1.e+32

  count = 0

  do n_iter=1, nt

  count = count + 1
  print *

  if (latr.gt.lat(1) .or. latr.lt.lat(nlat)        &
      .or. zr.lt.z(1) .or. zr.gt.z(nz)) then
    print *, 'RAY IS OUT OF THE DOMAIN.'
    endcond = 9
    EXIT
  end if

  if (lonr .gt. 360.d0) then
    lonr = lonr - 360.d0
  else if (lonr .lt. 0.d0) then
    lonr = lonr + 360.d0
  end if

  indx = int((lonr-lon(1))/dlon)+1
  indy = int((latr-lat(1))/dlat)+1
  indz = int((zr-z(1))/dz)+1
  deltax = r_ear*cos(latr*pi/180.)*(lonr-lon(indx))*pi/180.0
  deltay = r_ear*(latr-lat(indy))*pi/180.0
  deltaz = zr-z(indz)

  indxp1 = indx + 1
  indyp1 = indy + 1
  indzp1 = indz + 1
  if (indx .eq. nlon)  indxp1 = 1
  if (deltay .eq. 0.d0)  indyp1 = nlat
  if (deltaz .eq. 0.d0)  indzp1 = nz

  dx8 = r_ear*cos(latr*pi/180.0)*dlon*pi/180.0

  ! calculate u,v,N2,alpha,f0 at ray's position
  call interpol(u0(indx,indy,indz),      u0(indxp1,indy,indz),      &
                u0(indx,indyp1,indz),    u0(indxp1,indyp1,indz),    &
                u0(indx,indy,indzp1),    u0(indxp1,indy,indzp1),    &
                u0(indx,indyp1,indzp1),  u0(indxp1,indyp1,indzp1),  &
                dx8, dy8, dz8, deltax, deltay, deltaz, uu)

  call interpol(v0(indx,indy,indz),      v0(indxp1,indy,indz),      &
                v0(indx,indyp1,indz),    v0(indxp1,indyp1,indz),    &
                v0(indx,indy,indzp1),    v0(indxp1,indy,indzp1),    &
                v0(indx,indyp1,indzp1),  v0(indxp1,indyp1,indzp1),  &
                dx8, dy8, dz8, deltax, deltay, deltaz, vv)

  call interpol(n2(indx,indy,indz),      n2(indxp1,indy,indz),      &
                n2(indx,indyp1,indz),    n2(indxp1,indyp1,indz),    &
                n2(indx,indy,indzp1),    n2(indxp1,indy,indzp1),    &
                n2(indx,indyp1,indzp1),  n2(indxp1,indyp1,indzp1),  &
                dx8, dy8, dz8, deltax, deltay, deltaz, nn)

  call interpol(al2(indx,indy,indz),     al2(indxp1,indy,indz),     &
                al2(indx,indyp1,indz),   al2(indxp1,indyp1,indz),   &
                al2(indx,indy,indzp1),   al2(indxp1,indy,indzp1),   &
                al2(indx,indyp1,indzp1), al2(indxp1,indyp1,indzp1), &
                dx8, dy8, dz8, deltax, deltay, deltaz, aa)

  f0 = 2.d0 * ome_ear * sin(pi*latr/180.)

  hrho = dsqrt(0.25/aa)

  ! obtain 1st derivatives
  call grad_int(nlon,nlat,nz,dx8,dy8,dz8,u0,indx,indy,indz,  &
                deltax,deltay,deltaz,dudx,dudy,dudz)

  call grad_int(nlon,nlat,nz,dx8,dy8,dz8,v0,indx,indy,indz,  &
                deltax,deltay,deltaz,dvdx,dvdy,dvdz)

  call grad_int(nlon,nlat,nz,dx8,dy8,dz8,n2,indx,indy,indz,  &
                deltax,deltay,deltaz,dn2dx,dn2dy,dn2dz)

  call grad_int(nlon,nlat,nz,dx8,dy8,dz8,al2,indx,indy,indz, &
                deltax,deltay,deltaz,da2dx,da2dy,da2dz)

  ! Richardson number
  ri = nn / (dudz**2+dvdz**2)

!-----------------------------------------------------------------------
!  Calculate int. freq. and m2, and check the propagation condition
!-----------------------------------------------------------------------

  if (intopt.eq.1 .and. n_iter.eq.1) then
    ifreq = freq
    freq = ifreq + kwn * uu + lwn * vv
  end if

  ! calculate int. freq.
  ifreq = freq - kwn * uu - lwn * vv

  print*, '--------------------------------------------'
  print*, 'z      :', real(zr/1000.)
  print*, 'time   :', real(tsum/3600.)
  print*, '  f0     /    int. freq    /     N'
  print*, real(f0), real(ifreq), sqrt(real(nn))

  ! check C.L. (ifreq = f0)
  if (n_iter .eq. 1)  ifreq_last = ifreq
  if (abs(ifreq).le.abs(f0) .or. (ifreq*ifreq_last).le.0.d0) then
    print *, '------------------- WAVE MEET CRITICAL LEVEL'
    print *, 'u = ', real(uu)
    print *, 'v = ', real(vv)
    endcond = 2
    EXIT
  end if

  ! calculate m2, m
  mwn_dis = (kwn**2+lwn**2)*(nn-ifreq**2) / (ifreq**2-f0**2) - aa

  ! check the sign of m^2
  if (mwn_dis .gt. 0.d0) then

    m2sign = 1

    mwn_dis = sqrt(mwn_dis)
    if (ifreq .gt. 0.0)  mwn_dis = -mwn_dis    ! for upward propagating wave

!-----------------------------------------------------------------------
! (1) calculate m by the dispersion relation
    mwn = mwn_dis
!-----------------------------------------------------------------------
! (2) calculate m by the eq. of dm/dt
!   ! if (2) is chosen, EXIT below 'else' paragraph must be activated.
!   mwn = mwn_geq
!   if (n_iter .eq. 1)  mwn = mwn_dis
!   if (n_iter.ne.1 .and. (mwn*mwn_last).le.0.d0) then
!     print *, '------------------ REFLECTION LEVEL (m2 < 0)'
!     endcond = 1
!     print *, 'm = ', real(mwn_dis)
!     print *, 'u = ', real(uu)
!     print *, 'v = ', real(vv)
!     EXIT
!   end if
!-----------------------------------------------------------------------
    mwni = 0.

    print*, '  Lz'
    print*, real(2.*pi/mwn)

    ! check C.L. (Lz < 100 m)
    if (abs(2.*pi/mwn) .lt. 100.d0) then
      print *, '------------------- WAVE MEET CRITICAL LEVEL (Lz)'
      print *, 'u = ', real(uu)
      print *, 'v = ', real(vv)
      endcond = 2
      EXIT
    end if

  else

    m2sign = -1

    print *, '------------------ REFLECTION LEVEL (m2 < 0)'
    endcond = 1
    print *, 'm = ', real(mwn_dis)
    print *, 'u = ', real(uu)
    print *, 'v = ', real(vv)
    EXIT   ! if (2) is chosen, this must be activated.

    mwn  = 0.
    mwni = sqrt(abs(mwn_dis))

  end if


  if (m2sign .eq. 1) then

    ! parameter for WKB approximation
    if (mwn_last2 .ne. 0.d0) then
      delm1 = zr - dble(vars(n_iter-1,11))
      delm2 = dble(vars(n_iter-1,11) - vars(n_iter-2,11))

      rwkb = abs( (mwn-mwn_last2)/(delm2+delm1) ) / mwn_last**2
    end if

    if (n_iter .eq. 1)  mwn_geq = mwn_dis
    print*, '    m err.    /  R_WKB (at before step)'
    print*, real((mwn_geq-mwn_dis)/mwn_dis), real(rwkb)


    sig2 = kwn**2 + lwn**2 + mwn**2 + aa
    beta = 2.d0*ome_ear*cos(latr*pi/180.0)/r_ear

    cgx = uu + kwn*(nn-ifreq**2) / (ifreq*sig2)
    cgy = vv + lwn*(nn-ifreq**2) / (ifreq*sig2)
    cgz = -mwn*(ifreq**2-f0**2) / (ifreq*sig2)

!-----------------------------------------------------------------------
!  Obtain cgx,cgy,cgz,dk,dl,dm from ray-tracing equation
!-----------------------------------------------------------------------
    dt = -100. 

    print*, 'dt','     delx', '     dely'
    print*, real(dt), real(cgx*dt), real(cgy*dt)
    print*, ' cgz  :',real(cgz)

    delk = -( (dn2dx*(kwn**2+lwn**2) - da2dx*(ifreq**2-f0**2)) /     &
           (2.*ifreq*sig2) + kwn*dudx + lwn*dvdx ) * dt

    dell = -( (dn2dy*(kwn**2+lwn**2) - da2dy*(ifreq**2-f0**2) +      &
           2.*f0*beta*(mwn**2+aa)) / (2.*ifreq*sig2)      &
           + kwn*dudy + lwn*dvdy ) * dt

    delm = -( (dn2dz*(kwn**2+lwn**2) - da2dz*(ifreq**2-f0**2)) /     &
           (2.*ifreq*sig2) + kwn*dudz + lwn*dvdz ) * dt

!-----------------------------------------------------------------------
!  output
!-----------------------------------------------------------------------

    vars(n_iter,1)  = real(kwn)    ! k
    vars(n_iter,2)  = real(lwn)    ! l
    vars(n_iter,3)  = real(mwn)    ! Re{m}
    vars(n_iter,4)  = real(mwni)   ! Im{m}
    vars(n_iter,5)  = real(ifreq)  ! intrinsic freq.
    vars(n_iter,6)  = real(cgx)    ! cgx
    vars(n_iter,7)  = real(cgy)    ! cgy
    vars(n_iter,8)  = real(cgz)    ! cgz
    vars(n_iter,9)  = real(lonr)   ! lon.
    vars(n_iter,10) = real(latr)   ! lat.
    vars(n_iter,11) = real(zr)     ! height
    vars(n_iter,12) = real(tsum)   ! time
    vars(n_iter,13) = real(uu)     ! U
    vars(n_iter,14) = real(vv)     ! V
    vars(n_iter,15) = real(nn)     ! N^2
    vars(n_iter,16) = real(hrho)   ! H (density scale height)
    vars(n_iter,17) = real(ri)     ! Richardson no.
    vars(n_iter,18) = real(rwkb)   ! WKB validation term
    if (n_iter .ne. 1)  vars(n_iter-1,18) = real(rwkb)

  else     ! m2sign = -1

    vars(n_iter,1)  = real(kwn)    ! k
    vars(n_iter,2)  = real(lwn)    ! l
    vars(n_iter,3)  = real(mwn)    ! Re{m}
    vars(n_iter,4)  = real(mwni)   ! Im{m}
    vars(n_iter,5)  = real(ifreq)  ! intrinsic freq.
    vars(n_iter,6)  = 1.e+32       ! cgx
    vars(n_iter,7)  = 1.e+32       ! cgy
    vars(n_iter,8)  = 0.0          ! cgz
    vars(n_iter,9)  = 1.e+32       ! lon.
    vars(n_iter,10) = 1.e+32       ! lat.
    vars(n_iter,11) = real(zr)     ! height
    vars(n_iter,12) = 1.e+32       ! time
    vars(n_iter,13) = real(uu)     ! U
    vars(n_iter,14) = real(vv)     ! V
    vars(n_iter,15) = real(nn)     ! N^2
    vars(n_iter,16) = real(hrho)   ! H (density scale height)
    vars(n_iter,17) = real(ri)     ! Richardson no.
    vars(n_iter,18) = 1.e+32       ! WKB validation term

  end if


!  if ( dabs(tsum) .ge. abs(tlimit)*3600. ) then    ! time limit
!    print*, '---------------------------------- TIME OUT'
!    print*, 'u = ', real(uu)
!    print*, 'v = ', real(vv)
!    count = count + 1
!    endcond = 3
!    EXIT
!  end if


!-----------------------------------------------------------------------
!  Update longitude, latitude, r_ear, height, and WN(k,l,m) rapidly
!-----------------------------------------------------------------------

  if (m2sign .eq. 1) then
    kwn     = kwn + delk
    lwn     = lwn + dell
    mwn_geq = mwn + delm

    lonr = lonr + cgx*dt*180./pi/(r_ear*cos(latr*pi/180.0))
    latr = latr + cgy*dt*180./pi/r_ear
    zr   = zr + cgz*dt

    tsum = tsum + dt
  else
    zr = zr + delzs
  end if

  mwn_last2  = mwn_last
  mwn_last   = mwn
  ifreq_last = ifreq

  enddo  ! n_iter  : end DO-SENTENCE WITH TIME

  if (count .eq. 1) then
    print*, ' cannot propagate initially ...    -----------'
    count = 2
    if (endcond .eq. 1)  endcond = 0
    vars(1,9)  = real(lonr)
    vars(1,10) = real(latr)
    vars(1,11) = real(zr)
    vars(1,12) = real(tsum)
  end if

  print *
  print *
  print *

  if (endcond .eq. 0)  term = '0, no calculation performed'
  if (endcond .eq. 1)  term = '1, reflection level'
  if (endcond .eq. 2)  term = '2, critical level'
  if (endcond .eq. 3)  term = '3, time limitation'
  if (endcond .eq. 9)  term = '9, out of domain'

  nt = count - 1


  return

end subroutine ray_integ


subroutine grad_int(nx,ny,nz,dx,dy,dz,f,i0,j0,k0,delx,dely,delz,dfdx,dfdy,dfdz)

  use interpolate

  implicit none

  integer, intent(in)  ::  nx, ny, nz, i0, j0, k0
  real,    intent(in)  ::  f(nx,ny,nz)
  real*8,  intent(in)  ::  dx, dy, dz, delx, dely, delz
  real*8,  intent(out) ::  dfdx, dfdy, dfdz

  integer :: i,j,k, i1, j1, k1, ii
  real, dimension(i0:i0+1,j0:j0+1,k0:k0+1) :: dfdx01, dfdy01, dfdz01


  i1 = i0 + 1
  j1 = j0 + 1
  k1 = k0 + 1
  if (dely .eq. 0.d0)  j1 = j0
  if (delz .eq. 0.d0)  k1 = k0

  do k=k0, k1
  do j=j0, j1
  do i=i0, i1
    if (i.eq.1 .or. i.eq.nx+1) then
      dfdx01(i,j,k) = ( f(2,j,k)-f(nx,j,k) ) / real(2.*dx)
    else if (i .eq. nx) then
      dfdx01(i,j,k) = ( f(1,j,k)-f(nx-1,j,k) ) / real(2.*dx)
    else
      dfdx01(i,j,k) = ( f(i+1,j,k)-f(i-1,j,k) ) / real(2.*dx)
    end if
  enddo
  enddo    
  enddo    

  do k=k0, k1
  do i=i0, i1
    ii = i
    if (ii .eq. nx+1)  ii = 1
    do j=j0, j1
      if (j .eq. 1) then
        dfdy01(i,j,k) = 2.*(f(ii,2,k)-f(ii,1,k)) / real(dy)            &
                         - (f(ii,3,k)-f(ii,1,k)) / real(2.*dy)
      else if (j .eq. ny) then
        dfdy01(i,j,k) = 2.*(f(ii,ny,k)-f(ii,ny-1,k)) / real(dy)        &
                         - (f(ii,ny,k)-f(ii,ny-2,k)) / real(2.*dy)
      else
        dfdy01(i,j,k) = ( f(ii,j+1,k)-f(ii,j-1,k) ) / real(2.*dy)
      end if
    enddo
  enddo
  enddo

  do j=j0, j1
  do i=i0, i1
    ii = i
    if (ii .eq. nx+1)  ii = 1
    do k=k0, k1
      if (k .eq. 1) then
        dfdz01(i,j,k) = 2.*(f(ii,j,2)-f(ii,j,1)) / real(dz)            &
                         - (f(ii,j,3)-f(ii,j,1)) / real(2.*dz)
      else if (k .eq. nz) then
        dfdz01(i,j,k) = 2.*(f(ii,j,nz)-f(ii,j,nz-1)) / real(dz)        &
                         - (f(ii,j,nz)-f(ii,j,nz-2)) / real(2.*dz)
      else
        dfdz01(i,j,k) = ( f(ii,j,k+1)-f(ii,j,k-1) ) / real(2.*dz)
      end if
    end do
  end do
  end do

  call interpol(dfdx01(i0,j0,k0), dfdx01(i1,j0,k0),  &
                dfdx01(i0,j1,k0), dfdx01(i1,j1,k0),  &
                dfdx01(i0,j0,k1), dfdx01(i1,j0,k1),  &
                dfdx01(i0,j1,k1), dfdx01(i1,j1,k1),  &
                dx, dy, dz, delx, dely, delz, dfdx)

  call interpol(dfdy01(i0,j0,k0), dfdy01(i1,j0,k0),  &
                dfdy01(i0,j1,k0), dfdy01(i1,j1,k0),  &
                dfdy01(i0,j0,k1), dfdy01(i1,j0,k1),  &
                dfdy01(i0,j1,k1), dfdy01(i1,j1,k1),  &
                dx, dy, dz, delx, dely, delz, dfdy)

  call interpol(dfdz01(i0,j0,k0), dfdz01(i1,j0,k0),  &
                dfdz01(i0,j1,k0), dfdz01(i1,j1,k0),  &
                dfdz01(i0,j0,k1), dfdz01(i1,j0,k1),  &
                dfdz01(i0,j1,k1), dfdz01(i1,j1,k1),  &
                dx, dy, dz, delx, dely, delz, dfdz)


  return

end subroutine grad_int
END module raytrace


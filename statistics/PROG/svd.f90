program svd

  use svdpack
  use netcdfio

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nxop0 = 125, nyop0 = 106, nxoa0 = 80, nyoa0 = 177
  integer, parameter :: nxap0 = 125, nyap0 = 54
  integer, parameter :: nxop = 42, nyop = 21, nxoa = 27, nyoa = 35
  integer, parameter :: nxo = 320, nyo = 320
  integer, parameter :: nxap = 42, nyap = 18
  integer, parameter :: nxa = 320, nya = 160
  integer, parameter :: nt = 1200
  real,    parameter :: dt = 1./12.
  real,    parameter :: missvd = 1.e20, missvo = 1.e32
!--------------------------------------------------------------

  integer :: i,j,n, ncid, nf, ii, jj, tag, missn1, missn2, missn3
  real    :: pi

  real, dimension(nxop,nyop,nt) :: top
  real, dimension(nxoa,nyoa,nt) :: toa
  real, dimension(nxap,nyap,nt) :: uap

  real    :: lonop(nxop), latop(nyop), lonoa(nxoa), latoa(nyoa)
  real    :: lonap(nxap), latap(nyap)
  real    :: lono(nxo), lato(nyo), lona(nxa), lata(nya), t1(nt)
  real    :: lonop0(nxop0), latop0(nyop0), lonoa0(nxoa0), latoa0(nyoa0)
  real    :: lonap0(nxap0), latap0(nyap0)
  real    :: dsinphiop(nyop), dsinphioa(nyoa), normop(nxop,nyop), normoa(nxoa,nyoa)
  real    :: dsinphiap(nyap), normap(nxap,nyap)
  real    :: bndop(2,nyop0), bndoa(2,nyoa0), bndap(2,nyap0)
  real    :: avg1(nt), avg2(nt), tavg1, tavg2
  real, dimension(:,:,:), allocatable   :: toa1, toa2, top0, toa0, uap0

  integer :: tmpi, ntsf
  real*8  :: tmp1(nxa), tmp2(nya), tmp3(nxo), tmp4(nyo)
  real*8  :: tmp5(2,nya), tmp6(2,nyo)
  real    :: tmpr1, tmpr2

  character*128 :: ofname1, ofname2
  character*128 :: fdir(2)
  character*128 :: fn_to, fn_ua(5)


  pi = acos(-1.)

! read data ------------------------------------------------------------
  fdir(1) = '/export30/kyh/MIROC/pictl/'
  fdir(2) = '/export30/kyh/MIROC/CO2/'
  write(fn_to,'(a)') trim(fdir(1))//'tos_O1.nc'
  write(fn_ua(1),'(a)') trim(fdir(1))//'ua_A1_1_20.nc'
  write(fn_ua(2),'(a)') trim(fdir(1))//'ua_A1_21_40.nc'
  write(fn_ua(3),'(a)') trim(fdir(1))//'ua_A1_41_60.nc'
  write(fn_ua(4),'(a)') trim(fdir(1))//'ua_A1_61_80.nc'
  write(fn_ua(5),'(a)') trim(fdir(1))//'ua_A1_81_100.nc'

! axis
  call opennc(fn_ua(1),ncid)
  call dget1d(ncid,'lon',nxa,tmp1)
  call dget1d(ncid,'lat',nya,tmp2)
  call dget2d(ncid,'lat_bnds',2,nya,tmp5)
  call closenc(ncid)
  lona = real(tmp1)
  lata = real(tmp2)
  lonap0(:) = lona(108:232)
  latap0(:) = lata(54:107)
  bndap(:,:) = real(tmp5(:,54:107))

  call opennc(fn_to,ncid)
  call dget1d(ncid,'lon',nxo,tmp3)
  call dget1d(ncid,'lat',nyo,tmp4)
  call dget2d(ncid,'lat_bnds',2,nyo,tmp6)
  lono = real(tmp3)
  lato = real(tmp4)
  lonop0(:) = lono(108:232)
  latop0(:) = lato(108:213)
  lonoa0(:61) = lono(260:) - 360.
  lonoa0(62:) = lono(1:19)
  latoa0(:) = lato(90:266)
  bndop(:,:) = real(tmp6(:,108:213))
  bndoa(:,:) = real(tmp6(:,90:266))

  do n=1, nt
    t1(n) = real(n-1)/12.
  enddo

! variables
  allocate(top0(nxop0,nyop0,nt))                        ; top0 = 0.
  allocate(toa0(nxoa0,nyoa0,nt))                        ; toa0 = 0.
  allocate(toa1(61,nyoa0,nt))                           ; toa1 = 0.
  allocate(toa2(19,83,nt))                              ; toa2 = 0.
  call geta3d(ncid,'tos',108,nxop0,108,nyop0,1,nt,top0)
  call geta3d(ncid,'tos',260,61,90,nyoa0,1,nt,toa1)
  call geta3d(ncid,'tos',1,19,90,83,1,nt,toa2)
  call closenc(ncid)

  toa0 = 1.e20
  toa0(:61,:,:) = toa1(:,:,:)
  toa0(62:,:83,:) = toa2(:,:,:)
  deallocate(toa1)  ;  deallocate(toa2)

  missn1 = 0
  jj = 0
  do j=3, nyop0, 5
    jj = jj + 1
    ii = 0
    do i=2, nxop0, 3
      ii = ii + 1
      top(ii,jj,:) = top0(i,j,:)
      tag = 0
      do n=1, nt
        if (top0(i,j,n) .eq. missvd)  tag = 1
      enddo
      if (tag .eq. 1) then
        missn1 = missn1 + 1
        top(ii,jj,:) = missvo
      endif
    enddo
  enddo
  deallocate(top0)

  missn2 = 0
  jj = 0 
  do j=3, nyoa0, 5 
    jj = jj + 1 
    ii = 0
    do i=2, nxoa0, 3 
      ii = ii + 1
      toa(ii,jj,:) = toa0(i,j,:)
      tag = 0
      do n=1, nt 
        if (toa0(i,j,n) .eq. missvd)  tag = 1
      enddo
      if (tag .eq. 1) then
        missn2 = missn2 + 1
        toa(ii,jj,:) = missvo
      endif
    enddo
  enddo
  deallocate(toa0)

  allocate(uap0(nxap0,nyap0,nt))                        ; uap0 = 0.
  tmpi = 0
  do nf=1, 5
    call opennc(fn_ua(nf),ncid)
    call dilen(ncid,'time',ntsf)
    call geta4d(ncid,'ua',108,nxap0,54,nyap0,1,1,1,ntsf, &
                uap0(:,:,tmpi+1:tmpi+ntsf))
    call closenc(ncid)
    tmpi = tmpi + ntsf
  enddo

  missn3 = 0
  jj = 0
  do j=2, nyap0, 3
    jj = jj + 1
    ii = 0
    do i=2, nxap0, 3
      ii = ii + 1
      uap(ii,jj,:) = uap0(i,j,:)
      tag = 0
      do n=1, nt
        if (uap0(i,j,n) .eq. missvd)  tag = 1
      enddo
      if (tag .eq. 1) then
        missn3 = missn3 + 1
        uap(ii,jj,:) = missvo
      endif
    enddo
  enddo
  deallocate(uap0)

  ii = 0 
  do i=2, nxop0, 3
    ii = ii + 1 
    lonop(ii) = lonop0(i)
  enddo 
  ii = 0
  do i=2, nxoa0, 3
    ii = ii + 1
    lonoa(ii) = lonoa0(i)
  enddo
  ii = 0
  do i=2, nxap0, 3
    ii = ii + 1
    lonap(ii) = lonap0(i)
  enddo
  jj = 0
  do j=3, nyop0, 5
    jj = jj + 1
    latop(jj) = latop0(j)
    dsinphiop(jj) = sin(pi/180.*bndop(2,j+2))-sin(pi/180.*bndop(1,j-2))
  enddo
  jj = 0 
  do j=3, nyoa0, 5
    jj = jj + 1
    latoa(jj) = latoa0(j)
    dsinphioa(jj) = sin(pi/180.*bndoa(2,j+2))-sin(pi/180.*bndoa(1,j-2))
  enddo
  jj = 0
  do j=2, nyap0, 3
    jj = jj + 1
    latap(jj) = latap0(j)
    dsinphiap(jj) = sin(pi/180.*bndap(2,j+1))-sin(pi/180.*bndap(1,j-1))
  enddo


  do i=1, nxop
    normop(i,:) = dsinphiop(:) / dsinphiop(11)
  enddo
  do i=1, nxoa
    normoa(i,:) = dsinphioa(:) / dsinphioa(15)
  enddo
  do i=1, nxap
    normap(i,:) = dsinphiap(:) / dsinphiap(9)
  enddo


! (1) SSTA over the Pacific & Atlantic ocean, 1000 mb wind -------------
  avg1 = 0.
  tavg1 = 0.
  do n=1, nt
    tmpr2 = 0.
    do j=1, nyop
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxop
        if (top(i,j,n) .ne. missvo) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + top(i,j,n)
        end if
      enddo
      avg1(n) = avg1(n) + tmpr1 * dsinphiop(j)
      tmpr2 = tmpr2 + tmpi * dsinphiop(j)
    enddo
    avg1(n) = avg1(n) / tmpr2
    tavg1 = tavg1 + avg1(n)/nt
  enddo
  top = top - tavg1

  avg2 = 0.
  tavg2 = 0.
  do n=1, nt
    tmpr2 = 0.
    do j=1, nyoa
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxoa
        if (toa(i,j,n) .ne. missvo) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + toa(i,j,n)
        end if
      enddo
      avg2(n) = avg2(n) + tmpr1 * dsinphioa(j)
      tmpr2 = tmpr2 + tmpi * dsinphioa(j)
    enddo 
    avg2(n) = avg2(n) / tmpr2
    tavg2 = tavg2 + avg2(n)/nt
  enddo 
  toa = toa - tavg2

  avg1 = 0.
  tavg1 = 0.
  do n=1, nt
    tmpr2 = 0.
    do j=1, nyap
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxap
        if (uap(i,j,n) .ne. missvo) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + uap(i,j,n)
        end if
      enddo
      avg1(n) = avg1(n) + tmpr1 * dsinphiap(j)
      tmpr2 = tmpr2 + tmpi * dsinphiap(j)
    enddo 
    avg1(n) = avg1(n) / tmpr2
    tavg1 = tavg1 + avg1(n)/nt
  enddo
  uap = uap - tavg1

  do j=1, nyop
  do i=1, nxop
    if (top(i,j,1) .ne. missvo) then
      call detrend_season(nt,top(i,j,:))
    end if
  enddo
  enddo
  do j=1, nyoa
  do i=1, nxoa
    if (toa(i,j,1) .ne. missvo) then
      call detrend_season(nt,toa(i,j,:))
    end if
  enddo
  enddo
  do j=1, nyap
  do i=1, nxap
    if (uap(i,j,1) .ne. missvo) then
      call detrend_season(nt,uap(i,j,:))
    end if
  enddo
  enddo

! (2) SVD analysis -----------------------------------------------------
  print*, nxop*nyop-missn1, nxap*nyap-missn3
  print*, nxop*nyop-missn1, nxoa*nyoa-missn2

  ofname1 = '../svd_pac_uas'
  call svdnc(ofname1,nt,t1, &
             nxop,lonop,nyop,latop,top,missn1,missvo,1,normop, &
             nxap,lonap,nyap,latap,uap,missn3,missvo,1,normap, &
             0,10)

  ofname2 = '../svd_pac_atl'
  call svdnc(ofname2,nt,t1, &
             nxop,lonop,nyop,latop,top,missn1,missvo,1,normop, &
             nxoa,lonoa,nyoa,latoa,toa,missn2,missvo,1,normoa, &
             0,10)


  stop

end program


subroutine detrend_season(nt,x)

  integer, intent(in)    :: nt
  real,    intent(inout) :: x(nt)

  integer :: i,j,n, nn(12)
  real    :: ave(12)

  n   = 0
  nn  = 0
  ave = 0.
  do i=1, nt
    n = n + 1
    ave(n) = ave(n) + x(i)
    nn(n) = nn(n) + 1
    if (n .eq. 12)  n = 0
  enddo
  ave(:) = ave(:) / nn(:)

  n = 0
  do i=1, nt
    n = n + 1
    x(i) = x(i) - ave(n)
    if (n .eq. 12)  n = 0
  enddo

end subroutine detrend_season


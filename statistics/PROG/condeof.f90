program cond_eof_analysis

  use eofpack
  use netcdfio

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nxp0 = 45, nyp0 = 18, nxa0 = 80, nya0 = 106
  integer, parameter :: nxa = 27, nya = 35
  integer, parameter :: nxo = 320, nyo = 320
  integer, parameter :: nt = 1200
  real,    parameter :: dt = 1./12.
  real,    parameter :: missvd = 1.e20, missvo = 1.e32
!--------------------------------------------------------------
  integer :: i,j,n, ncid, ii, jj, tag, missn
  real    :: pi

  real, dimension(nxa,nya,nt) :: toa

  real    :: lona(nxa), lata(nya)
  real    :: lono(nxo), lato(nyo), t1(nt), inino(nt)
  real    :: lona0(nxa0), lata0(nya0)
  real    :: dsinphip(nyp0), dsinphia(nya), norma(nxa,nya)
  real    :: bndp(2,nyp0), bnda(2,nya0)
  real    :: avg1(nt), avg2(nt), tavg1, tavg2
  real, dimension(:,:,:), allocatable :: toa1, toa2, top0, toa0

  integer :: tmpi
  real*8  :: tmp1(nxa0), tmp2(nya0), tmp3(nxo), tmp4(nyo)
  real*8  :: tmp5(2,nya0), tmp6(2,nyo)
  real    :: tmpr1, tmpr2, tmpr3

  character*128 :: ofname1, ofname2
  character*128 :: fdir(2)
  character*128 :: fn_to


  pi = acos(-1.)

! read data ------------------------------------------------------------
  fdir(1) = '/export30/kyh/MIROC/pictl/'
  fdir(2) = '/export30/kyh/MIROC/CO2/'
  write(fn_to,'(a)') trim(fdir(1))//'tos_O1.nc'

! axis
  call opennc(fn_to,ncid)
  call dget1d(ncid,'lon',nxo,tmp3)
  call dget1d(ncid,'lat',nyo,tmp4)
  call dget2d(ncid,'lat_bnds',2,nyo,tmp6)
  lono = real(tmp3)
  lato = real(tmp4)
  lona0(:61) = lono(260:) - 360.
  lona0(62:) = lono(1:19)
  lata0(:) = lato(108:213)
  bndp(:,:) = real(tmp6(:,152:169))
  bnda(:,:) = real(tmp6(:,108:213))

  do n=1, nt
    t1(n) = real(n-1)/12.
  enddo

! variables
  allocate(top0(nxp0,nyp0,nt))                        ; top0 = 0.
  allocate(toa0(nxa0,nya0,nt))                        ; toa0 = 0.
  allocate(toa1(61,nya0,nt))                          ; toa1 = 0.
  allocate(toa2(19,65,nt))                            ; toa2 = 0.
  call geta3d(ncid,'tos',170,nxp0,152,nyp0,1,nt,top0)
  call geta3d(ncid,'tos',260,61,108,nya0,1,nt,toa1)
  call geta3d(ncid,'tos',1,19,108,65,1,nt,toa2)
  call closenc(ncid)
  toa0 = 1.e20
  toa0(:61,:,:) = toa1(:,:,:)
  toa0(62:,:65,:) = toa2(:,:,:)
  deallocate(toa1)  ;  deallocate(toa2)

  missn = 0
  jj = 0 
  do j=2, nya0, 3
    jj = jj + 1 
    ii = 0
    do i=2, nxa0, 3 
      ii = ii + 1
      toa(ii,jj,:) = toa0(i,j,:)
      tag = 0
      do n=1, nt 
        if (toa0(i,j,n) .eq. missvd)  tag = 1
      enddo
      if (tag .eq. 1) then
        missn = missn + 1
        toa(ii,jj,:) = missvo
      endif
    enddo
  enddo
  deallocate(toa0)

  dsinphip(:) = sin(pi/180.*bndp(2,:))-sin(pi/180.*bndp(1,:))

  ii = 0
  do i=2, nxa0, 3
    ii = ii + 1
    lona(ii) = lona0(i)
  enddo
  jj = 0 
  do j=2, nya0, 3
    jj = jj + 1
    lata(jj) = lata0(j)
    dsinphia(jj) = sin(pi/180.*bnda(2,j+1))-sin(pi/180.*bnda(1,j-1))
  enddo

  do i=1, nxa
    norma(i,:) = dsinphia(:) / dsinphia(18)
  enddo

! (1) Nino3.4 index & SSTA over Atlantic ocean -------------------------
  inino = 0.
  tmpr3 = 0.
  do n=1, nt
    tmpr2 = 0.
    do j=1, nyp0
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxp0
        if (top0(i,j,n) .ne. 1.e20) then
          tmpi = tmpi + 1 
          tmpr1 = tmpr1 + top0(i,j,n)
        end if
      enddo
      inino(n) = inino(n) + tmpr1 * dsinphip(j)
      tmpr2 = tmpr2 + tmpi * dsinphip(j)
    enddo
    inino(n) = inino(n) / tmpr2
    tmpr3 = tmpr3 + inino(n)/nt
  enddo
  inino(:) = inino(:) - tmpr3

  call detrend_season(nt,inino)

  avg1 = 0.
  tavg1 = 0.
  do n=1, nt
    tmpr2 = 0.
    do j=1, nya
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxa
        if (toa(i,j,n) .ne. missvo) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + toa(i,j,n)
        end if
      enddo
      avg1(n) = avg1(n) + tmpr1 * dsinphia(j)
      tmpr2 = tmpr2 + tmpi * dsinphia(j)
    enddo 
    avg1(n) = avg1(n) / tmpr2
    tavg1 = tavg1 + avg1(n)/nt
  enddo 
  toa = toa - tavg1

  do j=1, nya
  do i=1, nxa
    if (toa(i,j,1) .ne. missvo) then
      call detrend_season(nt,toa(i,j,:))
    end if
  enddo
  enddo

! (2) Conditional EOF analysis -----------------------------------------
  ofname1 = '../cond_eof' 
  call cond_eofnc(ofname1,nxa,lona,nya,lata,nt,t1,toa,missn,missvo,inino,&
                  30,1,norma)


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


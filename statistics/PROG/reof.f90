program reof_analysis_and_ssa

  use eofpack
  use ssapack
  use netcdfio

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nxp0 = 125, nyp0 = 106, nxa0 = 80, nya0 = 177
  integer, parameter :: nxp = 42, nyp = 21, nxa = 27, nya = 35
  integer, parameter :: nxo = 320, nyo = 320
  integer, parameter :: nt = 1200
  real,    parameter :: missvd = 1.e20, missvo = 1.e32
!--------------------------------------------------------------
  integer, parameter :: neof = 15, nch = 2, maxlag = 84

  integer :: i,j,n, ncid, ii, jj, tag, missn1, missn2
  real    :: pi

  real, dimension(nxp,nyp,nt) :: top
  real, dimension(nxa,nya,nt) :: toa
  real, dimension(nxp,nyp,neof) :: eofp
  real, dimension(nxa,nya,neof) :: eofa
  real, dimension(neof) :: ev100p, ev100a

  real    :: pc_pac(nt,3), pc_atl(nt,3), pcs(nt,nch), ch(nch)

  real    :: lonp(nxp), latp(nyp), lona(nxa), lata(nya)
  real    :: lono(nxo), lato(nyo), t1(nt)
  real    :: lonp0(nxp0), latp0(nyp0), lona0(nxa0), lata0(nya0)
  real    :: dsinphip(nyp), dsinphia(nya), normp(nxp,nyp), norma(nxa,nya)
  real    :: bndp(2,nyp0), bnda(2,nya0)
  real    :: avg1(nt), avg2(nt), tavg1, tavg2
  real, dimension(:,:,:), allocatable :: toa1, toa2, top0, toa0

  integer :: tmpi
  real*8  :: tmp1(nxa0), tmp2(nya0), tmp3(nxo), tmp4(nyo)
  real*8  :: tmp5(2,nya0), tmp6(2,nyo)
  real    :: tmpr1, tmpr2, tmp1d(nch)

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
  lonp0(:) = lono(108:232)
  latp0(:) = lato(108:213)
  lona0(:61) = lono(260:) - 360.
  lona0(62:) = lono(1:19)
  lata0(:) = lato(90:266)
  bndp(:,:) = real(tmp6(:,108:213))
  bnda(:,:) = real(tmp6(:,90:266))

  do n=1, nt
    t1(n) = real(n-1)/12.
  enddo

! variables
  allocate(top0(nxp0,nyp0,nt))                        ; top0 = 0.
  allocate(toa0(nxa0,nya0,nt))                        ; toa0 = 0.
  allocate(toa1(61,nya0,nt))                          ; toa1 = 0.
  allocate(toa2(19,83,nt))                            ; toa2 = 0.
  call geta3d(ncid,'tos',108,nxp0,108,nyp0,1,nt,top0)
  call geta3d(ncid,'tos',260,61,90,nya0,1,nt,toa1)
  call geta3d(ncid,'tos',1,19,90,83,1,nt,toa2)
  call closenc(ncid)
  toa0 = 1.e20
  toa0(:61,:,:) = toa1(:,:,:)
  toa0(62:,:83,:) = toa2(:,:,:)
  deallocate(toa1)  ;  deallocate(toa2)

  missn1 = 0
  jj = 0
  do j=3, nyp0, 5
    jj = jj + 1
    ii = 0
    do i=2, nxp0, 3
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
  do j=3, nya0, 5 
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
        missn2 = missn2 + 1
        toa(ii,jj,:) = missvo
      endif
    enddo
  enddo
  deallocate(toa0)

  ii = 0 
  do i=2, nxp0, 3
    ii = ii + 1 
    lonp(ii) = lonp0(i)
  enddo 
  ii = 0
  do i=2, nxa0, 3
    ii = ii + 1
    lona(ii) = lona0(i)
  enddo
  jj = 0
  do j=3, nyp0, 5
    jj = jj + 1
    latp(jj) = latp0(j)
    dsinphip(jj) = sin(pi/180.*bndp(2,j+2))-sin(pi/180.*bndp(1,j-2))
  enddo
  jj = 0 
  do j=3, nya0, 5
    jj = jj + 1
    lata(jj) = lata0(j)
    dsinphia(jj) = sin(pi/180.*bnda(2,j+2))-sin(pi/180.*bnda(1,j-2))
  enddo

  do i=1, nxp
    normp(i,:) = dsinphip(:) / dsinphip(11)
  enddo
  do i=1, nxa
    norma(i,:) = dsinphia(:) / dsinphia(15)
  enddo

! (1) SSTA over the Pacific & Atlantic ocean ---------------------------
  avg1 = 0.
  tavg1 = 0.
  do n=1, nt
    tmpr2 = 0.
    do j=1, nyp
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxp
        if (top(i,j,n) .ne. missvo) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + top(i,j,n)
        end if
      enddo
      avg1(n) = avg1(n) + tmpr1 * dsinphip(j)
      tmpr2 = tmpr2 + tmpi * dsinphip(j)
    enddo
    avg1(n) = avg1(n) / tmpr2
    tavg1 = tavg1 + avg1(n)/nt
  enddo
  top = top - tavg1

  avg2 = 0.
  tavg2 = 0.
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
      avg2(n) = avg2(n) + tmpr1 * dsinphia(j)
      tmpr2 = tmpr2 + tmpi * dsinphia(j)
    enddo 
    avg2(n) = avg2(n) / tmpr2
    tavg2 = tavg2 + avg2(n)/nt
  enddo 
  toa = toa - tavg2

  do j=1, nyp
  do i=1, nxp
    if (top(i,j,1) .ne. missvo) then
      call detrend_season(nt,top(i,j,:))
    end if
  enddo
  enddo
  do j=1, nya
  do i=1, nxa
    if (toa(i,j,1) .ne. missvo) then
      call detrend_season(nt,toa(i,j,:))
    end if
  enddo
  enddo

! (2) R-EOF analysis ---------------------------------------------------
!  ofname1 = '../eof_pac' 
!  call eofnc(ofname1,nxp,lonp,nyp,latp,nt,t1,top,missn1,missvo,30,1,normp)

!  ofname2 = '../eof_atl'
!  call eofnc(ofname2,nxa,lona,nya,lata,nt,t1,toa,missn2,missvo,30,1,norma)

  call opennc('../res/eof_pac.nc',ncid)
  call geta3d(ncid,'EOF',1,nxp,1,nyp,1,neof,eofp)
  call geta1d(ncid,'e_value_100',1,neof,ev100p)
  call closenc(ncid)
  call opennc('../res/eof_atl.nc',ncid)
  call geta3d(ncid,'EOF',1,nxa,1,nya,1,neof,eofa)
  call geta1d(ncid,'e_value_100',1,neof,ev100a)
  call closenc(ncid)

  ofname1 = '../reof_pac'
  call reofnc(ofname1,nxp,lonp,nyp,latp,nt,t1,top,missn1,missvo, &
              neof,eofp,ev100p,50,1,normp)

  ofname2 = '../reof_atl'
  call reofnc(ofname2,nxa,lona,nya,lata,nt,t1,toa,missn2,missvo, &
              neof,eofa,ev100a,50,1,norma)

! (3) MSSA -------------------------------------------------------------
  call opennc('../res/reof_pac.nc',ncid)
  call geta2d(ncid,'PC',1,nt,1,3,pc_pac)
  call closenc(ncid)
  call opennc('../res/reof_atl.nc',ncid)
  call geta2d(ncid,'PC',1,nt,1,3,pc_atl)
  call closenc(ncid)

  pcs(:,1) = pc_pac(:,1)
  pcs(:,2) = pc_atl(:,1)
  do i=1, nch
    ch(i) = real(i)
    tmp1d(i) = 1.
  enddo
  ofname1 = '../mssa11'
  call ssanc(ofname1,2,ch,1,(/1./),nt,t1,pcs,0,1.e32,maxlag,15,0,tmp1d)

  pcs(:,1) = pc_pac(:,3)
  pcs(:,2) = pc_atl(:,1)
  do i=1, nch
    ch(i) = real(i)
    tmp1d(i) = 1.
  enddo
  ofname1 = '../mssa31'
  call ssanc(ofname1,2,ch,1,(/1./),nt,t1,pcs,0,1.e32,maxlag,15,0,tmp1d)

  pcs(:,1) = pc_pac(:,1)
  pcs(:,2) = pc_atl(:,2) 
  do i=1, nch  
    ch(i) = real(i)  
    tmp1d(i) = 1.  
  enddo 
  ofname1 = '../mssa12' 
  call ssanc(ofname1,2,ch,1,(/1./),nt,t1,pcs,0,1.e32,maxlag,15,0,tmp1d)

  pcs(:,1) = pc_pac(:,3)
  pcs(:,2) = pc_atl(:,2) 
  do i=1, nch  
    ch(i) = real(i)  
    tmp1d(i) = 1.
  enddo 
  ofname1 = '../mssa32'
  call ssanc(ofname1,2,ch,1,(/1./),nt,t1,pcs,0,1.e32,maxlag,15,0,tmp1d)


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


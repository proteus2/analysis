program spectral_analysis_2b

  use netcdfio
  use specanal

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nxa = 320, nya = 160, nxo = 320, nyo = 320
  integer, parameter :: nt1all = 1200, nt2all = 960
  integer, parameter :: ts1 = 1, ts2 = 1, nt1 = nt1all, nt2 = nt2all
  integer, parameter :: nfu1 = 5        ! CTL (pre-ind.)
  integer, parameter :: nfu2 = 4        ! CO2
!--------------------------------------------------------------
  integer, parameter :: maxlag1 = nt1/2, maxlag2 = nt2/2
  real,    parameter :: dt = 1.
!--------------------------------------------------------------
  integer, parameter :: nps1 = maxlag1+1, nps2 = maxlag2+1
  real               :: freq1(nps1), freq2(nps2)

  integer :: i,j,n, nf, ncid, ii, jj
  real    :: pi

  real, dimension(nxo,nyo,nt1) :: to1
  real, dimension(nxo,nyo,nt2) :: to2

  real, dimension(nt1)         :: inino1, to1gave
  real, dimension(nt2)         :: inino2, to2gave
  real, dimension(nt1)         :: inino1de, to1gavede, tatl1de
  real, dimension(nt2)         :: inino2de, to2gavede, tatl2de
  real, dimension(nt1,12)      :: tatl1
  real, dimension(nt2,12)      :: tatl2

  real, dimension(nps1)        :: psin1, psto1, amp1, phs1, coh1
  real, dimension(nps2)        :: psin2, psto2, amp2, phs2, coh2
  real, dimension(nps1,12)     :: amp3, phs3, coh3
  real, dimension(nps2,12)     :: amp4, phs4, coh4

  real, dimension(nps1*2-1)    :: corr1, lag1
  real, dimension(nps2*2-1)    :: corr2, lag2
  real, dimension(nps1*2-1,12) :: corr3
  real, dimension(nps2*2-1,12) :: corr4

  real    :: lono(nxo), lato(nyo), t1(nt1), t2(nt2), lat_atl(12)
  real    :: dsinphio(nyo)
  integer :: istt(nyo), iend(nyo), jstt(12), jend(12)

  integer :: tmpi
  real*8  :: tmp1(nxo), tmp2(nyo)
  real*8  :: tmp3(2,nyo)
  real    :: tmpr1, tmpr2, tmpr3

  character*128 :: fdir(2)
  character*128 :: fn_to1
  character*128 :: fn_to2


  pi = acos(-1.)

! read data ------------------------------------------------------------
  fdir(1) = '/export30/kyh/MIROC/pictl/'
  fdir(2) = '/export30/kyh/MIROC/CO2/'
  write(fn_to1,'(a)') trim(fdir(1))//'tos_O1.nc'
  write(fn_to2,'(a)') trim(fdir(2))//'tos_O1.nc'

! axis
  call opennc(fn_to1,ncid)
  call dget1d(ncid,'lon',nxo,tmp1)
  call dget1d(ncid,'lat',nyo,tmp2)
  call dget2d(ncid,'lat_bnds',2,nyo,tmp3)
  call closenc(ncid)
  lono = real(tmp1)
  lato = real(tmp2)
  dsinphio(:) = real(dsin(pi/180.*tmp3(2,:))-dsin(pi/180.*tmp3(1,:)))

  do n=1, nt1
    t1(n) = real(n-1)/12.
  enddo
  do n=1, nt2
    t2(n) = real(n-1)/12.
  enddo

! variables
  call opennc(fn_to1,ncid)
  call geta3d(ncid,'tos',1,nxo,1,nyo,ts1,nt1,to1)
  call closenc(ncid)
  call opennc(fn_to2,ncid)
  call geta3d(ncid,'tos',1,nxo,1,nyo,ts2,nt2,to2)
  call closenc(ncid)


! (1) Nino3.4 index, global mean SST -----------------------------------
  inino1 = 0.
  tmpr3 = 0.
  do n=1, nt1
    tmpr2 = 0.
    do j=152, 169
      tmpi  = 0
      tmpr1 = 0.
      do i=170, 214
        if (to1(i,j,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to1(i,j,n)
        end if
      enddo
      inino1(n) = inino1(n) + tmpr1 * dsinphio(j)
      tmpr2 = tmpr2 + tmpi * dsinphio(j)
    enddo
    inino1(n) = inino1(n) / tmpr2
    tmpr3 = tmpr3 + inino1(n)/nt1
  enddo
  inino1(:) = inino1(:) - tmpr3

  inino2 = 0.
  tmpr3 = 0.
  do n=1, nt2
    tmpr2 = 0.
    do j=152, 169
      tmpi  = 0
      tmpr1 = 0.
      do i=170, 214
        if (to2(i,j,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to2(i,j,n)
        end if
      enddo
      inino2(n) = inino2(n) + tmpr1 * dsinphio(j)
      tmpr2 = tmpr2 + tmpi * dsinphio(j)
    enddo
    inino2(n) = inino2(n) / tmpr2
    tmpr3 = tmpr3 + inino2(n)/nt2
  enddo
  inino2(:) = inino2(:) - tmpr3

  to1gave = 0.
  do n=1, nt1
    tmpr2 = 0.
    do j=1, nyo
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxo
        if (to1(i,j,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to1(i,j,n)
        end if
      enddo
      to1gave(n) = to1gave(n) + tmpr1 * dsinphio(j)
      tmpr2 = tmpr2 + tmpi * dsinphio(j)
    enddo
    to1gave(n) = to1gave(n) / tmpr2
  enddo
  
  to2gave = 0.
  do n=1, nt2
    tmpr2 = 0.
    do j=1, nyo
      tmpi  = 0
      tmpr1 = 0.
      do i=1, nxo
        if (to2(i,j,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to2(i,j,n)
        end if
      enddo
      to2gave(n) = to2gave(n) + tmpr1 * dsinphio(j)
      tmpr2 = tmpr2 + tmpi * dsinphio(j)
    enddo 
    to2gave(n) = to2gave(n) / tmpr2
  enddo

  inino1de  = inino1
  to1gavede = to1gave
  inino2de  = inino2
  to2gavede = to2gave
  ! remove linear trend
  call detrend1(nt1,t1,inino1de)
  call detrend1(nt1,t1,to1gavede)
  call detrend1(nt2,t2,inino2de)
  call detrend1(nt2,t2,to2gavede)
  ! remove seasonal cycle
  call detrend_season(nt1,inino1de)
  call detrend_season(nt1,to1gavede)
  call detrend_season(nt2,inino2de)
  call detrend_season(nt2,to2gavede)

  call crosss(inino1de,to1gavede,dt,nt1,maxlag1,corr1, &
              freq1,psin1,psto1,amp1,phs1,coh1)
  call crosss(inino2de,to2gavede,dt,nt2,maxlag2,corr2, &
              freq2,psin2,psto2,amp2,phs2,coh2)

  do i=1, nps1*2-1
    lag1(i) = (i-1-maxlag1) * dt
  enddo
  do i=1, nps2*2-1
    lag2(i) = (i-1-maxlag2) * dt
  enddo


  call out1d2('../sa2b1pic.nc',1,(/'COR'/),(/corr1/),'lag',nps1*2-1,lag1, &
              3,(/'COH','PHASE','AMP'/),(/coh1,phs1,amp1/),               &
              'freq',nps1,freq1,'Cross-spectrum_pictl')
  call out1d2('../sa2b1co2.nc',1,(/'COR'/),(/corr2/),'lag',nps2*2-1,lag2, &
              3,(/'COH','PHASE','AMP'/),(/coh2,phs2,amp2/),               &
              'freq',nps2,freq2,'Cross-spectrum_CO2')

! (2) Basin averaged SST (Atlantic ocean) ------------------------------
  istt(54:63)   = 261 !
  istt(64:68)   = 259
  istt(69)      = 258
  istt(70:96)   = 259
  istt(97:178)  = 265
  istt(179:187) = 246
  istt(188:192) = 239
  istt(193:246) = 233
  istt(247:267) = 256
  iend(54:200)  = 20
  iend(201:237) = 314
  iend(238:245) = 1
  iend(246:267) = 9
  do j=54, 267
    if ( iend(j) .lt. istt(j) )  iend(j) = iend(j) + nxo
  enddo

  jstt = (/54,72,90,108,125,143,161,179,197,214,232,250/)
  jend = (/71,89,107,124,142,160,178,196,213,231,249,267/)

  do jj=1, 12
    lat_atl(jj) = (jj-6.5) * 10.
  enddo

  tatl1 = 0.
  do n=1, nt1
  do jj=1, 12
    tmpr2 = 0.
    do j=jstt(jj), jend(jj)
      tmpi  = 0
      tmpr1 = 0.
      do i=istt(j), iend(j)
        ii = i
        if (ii .gt. nxo)  ii = ii - nxo
        if (to1(ii,j,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to1(ii,j,n)
        end if
      enddo
      tatl1(n,jj) = tatl1(n,jj) + tmpr1 * dsinphio(j)
      tmpr2 = tmpr2 + tmpi * dsinphio(j)
    enddo
    tatl1(n,jj) = tatl1(n,jj) / tmpr2
  enddo
  enddo

  tatl2 = 0.
  do n=1, nt2
  do jj=1, 12
    tmpr2 = 0.
    do j=jstt(jj), jend(jj)
      tmpi  = 0
      tmpr1 = 0.
      do i=istt(j), iend(j)
        ii = i
        if (ii .gt. nxo)  ii = ii - nxo
        if (to2(ii,j,n) .ne. 1.e20) then
          tmpi = tmpi + 1
          tmpr1 = tmpr1 + to2(ii,j,n)
        end if
      enddo
      tatl2(n,jj) = tatl2(n,jj) + tmpr1 * dsinphio(j)
      tmpr2 = tmpr2 + tmpi * dsinphio(j)
    enddo
    tatl2(n,jj) = tatl2(n,jj) / tmpr2
  enddo
  enddo

  do jj=1, 12
    tatl1de(:) = tatl1(:,jj)
    call detrend1(nt1,t1,tatl1de)
    call detrend_season(nt1,tatl1de)
    call crosss(inino1de,tatl1de,dt,nt1,maxlag1,corr1, &
                freq1,psin1,psto1,amp1,phs1,coh1)
    corr3(:,jj) = corr1(:)
    coh3(:,jj) = coh1(:)
    amp3(:,jj) = amp1(:)
    phs3(:,jj) = phs1(:)
  enddo

  do jj=1, 12
    tatl2de(:) = tatl2(:,jj)
    call detrend1(nt2,t2,tatl2de)
    call detrend_season(nt2,tatl2de)
    call crosss(inino2de,tatl2de,dt,nt2,maxlag2,corr2, &
                freq2,psin2,psto2,amp2,phs2,coh2)
    corr4(:,jj) = corr2(:)
    coh4(:,jj) = coh2(:)
    amp4(:,jj) = amp2(:)
    phs4(:,jj) = phs2(:)
  enddo

  call out2d2('../sa2b2pic.nc',1,(/'COR'/),(/corr3/),       &
              'lag',nps1*2-1,lag1,'lat1',12,lat_atl,        &
              3,(/'COH','PHASE','AMP'/),(/coh3,phs3,amp3/), &
              'freq',nps1,freq1,'lat2',12,lat_atl,'Cross-spectrum_pictl')
  call out2d2('../sa2b2co2.nc',1,(/'COR'/),(/corr4/),       &
              'lag',nps2*2-1,lag2,'lat1',12,lat_atl,        &
              3,(/'COH','PHASE','AMP'/),(/coh4,phs4,amp4/), &
              'freq',nps2,freq2,'lat2',12,lat_atl,'Cross-spectrum_CO2')


end program


subroutine detrend1(nt, t, x)

  use regress

  integer,             intent(in)    :: nt
  real, dimension(nt), intent(inout) :: t, x

  integer :: i,j
  real    :: a(2), afunc(50), sig(nt), covar(12,12), chisq, trend(nt)

  a(:)   = 0.0         ! coeff. for ia=0
  sig(:) = 1.0 
  call lfit(t,x,sig,nt,a,(/1,1/),2,covar,12,chisq)
  trend = 0.
  do i=1, nt 
    call funcs(t(i),afunc,2)
    do j=1, 2 
      trend(i) = trend(i) + a(j)*afunc(j)
    end do 
    x(i) = x(i) - trend(i)
  end do

  return

end subroutine detrend1 


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

  return

end subroutine detrend_season


program nonlinear_regression

  use netcdfio

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nxo = 320, nyo = 320
  integer, parameter :: nt1all = 1200, nt2all = 960
  integer, parameter :: ts1 = 1, ts2 = 1, nt1 = nt1all, nt2 = nt2all
!--------------------------------------------------------------
  integer :: i,j,n, nf, ncid, ii, jj
  real    :: pi

  real, dimension(nxo,nyo,nt1) :: to1
  real, dimension(nxo,nyo,nt2) :: to2

  real, dimension(nt1)         :: inino1
  real, dimension(nt2)         :: inino2
  real, dimension(nt1,5)      :: tatl1
  real, dimension(nt2,5)      :: tatl2

  real    :: lono(nxo), lato(nyo), t1(nt1), t2(nt2), lat_atl(5)
  real    :: dsinphio(nyo)
  integer :: istt(nyo), iend(nyo), jstt(5), jend(5)

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

  call detrend1(nt1,t1,inino1)
  call detrend1(nt2,t2,inino2)
  call detrend_season(nt1,inino1)
  call detrend_season(nt2,inino2)

! (2) Basin averaged SST (Atlantic ocean) ------------------------------
  istt(116:178) = 265
  istt(179:187) = 246
  istt(188:192) = 239
  istt(193:205) = 233
  iend(116:200) = 20
  iend(201:205) = 314
  do j=116, 205
    if ( iend(j) .lt. istt(j) )  iend(j) = iend(j) + nxo
  enddo

  jstt = (/116,134,152,170,188/)
  jend = (/133,151,169,187,205/)

  do jj=1, 5
    lat_atl(jj) = (jj-3) * 10.
  enddo

  tatl1 = 0.
  do n=1, nt1
  do jj=1, 5
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
  do jj=1, 5
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

  do jj=1, 5
    call detrend1(nt1,t1,tatl1(:,jj))
    call detrend1(nt2,t2,tatl2(:,jj))
    call detrend_season(nt1,tatl1(:,jj))
    call detrend_season(nt2,tatl2(:,jj))
  enddo

  open(1,file='../nlreg_nino34_pic.asc')
  open(2,file='../nlreg_nino34_co2.asc')
  do n=1, nt1
    write(1,*) inino1(n)
  enddo
  do n=1, nt2
    write(2,*) inino2(n)
  enddo
  close(1)
  close(2)
  open(3,file='../nlreg_tatl5_pic.asc')
  open(4,file='../nlreg_tatl5_co2.asc')
  do jj=1, 5
    do n=1, nt1
      write(3,*) tatl1(n,jj)
    enddo
    do n=1, nt2
      write(4,*) tatl2(n,jj)
    enddo
  enddo
  close(3)
  close(4)


  stop

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


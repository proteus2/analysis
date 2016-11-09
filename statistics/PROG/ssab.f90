PROGRAM wavelet_analysis_b

  use ssapack
  use waveletpack
  use netcdfio

  implicit none

! Data --------------------------------------------------------
  integer, parameter :: nxo = 320, nyo = 320
  integer, parameter :: nt1all = 1200, nt2all = 960
  integer, parameter :: ts1 = 1, ts2 = 1, nt1 = nt1all, nt2 = nt2all
  real,    parameter :: dt = 1./12.
! SSA para. ---------------------------------------------------
  integer, parameter :: maxlag = 84, nout = 20
! Wavelet para. -----------------------------------------------
  integer, parameter :: ssub = 8
  real, parameter    :: s0 = dt !*2.
!--------------------------------------------------------------
  integer, parameter :: nscale1 = 1+int(log(nt1*dt/s0)/log(2.)*ssub+0.999)
  integer, parameter :: nscale2 = 1+int(log(nt2*dt/s0)/log(2.)*ssub+0.999)

  integer :: i,j,n, ncid

  real, dimension(nxo,nyo,nt1) :: to1
  real, dimension(nxo,nyo,nt2) :: to2
  real, dimension(nt1)         :: inino1, recon1
  real, dimension(nt2)         :: inino2, recon2
  real, dimension(nt1,nout)    :: rc1
  real, dimension(nt2,nout)    :: rc2

  real    :: lono(nxo), lato(nyo), t1(nt1), t2(nt2), dsinphio(nyo)

  integer :: tmpi
  real*8  :: tmp1(nxo), tmp2(nyo), tmp3(2,nyo)
  real    :: pi, tmpr1, tmpr2, tmpr3, tmp4(nt1), tmp5(nt2)

  character*128 :: fdir(2), fn_to1, fn_to2, ofname


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
    t1(n) = real(n-1)*dt
  enddo
  do n=1, nt2
    t2(n) = real(n-1)*dt
  enddo

! variables
  call opennc(fn_to1,ncid)
  call geta3d(ncid,'tos',1,nxo,1,nyo,ts1,nt1,to1)
  call closenc(ncid)
  call opennc(fn_to2,ncid)
  call geta3d(ncid,'tos',1,nxo,1,nyo,ts2,nt2,to2)
  call closenc(ncid)


! Nino3.4 index --------------------------------------------------------
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


  call ssanc('../ssabpic',nt1,dt,t1,inino1,1.e20,maxlag,nout)
  call ssanc('../ssabco2',nt2,dt,t2,inino2,1.e20,maxlag,nout)

  call opennc('../ssabpic_compo.nc',ncid)
  call get2d(ncid,'RC',nt1,nout,rc1)
  call closenc(ncid)
  call opennc('../ssabco2_compo.nc',ncid)
  call get2d(ncid,'RC',nt2,nout,rc2)
  call closenc(ncid)

  recon1 = 0.
  do i=1, 5
    recon1(:) = recon1(:) + rc1(:,i)
  enddo
  recon2 = 0. 
  do i=1, 5 
    recon2(:) = recon2(:) + rc2(:,i)
  enddo

  call wlncpack('../ssabwlpic',nt1,dt,t1,recon1,0,nscale1,s0,ssub,0,6.,1,&
           0,(/2.,7.9/),0,0,0.05,0)
  call wlncpack('../ssabwlco2',nt2,dt,t2,recon2,0,nscale2,s0,ssub,0,6.,1,&
           0,(/2.,7.9/),0,0,0.05,0)


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


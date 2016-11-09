program spectral_analysis_2a

  use netcdfio
  use specanal
  use random

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nx = 1000, maxlag = 500
  real,    parameter :: peri1 = 200., peri2 =  80., dx = 1., phase =135.
  real,    parameter :: amp1 = 2., amp2 = 3., gauss_sd = 0.
!--------------------------------------------------------------
  integer, parameter   :: nps = maxlag+1
  real, dimension(nx)  :: x1, x2, x3, t, wn
  real, dimension(nps) :: freq, ps11, ps12, ps21, ps22
  real, dimension(nps) :: csamp1, csphs1, coh1, csamp2, csphs2, coh2
  real, dimension(nps*2-1) :: corr1, corr2, lag
  integer :: i,j
  real    :: pi


  pi = acos(-1.)

  do i=1, nx
    t(i)  = real(i) * dx
    x1(i) = amp1*sin(2.*pi/peri1*t(i))
    x2(i) = amp2*cos(2.*pi/peri2*t(i))
    x3(i) = amp2*sin(2.*pi/peri1*t(i) + pi/180.*phase)
  enddo

  do i=1, nx
    wn(i) = random_normal()
  enddo
  wn = gauss_sd * wn

  x1 = x1 + wn
  x2 = x2 + wn

  do i=1, nps*2-1
    lag(i) = real(i)-maxlag-1
  enddo

  call crosss(x1,x2,dx,nx,maxlag,corr1,freq,ps11,ps12,csamp1,csphs1,coh1)
  call crosss(x1,x3,dx,nx,maxlag,corr2,freq,ps21,ps22,csamp2,csphs2,coh2)


  call out1d2('../sa2a1.nc',1,(/'COR'/),(/corr1/),'lag',nps*2-1,lag, &
              3,(/'COH','PHASE','AMP'/),(/coh1,csphs1,csamp1/), &
              'freq',nps,freq,'Idealized exp.')
  call out1d2('../sa2a2.nc',1,(/'COR'/),(/corr2/),'lag',nps*2-1,lag, &
              3,(/'COH','PHASE','AMP'/),(/coh2,csphs2,csamp2/), &
              'freq',nps,freq,'Idealized exp.')


end


 program k_zt

  use pwrspd
  use netcdfio

  implicit none

  integer, parameter :: nxp = 800
  integer, parameter :: nkk = nxp/2+1, k1 = 20
  integer, parameter :: eval = 1
  real,    parameter :: dx = 1., amp1 = 2., amp2 = 2.
  integer :: i
  real, dimension(nxp)            :: w, x, phs, w1, w2
  real, dimension(nkk)            :: kk, d
  real    :: pi
  character*128 :: fname


  pi = acos(-1.)

  w1 = 0.  ;  w2 = 0.
  do i=1, nxp
    phs(i) = 2.*pi * k1/real(nxp) * i
    w1(i) = amp1*sin(phs(i)) * exp(-(i-400.)**2/10000.)
    w2(i) = amp2*cos(phs(i)+pi/180.*30) * exp(-(i-350.)**2/5000.)
  enddo
  w = w1 !+ w2

  do i=1, nxp
    x(i) = i*dx
  enddo

  kk(:)=0.
  d(:)=0.
  call psd1d(nxp,0,x,dx,w,kk,d,0,eval)


  write(fname,'(a)') '../test.nc'
  call out1d2(trim(fname),1,(/'TS'/),w,'x',nxp,x, &
              1,(/'PSD'/),d,'k',nkk,kk,'PSD1D k')


  stop

 end

 program k_zt

  use fft
  use netcdfio

  implicit none

  integer, parameter :: nxp = 500
  integer, parameter :: nkk = nxp/2+1, k1 = 10
  integer, parameter :: eval = 1
  real,    parameter :: dx = 1., amp1 = 2., amp2 = 3.
  integer :: i
  real, dimension(nxp)            :: w, x, phs, w1, u1, u1t
  real, dimension(nkk)            :: kk, d
  real    :: pi, wei(nxp)
  double complex :: coef(nxp)
  character*128 :: fname


  pi = acos(-1.)

  do i=1, nxp
    phs(i) = 2.*pi * k1/real(nxp) * i
    w1(i) = amp1*cos(phs(i)) * exp(-(i-100.)**2/5000.)
    u1(i) = amp2*cos(phs(i)) * exp(-(i-100.)**2/5000.)
  enddo


  wei = 1.
  call fft1df(nxp,u1,coef)
  coef(:) = (0.d0,1.d0)*coef(:) * wei(:)
  coef(nxp/2+2:) = -coef(nxp/2+2:)
  call fft1db(nxp,coef,u1t)

!  cf(ij,:) = sf(ij,:) + iii*w2t(:)

  do i=1, nxp
    x(i) = i*dx
  enddo

  kk(:)=0.
  d(:)=0.
  call psd1d(nxp,0,x,dx,w,kk,d,0,eval)


  write(fname,'(a)') '../test.nc'
  call out1d2(trim(fname),2,(/'w','u'/),(/w1,u1t/),'x',nxp,x, &
              1,(/'PSD'/),d,'k',nkk,kk,'PSD1D k')


  stop

 end

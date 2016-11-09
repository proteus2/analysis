program spectral_analysis_1a

  use netcdfio
  use specanal
  use random

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nx = 1000, maxlag = 500
  real,    parameter :: peri1 = 300., peri2 =  50., dx = 1.
  real,    parameter :: amp1 = 2., amp2 = 3., gauss_sd = 2.
!--------------------------------------------------------------
  integer, parameter   :: nps = maxlag+1
  real, dimension(nx)  :: x1, x2, t, wn
  real, dimension(nps) :: freq
  real, dimension(nps) :: ps1, ps2
  integer :: i,j
  real    :: pi, sum1, sum2, sum3, sum4


  pi = acos(-1.)

  do i=1, nx
    t(i)  = real(i) * dx
    x1(i) = amp1*sin(2.*pi/peri1*t(i)) + amp2*cos(2.*pi/peri2*t(i))
    x2(i) = amp1*sin(2.*pi/peri1*t(i))*cos(2.*pi/peri2*t(i))
  enddo

  do i=1, nx
    wn(i) = random_normal()
  enddo
  wn = gauss_sd * wn

  x1 = x1 + wn
  x2 = x2 + wn

  call fps(x1,dx,nx,maxlag,freq,ps1,0.)
  call fps(x2,dx,nx,maxlag,freq,ps2,0.)

  sum1 = 0.  ;  sum2 = 0.  ;  sum3 = 0.  ;  sum4 = 0.
  do i=1, nx
    sum1 = sum1 + x1(i)**2/nx
    sum2 = sum2 + x2(i)**2/nx
  enddo
  do i=1, nps
    sum3 = sum3 + ps1(i)
    sum4 = sum4 + ps2(i)
  enddo
  sum3 = sum3 * 2.*pi*freq(2)
  sum4 = sum4 * 2.*pi*freq(2)

  print*, 'var   :', sum1, sum2
  print*, 'power :', sum3, sum4
  print*, 'raio  :', abs(sum3-sum1)/sum1, abs(sum4-sum2)/sum2

  call out1d2('../sa1a.nc',2,(/'X1','X2'/),(/x1,x2/),'t',nx,t, &
              2,(/'PS1','PS2'/),(/ps1,ps2/),'freq',nps,freq, &
              'Idealized exp.')


end


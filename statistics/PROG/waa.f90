program wavelet_analysis_a

  use waveletpack
  use random

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nx = 1000
  real,    parameter :: peri1 = 200., peri2 = 50., dx = 1.
  real,    parameter :: amp1 = 2., amp2 = 3., gauss_sd = 0.
! Wavelet para. -----------------------------------------------
  integer, parameter :: ssub = 8
  real, parameter    :: s0 = dx
!--------------------------------------------------------------
  integer, parameter :: nscale = 1+int(log(nx*dx/s0)/log(2.)*ssub+0.999)

  real, dimension(nx)  :: x1, x2, x3, t, wn
  integer :: i,j
  real    :: pi


  pi = acos(-1.)

  do i=1, nx
    t(i)  = real(i) * dx
    x1(i) = amp1*sin(2.*pi/peri1*t(i)) + amp2*sin(2.*pi/peri2*t(i))
    x2(i) = amp1*sin((2.*pi/peri1*t(i))**2)*sin(2.*pi/peri2*t(i))
  enddo

  do i=1, nx
    wn(i) = random_normal()
  enddo
  wn = gauss_sd * wn

  x1 = x1 + wn

  call wlncpack('../waa1',nx,dx,t,x1,0,nscale,s0,ssub,0,6.,1, &
           0,(/2.,7.9/),0,0,0.05,0)

  call wlncpack('../waa1',nx,dx,t,x2,0,nscale,s0,ssub,0,6.,0, &
           0,(/2.,7.9/),0,0,0.05,0)


end


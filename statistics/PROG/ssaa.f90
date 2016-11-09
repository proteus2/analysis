program ssa_a

  use ssapack
  use netcdfio
  use random

  implicit none

!--------------------------------------------------------------
  integer, parameter :: nx = 1000
  real,    parameter :: peri1 = 200., peri2 = 50., dx = 1.
  real,    parameter :: amp1 = 2., amp2 = 3., gauss_sd = 0.
! SSA para. ---------------------------------------------------
  integer, parameter :: maxlag = 200, nout = 5
!--------------------------------------------------------------
  real, dimension(nx)  :: x1, t, wn
  integer :: i,j
  real    :: pi


  pi = acos(-1.)

  do i=1, nx
    t(i)  = real(i) * dx
    x1(i) = amp1*sin(2.*pi/peri1*t(i)) + amp2*sin(2.*pi/peri2*t(i))
  enddo

  do i=1, nx
    wn(i) = random_normal()
  enddo
  wn = gauss_sd * wn

  x1 = x1 + wn

  call ssanc('../ssaa',nx,dx,t,x1,1.e32,maxlag,nout)


  stop

end


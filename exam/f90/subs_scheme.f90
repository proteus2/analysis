
! SUBROUTINES OF THREE SCHEMES BY WHICH THE DIFFUSION PROBLEM IS SOLVED.

SUBROUTINE diff_forward(nx,dx,dt,alpha,temp)

  implicit none

  integer,                intent(in)    ::  nx
  real,                   intent(in)    ::  dx, dt, alpha
  real,    dimension(nx), intent(inout) ::  temp

  integer             ::  i
  real                ::  coef
  real, dimension(nx) ::  dtemp

  coef = alpha*dt/dx/dx

  do i=2, nx-1
    dtemp(i) = coef*(temp(i-1)-2.*temp(i)+temp(i+1))
  enddo
  temp(2:nx-1) = temp(2:nx-1) + dtemp(2:nx-1)
  temp( 1) = 0.
  temp(nx) = 0.

  RETURN

END subroutine diff_forward


SUBROUTINE diff_backward(nx,dx,dt,alpha,temp)

  implicit none

  integer,             intent(in)    ::  nx
  real,                intent(in)    ::  dx, dt, alpha
  real, dimension(nx), intent(inout) ::  temp

  real                ::  coef
  real, dimension(nx) ::  a, b, c, d

  coef = alpha*dt/dx/dx

  a(:) = -coef
  b(:) = 1.+2.*coef
  c(:) = -coef
  d(:) = temp(:)

  b( 1) = 1.  ;  c( 1) = 0.
  b(nx) = 1.  ;  a(nx) = 0.

  call tridag(nx,a,b,c,d, temp)

  RETURN

END subroutine diff_backward


SUBROUTINE diff_crni(nx,dx,dt,alpha,temp)

  implicit none

  integer,             intent(in)    ::  nx
  real,                intent(in)    ::  dx, dt, alpha
  real, dimension(nx), intent(inout) ::  temp

  integer             ::  i
  real                ::  coef
  real, dimension(nx) ::  a, b, c, d

  coef = alpha*dt/dx/dx

  a(:) = -0.5*coef
  b(:) = 1.+coef
  c(:) = -0.5*coef
  do i=2, nx-1
    d(i) = (1.-coef)*temp(i) + 0.5*coef*(temp(i-1)+temp(i+1))
  enddo

  b( 1) = 1.  ;  c( 1) = 0.  ;  d( 1) = 0.
  b(nx) = 1.  ;  a(nx) = 0.  ;  d(nx) = 0.

  call tridag(nx,a,b,c,d, temp)

  RETURN

END subroutine diff_crni


SUBROUTINE tridag(n,a,b,c,d,x)

  implicit none

  integer,            intent(in)  ::  n
  real, dimension(n), intent(in)  ::  a, b, c, d
  real, dimension(n), intent(out) ::  x

  integer ::  i
  real    ::  bet, gam(n)

  bet = b(1)
  x(1) = d(1)/bet
  do i=2, n
    gam(i) = c(i-1)/bet
    bet = b(i)-a(i)*gam(i)
    if (bet == 0.)  pause ' tridag failed !!'
    x(i) = (d(i)-a(i)*x(i-1))/bet
  enddo
  do i=n-1, 1, -1
    x(i) = x(i) - gam(i+1)*x(i+1)
  enddo

  RETURN

END subroutine tridag


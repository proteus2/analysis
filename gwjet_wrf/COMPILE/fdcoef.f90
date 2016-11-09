SUBROUTINE fdcoef(m,n,z,x,coef)
!
! Calculate coefficients of finite difference formula for arbitrary spaced
!   grids, e.g., for the 2nd-order approximation (n = 2) of the 1st 
!   derivative (m = 1) of F :
!   dF/dx = sum( coef(0:n)*F(0:n) )
!
! Reference : Fornberg (1998, SIAM Rev.)
!
! Input parameters
!   z       location where approximations are to be accurate
!   x(0:n)  grid point locations
!   n       one less than total number of grid points
!   m       highest derivative for which weights are sought
! Output parameter
!   coef(0:n)  weights at grid locations x(0:n) for derivative of order m
!
  implicit none

  integer,              intent(in)  ::  m, n
  real,                 intent(in)  ::  z
  real, dimension(0:n), intent(in)  ::  x
  real, dimension(0:n), intent(out) ::  coef

  integer          ::  i,j,k, mn
  double precision ::  c1, c2, c3, c4, c5
  double precision, dimension(0:n,0:m) ::  c

  c1 = 1.d0
  c4 = dble(x(0) - z)
  c(:,:) = 0.d0
  c(0,0) = 1.d0
  do i=1, n
    mn = min(i,m)
    c2 = 1.d0
    c5 = c4
    c4 = dble(x(i) - z)
    do j=0, i-1
      c3 = dble(x(i) - x(j))
      c2 = c2*c3
      if (j == i-1) then
        do k=mn, 1, -1
          c(i,k) = c1*(k*c(i-1,k-1) - c5*c(i-1,k))/c2
        enddo
        c(i,0) = -c1*c5*c(i-1,0)/c2
      end if
      do k=mn, 1, -1
        c(j,k) = (c4*c(j,k) - k*c(j,k-1))/c3
      enddo
      c(j,0) = c4*c(j,0)/c3
    enddo
    c1 = c2
  enddo

  coef(:) = real(c(:,m))

END subroutine fdcoef


SUBROUTINE fdcoef_1d2o(x,coef)
!
! Calculate coefficients of finite difference formula for arbitrary spaced
!   grids: for the 2nd-order approximation of the 1st derivative of F at x(1)
!   dF/dx = sum( coef(0:2)*F(0:2) )
!
! Reference : Fornberg (1998, SIAM Rev.)
!
! Input parameters
!   x(0:2)  grid point locations
! Output parameter
!   coef(0:2)  weights at grid locations x(0:2)
!
  implicit none

  real, dimension(0:2), intent(in)  ::  x
  real, dimension(0:2), intent(out) ::  coef

  double precision ::  l, l1, l2

  l  = dble(x(2) - x(0))
  l1 = dble(x(1) - x(0))
  l2 = dble(x(1) - x(2))

  coef(:) = real( (/ l2/(l*l1),          &
                     1.d0/l1 + 1.d0/l2,  &
                     l1/(l*l2)*(-1.d0) /) )

END subroutine fdcoef_1d2o


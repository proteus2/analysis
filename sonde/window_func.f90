SUBROUTINE welch_window(n,w)

  implicit none

  integer,            intent(in)  ::  n
  real, dimension(n), intent(out) ::  w

  integer ::  i

  do i=1, n
    w(i) = 1. - (float(i*2-1-n)/(float(n-1)))**2
  enddo

END subroutine welch_window
 

SUBROUTINE poly_func(x,arr)

  use nrtype

  implicit none

  real(sp),               intent(in ) ::  x
  real(sp), dimension(:), intent(out) ::  arr

! DO NOT MODIFY THE ABOVE

  integer ::  im

  arr(1) = 1.
  do im=2, size(arr)
    arr(im) = arr(im-1)*x
  enddo

END subroutine poly_func
 

SUBROUTINE cub_spl_intp0(nxi,xi,yi,nxo,xo,yo)

  use nr,  only: spline, splint

  implicit none

  integer, intent(in)  :: nxi, nxo
  real,    intent(in)  :: xi(nxi), yi(nxi), xo(nxo)
  real,    intent(out) :: yo(nxo)

  integer :: i,j,k
  real    :: work(nxi)

  call spline(xi,yi,1.e32,1.e32,work)
  do k=1, nxo
    yo(k) = splint(xi,yi,work,xo(k))
  enddo

END subroutine cub_spl_intp0
 

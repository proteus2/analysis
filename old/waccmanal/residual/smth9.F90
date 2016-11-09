!2345
    subroutine smth9(nx    ,ny    ,p     ,q     ,src   )

    implicit none

    integer, intent(in) :: nx    ,ny
    real,    intent(in) :: p     ,q
    real, dimension(nx,ny), intent(inout)  :: src

    integer :: i     ,j
    real, dimension(nx,ny) :: tmp

    tmp(1:nx,1:ny) = 0.0

    do j=2,ny-1
      do i=2,nx-1
        tmp(i,j) = src(i,j) + (p/4.)*(src(i-1,j  )+src(i  ,j-1)+src(i+1,j  )+src(i  ,j+1)-4.*src(i,j)) +  &
                              (q/4.)*(src(i-1,j+1)+src(i-1,j-1)+src(i+1,j-1)+src(i+1,j+1)-4.*src(i,j))
      end do
    end do

    do j=2,ny-1
      do i=2,nx-1
        src(i,j) = tmp(i,j)
      end do
    end do
    
    return
    end subroutine smth9

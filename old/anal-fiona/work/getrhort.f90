 program get_standard_rho_sqrt

  implicit none

  integer :: i
  real    :: rho(472)
  character*128 :: temp

  open (1,file='../../ideal/work/id0/rhobar')
  open (10,file='rhoall')

  read (1,*) rho

  do i=1,472
    write (10,*) sqrt(rho(i)/rho(52))
  enddo
  close(1)
  close(10)

  end

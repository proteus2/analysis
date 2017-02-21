
! SUBROUTINES FOR I/O.

SUBROUTINE get_init(i_unit,filename,nx,f_init)

  implicit none

  integer,             intent(in)  ::  i_unit, nx
  character(len=*),    intent(in)  ::  filename
  real, dimension(nx), intent(out) ::  f_init

  open(i_unit, file = filename, status = 'old')
  read(i_unit,*) f_init
  close(i_unit)
  f_init( 1) = 0.
  f_init(nx) = 0.

  RETURN

END subroutine get_init


SUBROUTINE check_para(opt_method,alpha,lx,nx)

  implicit none

  integer, intent(in) ::  opt_method, nx
  real,    intent(in) ::  alpha, lx

  write(6,*)
  write(6,'(a)'        ) '--- INPUT PARAMETERS -----------------------'

  SELECT case ( opt_method )
    case ( 1 )
      write(6,'(a)'    ) ' Chosen by the forward Euler scheme'
    case ( 2 )
      write(6,'(a)'    ) ' Chosen by the backward Euler scheme'
    case ( 3 )
      write(6,'(a)'    ) ' Chosen by the Crank-Nicolson scheme'
  END select

  write(6,'(a,f10.2,a)') ' Diffusion coeff.     :', alpha, ' m^2/s'
  write(6,'(a,f10.2,a)') ' Length of the domain :', lx, ' m'
  write(6,'(a,i7.0)'   ) ' Number of grids      :', nx
  write(6,'(a)'        ) '--------------------------------------------'
  RETURN

END subroutine check_para


SUBROUTINE dump_field(opt_dump,i_unit,filehead,time,nx,x,temp)

  implicit none

  integer,             intent(in) ::  opt_dump, i_unit, nx
  real,                intent(in) ::  time
  real, dimension(nx), intent(in) ::  x, temp
  character(len=*),    intent(in) ::  filehead

  integer            ::  i
  character(len=128) ::  filename

  write(6,*)
  write(6,'(a,f10.2,a)') ' Dump time : ', time, ' s'

  SELECT case ( opt_dump )
    case ( 1 )
      do i=1, nx
        write(6,'(f10.2,a,f10.4,a)') x(i),' m : ', temp(i),' K'
      enddo
    case ( 2 )
      write(filename,'(a,i6.6)') trim(filehead)//'.', int(time)
      open(i_unit, file = filename)
      write(i_unit,'(a)') '     x (m)   temp (K) '
      write(i_unit,'(a)') ' ---------------------'
      do i=1, nx
        write(i_unit,'(f10.2,f10.4)') x(i), temp(i)
      enddo
      close(i_unit)
  END select

  RETURN

END subroutine dump_field


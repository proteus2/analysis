
  implicit none

!========================================================= PARAMETERS ==

  integer ::  opt_method, nx
  real    ::  alpha, lx, dx, dt, t_end

  namelist /PARA_EQN/ opt_method, nx, alpha, lx, dt, t_end

  integer            ::  opt_dump
  real               ::  dt_dump
  character(len=128) ::  file_in, file_out

  namelist /PARA_IO/ opt_dump, dt_dump, file_in, file_out


!====================================================== GENERAL VAR.S ==

  real, dimension(:), allocatable ::  x, temp
  real    ::  time
  integer ::  i, n, nt

!=======================================================================


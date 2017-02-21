
! MAIN PROGRAM FOR SOLVING THE HEAT EQUATION.

PROGRAM heat_main

  include 'decl.inc'

! READ INPUT PARAMETERS ================================================

  open(10, file = 'para.in', status = 'old')
  read(10, nml = PARA_EQN)
  read(10, nml = PARA_IO )
  close(10)
  call check_para(opt_method,alpha,lx,nx)

! READ THE INITIAL FIELD ===============================================

  allocate( temp(nx), x(nx) )
  dx = lx/float(nx-1)
  do i=1, nx
    x(i) = dx*float(i-1)
  enddo
  call get_init(20,file_in,nx, temp(:))

! INTEGRATION AND DUMP =================================================

  time = 0.
  nt = int(t_end/dt)

  NSTEP:  DO n=1, nt

    ! integration ------------------------------------------
    time = time + dt

    SELECT case ( opt_method )
      case ( 1 )
        call diff_forward(nx,dx,dt,alpha, temp)
      case ( 2 )
        call diff_backward(nx,dx,dt,alpha, temp)
      case ( 3 )
        call diff_crni(nx,dx,dt,alpha, temp)
      case default
        print*, ' Reset opt_method ( 1 - 3 ) !!'
        EXIT
    END select

    ! dump -------------------------------------------------
    if (mod(time,dt_dump) == 0)                                         &
       call dump_field(opt_dump,100,file_out,time,nx,x,temp)

  ENDDO  nstep

! END PROGRAM ==========================================================

  STOP

END program heat_main


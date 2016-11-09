PROGRAM test

  use integ

  implicit none

  integer ::  nx = 10
  real    ::  dx = 1., x0 = 4.0
  real    ::  h = -1., a = 30., b = -100.

  integer             ::  i
  real, dimension(nx) ::  x, f, intg1, intg2, intgr


  do i=1, nx
    x(i) = dx*(nx-i)
  enddo

  f = (a*x+b)*exp(x/h)
  print*
  print*, 'f(x)'
  write(6,'(5f15.6)') f

  intgr = h*((a*(x-h)+b)*exp(x/h)-(a*(x0-h)+b)*exp(x0/h))
  print*
  print*, 'integral (analytic)'
  write(6,'(5f15.6)') intgr

  call integ1d_exp((/nx,1,1,1/),1,x,f,h*0.9,x0,intg1)
  print*
  print*, 'integ1d_exp'
  write(6,'(5f15.6)') intg1
  print*
  print*, 'integ1d_exp err.'
  write(6,'(5f15.6)') abs((intg1-intgr)/intgr)

  call integ1d((/nx,1,1,1/),1,x,f,x0,intg2)
  print*
  print*, 'integ1d'
  write(6,'(5f15.6)') intg2
  print*
  print*, 'integ1d err.'
  write(6,'(5f15.6)') abs((intg2-intgr)/intgr)

  do i=nx, 2, -1
    intgr(i) = h*((a*(x(i)-h)+b)*exp(x(i)/h)-(a*(x(i-1)-h)+b)*exp(x(i-1)/h))
    call integ1d_exp((/2,1,1,1/),1,x(i-1:i),f(i-1:i),h*0.9,x(i-1),intg1(i-1:i))
    call integ1d((/2,1,1,1/),1,x(i-1:i),f(i-1:i),x(i-1),intg2(i-1:i))
  enddo

  print*
  print*
  print*, 'integral (analytic)'
  write(6,'(5f15.6)') intgr

  print*
  print*, 'integ1d_exp err.'
  write(6,'(5f15.6)') abs((intg1-intgr)/intgr)

  print*
  print*, 'integ1d err.'
  write(6,'(5f15.6)') abs((intg2-intgr)/intgr)

  print*
  print*, 'DIFF'
  write(6,'(5f15.6)') (abs(intg2-intgr)-abs(intg1-intgr))/abs(intgr)


  STOP

END program test


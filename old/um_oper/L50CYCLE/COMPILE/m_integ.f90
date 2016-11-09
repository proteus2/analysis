MODULE integ


  CONTAINS

!=======================================================================

SUBROUTINE integ1d(nd,id,x,f,x0,intg)

! integration from x0 to x

  implicit none

  integer,                                  intent(in) ::  nd(4), id
  real,                                     intent(in) ::  x0
  real, dimension(nd(id)),                  intent(in) ::  x
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  intg

  integer ::  i,j,k,n, i0
  integer ::  nd1, nd2, nd3, nd4, nx
  real    ::  f0, dx(nd(id)), dx05(nd(id)), d2i0, d2i0p1, temp


  nd1 = nd(1)  ;  nd2 = nd(2)  ;  nd3 = nd(3)  ;  nd4 = nd(4)
  nx  = nd(id)

  call intg_set(nx,x,x0, i0,dx)

  dx05(:) = 0.5*dx(:)

  d2i0   = abs(dx(i0  ))
  d2i0p1 = abs(dx(i0+1))
  temp   = d2i0 + d2i0p1
  d2i0   = d2i0   / temp
  d2i0p1 = d2i0p1 / temp


  SELECT case ( id )

  case ( 1 )

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
      f0 = d2i0p1*f(i0,j,k,n)+d2i0*f(i0+1,j,k,n)
      intg(i0:i0+1,j,k,n) = (f(i0:i0+1,j,k,n) + f0)*dx05(i0:i0+1)
      do i=i0+2, nx
        intg(i,j,k,n) = intg(i-1,j,k,n) + dx05(i)* &
                        (f(i-1,j,k,n)+f(i,j,k,n))
      enddo
      do i=i0-1, 1, -1
        intg(i,j,k,n) = intg(i+1,j,k,n) + dx05(i)* &
                        (f(i+1,j,k,n)+f(i,j,k,n))
      enddo
    enddo
    enddo
    enddo

  case ( 2 )

    do n=1, nd4
    do k=1, nd3
      do i=1, nd1
        f0 = d2i0p1*f(i,i0,k,n)+d2i0*f(i,i0+1,k,n)
        intg(i,i0:i0+1,k,n) = (f(i,i0:i0+1,k,n) + f0)*dx05(i0:i0+1)
      enddo
      do j=i0+2, nx
        intg(:,j,k,n) = intg(:,j-1,k,n) + dx05(j)* &
                        (f(:,j-1,k,n)+f(:,j,k,n))
      enddo
      do j=i0-1, 1, -1
        intg(:,j,k,n) = intg(:,j+1,k,n) + dx05(j)* &
                        (f(:,j+1,k,n)+f(:,j,k,n))
      enddo
    enddo
    enddo

  case ( 3 )

    do n=1, nd4
      do j=1, nd2
      do i=1, nd1
        f0 = d2i0p1*f(i,j,i0,n)+d2i0*f(i,j,i0+1,n)
        intg(i,j,i0:i0+1,n) = (f(i,j,i0:i0+1,n) + f0)*dx05(i0:i0+1)
      enddo
      enddo
      do k=i0+2, nx
        intg(:,:,k,n) = intg(:,:,k-1,n) + dx05(k)* &
                        (f(:,:,k-1,n)+f(:,:,k,n))
      enddo
      do k=i0-1, 1, -1
        intg(:,:,k,n) = intg(:,:,k+1,n) + dx05(k)* &
                        (f(:,:,k+1,n)+f(:,:,k,n))
      enddo
    enddo

  case ( 4 )

    do k=1, nd3
    do j=1, nd2
    do i=1, nd1
      f0 = d2i0p1*f(i,j,k,i0)+d2i0*f(i,j,k,i0+1)
      intg(i,j,k,i0:i0+1) = (f(i,j,k,i0:i0+1) + f0)*dx05(i0:i0+1)
    enddo
    enddo
    enddo
    do n=i0+2, nx
      intg(:,:,:,n) = intg(:,:,:,n-1) + dx05(n)* &
                      (f(:,:,:,n-1)+f(:,:,:,n))
    enddo
    do n=i0-1, 1, -1
      intg(:,:,:,n) = intg(:,:,:,n+1) + dx05(n)* &
                      (f(:,:,:,n+1)+f(:,:,:,n))
    enddo

  END select

  RETURN

END subroutine integ1d

!=======================================================================

SUBROUTINE integ1d_exp(nd,id,x,f,x_scale,x0,intg)

! integration from x0 to x
! for the case that y is assumed linear between each two adjacent grids,
!   when f = y*exp(x/x_scale)

  implicit none

  integer,                                  intent(in) ::  nd(4), id
  real,                                     intent(in) ::  x0, x_scale
  real, dimension(nd(id)),                  intent(in) ::  x
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  intg

  integer                                  ::  i,j,k,n, i0
  integer                                  ::  nd1, nd2, nd3, nd4, nx
  real                                     ::  f0, fdexp0, exp_xds0
  real                                     ::  d2i0, d2i0p1, temp
  real, dimension(nd(id))                  ::  dx, exp_xds
  real, dimension(nd(1),nd(2),nd(3),nd(4)) ::  fdexp
  logical                                  ::  l_x0i0, l_x0i0p1


  nd1 = nd(1)  ;  nd2 = nd(2)  ;  nd3 = nd(3)  ;  nd4 = nd(4)
  nx  = nd(id)

  call intg_set(nx,x,x0, i0,dx)

  d2i0   = abs(dx(i0  ))
  d2i0p1 = abs(dx(i0+1))
  temp   = d2i0 + d2i0p1
  d2i0   = d2i0   / temp
  d2i0p1 = d2i0p1 / temp

  exp_xds(:) = exp(x(:)/x_scale)
  exp_xds0 = exp(x0/x_scale)

  ! for numerical purpose (avoid NaN)
  if ( dx(i0) == 0. ) then
    l_x0i0 = .TRUE.
    dx(i0) = minval(abs(dx(:i0-1)))*1.e-2
  end if
  if ( dx(i0+1) == 0. ) then
    l_x0i0p1 = .TRUE.
    dx(i0+1) = minval(abs(dx(i0+2:)))*1.e-2
  end if


  SELECT case ( id )

  case ( 1 )

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
      fdexp(:,j,k,n) = f(:,j,k,n)/exp_xds(:)
    enddo
    enddo
    enddo

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
      fdexp0 = d2i0p1*fdexp(i0,j,k,n)+d2i0*fdexp(i0+1,j,k,n)
      f0 = fdexp0*exp_xds0
      intg(i0:i0+1,j,k,n) = x_scale*( (f(i0:i0+1,j,k,n)-f0) -        &
                          (fdexp(i0:i0+1,j,k,n)-fdexp0)/dx(i0:i0+1)* &
                          (exp_xds(i0:i0+1)-exp_xds0)*x_scale )
      do i=i0+2, nx
        intg(i,j,k,n) = intg(i-1,j,k,n) + x_scale*                 &
                        ( (f(i,j,k,n)-f(i-1,j,k,n)) -              &
                          (fdexp(i,j,k,n)-fdexp(i-1,j,k,n))/dx(i)* &
                          (exp_xds(i)-exp_xds(i-1))*x_scale )
      enddo
      do i=i0-1, 1, -1
        intg(i,j,k,n) = intg(i+1,j,k,n) + x_scale*                 &
                        ( (f(i,j,k,n)-f(i+1,j,k,n)) -              &
                          (fdexp(i,j,k,n)-fdexp(i+1,j,k,n))/dx(i)* &
                          (exp_xds(i)-exp_xds(i+1))*x_scale )
      enddo
    enddo
    enddo
    enddo

    if ( l_x0i0   )  intg(i0  ,:,:,:) = 0.
    if ( l_x0i0p1 )  intg(i0+1,:,:,:) = 0.

  case ( 2 )

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
      fdexp(:,j,k,n) = f(:,j,k,n)/exp_xds(j)
    enddo
    enddo
    enddo

    do n=1, nd4
    do k=1, nd3
      do i=1, nd1
        fdexp0 = d2i0p1*fdexp(i,i0,k,n)+d2i0*fdexp(i,i0+1,k,n)
        f0 = fdexp0*exp_xds0
        intg(i,i0:i0+1,k,n) = x_scale*( (f(i,i0:i0+1,k,n)-f0) -        &
                            (fdexp(i,i0:i0+1,k,n)-fdexp0)/dx(i0:i0+1)* &
                            (exp_xds(i0:i0+1)-exp_xds0)*x_scale )
      enddo
      do j=i0+2, nx
        intg(:,j,k,n) = intg(:,j-1,k,n) + x_scale*                 &
                        ( (f(:,j,k,n)-f(:,j-1,k,n)) -              &
                          (fdexp(:,j,k,n)-fdexp(:,j-1,k,n))/dx(j)* &
                          (exp_xds(j)-exp_xds(j-1))*x_scale )
      enddo
      do j=i0-1, 1, -1
        intg(:,j,k,n) = intg(:,j+1,k,n) + x_scale*                 &
                        ( (f(:,j,k,n)-f(:,j+1,k,n)) -              &
                          (fdexp(:,j,k,n)-fdexp(:,j+1,k,n))/dx(j)* &
                          (exp_xds(j)-exp_xds(j+1))*x_scale )
      enddo
    enddo
    enddo

    if ( l_x0i0   )  intg(:,i0  ,:,:) = 0.
    if ( l_x0i0p1 )  intg(:,i0+1,:,:) = 0.

  case ( 3 )

    do n=1, nd4
    do k=1, nd3
      fdexp(:,:,k,n) = f(:,:,k,n)/exp_xds(k)
    enddo
    enddo

    do n=1, nd4
      do j=1, nd2
      do i=1, nd1
        fdexp0 = d2i0p1*fdexp(i,j,i0,n)+d2i0*fdexp(i,j,i0+1,n)
        f0 = fdexp0*exp_xds0
        intg(i,j,i0:i0+1,n) = x_scale*( (f(i,j,i0:i0+1,n)-f0) -        &
                            (fdexp(i,j,i0:i0+1,n)-fdexp0)/dx(i0:i0+1)* &
                            (exp_xds(i0:i0+1)-exp_xds0)*x_scale )
      enddo
      enddo
      do k=i0+2, nx
        intg(:,:,k,n) = intg(:,:,k-1,n) + x_scale*                 &
                        ( (f(:,:,k,n)-f(:,:,k-1,n)) -              &
                          (fdexp(:,:,k,n)-fdexp(:,:,k-1,n))/dx(k)* &
                          (exp_xds(k)-exp_xds(k-1))*x_scale )
      enddo
      do k=i0-1, 1, -1
        intg(:,:,k,n) = intg(:,:,k+1,n) + x_scale*                 &
                        ( (f(:,:,k,n)-f(:,:,k+1,n)) -              &
                          (fdexp(:,:,k,n)-fdexp(:,:,k+1,n))/dx(k)* &
                          (exp_xds(k)-exp_xds(k+1))*x_scale )
      enddo
    enddo

    if ( l_x0i0   )  intg(:,:,i0  ,:) = 0.
    if ( l_x0i0p1 )  intg(:,:,i0+1,:) = 0.

  case ( 4 )

    do n=1, nd4
      fdexp(:,:,:,n) = f(:,:,:,n)/exp_xds(n)
    enddo

    do k=1, nd3
    do j=1, nd2
    do i=1, nd1
      fdexp0 = d2i0p1*fdexp(i,j,k,i0)+d2i0*fdexp(i,j,k,i0+1)
      f0 = fdexp0*exp_xds0
      intg(i,j,k,i0:i0+1) = x_scale*( (f(i,j,k,i0:i0+1)-f0) -        &
                          (fdexp(i,j,k,i0:i0+1)-fdexp0)/dx(i0:i0+1)* &
                          (exp_xds(i0:i0+1)-exp_xds0)*x_scale )
    enddo
    enddo
    enddo
    do n=i0+2, nx
      intg(:,:,:,n) = intg(:,:,:,n-1) + x_scale*                 &
                      ( (f(:,:,:,n)-f(:,:,:,n-1)) -              &
                        (fdexp(:,:,:,n)-fdexp(:,:,:,n-1))/dx(n)* &
                        (exp_xds(n)-exp_xds(n-1))*x_scale )
    enddo
    do n=i0-1, 1, -1
      intg(:,:,:,n) = intg(:,:,:,n+1) + x_scale*                 &
                      ( (f(:,:,:,n)-f(:,:,:,n+1)) -              &
                        (fdexp(:,:,:,n)-fdexp(:,:,:,n+1))/dx(n)* &
                        (exp_xds(n)-exp_xds(n+1))*x_scale )
    enddo

    if ( l_x0i0   )  intg(:,:,:,i0  ) = 0.
    if ( l_x0i0p1 )  intg(:,:,:,i0+1) = 0.

  END select


  RETURN

END subroutine integ1d_exp

!=======================================================================

SUBROUTINE intg_set(nx,x,x0,i0,dx)

  ! x(i0) <= x0 < x(i0+1), or x(i0) > x0 >= x(i0+1)

  implicit none

  integer,             intent(in) ::  nx
  real,                intent(in) ::  x0
  real, dimension(nx), intent(in) ::  x

  integer,             intent(out) ::  i0
  real, dimension(nx), intent(out) ::  dx

  integer ::  i


  if (x(1) < x(2))  then
    if ( x0 < x(1) .or. x0 > x(nx) ) then
      print*, 'Integration failed in INTG_SET: x0 out of the range.'
      STOP
    end if
    i0 = nx-1
    do i=2, nx
      if (x(i) > x0) then
        i0 = i - 1
        EXIT
      end if
    enddo
  else
    if ( x0 > x(1) .or. x0 < x(nx) ) then
      print*, 'Integration failed in INTG_SET: x0 out of the range.'
      STOP
    end if
    i0 = 1
    do i=nx-1, 1, -1
      if (x(i) > x0) then
        i0 = i
        EXIT
      end if
    enddo
  end if

  dx(i0:i0+1) = x(i0:i0+1) - x0
  do i=i0+2, nx
    dx(i) = x(i) - x(i-1)
  enddo
  do i=1, i0-1
    dx(i) = x(i) - x(i+1)
  enddo

  RETURN

END subroutine intg_set

!=======================================================================

END module integ


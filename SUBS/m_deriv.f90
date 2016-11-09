MODULE deriv

  implicit none

  CONTAINS
! consistency check for descending x, then fix m_temeq

!  l_bdy(2) :  0, for missing / 1, 1st-order differencing / 
!              2, periodic B.C. (2nd-order diff. on uniform grid)

!=======================================================================

SUBROUTINE deriv1d_uni(nd,id,dx,f,l_bdy,missvo,der)

  integer,                                  intent(in) ::  nd(4), id
  integer, dimension(2),                    intent(in) ::  l_bdy
  real,                                     intent(in) ::  dx, missvo
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer ::  i,j,k,n
  integer ::  nd1, nd2, nd3, nd4, nx
  real    ::  d2x


  nd1 = nd(1)  ;  nd2 = nd(2)  ;  nd3 = nd(3)  ;  nd4 = nd(4)
  nx  = nd(id)
  d2x = 2.*dx
  der(:,:,:,:) = missvo


  SELECT case ( id )

  case ( 1 )

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
      do i=2, nx-1
        der(i,j,k,n) = (f(i+1,j,k,n) - f(i-1,j,k,n)) / d2x
      enddo
    enddo
    enddo
    enddo

    if (l_bdy(1) == 1)  der(1 ,:,:,:) = (f(2 ,:,:,:) - f(1   ,:,:,:)) / dx
    if (l_bdy(2) == 1)  der(nx,:,:,:) = (f(nx,:,:,:) - f(nx-1,:,:,:)) / dx
    if (l_bdy(1) == 2)  der(1 ,:,:,:) = (f(2 ,:,:,:) - f(nx  ,:,:,:)) / d2x
    if (l_bdy(2) == 2)  der(nx,:,:,:) = (f(1 ,:,:,:) - f(nx-1,:,:,:)) / d2x

  case ( 2 )
      
    do n=1, nd4
    do k=1, nd3
      do j=2, nx-1
      do i=1, nd1
        der(i,j,k,n) = (f(i,j+1,k,n) - f(i,j-1,k,n)) / d2x
      enddo
      enddo
    enddo
    enddo

    if (l_bdy(1) == 1)  der(:,1 ,:,:) = (f(:,2 ,:,:) - f(:,1   ,:,:)) / dx
    if (l_bdy(2) == 1)  der(:,nx,:,:) = (f(:,nx,:,:) - f(:,nx-1,:,:)) / dx
    if (l_bdy(1) == 2)  der(:,1 ,:,:) = (f(:,2 ,:,:) - f(:,nx  ,:,:)) / d2x
    if (l_bdy(2) == 2)  der(:,nx,:,:) = (f(:,1 ,:,:) - f(:,nx-1,:,:)) / d2x

  case ( 3 )

    do n=1, nd4
      do k=2, nx-1
      do j=1, nd2
      do i=1, nd1
        der(i,j,k,n) = (f(i,j,k+1,n) - f(i,j,k-1,n)) / d2x
      enddo
      enddo
      enddo
    enddo

    if (l_bdy(1) == 1)  der(:,:,1 ,:) = (f(:,:,2 ,:) - f(:,:,1   ,:)) / dx
    if (l_bdy(2) == 1)  der(:,:,nx,:) = (f(:,:,nx,:) - f(:,:,nx-1,:)) / dx
    if (l_bdy(1) == 2)  der(:,:,1 ,:) = (f(:,:,2 ,:) - f(:,:,nx  ,:)) / d2x
    if (l_bdy(2) == 2)  der(:,:,nx,:) = (f(:,:,1 ,:) - f(:,:,nx-1,:)) / d2x

  case ( 4 )

    do n=2, nx-1
    do k=1, nd3
    do j=1, nd2
    do i=1, nd1
      der(i,j,k,n) = (f(i,j,k,n+1) - f(i,j,k,n-1)) / d2x
    enddo
    enddo
    enddo
    enddo

    if (l_bdy(1) == 1)  der(:,:,:,1 ) = (f(:,:,:,2 ) - f(:,:,:,1   )) / dx
    if (l_bdy(2) == 1)  der(:,:,:,nx) = (f(:,:,:,nx) - f(:,:,:,nx-1)) / dx
    if (l_bdy(1) == 2)  der(:,:,:,1 ) = (f(:,:,:,2 ) - f(:,:,:,nx  )) / d2x
    if (l_bdy(2) == 2)  der(:,:,:,nx) = (f(:,:,:,1 ) - f(:,:,:,nx-1)) / d2x

  END select

  RETURN

END subroutine deriv1d_uni

!=======================================================================

SUBROUTINE deriv1d(nd,id,x,f,l_bdy,missvo,der)

  integer,                                  intent(in) ::  nd(4), id
  integer, dimension(2),                    intent(in) ::  l_bdy
  real,                                     intent(in) ::  missvo
  real, dimension(nd(id)),                  intent(in) ::  x
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer ::  i,j,k,n
  integer ::  nd1, nd2, nd3, nd4, nx
  logical ::  l_uni
  real    ::  dcoef(3,nd(id)), dx1, dxn
  real*8  ::  dcoef8(3)


  dx1 = x(2 ) - x(1   )
  dxn = x(nx) - x(nx-1)

  l_uni = .TRUE.
  do i=2, nd(id)-1
    if ( abs((x(i+1)-x(i))/dx1-1.) > 1.e-3 )  l_uni = .FALSE.
  enddo
  if ( (.not. l_uni) .and. (l_bdy(1) == 2 .or. l_bdy(2) == 2) ) then
    print*, 'deriv1d: periodic B.C. on non-uniform grid is not allowed.'
    STOP
  end if

  IF ( l_uni ) then

  call deriv1d_uni(nd,id,x(2)-x(1),f,l_bdy,missvo, der)

  ELSE

  nd1 = nd(1)  ;  nd2 = nd(2)  ;  nd3 = nd(3)  ;  nd4 = nd(4)
  nx  = nd(id)
  der(:,:,:,:) = missvo

  do i=2, nx-1
    call fdcoef(2,3,dble(x(i)),dble(x(i-1:i+1)), dcoef8)
    dcoef(:,i) = dcoef8(:)
  enddo


  SELECT case ( id )

  case ( 1 )

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
      do i=2, nx-1
        der(i,j,k,n) = sum( dcoef(:,i)*f(i-1:i+1,j,k,n) )
      enddo
    enddo
    enddo
    enddo

    if (l_bdy(1) == 1)  der(1 ,:,:,:) = (f(2 ,:,:,:) - f(1   ,:,:,:)) / dx1
    if (l_bdy(2) == 1)  der(nx,:,:,:) = (f(nx,:,:,:) - f(nx-1,:,:,:)) / dxn

  case ( 2 )

    do n=1, nd4
    do k=1, nd3
      do j=2, nx-1
      do i=1, nd1
        der(i,j,k,n) = sum( dcoef(:,j)*f(i,j-1:j+1,k,n) )
      enddo
      enddo
    enddo
    enddo

    if (l_bdy(1) == 1)  der(:,1 ,:,:) = (f(:,2 ,:,:) - f(:,1   ,:,:)) / dx1
    if (l_bdy(2) == 1)  der(:,nx,:,:) = (f(:,nx,:,:) - f(:,nx-1,:,:)) / dxn

  case ( 3 )

    do n=1, nd4
      do k=2, nx-1
      do j=1, nd2
      do i=1, nd1
        der(i,j,k,n) = sum( dcoef(:,k)*f(i,j,k-1:k+1,n) )
      enddo
      enddo
      enddo
    enddo

    if (l_bdy(1) == 1)  der(:,:,1 ,:) = (f(:,:,2 ,:) - f(:,:,1   ,:)) / dx1
    if (l_bdy(2) == 1)  der(:,:,nx,:) = (f(:,:,nx,:) - f(:,:,nx-1,:)) / dxn

  case ( 4 )

    do n=2, nx-1
    do k=1, nd3
    do j=1, nd2
    do i=1, nd1
      der(i,j,k,n) = sum( dcoef(:,n)*f(i,j,k,n-1:n+1) )
    enddo
    enddo
    enddo
    enddo

    if (l_bdy(1) == 1)  der(:,:,:,1 ) = (f(:,:,:,2 ) - f(:,:,:,1   )) / dx1
    if (l_bdy(2) == 1)  der(:,:,:,nx) = (f(:,:,:,nx) - f(:,:,:,nx-1)) / dxn

  END select

  END if

  RETURN

END subroutine deriv1d

!=======================================================================

SUBROUTINE deriv1d_lon(nd,lon,lat,f,l_bdy,missvo,der)

  integer,                                  intent(in) ::  nd(4)
  integer, dimension(2),                    intent(in) ::  l_bdy
  real,                                     intent(in) ::  missvo
  real, dimension(nd(1)),                   intent(in) ::  lon
  real, dimension(nd(2)),                   intent(in) ::  lat
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer                ::  j
  real, dimension(nd(1)) ::  x

  include 'c_math.inc'
  include 'c_phys.inc'


  do j=1, nd(2)
    x(:) = lon(:)*deg2rad*r_earth*cos(lat(j)*deg2rad)
    call deriv1d((/nd(1),1,nd(3),nd(4)/),1,x,f(:,j,:,:),l_bdy,missvo,    &
                 der(:,j,:,:))
  enddo
  if (abs(lat(1    )) == 90.)  der(:,1    ,:,:) = missvo
  if (abs(lat(nd(2))) == 90.)  der(:,nd(2),:,:) = missvo

  RETURN

END subroutine deriv1d_lon

!=======================================================================

SUBROUTINE deriv1d_lat(nd,lat,f,l_bdy,missvo,der)

  integer,                                  intent(in) ::  nd(4)
  integer, dimension(2),                    intent(in) ::  l_bdy
  real,                                     intent(in) ::  missvo
  real, dimension(nd(2)),                   intent(in) ::  lat
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer ::  l_bdy2(2)
  real, dimension(nd(2)) ::  x

  include 'c_math.inc'
  include 'c_phys.inc'


  if ( l_bdy(1) == 2 .or. l_bdy(2) == 2 ) then
    print*, 'deriv1d_lat: B.C. should be checked.'  ;  STOP
  end if

  l_bdy2(:) = l_bdy(:)
  if (abs(lat(1    )) == 90.)  l_bdy2(1) = 0
  if (abs(lat(nd(2))) == 90.)  l_bdy2(2) = 0

  x(:) = lat(:)*deg2rad*r_earth
  call deriv1d(nd,2,x,f,l_bdy2,missvo, der)

  RETURN

END subroutine deriv1d_lat

!=======================================================================

SUBROUTINE deriv1d_exp(nd,id,x,f,l_exp,l_bdy,missvo,der)

  integer,                                  intent(in) ::  nd(4), id
  integer, dimension(2),                    intent(in) ::  l_bdy
  real,                                     intent(in) ::  missvo
  real, dimension(nd(id)),                  intent(in) ::  x
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f
  logical, dimension(2),                    intent(in) ::  l_exp

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer                                  ::  i,j,k,n
  real, dimension(nd(id))                  ::  xp
  real, dimension(nd(1),nd(2),nd(3),nd(4)) ::  fp, temp


  if ( l_bdy(1) == 2 .or. l_bdy(2) == 2 ) then
    print*, 'deriv1d_exp: B.C. should be checked.'  ;  STOP
  end if

  if ( l_exp(1) ) then
    if (x(1)*x(nd(id)) <= 0.) then
      print*, 'cannot use "deriv1d_exp"'  ;  STOP
    end if
    xp(:) = log(abs(x(:)))
  else
    xp(:) = x(:)
  end if

  if ( l_exp(2) ) then
    if (minval(f)*maxval(f) <= 0.) then
      print*, 'cannot use "deriv1d_exp"'  ;  STOP
    end if
    fp(:,:,:,:) = log(abs(f(:,:,:,:)))
  else
    fp(:,:,:,:) = f(:,:,:,:)
  end if

  call deriv1d(nd,id,xp,fp,l_bdy,missvo, temp)

  if ( l_exp(1) ) then
    SELECT case ( id )
    case ( 1 )
      do n=1, nd(4)
      do k=1, nd(3)
      do j=1, nd(2)
        temp(:,j,k,n) = temp(:,j,k,n)/x(:)
      enddo
      enddo
      enddo
    case ( 2 )
      do n=1, nd(4)
      do k=1, nd(3)
      do j=1, nd(2)
        temp(:,j,k,n) = temp(:,j,k,n)/x(j)
      enddo
      enddo
      enddo
    case ( 3 )
      do n=1, nd(4)
      do k=1, nd(3)
        temp(:,:,k,n) = temp(:,:,k,n)/x(k)
      enddo
      enddo
    case ( 4 )
      do n=1, nd(4)
        temp(:,:,:,n) = temp(:,:,:,n)/x(n)
      enddo
    END select
  end if

  if ( l_exp(2) ) then
    temp(:,:,:,:) = temp(:,:,:,:)*f(:,:,:,:)
  end if

  der(:,:,:,:) = temp(:,:,:,:)

  SELECT case ( id )
  case ( 1 )
    if (l_bdy(1) == 0)  der(1     ,:,:,:) = missvo
    if (l_bdy(2) == 0)  der(nd(id),:,:,:) = missvo
  case ( 2 )
    if (l_bdy(1) == 0)  der(:,1     ,:,:) = missvo
    if (l_bdy(2) == 0)  der(:,nd(id),:,:) = missvo
  case ( 3 )
    if (l_bdy(1) == 0)  der(:,:,1     ,:) = missvo
    if (l_bdy(2) == 0)  der(:,:,nd(id),:) = missvo
  case ( 4 )
    if (l_bdy(1) == 0)  der(:,:,:,1     ) = missvo
    if (l_bdy(2) == 0)  der(:,:,:,nd(id)) = missvo
  END select

  RETURN

END subroutine deriv1d_exp

!=======================================================================

SUBROUTINE div_lat(nd,lat,f,l_bdy,missvo,der)

  integer,                                  intent(in) ::  nd(4)
  integer, dimension(2),                    intent(in) ::  l_bdy
  real,                                     intent(in) ::  missvo
  real, dimension(nd(2)),                   intent(in) ::  lat
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer ::  j
  integer ::  ny, l_bdy2(2)
  real    ::  coslat(nd(2)), coslat2(2)
  real    ::  fc(nd(1),nd(2),nd(3),nd(4))

  include 'c_math.inc'


  if ( l_bdy(1) == 2 .or. l_bdy(2) == 2 ) then
    print*, 'div_lat: B.C. should be checked.'  ;  STOP
  end if

  coslat(:) = cos(lat(:)*deg2rad)

  l_bdy2(:) = l_bdy(:)

  ny = nd(2)
  if (abs(lat(1)) == 90.) then
    coslat(1) = 0.
    l_bdy2(1) = 0
  end if
  if (abs(lat(ny)) == 90.) then
    coslat(ny) = 0.
    l_bdy2(2 ) = 0
  end if
  coslat2(1) = cos(0.5*(lat(1 )+lat(2   ))*deg2rad)
  coslat2(2) = cos(0.5*(lat(ny)+lat(ny-1))*deg2rad)

  do j=1, nd(2)
    fc(:,j,:,:) = f(:,j,:,:)*coslat(j)
  enddo

  call deriv1d_lat(nd,lat,fc,l_bdy2,missvo, der)

  do j=2, nd(2)-1
    der(:,j,:,:) = der(:,j,:,:)/coslat(j)
  enddo
  if (l_bdy2(1) == 1)  der(:,1 ,:,:) = der(:,1 ,:,:)/coslat2(1)
  if (l_bdy2(2) == 1)  der(:,ny,:,:) = der(:,ny,:,:)/coslat2(2)

  RETURN

END subroutine div_lat

!=======================================================================

END module deriv

!=======================================================================

SUBROUTINE fdcoef(mord,nord,x0,grid,coef)

      implicit none
      save

!..This routine implements simple recursions for calculating the weights
!..of finite difference formulas for any order of derivative and any order
!..of accuracy on one-dimensional grids with arbitrary spacing.

!..from Bengt Fornberg's article
!..Generation of finite difference formulas on arbitrary spaced grids.
!..Math. Comp., 51(184):699-706, 1988.


!..input:
!..mord       = the order of the derivative
!..nord       = order of accuracy n
!..x0         = point at which to evaluate the coefficients
!..grid(nord) = array containing the grid

!..output:
!..coef(nord) = coefficients of the finite difference formula.


!..declare the pass
      integer          mord,nord
      double precision x0,grid(nord),coef(nord)


!..local variables
      integer          nu,nn,mm,nmmin,mmax,nmax
      parameter        (mmax=8, nmax=10)
      double precision weight(mmax,nmax,nmax),c1,c2,c3,c4,pfac


!..zero the weight array
      do nu=1,nord
       do nn=1,nord
        do mm=1,mord
         weight(mm,nn,nu) = 0.0d0
        enddo
       enddo
      enddo

      weight(1,1,1) = 1.0d0
      c1            = 1.0d0
      nmmin         = min(nord,mord)
      do nn = 2,nord
       c2 = 1.0d0
       do nu=1,nn-1
        c3 = grid(nn) - grid(nu)
        c2 = c2 * c3
        c4 = 1.0d0/c3
        pfac = grid(nn) - x0
        weight(1,nn,nu) = c4 * ( pfac * weight(1,nn-1,nu) )
        do mm=2,nmmin
         weight(mm,nn,nu) = c4 * ( pfac * weight(mm,nn-1,nu)            &
                            - dfloat(mm-1)*weight(mm-1,nn-1,nu) )
        enddo
       enddo
       pfac = (grid(nn-1) - x0)
       weight(1,nn,nn) = c1/c2 * (-pfac*weight(1,nn-1,nn-1))
       c4 = c1/c2
       do mm=2,nmmin
        weight(mm,nn,nn) = c4 * (dfloat(mm-1)*weight(mm-1,nn-1,nn-1) -  &
                                  pfac*weight(mm,nn-1,nn-1))
       enddo
       c1 = c2
      enddo

!..load the coefficients
      do nu = 1,nord
       coef(nu) = weight(mord,nord,nu)
      enddo

      return

END subroutine fdcoef


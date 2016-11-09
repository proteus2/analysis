MODULE deriv


  CONTAINS
! consistency check for descending x, then fix m_temeq

!=======================================================================

SUBROUTINE deriv1d_uni(nd,id,dx,f,l_bdy,missvo,der)

  implicit none

  integer,                                  intent(in) ::  nd(4), id
  real,                                     intent(in) ::  dx, missvo
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f
  logical, dimension(2),                    intent(in) ::  l_bdy

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

    if ( l_bdy(1) )  der(1 ,:,:,:) = (f(2 ,:,:,:) - f(1   ,:,:,:)) / dx
    if ( l_bdy(2) )  der(nx,:,:,:) = (f(nx,:,:,:) - f(nx-1,:,:,:)) / dx

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

    if ( l_bdy(1) )  der(:,1 ,:,:) = (f(:,2 ,:,:) - f(:,1   ,:,:)) / dx
    if ( l_bdy(2) )  der(:,nx,:,:) = (f(:,nx,:,:) - f(:,nx-1,:,:)) / dx

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

    if ( l_bdy(1) )  der(:,:,1 ,:) = (f(:,:,2 ,:) - f(:,:,1   ,:)) / dx
    if ( l_bdy(2) )  der(:,:,nx,:) = (f(:,:,nx,:) - f(:,:,nx-1,:)) / dx

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

    if ( l_bdy(1) )  der(:,:,:,1 ) = (f(:,:,:,2 ) - f(:,:,:,1   )) / dx
    if ( l_bdy(2) )  der(:,:,:,nx) = (f(:,:,:,nx) - f(:,:,:,nx-1)) / dx

  END select

  RETURN

END subroutine deriv1d_uni

!=======================================================================

SUBROUTINE deriv1d(nd,id,l_uni,x,f,l_bdy,missvo,der)

  implicit none

  integer,                                  intent(in) ::  nd(4), id
  real,                                     intent(in) ::  missvo
  real, dimension(nd(id)),                  intent(in) ::  x
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f
  logical,                                  intent(in) ::  l_uni
  logical, dimension(2),                    intent(in) ::  l_bdy

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer ::  i,j,k,n
  integer ::  nd1, nd2, nd3, nd4, nx
  real    ::  dcoef(3,nd(id)), dx1, dxn
  real*8  ::  dcoef8(3)


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
  dx1 = x(2 ) - x(1   )
  dxn = x(nx) - x(nx-1)


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

    if ( l_bdy(1) )  der(1 ,:,:,:) = (f(2 ,:,:,:) - f(1   ,:,:,:)) / dx1
    if ( l_bdy(2) )  der(nx,:,:,:) = (f(nx,:,:,:) - f(nx-1,:,:,:)) / dxn

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

    if ( l_bdy(1) )  der(:,1 ,:,:) = (f(:,2 ,:,:) - f(:,1   ,:,:)) / dx1
    if ( l_bdy(2) )  der(:,nx,:,:) = (f(:,nx,:,:) - f(:,nx-1,:,:)) / dxn

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

    if ( l_bdy(1) )  der(:,:,1 ,:) = (f(:,:,2 ,:) - f(:,:,1   ,:)) / dx1
    if ( l_bdy(2) )  der(:,:,nx,:) = (f(:,:,nx,:) - f(:,:,nx-1,:)) / dxn

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

    if ( l_bdy(1) )  der(:,:,:,1 ) = (f(:,:,:,2 ) - f(:,:,:,1   )) / dx1
    if ( l_bdy(2) )  der(:,:,:,nx) = (f(:,:,:,nx) - f(:,:,:,nx-1)) / dxn

  END select

  END if

  RETURN

END subroutine deriv1d

!=======================================================================

SUBROUTINE deriv1d_exp(nd,id,l_uni,x,f,f_scale,l_bdy,missvo,der)

  implicit none

  integer,                                  intent(in) ::  nd(4), id
  real,                                     intent(in) ::  f_scale
  real,                                     intent(in) ::  missvo
  real, dimension(nd(id)),                  intent(in) ::  x
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  f
  logical,                                  intent(in) ::  l_uni
  logical, dimension(2),                    intent(in) ::  l_bdy

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(out) ::  der

  integer ::  i,j,k,n
  real    ::  exp_xds(nd(id))
  real, dimension(nd(1),nd(2),nd(3),nd(4)) ::  temp1, temp2


  exp_xds(:) = exp(x(:)/f_scale)

  SELECT case ( id )
  case ( 1 )
    do n=1, nd(4)
    do k=1, nd(3)
    do j=1, nd(2)
      der(:,j,k,n) = exp_xds(:)
    enddo
    enddo
    enddo
  case ( 2 )
    do n=1, nd(4)
    do k=1, nd(3)
    do j=1, nd(2)
      der(:,j,k,n) = exp_xds(j)
    enddo
    enddo
    enddo
  case ( 3 )
    do n=1, nd(4)
    do k=1, nd(3)
      der(:,:,k,n) = exp_xds(k)
    enddo
    enddo
  case ( 4 )
    do n=1, nd(4)
      der(:,:,:,n) = exp_xds(n)
    enddo
  END select

  temp1(:,:,:,:) = f(:,:,:,:)/der(:,:,:,:)

  call deriv1d(nd,id,l_uni,x,temp1,l_bdy,missvo, temp2)

  der(:,:,:,:) = der(:,:,:,:)*(temp2(:,:,:,:)+temp1(:,:,:,:)/f_scale)

  SELECT case ( id )
  case ( 1 )
    if ( .not. l_bdy(1) )  der(1     ,:,:,:) = missvo
    if ( .not. l_bdy(2) )  der(nd(id),:,:,:) = missvo
  case ( 2 )
    if ( .not. l_bdy(1) )  der(:,1     ,:,:) = missvo
    if ( .not. l_bdy(2) )  der(:,nd(id),:,:) = missvo
  case ( 3 )
    if ( .not. l_bdy(1) )  der(:,:,1     ,:) = missvo
    if ( .not. l_bdy(2) )  der(:,:,nd(id),:) = missvo
  case ( 4 )
    if ( .not. l_bdy(1) )  der(:,:,:,1     ) = missvo
    if ( .not. l_bdy(2) )  der(:,:,:,nd(id)) = missvo
  END select

  RETURN

END subroutine deriv1d_exp

!=======================================================================

SUBROUTINE div_lat_uni(nx,ny,nz,nt,lat,f,l_bdy,missvo,divy)

! boundary condition is reset, if it has poles.

  implicit none

  integer,                         intent(in) ::  nx, ny, nz, nt
  real,                            intent(in) ::  missvo
  real,    dimension(ny),          intent(in) ::  lat
  real,    dimension(nx,ny,nz,nt), intent(in) ::  f
  logical, dimension(2),           intent(in) ::  l_bdy

  real, dimension(nx,ny,nz,nt), intent(out) ::  divy

  integer                      ::  j,k,n, iy1, iy2
  real                         ::  dy, cosphi2(2)
  real, dimension(ny)          ::  phi, tempy
  real, dimension(nx,ny)       ::  cosphi
  real, dimension(nx,ny,nz,nt) ::  temp
  logical                      ::  l_bdy2(2)

  include 'c_math.inc'
  include 'c_phys.inc'


  l_bdy2 = l_bdy

  tempy(:) = cos(lat(:)*deg2rad)
  if (abs(lat(1)) == 90.) then
    tempy(1) = 0.
    l_bdy2(1) = .FALSE.
  end if
  if (abs(lat(ny)) == 90.) then
    tempy(ny) = 0.
    l_bdy2(2) = .FALSE.
  end if
  do j=1, ny
    cosphi(:,j) = tempy(j)
  enddo
  cosphi2(1) = cos(0.5*(lat(1   )+lat(2 ))*deg2rad)
  cosphi2(2) = cos(0.5*(lat(ny-1)+lat(ny))*deg2rad)

  do n=1, nt
  do k=1, nz
    temp(:,:,k,n) = f(:,:,k,n)*cosphi(:,:)
  enddo
  enddo

  dy = r_earth*(lat(2)-lat(1))*deg2rad
  call deriv1d_uni((/nx,ny,nz,nt/),2,dy,temp,l_bdy2,missvo, divy)

  do n=1, nt
  do k=1, nz
    divy(:,2:ny-1,k,n) = divy(:,2:ny-1,k,n) / cosphi(:,2:ny-1)
  enddo
  enddo
  if ( l_bdy2(1) )  divy(:,1 ,k,n) = divy(:,1 ,k,n) / cosphi2(1)
  if ( l_bdy2(2) )  divy(:,ny,k,n) = divy(:,ny,k,n) / cosphi2(2)

  RETURN

END subroutine div_lat_uni

!=======================================================================

SUBROUTINE div_lat(nx,ny,nz,nt,l_uni,lat,f,l_bdy,missvo,divy)

! boundary condition is reset, if it has poles.

  implicit none

  integer,                         intent(in) ::  nx, ny, nz, nt
  real,                            intent(in) ::  missvo
  real,    dimension(ny),          intent(in) ::  lat
  real,    dimension(nx,ny,nz,nt), intent(in) ::  f
  logical,                         intent(in) ::  l_uni
  logical, dimension(2),           intent(in) ::  l_bdy

  real, dimension(nx,ny,nz,nt), intent(out) ::  divy

  integer                      ::  j,k,n, iy1, iy2
  real                         ::  cosphi2(2)
  real, dimension(ny)          ::  y, tempy
  real, dimension(nx,ny)       ::  cosphi
  real, dimension(nx,ny,nz,nt) ::  temp
  logical                      ::  l_bdy2(2)

  include 'c_math.inc'
  include 'c_phys.inc'


  IF ( l_uni ) then

  call div_lat_uni(nx,ny,nz,nt,lat,f,l_bdy,missvo,divy)

  ELSE

  l_bdy2 = l_bdy

  tempy(:) = lat(:)*deg2rad
  y(:) = r_earth*tempy(:)

  tempy(:) = cos(tempy(:))
  if (abs(lat(1)) == 90.) then
    tempy(1) = 0.
    l_bdy2(1) = .FALSE.
  end if
  if (abs(lat(ny)) == 90.) then
    tempy(ny) = 0.
    l_bdy2(2) = .FALSE.
  end if
  do j=1, ny
    cosphi(:,j) = tempy(j)
  enddo
  cosphi2(1) = cos(0.5*(lat(1   )+lat(2 ))*deg2rad)
  cosphi2(2) = cos(0.5*(lat(ny-1)+lat(ny))*deg2rad)

  do n=1, nt
  do k=1, nz
    temp(:,:,k,n) = f(:,:,k,n)*cosphi(:,:)
  enddo
  enddo

  call deriv1d((/nx,ny,nz,nt/),2,.FALSE.,y,temp,l_bdy2,missvo, divy)

  do n=1, nt
  do k=1, nz
    divy(:,2:ny-1,k,n) = divy(:,2:ny-1,k,n) / cosphi(:,2:ny-1)
  enddo
  enddo
  if ( l_bdy2(1) )  divy(:,1 ,k,n) = divy(:,1 ,k,n) / cosphi2(1)
  if ( l_bdy2(2) )  divy(:,ny,k,n) = divy(:,ny,k,n) / cosphi2(2)

  END if

  RETURN

END subroutine div_lat

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

!=======================================================================

END module deriv


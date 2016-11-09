MODULE utiletc

  implicit none

!  Include files used here
!    c_math


  CONTAINS

!=======================================================================

SUBROUTINE reorder(nd,id,vari,varo)

  integer, intent(in) ::  nd(4), id(4)

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(inout) ::  vari
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(inout) ::  varo

  real, dimension(nd(1),nd(2),nd(3),nd(4)) ::  temp

  integer ::  n

  temp(:,:,:,:) = vari(:,:,:,:)
  do n=1, 4
    if (id(n) /= 0)  call reverse_dim(nd,n,temp, temp)
  enddo
  varo(:,:,:,:) = temp(:,:,:,:)

END SUBROUTINE reorder

!=======================================================================

SUBROUTINE reverse_dim(nd,id,vari, varo)

  integer, intent(in) ::  nd(4), id

  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(inout) ::  vari
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(inout) ::  varo

  real, dimension(nd(1),nd(2),nd(3),nd(4)) ::  temp
  integer ::  i,j,k,n
  integer ::  nd1, nd2, nd3, nd4

  nd1 = nd(1)  ;  nd2 = nd(2)  ;  nd3 = nd(3)  ;  nd4 = nd(4)
  temp(:,:,:,:) = vari(:,:,:,:)

  SELECT case ( id )

  case ( 1 )

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
    do i=1, nd1
      varo(i,j,k,n) = temp(nd1+1-i,j,k,n)
    enddo
    enddo
    enddo
    enddo

  case ( 2 )

    do n=1, nd4
    do k=1, nd3
    do j=1, nd2
      varo(:,j,k,n) = temp(:,nd2+1-j,k,n)
    enddo
    enddo
    enddo

  case ( 3 )

    do n=1, nd4
    do k=1, nd3
      varo(:,:,k,n) = temp(:,:,nd3+1-k,n)
    enddo
    enddo

  case ( 4 )

    do n=1, nd4
      varo(:,:,:,n) = temp(:,:,:,nd4+1-n)
    enddo

  END select

END subroutine reverse_dim

!=======================================================================

SUBROUTINE f_coriolis(nlat,lat,f)

  integer,               intent(in) ::  nlat
  real, dimension(nlat), intent(in) ::  lat

  real, dimension(nlat), intent(out) ::  f

  integer ::  j
  real    ::  ome_earth

  include 'c_math.inc'

  ome_earth = 2.*pi/86400.

  f(:) = 2.*ome_earth*sin(lat(:)*deg2rad)

END subroutine f_coriolis

!=======================================================================

SUBROUTINE coslat(nlat,lat,cosphi)

  integer,               intent(in) ::  nlat
  real, dimension(nlat), intent(in) ::  lat

  real, dimension(nlat), intent(out) ::  cosphi

  include 'c_math.inc'

  cosphi(:) = cos(lat(:)*deg2rad)
  if (abs(lat(1   )) == 90.)  cosphi(1   ) = 0.
  if (abs(lat(nlat)) == 90.)  cosphi(nlat) = 0.

END subroutine coslat

!=======================================================================

SUBROUTINE tridag(n,a,b,c,d)

  integer,            intent(in)    ::  n
  real, dimension(n), intent(in)    ::  a, b, c
  real, dimension(n), intent(inout) ::  d

  integer ::  i
  real    ::  tem1, tem2(n-1)

  d(1) = d(1)/b(1)
  tem1 = b(1)

  do i=2, n
    tem2(i) = c(i-1)/tem1
    tem1 = b(i)-a(i)*tem2(i)
    d(i) = (d(i)-a(i)*d(i-1))/tem1
  enddo

  do i=n-1, 1, -1
    d(i) = d(i)-tem2(i+1)*d(i+1)
  enddo

END subroutine tridag

!=======================================================================

SUBROUTINE check_ex(fname,existence)

  character(len=*), intent(in)  ::  fname
  logical,          intent(out) ::  existence

  inquire(file=trim(fname), exist=existence)
  if (.not. existence)  print*, '    ',trim(fname),' not found.'

END subroutine check_ex

!=======================================================================

END module utiletc


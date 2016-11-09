! to be modified for ununiform grid (now, zonal_avg only)
MODULE avg


  CONTAINS

! zonal_avg
! merid_avg
! zonal_tempo_avg
! corr_zonal_avg
! div_yz


! Routines =====================================================================


! avg_d1, 2, 3, 4 - dim2 : if additional one appended (lat), weighting
SUBROUTINE avg_d1(nd,vari,missv,avg1)

  implicit none

  integer,                                  intent(in) ::  nd(4)
  real,                                     intent(in) ::  missv
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  vari

  real, dimension(nd(2),nd(3),nd(4)), intent(out) ::  avg1


  integer ::  i,j,k,n
  integer ::  num


  if (missv == 1.) then

    do n=1, nd(4)
    do k=1, nd(3)
    do j=1, nd(2)
      avg1(j,k,n) = sum(vari(:,j,k,n))/nd(1)
    enddo
    enddo
    enddo
    
  else

    avg1(:,:,:) = 0.
    do n=1, nd(4)
    do k=1, nd(3)
    do j=1, nd(2)
      num = 0
      do i=1, nd(1)
        if (vari(i,j,k,n) /= missv) then
          num = num + 1
          avg1(j,k,n) = avg1(j,k,n) + vari(i,j,k,n)
        end if
      enddo
      if (num /= 0) then
        avg1(j,k,n) = avg1(j,k,n) / float(num)
      else
        avg1(j,k,n) = missv
      end if
    enddo
    enddo
    enddo

  end if

  RETURN

END subroutine avg_d1


SUBROUTINE avg_d2(nd,vari,missv,avg1)

  implicit none

  integer,                                  intent(in) ::  nd(4)
  real,                                     intent(in) ::  missv
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  vari

  real, dimension(nd(1),nd(3),nd(4)), intent(out) ::  avg1


  integer ::  i,j,k,n
  integer ::  num


  if (missv == 1.) then

    do n=1, nd(4)
    do k=1, nd(3)
    do i=1 ,nd(1)
      avg1(i,k,n) = sum(vari(i,:,k,n))/nd(2)
    enddo
    enddo
    enddo

  else

    avg1(:,:,:) = 0.
    do n=1, nd(4)
    do k=1, nd(3)
    do i=1, nd(1)
      num = 0
      do j=1, nd(2)
        if (vari(i,j,k,n) /= missv) then
          num = num + 1
          avg1(i,k,n) = avg1(i,k,n) + vari(i,j,k,n)
        end if
      enddo
      if (num /= 0) then
        avg1(i,k,n) = avg1(i,k,n) / float(num)
      else
        avg1(i,k,n) = missv
      end if
    enddo
    enddo
    enddo

  end if

  RETURN

END subroutine avg_d2


SUBROUTINE avg_d3(nd,vari,missv,avg1)

  implicit none

  integer,                                  intent(in) ::  nd(4)
  real,                                     intent(in) ::  missv
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  vari

  real, dimension(nd(1),nd(2),nd(4)), intent(out) ::  avg1


  integer ::  i,j,k,n
  integer ::  num


  if (missv == 1.) then

    do n=1, nd(4)
    do j=1, nd(2)
    do i=1 ,nd(1)
      avg1(i,j,n) = sum(vari(i,j,:,n))/nd(3)
    enddo
    enddo
    enddo

  else

    avg1(:,:,:) = 0.
    do n=1, nd(4)
    do j=1, nd(2)
    do i=1, nd(1)
      num = 0
      do k=1, nd(3)
        if (vari(i,j,k,n) /= missv) then
          num = num + 1
          avg1(i,j,n) = avg1(i,j,n) + vari(i,j,k,n)
        end if
      enddo
      if (num /= 0) then
        avg1(i,j,n) = avg1(i,j,n) / float(num)
      else
        avg1(i,j,n) = missv
      end if
    enddo
    enddo
    enddo

  end if

  RETURN

END subroutine avg_d3


SUBROUTINE avg_d4(nd,vari,missv,avg1)

  implicit none

  integer,                                  intent(in) ::  nd(4)
  real,                                     intent(in) ::  missv
  real, dimension(nd(1),nd(2),nd(3),nd(4)), intent(in) ::  vari

  real, dimension(nd(1),nd(2),nd(3)), intent(out) ::  avg1


  integer ::  i,j,k,n
  integer ::  num


  if (missv == 1.) then

    do k=1, nd(3)
    do j=1, nd(2)
    do i=1 ,nd(1)
      avg1(i,j,k) = sum(vari(i,j,k,:))/nd(4)
    enddo
    enddo
    enddo

  else

    avg1(:,:,:) = 0.
    do k=1, nd(3)
    do j=1, nd(2)
    do i=1, nd(1)
      num = 0
      do n=1, nd(4)
        if (vari(i,j,k,n) /= missv) then
          num = num + 1
          avg1(i,j,k) = avg1(i,j,k) + vari(i,j,k,n)
        end if
      enddo
      if (num /= 0) then
        avg1(i,j,k) = avg1(i,j,k) / float(num)
      else
        avg1(i,j,k) = missv
      end if
    enddo
    enddo
    enddo

  end if

  RETURN

END subroutine avg_d4


SUBROUTINE zmflx(nx,ny,nz,nt,var1,var2,missv,mflx)

  implicit none

  integer,                      intent(in) ::  nx, ny, nz, nt
  real,                         intent(in) ::  missv
  real, dimension(nx,ny,nz,nt), intent(in) ::  var1, var2

  real, dimension(ny,nz,nt), intent(out) ::  mflx

  integer ::  i,j,k,n
  integer ::  num
  real    ::  temp(nx,ny,nz,nt)


  if (missv == 1.) then

    temp(:,:,:,:) = var1(:,:,:,:)*var2(:,:,:,:)
    call avg_d1((/nx,ny,nz,nt/),temp,1., mflx)

  else

    mflx(:,:,:) = 0.
    do n=1, nt
    do k=1, nz
    do j=1, ny
      num = 0
      do i=1, nx
        if ( var1(i,j,k,n) /= missv .and. var2(i,j,k,n) /= missv ) then
          num = num + 1
          mflx(j,k,n) = mflx(j,k,n) + var1(i,j,k,n)*var2(i,j,k,n)
        end if
      enddo
      if (num /= 0) then
        mflx(j,k,n) = mflx(j,k,n) / float(num)
      else
        mflx(j,k,n) = missv
      end if
    enddo
    enddo
    enddo

  end if

  RETURN

END subroutine zmflx


SUBROUTINE zm_prt(nx,ny,nz,nt,var,missv,zm,prt)

  implicit none

  integer,                      intent(in) ::  nx, ny, nz, nt
  real,                         intent(in) ::  missv
  real, dimension(nx,ny,nz,nt), intent(in) ::  var

  real, dimension(ny,nz,nt),    intent(out) ::  zm
  real, dimension(nx,ny,nz,nt), intent(out) ::  prt

  integer ::  j,k,n


  call avg_d1((/nx,ny,nz,nt/),var,missv, zm)

  do n=1, nt
  do k=1, nz
  do j=1, ny
    prt(:,j,k,n) = var(:,j,k,n) - zm(j,k,n)
  enddo
  enddo
  enddo

  RETURN

END subroutine zm_prt


SUBROUTINE merid_avg(nx,ny,nz,nt,missv,lat,vari,varo)

  implicit none

  integer,                         intent(in)  ::  nx, ny, nz, nt
  real,                            intent(in)  ::  missv
  real,    dimension(ny),          intent(in)  ::  lat
  real,    dimension(nx,ny,nz,nt), intent(in)  ::  vari
  real,    dimension(nx,nz,nt),    intent(out) ::  varo

  integer ::  i,j,k,n
  real    ::  cosphi(ny), sumcos(nx)

  include 'c_math.inc'


  do j=1, ny
    cosphi(j) = cos(lat(j)/180.*pi)
  enddo

  varo(:,:,:) = 0.
  do n=1, nt
  do k=1, nz
    sumcos(:) = 0.
    do j=1, ny
    do i=1, nx
      if (vari(i,j,k,n) /= missv) then
        sumcos(i) = sumcos(i) + cosphi(j)
        varo(i,k,n) = varo(i,k,n) + vari(i,j,k,n)*cosphi(j)
      end if
    enddo
    enddo
    do i=1, nx
      if (sumcos(i) /= 0.) then
        varo(i,k,n) = varo(i,k,n) / sumcos(i)
      else
        varo(i,k,n) = missv
      end if
    enddo
  enddo
  enddo


  RETURN

END subroutine merid_avg


SUBROUTINE corr_zonal_avg(nx,ny,nz,nt,missv,vari1,vari2, varo)

  implicit none

  integer,                         intent(in)  ::  nx, ny, nz, nt
  real,                            intent(in)  ::  missv
  real,    dimension(nx,ny,nz,nt), intent(in)  ::  vari1, vari2
  real,    dimension(ny,nz,nt),    intent(out) ::  varo

  integer ::  i,j,k,n
  real    ::  temp(nx,ny,nz,nt)


  do n=1, nt
  do k=1, nz
  do j=1, ny
  do i=1, nx
    if ( vari1(i,j,k,n) /= missv .and. vari2(i,j,k,n) /= missv ) then
      temp(i,j,k,n) = vari1(i,j,k,n)*vari2(i,j,k,n)
    else
      temp(i,j,k,n) = missv
    end if
  enddo
  enddo
  enddo
  enddo

!  call zonal_avg(nx,ny,nz,nt,missv,temp,varo)


  RETURN

END subroutine corr_zonal_avg


SUBROUTINE div_yz(ny,nz,nt,var_y,var_z,lat,zh,missv, div)

  implicit none

  integer,                   intent(in)  ::  ny, nz, nt
  real,                      intent(in)  ::  missv
  real, dimension(ny,nz,nt), intent(in)  ::  var_y, var_z
  real, dimension(ny),       intent(in)  ::  lat
  real, dimension(nz),       intent(in)  ::  zh
  real, dimension(ny,nz,nt), intent(out) ::  div

  integer ::  j,k,n
  real    ::  cosphi(ny), temp(ny,nz), grady(ny,nz), gradz(ny,nz)

  include 'c_math.inc'


  do j=1, ny
    cosphi(j) = cos(lat(j)*deg2rad)
  enddo

  do n=1, nt

    do k=1, nz
    do j=1, ny
      if (var_y(j,k,n) /= missv) then
        temp(j,k) = var_y(j,k,n)*cosphi(j)
      else
        temp(j,k) = missv
      end if
    enddo
    enddo
!    call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)
    grady( 1,:) = missv
    grady(ny,:) = missv
    do k=1, nz
    do j=2, ny-1
      if ( grady(j,k) /= missv )  grady(j,k) = grady(j,k)/cosphi(j)
    enddo
    enddo

!    call gradz_2nd(1,ny,nz,1,var_z(:,:,n),zh,missv, gradz)

    do k=1, nz
      do j=2, ny-1
        if ( grady(j,k) /= missv .and. gradz(j,k) /= missv ) then
          div(j,k,n) = grady(j,k) + gradz(j,k)
        else
          div(j,k,n) = missv
        end if
      enddo
      div( 1,k,n) = missv
      div(ny,k,n) = missv
    enddo

  enddo


  RETURN

END subroutine div_yz

! End routines =================================================================


END module avg


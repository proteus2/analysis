MODULE UM_anal

  implicit none

  real, parameter ::  g = 9.80665, rd = 287.05, cp = 1005.
  real, parameter ::  kappa = rd/cp
  real, parameter ::  r_earth = 6371229.
  real, parameter ::  pi = 3.14159265358979323846
  real, parameter ::  ome_earth = 2.*pi/86400.
  real, parameter ::  deg2rad = pi/180.

  contains

! zonal_avg
! tempo_avg
! zonal_tempo_avg
! corr_zonal_avg
! int_z2p
! gradx_2nd
! grady_2nd
! gradz_2nd
! u2rho
! v2rho
! t2rho
! div_yz


! Routines =====================================================================


SUBROUTINE zonal_avg(nx,ny,nz,nt,missv,vari,varo)

  implicit none

  integer,                         intent(in)  ::  nx, ny, nz, nt
  real,                            intent(in)  ::  missv
  real,    dimension(nx,ny,nz,nt), intent(in)  ::  vari
  real,    dimension(ny,nz,nt),    intent(out) ::  varo

  integer ::  i,j,k,n
  integer ::  num


  varo(:,:,:) = 0.
  do n=1, nt
  do k=1, nz
  do j=1, ny
    num = 0
    do i=1, nx
      if (vari(i,j,k,n) /= missv) then
        num = num + 1
        varo(j,k,n) = varo(j,k,n) + vari(i,j,k,n)
      end if
    enddo
    if (num /= 0) then
      varo(j,k,n) = varo(j,k,n) / float(num)
    else
      varo(j,k,n) = missv
    end if
  enddo
  enddo
  enddo


  RETURN

END subroutine zonal_avg


SUBROUTINE merid_avg(nx,ny,nz,nt,missv,lat,vari,varo)

  implicit none

  integer,                         intent(in)  ::  nx, ny, nz, nt
  real,                            intent(in)  ::  missv
  real,    dimension(ny),          intent(in)  ::  lat
  real,    dimension(nx,ny,nz,nt), intent(in)  ::  vari
  real,    dimension(nx,nz,nt),    intent(out) ::  varo

  integer ::  i,j,k,n
  real    ::  cosphi(ny), sumcos(nx)
  real    ::  pi


  pi = acos(-1.)
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


SUBROUTINE tempo_avg(nx,ny,nz,nt,missv,vari,varo)

  implicit none

  integer,                         intent(in)  ::  nx, ny, nz, nt
  real,                            intent(in)  ::  missv
  real,    dimension(nx,ny,nz,nt), intent(in)  ::  vari
  real,    dimension(nx,ny,nz),    intent(out) ::  varo

  integer ::  i,j,k,n
  integer ::  num(nx,ny,nz)


  varo(:,:,:) = 0.
  num (:,:,:) = 0.
  do n=1, nt
  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (vari(i,j,k,n) /= missv) then
      num(i,j,k) = num(i,j,k) + 1
      varo(i,j,k) = varo(i,j,k) + vari(i,j,k,n)
    end if
  enddo
  enddo
  enddo
  enddo

  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (num(i,j,k) /= 0) then
      varo(i,j,k) = varo(i,j,k) / num(i,j,k)
    else
      varo(i,j,k) = missv
    end if
  enddo
  enddo
  enddo


  RETURN

END subroutine tempo_avg


SUBROUTINE zonal_tempo_avg(nx,ny,nz,nt,missv,vari,varo)

  implicit none

  integer,                         intent(in)  ::  nx, ny, nz, nt
  real,                            intent(in)  ::  missv
  real,    dimension(nx,ny,nz,nt), intent(in)  ::  vari
  real,    dimension(ny,nz),       intent(out) ::  varo

  real ::  zavg(ny,nz,nt)


  call zonal_avg(nx,ny,nz,nt,missv,vari,zavg)

  call tempo_avg(1,ny,nz,nt,missv,zavg,varo)


  RETURN

END subroutine zonal_tempo_avg


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

  call zonal_avg(nx,ny,nz,nt,missv,temp,varo)


  RETURN

END subroutine corr_zonal_avg


SUBROUTINE int_z2p(nx,ny,nz,nt,k1,varz,pz,np,p,missv, varp)

  implicit none

  integer,                         intent(in)  ::  nx, ny, nz, nt, k1, np
  real,                            intent(in)  ::  missv
  real, dimension(np),             intent(in)  ::  p
  real, dimension(nx,ny,k1:nz,nt), intent(in)  ::  varz
  real, dimension(nx,ny, 0:nz,nt), intent(in)  ::  pz
  real, dimension(nx,ny,   np,nt), intent(out) ::  varp

  integer ::  i,j,k,n, kz


  do n=1, nt
  do k=1, np
  do j=1, ny
  do i=1, nx

    ! if p(k) is out of the atmosphere-model domain
    if ( p(k) > pz(i,j, 0,n) .or. p(k) < pz(i,j,nz,n) ) then
      varp(i,j,k,n) = missv
      CYCLE
    end if

!    if (p(k) >= pz(i,j,k1,n)) then
!      varp(i,j,k,n) = varz(i,j,k1,n)
!      CYCLE
!    end if

    do kz=k1+1, nz
      ! varp is interpolated (or, sometimes extrapolated near sfc)
      if (p(k) >= pz(i,j,kz,n)) then
        varp(i,j,k,n) = ( varz(i,j,kz-1,n)*log(pz(i,j,kz,n)/p(k)) +   &
                          varz(i,j,kz,  n)*log(p(k)/pz(i,j,kz-1,n)) ) &
                        / log(pz(i,j,kz,n)/pz(i,j,kz-1,n))
        EXIT
      end if
    enddo

  enddo
  enddo
  enddo
  enddo


  RETURN

END subroutine int_z2p


SUBROUTINE gradx_2nd(nx,ny,nz,nt,var,lon,lat,missv, gradx)

  implicit none

  integer,                      intent(in)  ::  nx, ny, nz, nt
  real,                         intent(in)  ::  missv
  real, dimension(nx,ny,nz,nt), intent(in)  ::  var
  real, dimension(nx)         , intent(in)  ::  lon
  real, dimension(ny)         , intent(in)  ::  lat
  real, dimension(nx,ny,nz,nt), intent(out) ::  gradx

  integer ::  i,j,k,n
  real    ::  c_2dx(ny)
  logical ::  ltag(nx,ny,nz,nt)


  do j=1, ny
    c_2dx(j) = (lon(3)-lon(1))*deg2rad*r_earth*cos(lat(j)*deg2rad)
  enddo

  do n=1, nt
  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (var(i,j,k,n) /= missv) then
      ltag(i,j,k,n) = .TRUE.
    else
      ltag(i,j,k,n) = .FALSE.
    end if
  enddo
  enddo
  enddo
  enddo

  do n=1, nt
  do k=1, nz
  do j=1, ny
    do i=2, nx-1
      if ( ltag(i-1,j,k,n) .and. ltag(i+1,j,k,n) ) then
        gradx(i,j,k,n) = (var(i+1,j,k,n)-var(i-1,j,k,n))/c_2dx(j)
      else
        gradx(i,j,k,n) = missv
      end if
    enddo
    if ( ltag(nx,j,k,n) .and. ltag(2,j,k,n) ) then
      gradx(1,j,k,n) = (var(2,j,k,n)-var(nx  ,j,k,n))/c_2dx(j)
    else
      gradx(1,j,k,n) = missv
    endif
    if ( ltag(nx-1,j,k,n) .and. ltag(1,j,k,n) ) then
      gradx(nx,j,k,n) = (var(1,j,k,n)-var(nx-1,j,k,n))/c_2dx(j)
    else
      gradx(nx,j,k,n) = missv
    end if
  enddo
  enddo
  enddo


  RETURN

END subroutine gradx_2nd


SUBROUTINE grady_2nd(nx,ny,nz,nt,var,lat,missv, grady)

  implicit none

  integer,                      intent(in)  ::  nx, ny, nz, nt
  real,                         intent(in)  ::  missv
  real, dimension(nx,ny,nz,nt), intent(in)  ::  var
  real, dimension(ny)         , intent(in)  ::  lat
  real, dimension(nx,ny,nz,nt), intent(out) ::  grady

  integer ::  i,j,k,n
  real    ::  c_2dy
  logical ::  ltag(nx,ny,nz,nt)


  c_2dy = (lat(3)-lat(1))*deg2rad*r_earth


  IF (missv /= 0.) THEN

  do n=1, nt
  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (var(i,j,k,n) /= missv) then
      ltag(i,j,k,n) = .TRUE.
    else
      ltag(i,j,k,n) = .FALSE.
    end if
  enddo
  enddo
  enddo
  enddo

  grady(:,:,:,:) = missv
  do n=1, nt
  do k=1, nz
    do i=1, nx
      if ( ltag(i,1,k,n) ) then
        if ( ltag(i,3,k,n) )  &
           grady(i,2,k,n) = (var(i,3,k,n)-var(i,1,k,n))/c_2dy
      else
        if ( ltag(i,3,k,n) .and. ltag(i,2,k,n) )  &
           grady(i,2,k,n) = (var(i,3,k,n)-var(i,2,k,n))/(c_2dy/2.)
      end if
      if ( ltag(i,ny,k,n) ) then
        if ( ltag(i,ny-2,k,n) )  &
           grady(i,ny-1,k,n) = (var(i,ny,k,n)-var(i,ny-2,k,n))/c_2dy
      else
        if ( ltag(i,ny-2,k,n) .and. ltag(i,ny-1,k,n) )  &
           grady(i,ny-1,k,n) = (var(i,ny-1,k,n)-var(i,ny-2,k,n))/(c_2dy/2.)
      end if
    enddo
    do j=3, ny-2
    do i=1, nx
      if ( ltag(i,j-1,k,n) .and. ltag(i,j+1,k,n) )  &
         grady(i,j,k,n) = (var(i,j+1,k,n)-var(i,j-1,k,n))/c_2dy
    enddo
    enddo
  enddo
  enddo

  ELSE

  do n=1, nt
  do k=1, nz
    do j=2, ny-1
    do i=1, nx
      grady(i,j,k,n) = (var(i,j+1,k,n)-var(i,j-1,k,n))/c_2dy
    enddo
    enddo
    grady(:, 1,k,n) = 0.
    grady(:,ny,k,n) = 0.
  enddo
  enddo

  END IF  ! missv /= 0.


  RETURN

END subroutine grady_2nd


SUBROUTINE gradz_2nd(nx,ny,nz,nt,var,z,missv, gradz)

  implicit none

  integer,                      intent(in)  ::  nx, ny, nz, nt
  real,                         intent(in)  ::  missv
  real, dimension(nx,ny,nz,nt), intent(in)  ::  var
  real, dimension(nz)         , intent(in)  ::  z
  real, dimension(nx,ny,nz,nt), intent(out) ::  gradz

  integer ::  i,j,k,n
  integer ::  itag1(nx,ny,nt), itag2(nx,ny,nt)
  real    ::  c_2dz(nz)
  logical ::  ltag(nx,ny,nz,nt)


  do k=2, nz-1
    c_2dz(k) = z(k+1) - z(k-1)
  enddo
  c_2dz( 1) = z( 2) - z(  1 )
  c_2dz(nz) = z(nz) - z(nz-1)


  IF (missv /= 0.) THEN

  do n=1, nt
  do k=1, nz
  do j=1, ny
  do i=1, nx
    if (var(i,j,k,n) /= missv) then
      ltag(i,j,k,n) = .TRUE.
    else
      ltag(i,j,k,n) = .FALSE.
    end if
  enddo
  enddo
  enddo
  enddo

  itag1(:,:,:) = -999
  itag2(:,:,:) = -999
  do n=1, nt
  do j=1, ny
  do i=1, nx
    do k=1, nz
      if ( ltag(i,j,k,n) ) then
        itag1(i,j,n) = k
        EXIT
      end if
    enddo
    do k=nz, 1, -1
      if ( ltag(i,j,k,n) ) then
        itag2(i,j,n) = k
        EXIT
      end if
    enddo
  enddo
  enddo
  enddo

  gradz(:,:,:,:) = missv
  do n=1, nt
  do j=1, ny
  do i=1, nx
    do k=itag1(i,j,n)+1, itag2(i,j,n)-1
      gradz(i,j,k,n) = (var(i,j,k+1,n)-var(i,j,k-1,n))/c_2dz(k)
    enddo
    if (itag1(i,j,n) /= itag2(i,j,n)) then
      k = itag1(i,j,n)
      gradz(i,j,k,n) = (var(i,j,k+1,n)-var(i,j,k,n))/(z(k+1)-z(k))
      k = itag2(i,j,n)
      gradz(i,j,k,n) = (var(i,j,k,n)-var(i,j,k-1,n))/(z(k)-z(k-1))
    end if
  enddo
  enddo
  enddo

  ELSE

  do n=1, nt
    do k=2, nz-1
    do j=1, ny
    do i=1, nx
      gradz(i,j,k,n) = (var(i,j,k+1,n)-var(i,j,k-1,n))/c_2dz(k)
    enddo
    enddo
    enddo
    do j=1, ny
    do i=1, nx
      gradz(i,j, 1,n) = (var(i,j, 2,n)-var(i,j,  1 ,n))/c_2dz( 1)
      gradz(i,j,nz,n) = (var(i,j,nz,n)-var(i,j,nz-1,n))/c_2dz(nz)
    enddo
    enddo
  enddo

  END IF  ! missv /= 0.


  RETURN

END subroutine gradz_2nd


SUBROUTINE u2rho(nx,ny,nz,nt,u_in,missv, u_out)

  implicit none

  integer,                      intent(in)  ::  nx, ny, nz, nt
  real,                         intent(in)  ::  missv
  real, dimension(nx,ny,nz,nt), intent(in)  ::  u_in
  real, dimension(nx,ny,nz,nt), intent(out) ::  u_out

  integer ::  i,j,k,n
  logical ::  ltag(nx,ny,nz)


  do n=1, nt

    do k=1, nz
    do j=1, ny
    do i=1, nx
      if (u_in(i,j,k,n) /= missv) then
        ltag(i,j,k) = .TRUE.
      else
        ltag(i,j,k) = .FALSE.
      end if
    enddo
    enddo
    enddo

    do k=1, nz
    do j=1, ny
      if ( ltag(nx,j,k) .and. ltag(1,j,k) ) then
        u_out(1,j,k,n) = 0.5*(u_in(nx,j,k,n)+u_in(1,j,k,n))
      else
        u_out(1,j,k,n) = missv
      end if
      do i=2, nx
        if ( ltag(i-1,j,k) .and. ltag(i,j,k) ) then
          u_out(i,j,k,n) = 0.5*(u_in(i-1,j,k,n)+u_in(i,j,k,n))
        else
          u_out(i,j,k,n) = missv
        end if
      enddo
    enddo
    enddo

  enddo


  RETURN

END subroutine u2rho


SUBROUTINE v2rho(nx,ny,nz,nt,v_in,missv, v_out)

  implicit none

  integer,                      intent(in)    ::  nx, ny, nz, nt
  real,                         intent(in)    ::  missv
  real, dimension(nx,ny-1,nz,nt), intent(in)  ::  v_in
  real, dimension(nx,ny  ,nz,nt), intent(out) ::  v_out

  integer ::  i,j,k,n
  logical ::  ltag(nx,ny-1,nz)


  do n=1, nt

    do k=1, nz
    do j=1, ny-1
    do i=1, nx
      if (v_in(i,j,k,n) /= missv) then
        ltag(i,j,k) = .TRUE.
      else
        ltag(i,j,k) = .FALSE.
      end if
    enddo
    enddo
    enddo

    do k=1, nz
      do j=2, ny-1
      do i=1, nx
        if ( ltag(i,j-1,k) .and. ltag(i,j,k) ) then
          v_out(i,j,k,n) = 0.5*(v_in(i,j-1,k,n)+v_in(i,j,k,n))
        else
          v_out(i,j,k,n) = missv
        end if
      enddo
      enddo
      v_out(:, 1,k,n) = 0.
      v_out(:,ny,k,n) = 0.
    enddo

  enddo


  RETURN

END subroutine v2rho


SUBROUTINE t2rho(nx,ny,nz,nt,t_in,missv, t_out)

  implicit none

  integer,                      intent(in)  ::  nx, ny, nz, nt
  real,                         intent(in)  ::  missv
  real, dimension(nx,ny,nz,nt), intent(in)  ::  t_in
  real, dimension(nx,ny,nz,nt), intent(out) ::  t_out

  integer ::  i,j,k,n
  logical ::  ltag(nx,ny,nz)


  do n=1, nt

    do k=1, nz
    do j=1, ny
    do i=1, nx
      if (t_in(i,j,k,n) /= missv) then
        ltag(i,j,k) = .TRUE.
      else
        ltag(i,j,k) = .FALSE.
      end if
    enddo
    enddo
    enddo

    t_out(:,:,1,n) = missv
    do k=2, nz
    do j=1, ny
    do i=1, nx
      if ( ltag(i,j,k-1) .and. ltag(i,j,k) ) then
        t_out(i,j,k,n) = 0.5*(t_in(i,j,k-1,n)+t_in(i,j,k,n))
      else
        t_out(i,j,k,n) = missv
      end if
    enddo
    enddo
    enddo

  enddo


  RETURN

END subroutine t2rho


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
    call grady_2nd(1,ny,nz,1,temp,lat,missv, grady)
    grady( 1,:) = missv
    grady(ny,:) = missv
    do k=1, nz
    do j=2, ny-1
      if ( grady(j,k) /= missv )  grady(j,k) = grady(j,k)/cosphi(j)
    enddo
    enddo

    call gradz_2nd(1,ny,nz,1,var_z(:,:,n),zh,missv, gradz)

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


END module UM_anal


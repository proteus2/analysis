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
! hintp_s
! vintp_p
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


  IF (missv == 0.) then

  do n=1, nt
  do k=1, nz
  do j=1, ny
    varo(j,k,n) = sum(vari(:,j,k,n))/nx
  enddo
  enddo
  enddo

  ELSE

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

  END if


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


  if (missv == 0.) then

    temp(:,:,:,:) = vari1(:,:,:,:)*vari2(:,:,:,:)

  else

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

  end if

  call zonal_avg(nx,ny,nz,nt,missv,temp,varo)


  RETURN

END subroutine corr_zonal_avg


SUBROUTINE hintp_s(nx,ny,nz,lon1,var, var_o)

  implicit none

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  lon1
  real, dimension(nx,ny,nz), intent(in) ::  var

  real, dimension(nx,ny/2*2+1,nz), intent(out) ::  var_o

  integer                   ::  i,j,k
  real, dimension(nx,ny,nz) ::  temp


  if ( lon1 == 0. .and. ny/2*2 /= ny ) then
    var_o(:,:,:) = var(:,:,:)
    RETURN
  end if

  if (lon1 /= 0.) then
    do k=1, nz
    do j=1, ny
    do i=2, nx
      var_o(i,j,k) = 0.5*(var(i-1,j,k)+var(i,j,k))
    enddo
    enddo
    enddo
    var_o(1,1:ny,:) = 0.5*(var(nx,1:ny,:)+var(1,1:ny,:))
  else
    var_o(:,1:ny,:) = var(:,1:ny,:)
  end if

  if (ny/2*2 == ny) then
    temp(:,:,:) = var_o(:,1:ny,:)
    do k=1, nz
    do j=2, ny
      var_o(:,j,k) = 0.5*(temp(:,j-1,k)+temp(:,j,k))
    enddo
    enddo
    var_o(:,1   ,:) = 0.
    var_o(:,ny+1,:) = 0.
  end if


  RETURN

END subroutine hintp_s


SUBROUTINE vintp_p(nx,ny,nz,exner,var,np,exner_o,fill_extp, var_o)

  implicit none

  integer,                   intent(in) ::  nx, ny, nz, np
  real,                      intent(in) ::  fill_extp
  real, dimension(np),       intent(in) ::  exner_o
  real, dimension(nx,ny,nz), intent(in) ::  var, exner

  real, dimension(nx,ny,np), intent(out) ::  var_o

  integer ::  i,j,k,kp
  integer ::  last, k0(nx,ny)


  last = 1


  N_LEV:   DO kp=1, np


  do j=1, ny
  do i=1, nx
    if (exner_o(kp) <= exner(i,j,last) ) then
      ! searching up
      k0(i,j) = -1  ! top level
      do k=last+1, nz
        if (exner_o(kp) > exner(i,j,k)) then
          k0(i,j) = k - 1
          EXIT
        end if
      enddo
    else
      ! searching down
      k0(i,j) = 0  ! bottom level
      do k=last-1, 1, -1
        if (exner_o(kp) <= exner(i,j,k)) then
          k0(i,j) = k
          EXIT
        end if
      enddo

    end if
    last = max(k0(i,j), 1)
  enddo
  enddo

  IF (fill_extp == 1.) then

  do j=1, ny
  do i=1, nx
    if (k0(i,j) == -1) then
      var_o(i,j,kp) = var(i,j,nz)
    else if (k0(i,j) == 0) then
      var_o(i,j,kp) = var(i,j,1)
    else
      var_o(i,j,kp) = ( (exner_o(kp)-exner(i,j,k0(i,j)))*           &
                        var(i,j,k0(i,j)+1) -                        &
                        (exner_o(kp)-exner(i,j,k0(i,j)+1))*         &
                        var(i,j,k0(i,j)) ) /                        &
                      ( exner(i,j,k0(i,j)+1) - exner(i,j,k0(i,j)) )
    end if
  enddo
  enddo

  ELSE

  do j=1, ny
  do i=1, nx
    if (k0(i,j) == -1) then
      var_o(i,j,kp) = fill_extp
    else if (k0(i,j) == 0) then
      var_o(i,j,kp) = fill_extp
    else
      var_o(i,j,kp) = ( (exner_o(kp)-exner(i,j,k0(i,j)))*           &
                        var(i,j,k0(i,j)+1) -                        &
                        (exner_o(kp)-exner(i,j,k0(i,j)+1))*         &
                        var(i,j,k0(i,j)) ) /                        &
                      ( exner(i,j,k0(i,j)+1) - exner(i,j,k0(i,j)) )
    end if
  enddo
  enddo

  END IF  ! fill_extp == 1.


  ENDDO  N_LEV


  RETURN

END subroutine vintp_p


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


  IF (missv == 0.) THEN

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

  ELSE

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

  END if  ! missv == 0.


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


  IF (missv == 0.) THEN

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

  ELSE

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

  END if  ! missv == 0.


  RETURN

END subroutine gradz_2nd


SUBROUTINE gradz_2nd_irr(nx,ny,nz,nt,var,z, gradz)

  implicit none

  integer,                      intent(in)  ::  nx, ny, nz, nt
  real, dimension(nx,ny,nz,nt), intent(in)  ::  var
  real, dimension(nz),          intent(in)  ::  z
  real, dimension(nx,ny,nz,nt), intent(out) ::  gradz

  integer ::  i,j,k,n
  real    ::  dcoef(3,nz), dz1, dzn
  real*8  ::  dcoef8(3)


  do k=2, nz-1
    call fdcoef(2,3,dble(z(k)),dble(z(k-1:k+1)), dcoef8)
    dcoef(:,k) = dcoef8(:)
  enddo
  dz1 = z(2 ) - z(1   )
  dzn = z(nz) - z(nz-1)

  do n=1, nt
  do k=2, nz-1
  do j=1, ny
  do i=1, nx
    gradz(i,j,k,n) = sum( dcoef(:,k)*var(i,j,k-1:k+1,n) )
  enddo
  enddo
  enddo
  enddo

  gradz(:,:,1 ,:) = (var(:,:,2 ,:) - var(:,:,1   ,:)) / dz1
  gradz(:,:,nz,:) = (var(:,:,nz,:) - var(:,:,nz-1,:)) / dzn


  RETURN

END subroutine gradz_2nd_irr


SUBROUTINE grad2z_2nd_irr(nx,ny,nz,nt,var,z,missv, gradz)

  implicit none

  integer,                      intent(in)  ::  nx, ny, nz, nt
  real,                         intent(in)  ::  missv
  real, dimension(nx,ny,nz,nt), intent(in)  ::  var
  real, dimension(nz),          intent(in)  ::  z
  real, dimension(nx,ny,nz,nt), intent(out) ::  gradz

  integer ::  i,j,k,n
  real    ::  dcoef(3,nz)
  real*8  ::  dcoef8(3)


  do k=2, nz-1
    call fdcoef(3,3,dble(z(k)),dble(z(k-1:k+1)), dcoef8)
    dcoef(:,k) = dcoef8(:)
  enddo

  do n=1, nt
  do k=2, nz-1
  do j=1, ny
  do i=1, nx
    gradz(i,j,k,n) = sum( dcoef(:,k)*var(i,j,k-1:k+1,n) )
  enddo
  enddo
  enddo
  enddo

  gradz(:,:,1 ,:) = missv
  gradz(:,:,nz,:) = missv


  RETURN

END subroutine grad2z_2nd_irr


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


! End routines =================================================================


END module UM_anal


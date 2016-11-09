PROGRAM W_continuity_eqn

  use hadgem
  use netio

  implicit none

  integer, parameter ::  nv = 1
  real ::  z_bot

  namelist /ANALCASE/ EXPNAME, YYYY, MM, DD, HH, REFDATE
  namelist /PARAM/ Z_BOT
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, VAR_I_NAME, FILE_O

  integer ::  iz, ihour, i_time, i_time_last
  integer ::  k
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  var5d
  real, dimension(:,:,:),     allocatable ::  u, v, rho, w

  type(vset), dimension(nv) ::  set

  ovarname(1) = 'dz_dt'

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  year = yyyy
  mon  = mm(1)

  call initialize

  i_time = 0

  L_DATE:  DO date=dd(1)-nday_i, dd(1)-1
  !---------------------------------------------------------------------
  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
  i_time = i_time + 1

  day_from_ref = get_dayfromref(year,mon,date,hour)

  ex0 = .TRUE.
  do iv_i=1, 3
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  ! get variable
  allocate( u(nx,ny,nz), v(nx,ny,nz), rho(nx,ny,nz) )

  call get_3var

  ! calculate w
  call w_cont(nx,ny,nz,lat,ht,ht_th(nz),u,v,rho, w)

  w(:,1 ,:) = spread(sum(w(:,2   ,:), dim=1)/float(nx),1,nx)
  w(:,ny,:) = spread(sum(w(:,ny-1,:), dim=1)/float(nx),1,nx)

  deallocate( u, v, rho )

  var5d(:,:,:,i_time,1) = w(:,:,:)

  t(i_time) = day_from_ref

  hour = hour + 24/nhour
  !---------------------------------------------------------------------
  ENDDO  L_HOUR

  !---------------------------------------------------------------------
  ENDDO  L_DATE

  nt = i_time

  nd1a = NX
  nd2a = NY
  nd3a = NZ
  nd4a = NT

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = var5d(:,:,:,:,iv)
  enddo

  if (z_bot /= 0.) then
    iz = nd3a + 1
    do k=2, nd3a
      if (ht(k) > z_bot*1.e3) then  ;  iz = k - 1  ;  EXIT  ;  end if
    enddo
    do iv=1, nv
      set(iv)%var_out(:,:,:iz-1,:) = 1.e32
    enddo
  end if


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'W*RHO retrieved from continuity eqn.')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nhour = hh(2)

  ndate = 100  ;  date = dd(1)-nday_i  ;  hour = hh(1)  ! for get_ifilename
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 3  ! rho
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  nt = nday_i*nhour

  allocate( var5d(nx,ny,nz,nt,nv), t(nt) )
  var5d(:,:,:,:,:) = 0.  ;  t(:) = 0.
  allocate( w(nx,ny,nz) )

  END subroutine initialize

  SUBROUTINE get_3var

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz   ! for get_ivar3d

  ! read 3 var.s
  print*, trim(file_i(1))
  iv_i = 1  ;  u  (:,:,:) = get_ivar3d()
  iv_i = 2  ;  v  (:,:,:) = get_ivar3d()
  print*, trim(file_i(3))
  iv_i = 3  ;  rho(:,:,:) = get_ivar3d()
  print*, 'time index :', it_i(2:3)

  ! corrections - error in u(:,:,1) and rho(:,:,1) at 00 UTC
  u(:,:,1) = u(:,:,2)*0.5
  rho(:,:,1) = exp(log(rho(:,:,2))*2.-log(rho(:,:,3)))

  END subroutine get_3var

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'longitude','latitude','hybrid_ht','t'/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lon
  set(iv)%axis2 = lat
  set(iv)%axis3 = ht
  set(iv)%axis4 = t
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var5d, w )
  deallocate( lon, lat, ht, ht_th, t )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program W_continuity_eqn


SUBROUTINE w_cont(nx,ny,nz,lat,z,zt,u,v,rho, w)

  use nr, only: spline

  implicit none

  integer,                   intent(in) ::  nx, ny, nz
  real,                      intent(in) ::  zt
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  z
  real, dimension(nx,ny,nz), intent(in) ::  u, v, rho

  real, dimension(nx,ny,nz), intent(out) ::  w

  integer                   ::  i,j,k
  real                      ::  d2lon
  real, dimension(ny)       ::  d2lat, coslat
  real, dimension(nz)       ::  rz, rz2
  real, dimension(nz*2)     ::  z2, divh2, divhpp2
  real, dimension(nx,ny,nz) ::  mu, mv, divh

  real, parameter ::  a_earth = 6371229.
  real, parameter ::  pi = 3.14159265358979323846
  real, parameter ::  deg2rad = pi/180.

  d2lon = 2.*(2.*pi)/float(nx)
  d2lat(2:ny-1) = (lat(3:ny) - lat(1:ny-2))*deg2rad
  coslat(:) = cos(lat(:)*deg2rad)
  if (abs(lat(1 )) == 90.)  coslat(1 ) = 0.
  if (abs(lat(ny)) == 90.)  coslat(ny) = 0.
  rz(:) = a_earth + z(:)
  rz2(:) = rz(:)*rz(:)

  mu(:,:,:) = rho(:,:,:)*u(:,:,:)
  mv(:,:,:) = rho(:,:,:)*v(:,:,:)
  do j=1, ny
    mv(:,j,:) = mv(:,j,:)*coslat(j)
  enddo
  do k=1, nz
    mu(:,:,k) = mu(:,:,k)*rz(k)
    mv(:,:,k) = mv(:,:,k)*rz(k)
  enddo

  divh(2:nx-1,:,:) = (mu(3:nx,:,:) - mu(1:nx-2,:,:))/d2lon
  divh(1     ,:,:) = (mu(2   ,:,:) - mu(nx    ,:,:))/d2lon
  divh(nx    ,:,:) = (mu(1   ,:,:) - mu(nx-1  ,:,:))/d2lon
  do k=1, nz
  do j=2, ny-1
    divh(:,j,k) = (divh(:,j,k) + (mv(:,j+1,k) - mv(:,j-1,k))/d2lat(j)) &
                  / coslat(j)
  enddo
  enddo

  do j=2, ny-1
  do i=1, nx
    z2(1   :nz  ) = z(:)
    z2(nz+1:nz*2) = zt*2. - z(nz:1:-1)
    divh2(1   :nz  ) = divh(i,j,  :    )
    divh2(nz+1:nz*2) = divh(i,j,nz:1:-1)
    call spline(z2,divh2,1.e32,1.e32,divhpp2)
    do k=1, nz
      call intg_cubicsp(nz*2,z2,divh2,divhpp2,1,z2(k),z2(k+1), w(i,j,k))
    enddo
  enddo
  enddo
  w(:,:,nz) = 0.5*w(:,:,nz)
  do k=nz-1, 1, -1
    w(:,:,k) = w(:,:,k+1) + w(:,:,k)
  enddo

  do k=1, nz
    w(:,:,k) = w(:,:,k)/rz2(k)
  enddo
  w(:,:,:) = w(:,:,:)/rho(:,:,:)

  w(:,1 ,:) = 1.e32
  w(:,ny,:) = 1.e32

END subroutine w_cont

SUBROUTINE intg_cubicsp(na,xa,ya,y2a,nf,xi,xf,f)

  implicit none

  integer,             intent(in) ::  na, nf
  real,                intent(in) ::  xi
  real, dimension(na), intent(in) ::  xa, ya, y2a
  real, dimension(nf), intent(in) ::  xf

  real, dimension(nf), intent(out) ::  f

  integer ::  i,j, j0i
  integer, dimension(0:nf) ::  j0
  real,    dimension(na)   ::  d, c1a, c2a, fa
  real,    dimension(0:nf) ::  xif, aa, bb, c1, c2, fsp

  if ( (xf(1) - xi)*(xf(nf) - xf(1)) < 0. ) then
    print*, 'INTG_CUBICSP: XI, ..., XF(NF) should be monotonic.' ; STOP
  end if
  if ( (xa(na) - xa(1))*(xf(nf) - xi) < 0. ) then
    print*, 'INTG_CUBICSP: Direction of XA and XF should be same.' ; STOP
  end if
  if ( (xa(1) - xi)*(xa(1) - xa(2)) < 0. .or.                            &
       (xa(na) - xf(nf))*(xa(na) - xa(na-1)) < 0. ) then
    print*, 'INTG_CUBICSP: XI and XF should be in [XA(1),XA(NA)].' ; STOP
  end if

  xif(0) = xi
  xif(1:nf) = xf(:)
  do i=0, nf
    do j=na, 2, -1
      if ( (xif(i)-xa(j))*(xif(i)-xa(j-1)) <= 0. ) then
        j0(i) = j - 1  ;  EXIT  ! take min(j, na-1), when xif(i) = xa(j)
      end if
    enddo
  enddo

  do j=j0(0), j0(nf)
    d(j) = xa(j+1) - xa(j)
  enddo

  do j=j0(0), j0(nf)-1
    c1a(j) = 0.5*d(j)
    c2a(j) = d(j)*d(j)*d(j)/12.
  enddo
  fa(j0(0)) = 0.
  do j=j0(0), j0(nf)-1
    fa(j+1) = fa(j) + c1a(j)*(ya (j)+ya (j+1)) +                         &
                      c2a(j)*(y2a(j)+y2a(j+1))*(-0.5)
  enddo

  do i=0, nf
    aa(i) = (xa(j0(i)+1) - xif(i))/d(j0(i))
    bb(i) = (xif(i) - xa(j0(i))  )/d(j0(i))
    c1(i) = 0.5*d(j0(i))
    c2(i) = d(j0(i))**3/12.
  enddo
  aa(:) = aa(:)*aa(:)  ;  bb(:) = bb(:)*bb(:)
  if ( sum(bb(:)) == 0. ) then
    fsp(:) = 0.
  else
    do i=0, nf
      fsp(i) = c1(i)*( (1.-aa(i))*ya(j0(i)) + bb(i)*ya(j0(i)+1) ) +      &
               c2(i)*( ((1.-0.5*aa(i))*aa(i)-0.5)*y2a(j0(i)) +           &
                       (0.5*bb(i)-1.)*bb(i)*y2a(j0(i)+1) )
    enddo
  end if

  f(:) = fa(j0(1:nf)) + fsp(1:nf) - fsp(0)

END subroutine intg_cubicsp


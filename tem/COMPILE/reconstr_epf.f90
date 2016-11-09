PROGRAM RC_EPF_RA

  use reanal
  use netio

  implicit none

  integer, parameter ::  nv = 4, nrc = 10

  integer ::  opt_30d, k_largest, period_smallest, nmon_patch

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ K_LARGEST, PERIOD_SMALLEST, LAT_RNG, P_RNG, NMON_PATCH
  namelist /FILEIO/ NT_F4, P_PREDEF, FILE_I_HEAD, FILE_I_FORM,           &
                    FILE_I_XXXX, VAR_I_NAME, FILE_I_HEAD2, FILE_I_FORM2, &
                    FILE_I_XXXX2, VAR_I2, VAR_I_NAME2, FILE_O

  integer ::  iz, imon, ihour, i_time
  integer ::  nk, nome
  integer ::  iy1(2), iz1(2), iy2(2), iz2(2), ny2, nz2, nta
  integer ::  j,k,n, tmpi
  real    ::  p_predef1, wgt
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  var4pn
  real, dimension(:,:,:,:), allocatable ::  varo, var4d
  real, dimension(:,:,:),   allocatable ::  kr, o_hat
  real, dimension(:,:),     allocatable ::  um
  real, dimension(:),       allocatable ::  lat0, p0a, kwn, ome, or, ff

  type(vset), dimension(nv) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  do iv=1, nv
    ovarname(iv) = trim(var_i_name(iv))
  enddo

! GET AXES AND INITIALIZE ARRAYS

  call initialize

  ! get mean U
  call switch_para_in
  allocate( um(ny2,nz2) )
  um(:,:) = 0.

  deallocate( lon, lat, p, dim4, t2pt )
  iv_i = 1
  file_i(iv_i) = get_ifilename()
  p_predef(1) = p_predef1
  call getdim(file_i(iv_i),var_i_name(iv_i))
  call get_iouter(lat    ,lat_rng    , iy1)
  call get_iouter(p*(-1.),p_rng*(-1.), iz1)
  if ( any(lat(iy1(1):iy1(2)) /= lat0) .or. &
       any(p  (iz1(1):iz1(2)) /= p0a ) ) then
    print*, 'Check lat0, p0'  ;  STOP
  end if
  if ( l_rev(2) ) then
    iy1(:) = ny + 1 - iy1(:)
    tmpi = iy1(1)
    iy1(1) = iy1(2)  ;  iy1(2) = tmpi
  end if
  if ( l_rev(3) ) then
    iz1(:) = nz + 1 - iz1(:)
    tmpi = iz1(1)
    iz1(1) = iz1(2)  ;  iz1(2) = tmpi
  end if

  i_time = 0

  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------

  ndate = get_ndate()

  do date=1, ndate
    hour = hh(1)
    do ihour=1, nhour

      i_time = i_time + 1

      call get_1var

      hour = hour + 24/nhour

    enddo
  enddo

  mon = mon + 1
  if (mon == 13) then
    year = year + 1  ;  mon = 1
  end if

  !---------------------------------------------------------------------
  ENDDO  L_MON

  nt = i_time

  um(:,:) = um(:,:)/float(nt)

  ! rewind
  year = yyyy
  mon  = mm(1)

  call switch_para_in

  allocate( var4d    (nk,-nome+1:nome,ny2,nv)   )
  allocate( var4pn   (nk,-nome+1:nome,ny2,nv,2) )
  allocate( o_hat    (nk,-nome+1:nome,ny2)      )

  L_LEV:  DO k=1, nz2

  call get_4var

  print*, ' Calculating...', k, '/', nz2
  call reconstr
  print*, ' .'

  ENDDO  L_LEV

  deallocate( var4d, var4pn, o_hat )

  nd1a = NY2
  nd2a = NZ2
  nd3a = NRC
  nd4a = 1

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,1) = varo(:,:,iv,:)
  enddo

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'Reconstructed EP flux')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  real, parameter ::  a_earth = 6371229.
  real, parameter ::  pi = 3.14159265358979323846
  real, parameter ::  ome_earth = 7.292116e-5
  real, parameter ::  deg2rad = pi/180.

  year = yyyy
  mon  = mm(1)

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()

  iv_i = 1
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  p_predef1 = p_predef(1)
  p_predef(1) = -999.
  call getdim(file_i(iv_i),var_i_name(iv_i))
  l_rev(3) = ( .not. l_rev(3) )
  p(:) = p(nz:1:-1)

  call get_iouter(p         ,lat_rng    , iy2)
  call get_iouter(dim4*(-1.),p_rng*(-1.), iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  allocate( lat0(ny2), p0a(nz2) )
  lat0(:) = p(iy2(1):iy2(2))  ;  p0a(:) = dim4(iz2(1):iz2(2))

  nta = ny
  if ( nta /= (nmon+2*nmon_patch)*30*nhour .or. opt_30d == 0 ) then
    print*, 'not yet programmed (for opt_30d = 0)'  ;  STOP
  end if

  nk = nx - 1
  if (k_largest /= -999)  nk = k_largest

  nome = nta/2
  if (period_smallest /= -999)  nome = nta/(period_smallest*nhour)

  allocate( kwn(nk), ome(-nome+1:nome) )
  kwn(:) = lon(2:nk+1)
  ome(0:nome) = lat(1:nome+1)
  ome(:-1   ) = ome(nome-1:1:-1)*(-1.)

  allocate( kr(nk,ny2,nz2), or(-nome+1:nome) )
  do k=1, nz2
  do j=1, ny2
    kr(:,j,k) = kwn(:)/(a_earth*cos(lat0(j)*deg2rad))
  enddo
  enddo
  or(:) = (-1.)*ome(:)*(2.*pi/(float(nmon+2*nmon_patch)*30.*86400.))

  allocate( ff(ny2) )
  do j=1, ny2
    ff(j) = (2.*ome_earth*sin(lat0(j)*deg2rad))**2
  enddo

  allocate( varo(ny2,nz2,nv,nrc) )

  END subroutine initialize

  SUBROUTINE get_1var

  real, dimension(nx,ny2,nz2) ::  var3d

  iv_i = 1
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if

  print*, trim(file_i(1))
  iv_i = 1  ;  var3d(:,:,:) = get_ivara3d(1,nx,iy1(1),ny2,iz1(1),nz2)
  um(:,:) = um(:,:) + sum(var3d(:,:,:), dim=1)/float(nx)
  print*, 'time index :', it_i(1)

  END subroutine get_1var

  SUBROUTINE get_4var

  ex0 = .TRUE.
  do iv_i=1, nv
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
    print*, trim(file_i(iv_i))
  enddo
  if (.not. ex0)  STOP

  ! read 4 var.s
  l_rev(:) = .False.
  do iv_i=1, nv
    var4d(:,0:nome,:,iv_i:iv_i) = get_ivara4d(2,nk,1         ,nome+1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
    var4d(:, :-1  ,:,iv_i:iv_i) = get_ivara4d(2,nk,nta-nome+2,nome-1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
  enddo
  var4d(:,:,:,:) = var4d(:,:,:,:)*2.
  var4d(nk,:   ,:,:) = 0.
  var4d(: ,nome,:,:) = 0.

  END subroutine get_4var

  SUBROUTINE reconstr

  integer ::  k1, k2, o1, o2

  do j=1, ny2
  do n=-nome+1,nome
    o_hat(:,n,j) = or(n) - um(j,k)*kr(:,j,k)
  enddo
  enddo

  var4pn(:,:,:,:,1) = var4d(:,:,:,:)
  var4pn(:,:,:,:,2) = var4d(:,:,:,:)
  where ( o_hat > 0. )
    var4pn(:,:,:,1,2) = 0.
    var4pn(:,:,:,2,2) = 0.
    var4pn(:,:,:,3,2) = 0.
    var4pn(:,:,:,4,2) = 0.
  else where
    var4pn(:,:,:,1,1) = 0.
    var4pn(:,:,:,2,1) = 0.
    var4pn(:,:,:,3,1) = 0.
    var4pn(:,:,:,4,1) = 0.
  end where

! total eastward/westward waves
  varo(:,k,:,1) = sum(sum(var4pn(:,:,:,:,1), dim=1), dim=1)
  varo(:,k,:,2) = sum(sum(var4pn(:,:,:,:,2), dim=1), dim=1)

! diurnal migrating tide (positive ome : westward)
  k1 = 1  ;  k2 = 1
  o1 = (nmon+2*nmon_patch)*30 - 1  ;  o2 = (nmon+2*nmon_patch)*30 + 1
  varo(:,k,:,3) = sum(sum(var4pn(k1:k2,o1:o2,:,:,2), dim=1), dim=1)
  var4pn(k1:k2,o1:o2,:,:,2) = 0.

  if ( nome > (nmon+2*nmon_patch)*60 ) then
    ! semi-diurnal migrating tide
    k1 = 2  ;  k2 = 2
    o1 = (nmon+2*nmon_patch)*60 - 1  ;  o2 = (nmon+2*nmon_patch)*60 + 1
    varo(:,k,:,4) = sum(sum(var4pn(k1:k2,o1:o2,:,:,2), dim=1), dim=1)
    var4pn(k1:k2,o1:o2,:,:,2) = 0.
  end if

! eastward/westward large-scale waves ( k = 1-6, period > 4/3 days )
  k1 = 1  ;  k2 = 6
  o1 = -(nmon+2*nmon_patch)*22 + 1  ;  o2 = (nmon+2*nmon_patch)*22 - 1
  varo(:,k,:,5) = sum(sum(var4pn(k1:k2,o1:o2,:,:,1), dim=1), dim=1)
  varo(:,k,:,6) = sum(sum(var4pn(k1:k2,o1:o2,:,:,2), dim=1), dim=1)
  var4pn(k1:k2,o1:o2,:,:,:) = 0.

! eastward/westward large-scale waves II ( k = 7-20, period > 4/3 days )
  k1 = 7  ;  k2 = 20
  o1 = -(nmon+2*nmon_patch)*22 + 1  ;  o2 = (nmon+2*nmon_patch)*22 - 1
  varo(:,k,:,7) = sum(sum(var4pn(k1:k2,o1:o2,:,:,1), dim=1), dim=1)
  varo(:,k,:,8) = sum(sum(var4pn(k1:k2,o1:o2,:,:,2), dim=1), dim=1)
  var4pn(k1:k2,o1:o2,:,:,:) = 0.

! GWs
  do iv=1, nv
    varo(:,k,:,9 ) = sum(sum(var4pn(:,:,:,:,1), dim=1), dim=1)
    varo(:,k,:,10) = sum(sum(var4pn(:,:,:,:,2), dim=1), dim=1)
  end do

  do iv=1, nv
    if (var4d(3,3,1    ,iv) == 2.e32)  varo(1    ,k,iv,:) = 1.e32
    if (var4d(3,3,ny2  ,iv) == 2.e32)  varo(ny2  ,k,iv,:) = 1.e32
    if (var4d(3,3,2    ,iv) == 2.e32)  varo(2    ,k,iv,:) = 1.e32
    if (var4d(3,3,ny2-1,iv) == 2.e32)  varo(ny2-1,k,iv,:) = 1.e32
  end do
  if ( k <= 2 .or. k >= nz2-1 ) then
    do iv=1, nv
      if (var4d(3,3,3,iv) == 2.e32)  varo(:,k,iv,:) = 1.e32
    enddo
  end if

  END subroutine reconstr

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lat   ','p ','wg',' '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat0
  set(iv)%axis2 = p0a
  set(iv)%axis3 = (/1,2,3,4,5,6,7,8,9,10/)*1.
  set(iv)%axis4 = -999.

  END subroutine setdim

  SUBROUTINE finalize

  deallocate( varo, um )
  deallocate( lon, lat, p, lat0, p0a, ff )
  deallocate( kwn, ome, kr, or )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program RC_EPF_RA


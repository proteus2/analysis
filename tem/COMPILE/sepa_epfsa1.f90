PROGRAM RC_EPFSA_RA

  use reanal
  use netio

  implicit none

  integer, parameter ::  nv = 8, nrc = 22

  integer ::  opt_30d, k_largest, period_smallest, nmon_patch

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ K_LARGEST, PERIOD_SMALLEST, LAT_RNG, P_RNG, NMON_PATCH
  namelist /FILEIO/ NT_F4, P_PREDEF, FILE_I_HEAD, FILE_I_FORM,           &
                    FILE_I_XXXX, VAR_I_NAME, FILE_I_HEAD2, FILE_I_FORM2, &
                    FILE_I_XXXX2, VAR_I2, VAR_I_NAME2, FILE_O, NL_AUX

  integer ::  iz, imon, ihour, i_time
  integer ::  nk, nome
  integer ::  iy1(2), iz1(2), iy2(2), iz2(2), ny2, nz2, nta
  integer ::  i,j,k,n, tmpi
  real    ::  p_predef1, wgt
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:), allocatable ::  varo, var4d, var4d2, vtmp
  real, dimension(:,:,:),   allocatable ::  kr
  real, dimension(:,:),     allocatable ::  um, okbi, lzrgw, or2d
  real, dimension(:,:),     allocatable ::  tmp1, tmp2, tmp3, tmp4
  real, dimension(:),       allocatable ::  lat0, p0a, kwn, ome, or, kr0

  type(vset), dimension(nv) ::  set

  real, parameter ::  nbv = 2.5e-2
  real, parameter ::  lz0_rgw = 2.4e3

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  do iv=1, nv
    ovarname(iv) = trim(var_i_name(iv)(1:len_trim(var_i_name(iv))-2))
  enddo

! GET AXES AND INITIALIZE ARRAYS

  call initialize

!U  ! get mean U
!U  call switch_para_in
!U  allocate( um(ny2,nz2) )
!U  um(:,:) = 0.
!U
!U  deallocate( lon, lat, p, dim4, t2pt )
!U  iv_i = 1
!U  file_i(iv_i) = get_ifilename()
!U  p_predef(1) = p_predef1
!U  call getdim(file_i(iv_i),var_i_name(iv_i))
!U  call get_iouter(lat    ,lat_rng    , iy1)
!U  call get_iouter(p*(-1.),p_rng*(-1.), iz1)
!U  if ( any(lat(iy1(1):iy1(2)) /= lat0) .or. &
!U       any(p  (iz1(1):iz1(2)) /= p0a ) ) then
!U    print*, 'Check lat0, p0'  ;  STOP
!U  end if
!U  if ( l_rev(2) ) then
!U    iy1(:) = ny + 1 - iy1(:)
!U    tmpi = iy1(1)
!U    iy1(1) = iy1(2)  ;  iy1(2) = tmpi
!U  end if
!U  if ( l_rev(3) ) then
!U    iz1(:) = nz + 1 - iz1(:)
!U    tmpi = iz1(1)
!U    iz1(1) = iz1(2)  ;  iz1(2) = tmpi
!U  end if
!U
!U  i_time = 0
!U
!U  L_MON:  DO imon=1, nmon
!U  !---------------------------------------------------------------------
!U
!U  ndate = get_ndate()
!U
!U  do date=1, ndate
!U    hour = hh(1)
!U    do ihour=1, nhour
!U
!U      i_time = i_time + 1
!U
!U      call get_1var
!U
!U      hour = hour + 24/nhour
!U
!U    enddo
!U  enddo
!U
!U  mon = mon + 1
!U  if (mon == 13) then
!U    year = year + 1  ;  mon = 1
!U  end if
!U
!U  !---------------------------------------------------------------------
!U  ENDDO  L_MON
!U
!U  nt = i_time
!U
!U  um(:,:) = um(:,:)/float(nt)
!U
!U
!U  ! rewind
!U  year = yyyy
!U  mon  = mm(1)
!U
!U  call switch_para_in

  allocate( var4d    (nk,-nome+1:nome,ny2,nv*2)   )
  allocate( var4d2   (nk,-nome+1:nome,ny2,nv  )   )
  allocate( vtmp     (nk,-nome+1:nome,ny2,4) )
  allocate( okbi     (nk,-nome+1:nome) )
  allocate( lzrgw    (nk,-nome+1:nome) )
  allocate( tmp1(ny2,nv*2), tmp2(ny2,nv*2) )

  L_LEV:  DO k=1, nz2

  call get_4var

  print*, ' Calculating...', k, '/', nz2
  call sepa_sum
  print*, ' .'

  ENDDO  L_LEV

  deallocate( tmp1, tmp2 )
  deallocate( var4d, vtmp, okbi, lzrgw )

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
  real, parameter ::  deg2rad = pi/180.

  year = yyyy
  mon  = mm(1)

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()

  iv_i = nv*2
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  p_predef1 = p_predef(1)
  p_predef(1) = -999.
  call getdim(file_i(iv_i),trim(var_i_name(iv_i)))
  l_rev(3) = ( .not. l_rev(3) )
  p(:) = p(nz:1:-1)

  call get_iouter(p         ,lat_rng    , iy2)
  call get_iouter(dim4*(-1.),p_rng*(-1.), iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  allocate( lat0(ny2), p0a(nz2) )
  lat0(:) = p(iy2(1):iy2(2))  ;  p0a(:) = dim4(iz2(1):iz2(2))

  nta = ny
  if ( opt_30d == 0 .and. nmon /= 1 ) then
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

  allocate( kr(nk,ny2,nz2), or(-nome+1:nome), or2d(nk,-nome+1:nome) )
  allocate( kr0(nk) )
  do k=1, nz2
  do j=1, ny2
    kr(:,j,k) = kwn(:)/(a_earth*cos(lat0(j)*deg2rad))
  enddo
  enddo
  or(:) = (-1.)*ome(:)*(2.*pi/(float(nmon+2*nmon_patch)*30.*86400.))
  or2d(:,:) = spread(or,1,nk)
  kr0(:) = kwn(:)/a_earth

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
  do iv_i=1, nv*4
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
!    print*, trim(file_i(iv_i))
  enddo
  if (.not. ex0)  STOP

  ! read 16 var.s
  l_rev(:) = .False.
  do iv_i=1, nv*2
    var4d(:,0:nome,:,iv_i:iv_i) = get_ivara4d(2,nk,1         ,nome+1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
    var4d(:, :-1  ,:,iv_i:iv_i) = get_ivara4d(2,nk,nta-nome+2,nome-1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
  enddo
  var4d(:,:,:,:) = var4d(:,:,:,:)*2.
  var4d(nk,:   ,:,:) = 0.
  var4d(: ,nome,:,:) = 0.

  ! read 16 var.s
  do iv_i=nv*2+1, nv*3
    iv = iv_i - nv*2
    var4d2(:,0:nome,:,iv:iv) = get_ivara4d(2,nk,1         ,nome+1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
    var4d2(:, :-1  ,:,iv:iv) = get_ivara4d(2,nk,nta-nome+2,nome-1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
  enddo
  do iv_i=nv*3+1, nv*4
    iv = iv_i - nv*3
    var4d2(:,0:nome,:,iv:iv) = var4d2(:,0:nome,:,iv:iv) +                &
       get_ivara4d(2,nk,1         ,nome+1,iy2(1),ny2,iz2(1)-1+k,1)
    var4d2(:, :-1  ,:,iv:iv) = var4d2(:, :-1  ,:,iv:iv) +                &
       get_ivara4d(2,nk,nta-nome+2,nome-1,iy2(1),ny2,iz2(1)-1+k,1)
  enddo
  var4d2(:,:,:,:) = var4d2(:,:,:,:)*2.
  var4d2(nk,:   ,:,:) = 0.
  var4d2(: ,nome,:,:) = 0.

  END subroutine get_4var

  SUBROUTINE sepa_sum

  integer ::  k1, k2, o1a, o2a, o1r, o2r, o1m, o2m

  real, parameter ::  pi = 3.14159265358979323846
  real, parameter ::  beta = 2.*(7.292116e-5)/6371229.

  okbi(:,:) = or2d(:,:)*spread(kr0(:),2,nome*2)/beta

  lzrgw = -999.e3
  where ( okbi > -1. )
    lzrgw(:,:) = 2.*pi*or2d(:,:)**2/(nbv*beta*(1.+okbi(:,:)))
  end where

! low-pass filter for large-scale waves ( k = 1-20, period > 4/3 days )
! This should be same with sepa_epf.f90
  k1 = 1  ;  k2 = 20
  o1a = -(nmon+2*nmon_patch)*22 + 1  ;  o2a = (nmon+2*nmon_patch)*22 - 1
  var4d(k1:k2,:o1a-1,:,:) = 0.  ;  var4d2(k1:k2,:o1a-1,:,:) = 0.
  var4d(k1:k2,o2a+1:,:,:) = 0.  ;  var4d2(k1:k2,o2a+1:,:,:) = 0.
  if ( k2 < nk ) then
    var4d(k2+1:,:,:,:) = 0.  ;  var4d2(k2+1:,:,:,:) = 0.
  end if

! low-freq. IGW(7,8), RW(4) ( period <, >= 2.5 days ) step 1
  o1r = -(nmon+2*nmon_patch)*12  ;  o2r = (nmon+2*nmon_patch)*12

  ! IGW: symm + antisymm
  varo(:,k,:,7) = sum(sum(var4d(:,:o1r-1,:,1:8 ), dim=1), dim=1) + &
                  sum(sum(var4d(:,:o1r-1,:,9:16), dim=1), dim=1)
  varo(:,k,:,8) = sum(sum(var4d(:,o2r+1:,:,1:8 ), dim=1), dim=1) + &
                  sum(sum(var4d(:,o2r+1:,:,9:16), dim=1), dim=1)
  ! IGW: symm-antisymm cross-correl.
  varo(:,k,:,17) = sum(sum(var4d2(:,:o1r-1,:,:), dim=1), dim=1)
  varo(:,k,:,18) = sum(sum(var4d2(:,o2r+1:,:,:), dim=1), dim=1)
  ! RW: symm + antisymm
  varo(:,k,:,1) = sum(sum(var4d(:,o1r:-1,:,1:8 ), dim=1), dim=1) + &
                  sum(sum(var4d(:,o1r:-1,:,9:16), dim=1), dim=1)
  varo(:,k,:,2) = sum(sum(var4d(:,0 :o2r,:,1:8 ), dim=1), dim=1) + &
                  sum(sum(var4d(:,0 :o2r,:,9:16), dim=1), dim=1)
  varo(:,k,:,4) = varo(:,k,:,1) + varo(:,k,:,2)
  ! RW: symm-antisymm cross-correl.
  varo(:,k,:,11) = sum(sum(var4d2(:,o1r:-1,:,:), dim=1), dim=1)
  varo(:,k,:,12) = sum(sum(var4d2(:,0 :o2r,:,:), dim=1), dim=1)
  varo(:,k,:,14) = varo(:,k,:,11) + varo(:,k,:,12)

! total large-scale waves
  varo(:,k,:,1) = varo(:,k,:,1) + varo(:,k,:,7)
  varo(:,k,:,2) = varo(:,k,:,2) + varo(:,k,:,8)
  varo(:,k,:,11) = varo(:,k,:,11) + varo(:,k,:,17)
  varo(:,k,:,12) = varo(:,k,:,12) + varo(:,k,:,18)

! Fs_H, Fs_M, Fa_H, Fa_M
  vtmp = 0.
  do i=-2, 2
  do j=4, ny2-3
    vtmp(:,:,j,1) = vtmp(:,:,j,1) + var4d(:,:,j+i,2 )
    vtmp(:,:,j,2) = vtmp(:,:,j,2) + var4d(:,:,j+i,6 )
    vtmp(:,:,j,3) = vtmp(:,:,j,3) + var4d(:,:,j+i,10)
    vtmp(:,:,j,4) = vtmp(:,:,j,4) + var4d(:,:,j+i,14)
  enddo
  enddo
  vtmp(:,:,:,:) = vtmp(:,:,:,:)*0.2
  vtmp(:,:,:,1) = vtmp(:,:,:,1) - vtmp(:,:,:,2)
  vtmp(:,:,:,3) = vtmp(:,:,:,3) - vtmp(:,:,:,4)

! KW(3)
  do n=-nome+1, -1
  do i=1, nk
    ! |F_uw/F_vT| < 1 : Filtering from KW
    do j=ny2/2+1, ny2-3
      if ( abs(vtmp(i,n,j,2)/vtmp(i,n,j,1)) < 1. .or. &
           abs(vtmp(i,n,ny2+1-j,2)/vtmp(i,n,ny2+1-j,1)) < 1. ) then
        var4d(i,n,j:,1:8) = 0.
        var4d(i,n,:ny2+1-j,1:8) = 0.
        EXIT
      end if
    enddo
  enddo
  enddo
  var4d(:,-nome+1:-1,1:3      ,1:8) = 0.
  var4d(:,-nome+1:-1,ny2-2:ny2,1:8) = 0.

  tmp1(:,:) = sum(sum(var4d(:,o1r:-1,:,1:8), dim=1), dim=1)
  tmp2(:,:) = sum(sum(var4d(:,:o1r-1,:,1:8), dim=1), dim=1)
  varo(:,k,:,3) = tmp1(:,:) + tmp2(:,:)

  ! RW(4), IGW(7) step 2
  varo(:,k,:,4) = varo(:,k,:,4) - tmp1(:,:)
  varo(:,k,:,7) = varo(:,k,:,7) - tmp2(:,:)

! MRGW(5,6) ( period >= 2 days )
  o1m = -(nmon+2*nmon_patch)*15  ;  o2m = (nmon+2*nmon_patch)*15

  do n=o1m, o2m
  do i=1, nk
    if (lzrgw(i,n) <= lz0_rgw) then
      var4d(i,n,:,9:16) = 0.
!    else
!      ! F_uw/F_vT > 0 : Filtering from MRGW
!      do j=ny2/2+1, ny2-3
!        if ( vtmp(i,n,j,4)/vtmp(i,n,j,3) > 0. .or. &
!             vtmp(i,n,ny2+1-j,4)/vtmp(i,n,ny2+1-j,3) > 0. ) then
!          var4d(i,n,j:,9:16) = 0.
!          var4d(i,n,:ny2+1-j,9:16) = 0.
!          EXIT
!        end if
!      enddo
    end if
  enddo
  enddo
  var4d(:,o1m:o2m,1:3      ,9:16) = 0.
  var4d(:,o1m:o2m,ny2-2:ny2,9:16) = 0.

  varo(:,k,:,5) = sum(sum(var4d(:,o1r:-1,:,9:16), dim=1), dim=1)
  varo(:,k,:,6) = sum(sum(var4d(:,0:o2r ,:,9:16), dim=1), dim=1)
  ! o1r > o1m ; o2r < o2m

  ! RW(4) step 3
  varo(:,k,:,4) = varo(:,k,:,4) - (varo(:,k,:,5) + varo(:,k,:,6))

  tmp1(:,:) = sum(sum(var4d(:,o1m:o1r-1,:,9:16), dim=1), dim=1)
  tmp2(:,:) = sum(sum(var4d(:,o2r+1:o2m,:,9:16), dim=1), dim=1)

  varo(:,k,:,5) = varo(:,k,:,5) + tmp1(:,:)
  varo(:,k,:,6) = varo(:,k,:,6) + tmp2(:,:)

  ! IGW(7,8) step 3
  varo(:,k,:,7) = varo(:,k,:,7) - tmp1(:,:)
  varo(:,k,:,8) = varo(:,k,:,8) - tmp2(:,:)

! temporary : (S+A) + symm-antisymm cross-correl.
  ! total large-scale waves
  varo(:,k,:,9 ) = varo(:,k,:,1) + varo(:,k,:,11)
  varo(:,k,:,10) = varo(:,k,:,2) + varo(:,k,:,12)
  ! E+W
  varo(:,k,:,19) = varo(:,k,:,9) + varo(:,k,:,10)
  ! RW
  varo(:,k,:,20) = varo(:,k,:,4) + varo(:,k,:,14)
  ! IGW
  varo(:,k,:,21) = varo(:,k,:,7) + varo(:,k,:,17)
  varo(:,k,:,22) = varo(:,k,:,8) + varo(:,k,:,18)

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

  END subroutine sepa_sum

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
  set(iv)%axis3 = (/1,2,3,4,5,6,7,8,9,10,                                &
                    11,12,13,14,15,16,17,18,19,20,21,22/)*1.
  set(iv)%axis4 = -999.

  END subroutine setdim

  SUBROUTINE finalize

!U  deallocate( um )
  deallocate( varo )
  deallocate( lon, lat, p, lat0, p0a )
  deallocate( kwn, ome, kr, or, or2d, kr0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program RC_EPFSA_RA


PROGRAM RC_EPFSA_UM_z

  use hadgem
  use netio

  implicit none

  integer, parameter ::  nv = 8, nrc = 11, nrco = 6

  integer ::  k_largest, period_smallest, nmon_patch

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ K_LARGEST, PERIOD_SMALLEST, LAT_RNG, Z_RNG, NMON_PATCH
  namelist /FILEIO/ DAY1, NDAY_I, FID, FILE_I_HEAD, FILE_I_FORM,         &
                    FILE_I_XXXX, VAR_I_NAME, FILE_I_HEAD2, FILE_I_FORM2, &
                    FILE_I_XXXX2, VAR_I_NAME2, FILE_O, NL_AUX

  integer ::  iz, imon, ihour, i_time
  integer ::  nk, nome, rcs(nrco), orcs(nrco)
  integer ::  iy1(2), iz1(2), iy2(2), iz2(2), ny2, nz2, nta
  integer ::  i,j,k,n, irc
  real    ::  wgt
  character(len=32), dimension(nv) ::  ovarname
  character(len=128) ::  file_o_rc, title_o_rc

  real, dimension(:,:,:,:,:,:), allocatable ::  var6d
  real, dimension(:,:,:,:), allocatable ::  var4d, var4d2, vtmp
  real, dimension(:,:,:),   allocatable ::  kr
  real, dimension(:,:),     allocatable ::  um, okbi, lzrgw, or2d
  real, dimension(:,:),     allocatable ::  tmp1, tmp2, tmp3, tmp4
  real, dimension(:),       allocatable ::  lat0, ht0, kwn, ome, or, kr0

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
!U  deallocate( lon, lat, ht, ht_th, dim4 )
!U  iv_i = 1
!U  file_i(iv_i) = get_ifilename()
!U  call getdim(file_i(iv_i),var_i_name(iv_i))
!U  call get_iouter(lat    ,lat_rng, iy1)
!U  call get_iouter(ht/1.e3,z_rng  , iz1)
!U  if ( any(lat(iy1(1):iy1(2)) /= lat0) .or. &
!U       any(ht (iz1(1):iz1(2)) /= ht0 ) ) then
!U    print*, 'Check lat0, ht0'  ;  STOP
!U  end if
!U
!U  do imon=1, nmon
!U    iv_i = 1
!U    file_i(iv_i) = get_ifilename()  ;  print*, trim(file_i(1))
!U    um(:,:) = um(:,:) + reshape(                                         &
!U       sum(get_ivara4d(1,nx,iy1(1),ny2,iz1(1),nz2,1,1),dim=1)/float(nx), &
!U       (/ny2,nz2/) )
!U    mon = mon + 1
!U    if (mon == 13) then
!U      year = year + 1  ;  mon = 1
!U    end if
!U  enddo
!U  um(:,:) = um(:,:)/float(nmon)
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
  call sepa_spl
  print*, ' .'

  ENDDO  L_LEV

  deallocate( tmp1, tmp2 )
  deallocate( var4d, vtmp, okbi, lzrgw )

  nd1a = NK+1
  nd2a = NTA
  nd3a = NY2
  nd4a = NZ2

  orcs = (/3,4,5,7,14,17/)
  rcs = (/2,3,4,5,7 ,8 /)  ! 3,4,(5,6),(7,8),14,(17,18)

  do irc=1, nrco

    do iv=1, nv
      call setdim
      if ( .not. allocated(set(iv)%var_out) )  &
         allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
      set(iv)%var_out(:,:,:,:) = 1.e32

      set(iv)%var_out(2:,         1:nome+1,:,:) = var6d(:,0:nome,:,:,iv,rcs(irc))
      set(iv)%var_out(2:,nta-nome+2:nta   ,:,:) = var6d(:, :-1  ,:,:,iv,rcs(irc))
    enddo

! DUMP

    write(file_o_rc,'(a,i2.2,a)') trim(file_o)//'_', orcs(irc), '.nc'
    write(title_o_rc,'(a,i2.1)') 'Reconstructed EP flux, ', orcs(irc)

    write(6,*) trim(file_o_rc)

    call outnc(trim(file_o_rc),nv,set,trim(title_o_rc))

  enddo

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
  call getdim(file_i(iv_i),trim(var_i_name(iv_i)))

  call get_iouter(ht       ,lat_rng, iy2)
  call get_iouter(dim4/1.e3,z_rng  , iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  allocate( lat0(ny2), ht0(nz2) )
  lat0(:) = ht(iy2(1):iy2(2))  ;  ht0(:) = dim4(iz2(1):iz2(2))

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
    kr(:,j,k) = kwn(:)/((a_earth+ht0(k))*cos(lat0(j)*deg2rad))
  enddo
  enddo
  or(:) = (-1.)*ome(:)*(2.*pi/(float(nmon+2*nmon_patch)*30.*86400.))
  or2d(:,:) = spread(or,1,nk)
  kr0(:) = kwn(:)/a_earth

  allocate( var6d(nk,-nome+1:nome,ny2,nz2,nv,nrc) )
  var6d(:,:,:,:,:,:) = 0.

  END subroutine initialize

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
  do iv_i=1, nv*2
    var4d(:,0:nome,:,iv_i:iv_i) = get_ivara4d(2,nk,1         ,nome+1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
    var4d(:, :-1  ,:,iv_i:iv_i) = get_ivara4d(2,nk,nta-nome+2,nome-1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
  enddo
!o  var4d(:,:,:,:) = var4d(:,:,:,:)*2.
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
!o  var4d2(:,:,:,:) = var4d2(:,:,:,:)*2.
  var4d2(nk,:   ,:,:) = 0.
  var4d2(: ,nome,:,:) = 0.

  END subroutine get_4var

  SUBROUTINE sepa_spl

  integer ::  k1, k2, o1a, o2a, o1r, o2r, o1m, o2m, o1m0, o2m0

  real, parameter ::  pi = 3.14159265358979323846
  real, parameter ::  beta = 2.*(7.292116e-5)/6371229.

  okbi(:,:) = or2d(:,:)*spread(kr0(:),2,nome*2)/beta

!  lzrgw = 999.e3
!  where ( okbi > -1. )
!    lzrgw(:,:) = 2.*pi*or2d(:,:)**2/(nbv*beta*(1.+okbi(:,:)))
!  end where

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

!o  ! IGW: symm + antisymm
!o  varo(:,k,:,7) = sum(sum(var4d(:,:o1r-1,:,1:8 ), dim=1), dim=1) + &
!o                  sum(sum(var4d(:,:o1r-1,:,9:16), dim=1), dim=1)
!o  varo(:,k,:,8) = sum(sum(var4d(:,o2r+1:,:,1:8 ), dim=1), dim=1) + &
!o                  sum(sum(var4d(:,o2r+1:,:,9:16), dim=1), dim=1)
!o  ! IGW: symm-antisymm cross-correl.
!o  varo(:,k,:,17) = sum(sum(var4d2(:,:o1r-1,:,:), dim=1), dim=1)
!o  varo(:,k,:,18) = sum(sum(var4d2(:,o2r+1:,:,:), dim=1), dim=1)
!o  ! RW: symm + antisymm
!o  varo(:,k,:,1) = sum(sum(var4d(:,o1r:-1,:,1:8 ), dim=1), dim=1) + &
!o                  sum(sum(var4d(:,o1r:-1,:,9:16), dim=1), dim=1)
!o  varo(:,k,:,2) = sum(sum(var4d(:,0 :o2r,:,1:8 ), dim=1), dim=1) + &
!o                  sum(sum(var4d(:,0 :o2r,:,9:16), dim=1), dim=1)
!o  varo(:,k,:,4) = varo(:,k,:,1) + varo(:,k,:,2)
!o  ! RW: symm-antisymm cross-correl.
!o  varo(:,k,:,11) = sum(sum(var4d2(:,o1r:-1,:,:), dim=1), dim=1)
!o  varo(:,k,:,12) = sum(sum(var4d2(:,0 :o2r,:,:), dim=1), dim=1)
!o  varo(:,k,:,14) = varo(:,k,:,11) + varo(:,k,:,12)

  ! IGW: symm + antisymm ( 7,8 --> 5 )
  var6d(:,:o1r-1,:,k,:,5) = var4d(:,:o1r-1,:,1:8 ) + &
                            var4d(:,:o1r-1,:,9:16)
  var6d(:,o2r+1:,:,k,:,5) = var4d(:,o2r+1:,:,1:8 ) + &
                            var4d(:,o2r+1:,:,9:16)
  ! IGW: symm-antisymm cross-correl. ( 17,18 --> 8 )
  var6d(:,:o1r-1,:,k,:,8) = var4d2(:,:o1r-1,:,:)
  var6d(:,o2r+1:,:,k,:,8) = var4d2(:,o2r+1:,:,:)
  ! RW: symm + antisymm ( 4 --> 3 )
  var6d(:,o1r:o2r,:,k,:,3) = var4d(:,o1r:o2r,:,1:8 ) + &
                             var4d(:,o1r:o2r,:,9:16)
  ! RW: symm-antisymm cross-correl. ( 14 --> 7 )
  var6d(:,o1r:o2r,:,k,:,7) = var4d2(:,o1r:o2r,:,:)


! total large-scale waves

!o  varo(:,k,:,1) = varo(:,k,:,1) + varo(:,k,:,7)
!o  varo(:,k,:,2) = varo(:,k,:,2) + varo(:,k,:,8)
!o  varo(:,k,:,11) = varo(:,k,:,11) + varo(:,k,:,17)
!o  varo(:,k,:,12) = varo(:,k,:,12) + varo(:,k,:,18)

  ! ( 1,2 --> 1  /  11,12 --> 6 )
  var6d(:,:,:,k,:,1) = var4d(:,:,:,1:8) + var4d(:,:,:,9:16)
  var6d(:,:,:,k,:,6) = var4d2(:,:,:,:)


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
!o  varo(:,k,:,3) = tmp1(:,:) + tmp2(:,:)

  ! ( 3 --> 2 )
  var6d(:,:-1,:,k,:,2) = var4d(:,:-1,:,1:8)

  ! RW(4), IGW(7) step 2
!o  varo(:,k,:,4) = varo(:,k,:,4) - tmp1(:,:)
!o  varo(:,k,:,7) = varo(:,k,:,7) - tmp2(:,:)

  ! ( 4 --> 3  /  7,8 --> 5 )
  var6d(:,o1r:-1,:,k,:,3) = var6d(:,o1r:-1,:,k,:,3) - var4d(:,o1r:-1,:,1:8)
  var6d(:,:o1r-1,:,k,:,5) = var6d(:,:o1r-1,:,k,:,5) - var4d(:,:o1r-1,:,1:8)


! MRGW(5,6) ( 2 <= period <= 10 days )
  o1m  = -(nmon+2*nmon_patch)*15  ;  o2m  = (nmon+2*nmon_patch)*15
  o1m0 = -(nmon+2*nmon_patch)*3   ;  o2m0 = (nmon+2*nmon_patch)*3

  var4d(:,o1m0+1:o2m0-1,:,9:16) = 0.
  do n=o1m, o2m
    if ( n > o1m0 .and. n < o2m0 )  CYCLE
  do i=1, nk
!    if (lzrgw(i,n) <= lz0_rgw) then
!      var4d(i,n,:,9:16) = 0.
!    else
      ! F_uw/F_vT > 0 : Filtering from MRGW
      do j=ny2/2+1, ny2-3
        if ( vtmp(i,n,j,4)/vtmp(i,n,j,3) > 0. .or. &
             vtmp(i,n,ny2+1-j,4)/vtmp(i,n,ny2+1-j,3) > 0. ) then
          var4d(i,n,j:,9:16) = 0.
          var4d(i,n,:ny2+1-j,9:16) = 0.
          EXIT
        end if
      enddo
!    end if
  enddo
  enddo
  var4d(:,o1m:o2m,1:3      ,9:16) = 0.
  var4d(:,o1m:o2m,ny2-2:ny2,9:16) = 0.

!o  varo(:,k,:,5) = sum(sum(var4d(:,o1r:-1,:,9:16), dim=1), dim=1)
!o  varo(:,k,:,6) = sum(sum(var4d(:,0:o2r ,:,9:16), dim=1), dim=1)

  ! ( 5,6 --> 4 )
  var6d(:,o1r:o2r,:,k,:,4) = var4d(:,o1r:o2r,:,9:16)
  ! o1r > o1m ; o2r < o2m

  ! RW(4) step 3
!o  varo(:,k,:,4) = varo(:,k,:,4) - (varo(:,k,:,5) + varo(:,k,:,6))

  ! ( 4 --> 3 )
  var6d(:,o1r:o2r,:,k,:,3) = var6d(:,o1r:o2r,:,k,:,3) - var4d(:,o1r:o2r,:,9:16)

  tmp1(:,:) = sum(sum(var4d(:,o1m:o1r-1,:,9:16), dim=1), dim=1)
  tmp2(:,:) = sum(sum(var4d(:,o2r+1:o2m,:,9:16), dim=1), dim=1)

!o  varo(:,k,:,5) = varo(:,k,:,5) + tmp1(:,:)
!o  varo(:,k,:,6) = varo(:,k,:,6) + tmp2(:,:)

  ! ( 5,6 --> 4 )
  var6d(:,o1m:o1r-1,:,k,:,4) = var6d(:,o1m:o1r-1,:,k,:,4) +  &
                               var4d(:,o1m:o1r-1,:,9:16)
  var6d(:,o2r+1:o2m,:,k,:,4) = var6d(:,o2r+1:o2m,:,k,:,4) +  &
                               var4d(:,o2r+1:o2m,:,9:16)

  ! IGW(7,8) step 3
!o  varo(:,k,:,7) = varo(:,k,:,7) - tmp1(:,:)
!o  varo(:,k,:,8) = varo(:,k,:,8) - tmp2(:,:)

  ! ( 7,8 --> 5 )
  var6d(:,o1m:o1r-1,:,k,:,5) = var6d(:,o1m:o1r-1,:,k,:,5) -  &
                               var4d(:,o1m:o1r-1,:,9:16)
  var6d(:,o2r+1:o2m,:,k,:,5) = var6d(:,o2r+1:o2m,:,k,:,5) -  &
                               var4d(:,o2r+1:o2m,:,9:16)


! temporary : (S+A) + symm-antisymm cross-correl.

!o  ! total large-scale waves
!o  varo(:,k,:,9 ) = varo(:,k,:,1) + varo(:,k,:,11)
!o  varo(:,k,:,10) = varo(:,k,:,2) + varo(:,k,:,12)
!o  ! E+W
!o  varo(:,k,:,19) = varo(:,k,:,9) + varo(:,k,:,10)
!o  ! RW
!o  varo(:,k,:,20) = varo(:,k,:,4) + varo(:,k,:,14)
!o  ! IGW
!o  varo(:,k,:,21) = varo(:,k,:,7) + varo(:,k,:,17)
!o  varo(:,k,:,22) = varo(:,k,:,8) + varo(:,k,:,18)

  ! total large-scale waves ( 1,2 --> 1  /  11,12 --> 6  /  9,10 (19) --> 9 )
  var6d(:,:,:,k,:,9 ) = var6d(:,:,:,k,:,1) + var6d(:,:,:,k,:,6)
  ! RW ( 4 --> 3  /  14 --> 7  /  20 --> 10 )
  var6d(:,:,:,k,:,10) = var6d(:,:,:,k,:,3) + var6d(:,:,:,k,:,7)
  ! IGW ( 7,8 --> 5  /  17,18 --> 8  /  21,22 --> 11 )
  var6d(:,:,:,k,:,11) = var6d(:,:,:,k,:,5) + var6d(:,:,:,k,:,8)

!yh+ test
!varo(:,k,:,:) = sum(sum(var6d(:,:,:,k,:,:), dim=1), dim=1) - varo(:,k,:,:)
!yh-

  do iv=1, nv
    if (var4d(3,3,1    ,iv) == 2.e32) then
      var6d(:,:,1,k,iv,:) = 1.e32
    end if
    if (var4d(3,3,ny2  ,iv) == 2.e32) then
      var6d(:,:,ny2,k,iv,:) = 1.e32
    end if
    if (var4d(3,3,2    ,iv) == 2.e32) then
      var6d(:,:,2,k,iv,:) = 1.e32
    end if
    if (var4d(3,3,ny2-1,iv) == 2.e32) then
      var6d(:,:,ny2-1,k,iv,:) = 1.e32
    end if
  end do
  if ( k <= 2 .or. k >= nz2-1 ) then
    do iv=1, nv
      if (var4d(3,3,3,iv) == 2.e32) then
        var6d(:,:,:,k,iv,:) = 1.e32
      end if
    enddo
  end if

  END subroutine sepa_spl

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'k_wn  ','ome_fr','lat   ','z     '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  if ( .not. allocated(set(iv)%axis1) ) then
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
  end if
  set(iv)%axis1(:) = lon(1:nk+1)
  set(iv)%axis2(:) = lat(1:nta)
  set(iv)%axis3 = lat0
  set(iv)%axis4 = ht0

  END subroutine setdim

  SUBROUTINE finalize

!U  deallocate( um )
  deallocate( var6d )
  deallocate( lon, lat, ht, lat0, ht0 )
  deallocate( kwn, ome, kr, or, or2d, kr0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program RC_EPFSA_UM_z


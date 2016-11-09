PROGRAM RC_VARSA_UM_z

  use reanal
  use netio

  implicit none

  integer, parameter ::  nw = 6, nvi = 4
  integer, parameter ::  nv = nvi*nw*2, nve = 8, nrc = 22

  integer ::  opt_30d, k_largest, period_smallest, nmon_patch

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ K_LARGEST, PERIOD_SMALLEST, LAT_RNG, P_RNG, NMON_PATCH
  namelist /FILEIO/ NT_F4, P_PREDEF, FILE_I_HEAD, FILE_I_FORM,           &
                    FILE_I_XXXX, VAR_I_NAME, FILE_I_HEAD2, FILE_I_FORM2, &
                    FILE_I_XXXX2, VAR_I2, VAR_I_NAME2, FILE_O, NL_AUX

  integer ::  iz, imon, ihour, i_time
  integer ::  nk, nome
  integer ::  iy1(2), iz1(2), iy2(2), iz2(2), ny2, nz2, nta
  integer ::  i,j,k,n, i_sa, i_ri, ii, i_wv, ivi
  real    ::  p_predef1, wgt
  character(len=1),  dimension(2)   ::  sa, ri
  character(len=5),  dimension(nw)  ::  wv
  character(len=32), dimension(nvi) ::  var_i_tmp
  character(len=32), dimension(nv)  ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  var5d, var5dw, varo
  real, dimension(:,:,:,:), allocatable ::  var4d, var4d2, vtmp
  real, dimension(:,:,:),   allocatable ::  kr
  real, dimension(:,:),     allocatable ::  um, okbi, lzrgw, or2d
  real, dimension(:,:),     allocatable ::  tmp1, tmp2, tmp3, tmp4
  real, dimension(:),       allocatable ::  lat0, p0a, kwn, ome, or, kr0
  real, dimension(:),       allocatable ::  kwn0, ome0

  type(vset), dimension(nv) ::  set

  real, parameter ::  nbv = 2.5e-2
  real, parameter ::  lz0_rgw = 2.4e3

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  call initialize

  iv = 0
  do ivi=1, nvi
  do i_wv=1, 6
  do i_ri=1, 2
    iv = iv + 1
    ovarname(iv) = 'fc'//ri(i_ri)//'_'//trim(var_i_name2(ivi))//'_'//trim(wv(i_wv))
  enddo
  enddo
  enddo

!U  ! get mean U
!U  call switch_para_in
!U  allocate( um(ny2,nz2) )
!U  um(:,:) = 0.
!U
!U  deallocate( lon, lat, p, dim4 )
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

  allocate( var5d    (nk,-nome+1:nome,ny,nvi*2,2)  )
  allocate( var5dw   (nk,-nome+1:nome,ny,nvi*2,nw) )
  allocate( var4d    (nk,-nome+1:nome,ny2,nve*2)   )
  allocate( var4d2   (nk,-nome+1:nome,ny2,nve  )   )
  allocate( vtmp     (nk,-nome+1:nome,ny2,4) )
  allocate( okbi     (nk,-nome+1:nome) )
  allocate( lzrgw    (nk,-nome+1:nome) )
  allocate( tmp1(ny2,nve*2), tmp2(ny2,nve*2) )

  L_LEV:  DO k=1, nz2

  call get_8var

  print*, ' Calculating...', k, '/', nz2
  call reconstr
  print*, ' .'

  iv = 0
  do ivi=1, nvi
  do i_wv=1, 6
  do i_ri=1, 2
    iv = iv + 1
    ii = i_ri + (ivi-1)*2
    varo(1:nk,         1:nome+1,:,k,iv) = var5dw(:,0:nome,:,ii,i_wv)
    varo(1:nk,nta-nome+2:nta   ,:,k,iv) = var5dw(:, :-1  ,:,ii,i_wv)
  enddo
  enddo
  enddo

  ENDDO  L_LEV

  deallocate( tmp1, tmp2 )
  deallocate( var5d, var4d, vtmp, okbi, lzrgw )

  nd1a = NK + 1
  nd2a = NTA
  nd3a = NY
  nd4a = NZ2

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = varo(:,:,:,:,iv)
  enddo

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'Equatorial wave modes')

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

  iv_i = nve*2
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
  allocate( kwn0(0:nk), ome0(nta) )
  kwn0(:) = lon(1:nk+1)
  ome0(:) = lat(:)
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

  sa = (/'s','a'/)
  ri = (/'r','i'/)
  wv = (/'k    ','mrg','r_s','r_a','ig_s','ig_a'/)

  ny = (ny2+1)/2

  allocate( varo(0:nk,nta,ny,nz2,nv) )
  varo = 0.

  END subroutine initialize

  SUBROUTINE get_8var

  ex0 = .TRUE.
  do iv_i=1, nve*4
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
!    print*, trim(file_i(iv_i))
  enddo
  if (.not. ex0)  STOP

  ! read 16 var.s
  do iv_i=1, nve*2
    var4d(:,0:nome,:,iv_i:iv_i) = get_ivara4d(2,nk,1         ,nome+1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
    var4d(:, :-1  ,:,iv_i:iv_i) = get_ivara4d(2,nk,nta-nome+2,nome-1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
  enddo
  var4d(:,:,:,:) = var4d(:,:,:,:)*2.
  var4d(nk,:   ,:,:) = 0.
  var4d(: ,nome,:,:) = 0.

  ! read 16 var.s
  do iv_i=nve*2+1, nve*3
    iv = iv_i - nve*2
    var4d2(:,0:nome,:,iv:iv) = get_ivara4d(2,nk,1         ,nome+1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
    var4d2(:, :-1  ,:,iv:iv) = get_ivara4d(2,nk,nta-nome+2,nome-1, &
                   iy2(1),ny2,iz2(1)-1+k,1)
  enddo
  do iv_i=nve*3+1, nve*4
    iv = iv_i - nve*3
    var4d2(:,0:nome,:,iv:iv) = var4d2(:,0:nome,:,iv:iv) +                &
       get_ivara4d(2,nk,1         ,nome+1,iy2(1),ny2,iz2(1)-1+k,1)
    var4d2(:, :-1  ,:,iv:iv) = var4d2(:, :-1  ,:,iv:iv) +                &
       get_ivara4d(2,nk,nta-nome+2,nome-1,iy2(1),ny2,iz2(1)-1+k,1)
  enddo
  var4d2(:,:,:,:) = var4d2(:,:,:,:)*2.
  var4d2(nk,:   ,:,:) = 0.
  var4d2(: ,nome,:,:) = 0.

  call switch_para_in

  ex0 = .TRUE.
  do iv_i=1, nvi
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
!    print*, trim(file_i(iv_i))
  enddo
  if (.not. ex0)  STOP

  ! read 16 var.s
  var_i_tmp(:) = var_i_name(1:nvi)
  do iv_i=1, nvi
  do i_sa=1, 2
  do i_ri=1, 2
    ii = i_ri + (iv_i-1)*2
    var_i_name(iv_i) = 'fc'//ri(i_ri)//'_'//trim(var_i_tmp(iv_i))//'_'//sa(i_sa)
    if (iv_i == 3)  var_i_name(iv_i) =                                 &
                     'fc'//ri(i_ri)//'_'//trim(var_i_tmp(iv_i))//'_'//sa(3-i_sa)
    var5d(:,0:nome,:,ii:ii,i_sa) = get_ivara4d(2,nk,1         ,nome+1, &
                   1,ny,iz2(1)-1+k,1)
    var5d(:, :-1  ,:,ii:ii,i_sa) = get_ivara4d(2,nk,nta-nome+2,nome-1, &
                   1,ny,iz2(1)-1+k,1)
  enddo
  enddo
  enddo
  var_i_name(1:nvi) = var_i_tmp(:)
  var5d(nk,:   ,:,:,:) = 0.
  var5d(: ,nome,:,:,:) = 0.

  call switch_para_in

  END subroutine get_8var

  SUBROUTINE reconstr

  integer ::  k1, k2, o1a, o2a, o1r, o2r, o1m, o2m, o1m0, o2m0

  real, parameter ::  pi = 3.14159265358979323846
  real, parameter ::  beta = 2.*(7.292116e-5)/6371229.

  okbi(:,:) = or2d(:,:)*spread(kr0(:),2,nome*2)/beta

!  lzrgw = 999.e3
!  where ( okbi > -1. )
!    lzrgw(:,:) = 2.*pi*or2d(:,:)**2/(nbv*beta*(1.+okbi(:,:)))
!  end where

! low-pass filter for large-scale waves ( k = 1-20, period > 4/3 days )
! This should be same with reconstr_epf.f90
  k1 = 1  ;  k2 = 20
  o1a = -(nmon+2*nmon_patch)*22 + 1  ;  o2a = (nmon+2*nmon_patch)*22 - 1
  var4d(k1:k2,:o1a-1,:,:) = 0.  ;  var4d2(k1:k2,:o1a-1,:,:) = 0.
  var4d(k1:k2,o2a+1:,:,:) = 0.  ;  var4d2(k1:k2,o2a+1:,:,:) = 0.
  var5d(k1:k2,:o1a-1,:,:,:) = 0.
  var5d(k1:k2,o2a+1:,:,:,:) = 0.
  if ( k2 < nk ) then
    var4d(k2+1:,:,:,:) = 0.  ;  var4d2(k2+1:,:,:,:) = 0.
    var5d(k2+1:,:,:,:,:) = 0.
  end if

  var5dw(:,:,:,:,1) = var5d(:,:,:,:,1)  ! K
  var5dw(:,:,:,:,2) = var5d(:,:,:,:,2)  ! MRG
  var5dw(:,:,:,:,3) = var5d(:,:,:,:,1)  ! R_S
  var5dw(:,:,:,:,4) = var5d(:,:,:,:,2)  ! R_A
  var5dw(:,:,:,:,5) = var5d(:,:,:,:,1)  ! IG_S
  var5dw(:,:,:,:,6) = var5d(:,:,:,:,2)  ! IG_A

! low-freq. IGW(7,8), RW(4) ( period <, >= 2.5 days ) step 1
  o1r = -(nmon+2*nmon_patch)*12  ;  o2r = (nmon+2*nmon_patch)*12

  ! IGW: symm + antisymm
!!  varo(:,k,:,7) = sum(sum(var4d(:,:o1r-1,:,1:8 ), dim=1), dim=1) + &
!!                  sum(sum(var4d(:,:o1r-1,:,9:16), dim=1), dim=1)
!!  varo(:,k,:,8) = sum(sum(var4d(:,o2r+1:,:,1:8 ), dim=1), dim=1) + &
!!                  sum(sum(var4d(:,o2r+1:,:,9:16), dim=1), dim=1)
  var5dw(:,o1r:o2r,:,:,5:6) = 0.
  ! RW: symm + antisymm
!!  varo(:,k,:,1) = sum(sum(var4d(:,o1r:-1,:,1:8 ), dim=1), dim=1) + &
!!                  sum(sum(var4d(:,o1r:-1,:,9:16), dim=1), dim=1)
!!  varo(:,k,:,2) = sum(sum(var4d(:,0 :o2r,:,1:8 ), dim=1), dim=1) + &
!!                  sum(sum(var4d(:,0 :o2r,:,9:16), dim=1), dim=1)
!!  varo(:,k,:,4) = varo(:,k,:,1) + varo(:,k,:,2)
  var5dw(:,:o1r-1,:,:,3:4) = 0.
  var5dw(:,o2r+1:,:,:,3:4) = 0.

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
!!        var4d(i,n,j:,1:8) = 0.
!!        var4d(i,n,:ny2+1-j,1:8) = 0.
        var5dw(i,n,j-ny2/2:,:,1) = 0.
        EXIT
      end if
    enddo
  enddo
  enddo
!!  var4d(:,-nome+1:-1,1:3      ,1:8) = 0.
!!  var4d(:,-nome+1:-1,ny2-2:ny2,1:8) = 0.
  var5dw(:,-nome+1:-1,ny-2:ny,:,1) = 0.

!!  tmp1(:,:) = sum(sum(var4d(:,o1r:-1,:,1:8), dim=1), dim=1)
!!  tmp2(:,:) = sum(sum(var4d(:,:o1r-1,:,1:8), dim=1), dim=1)
!!  varo(:,k,:,3) = tmp1(:,:) + tmp2(:,:)
  var5dw(:,0:,:,:,1) = 0.

  ! RW(4), IGW(7) step 2
!!  varo(:,k,:,4) = varo(:,k,:,4) - tmp1(:,:)
!!  varo(:,k,:,7) = varo(:,k,:,7) - tmp2(:,:)
if (nvi /= 4) then
  print*, 'modify the code'  ;  STOP
end if
  where ( var5dw(:,:,:,1,1) /= 0. )
    var5dw(:,:,:,1,3) = 0.  ;  var5dw(:,:,:,2,3) = 0.
    var5dw(:,:,:,3,3) = 0.  ;  var5dw(:,:,:,4,3) = 0.
    var5dw(:,:,:,5,3) = 0.  ;  var5dw(:,:,:,6,3) = 0.
    var5dw(:,:,:,7,3) = 0.  ;  var5dw(:,:,:,8,3) = 0.
    var5dw(:,:,:,1,5) = 0.  ;  var5dw(:,:,:,2,5) = 0.
    var5dw(:,:,:,3,5) = 0.  ;  var5dw(:,:,:,4,5) = 0.
    var5dw(:,:,:,5,5) = 0.  ;  var5dw(:,:,:,6,5) = 0.
    var5dw(:,:,:,7,5) = 0.  ;  var5dw(:,:,:,8,5) = 0.
  end where

! MRGW(5,6) ( 2 <= period <= 10 days )
  o1m  = -(nmon+2*nmon_patch)*15  ;  o2m  = (nmon+2*nmon_patch)*15
  o1m0 = -(nmon+2*nmon_patch)*3   ;  o2m0 = (nmon+2*nmon_patch)*3

!!  var4d(:,o1m0+1:o2m0-1,:,9:16) = 0.
  var5dw(:,o1m0+1:o2m0-1,:,:,2) = 0.
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
!!          var4d(i,n,j:,9:16) = 0.
!!          var4d(i,n,:ny2+1-j,9:16) = 0.
          var5dw(i,n,j-ny2/2:,:,2) = 0.
          EXIT
        end if
      enddo
!    end if
  enddo
  enddo
!!  var4d(:,o1m:o2m,1:3      ,9:16) = 0.
!!  var4d(:,o1m:o2m,ny2-2:ny2,9:16) = 0.
  var5dw(:,o1m:o2m,ny-2:ny,:,2) = 0.

!!  varo(:,k,:,5) = sum(sum(var4d(:,o1r:-1,:,9:16), dim=1), dim=1)
!!  varo(:,k,:,6) = sum(sum(var4d(:,0:o2r ,:,9:16), dim=1), dim=1)
!!  ! o1r > o1m ; o2r < o2m

  ! RW(4) step 3
!!  varo(:,k,:,4) = varo(:,k,:,4) - (varo(:,k,:,5) + varo(:,k,:,6))

!!  tmp1(:,:) = sum(sum(var4d(:,o1m:o1r-1,:,9:16), dim=1), dim=1)
!!  tmp2(:,:) = sum(sum(var4d(:,o2r+1:o2m,:,9:16), dim=1), dim=1)

!!  varo(:,k,:,5) = varo(:,k,:,5) + tmp1(:,:)
!!  varo(:,k,:,6) = varo(:,k,:,6) + tmp2(:,:)
  var5dw(:,:o1m-1,:,:,2) = 0.
  var5dw(:,o2m+1:,:,:,2) = 0.
 
  ! IGW(7,8) step 3
!!  varo(:,k,:,7) = varo(:,k,:,7) - tmp1(:,:)
!!  varo(:,k,:,8) = varo(:,k,:,8) - tmp2(:,:)

if (nvi /= 4) then
  print*, 'modify the code'  ;  STOP
end if
  where ( var5dw(:,:,:,1,2) /= 0. )
    var5dw(:,:,:,1,4) = 0.  ;  var5dw(:,:,:,2,4) = 0.
    var5dw(:,:,:,3,4) = 0.  ;  var5dw(:,:,:,4,4) = 0.
    var5dw(:,:,:,5,4) = 0.  ;  var5dw(:,:,:,6,4) = 0.
    var5dw(:,:,:,7,4) = 0.  ;  var5dw(:,:,:,8,4) = 0.
    var5dw(:,:,:,1,6) = 0.  ;  var5dw(:,:,:,2,6) = 0.
    var5dw(:,:,:,3,6) = 0.  ;  var5dw(:,:,:,4,6) = 0.
    var5dw(:,:,:,5,6) = 0.  ;  var5dw(:,:,:,6,6) = 0.
    var5dw(:,:,:,7,6) = 0.  ;  var5dw(:,:,:,8,6) = 0.
  end where

  var5dw(:,:,ny-2:ny,:,:) = 0.

  END subroutine reconstr

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'k_wn   ','ome_fr','lat','p'/)
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = kwn0
  set(iv)%axis2 = ome0
  set(iv)%axis3 = lat0(ny2/2+1:ny2)
  set(iv)%axis4 = p0a

  END subroutine setdim

  SUBROUTINE finalize

!U  deallocate( um )
  deallocate( varo )
  deallocate( lon, lat, p, lat0, p0a )
  deallocate( kwn, ome, kr, or, or2d, kr0 )
  deallocate( kwn0, ome0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program RC_VARSA_UM_z


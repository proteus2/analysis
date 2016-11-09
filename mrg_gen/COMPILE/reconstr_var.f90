PROGRAM RC_VAR_UM_z

  use fft
  use hadgem
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 12

  integer ::  k_largest, period_smallest, nmon_patch

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D, NRES
  namelist /PARAM/ VAR_NAME, K_LARGEST, PERIOD_SMALLEST, LAT_RNG, Z_RNG, &
                   NMON_PATCH
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, FILE_I_HEAD2, FILE_I_FORM2,             &
                    FILE_I_XXXX2, VAR_I_NAME2, FILE_O

  integer ::  imon, idir, l
  integer ::  iy1(2), iz1(2), iy2(2), iz2(2), nx2, ny2, nz2, nta,        &
              nt0, tmpi, nk, nome
  real    ::  dt, inv_pcfsum(0:1)
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  varo
  real, dimension(:,:,:,:),   allocatable ::  var4d, rcs, arcs
  real, dimension(:,:,:),     allocatable ::  kr, o_hat
  real, dimension(:,:),       allocatable ::  pcf, um
  real, dimension(:),         allocatable ::  lat0, ht0, lon0, time0,    &
                                              lat0s, tmp1d, tmp1d2,      &
                                              kwn, ome, or
  complex, dimension(:,:,:,:), allocatable ::  fc4d
  complex, dimension(:,:,:),   allocatable ::  fc3d

  type(vset), dimension(nv) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  ovarname(1 ) = trim(var_name(1))//'_e_s  '
  ovarname(2 ) = trim(var_name(1))//'_w_s  '
  ovarname(3 ) = trim(var_name(1))//'_e_a  '
  ovarname(4 ) = trim(var_name(1))//'_w_a  '
  ovarname(5 ) = trim(var_name(1))//'_e_s0 '
  ovarname(6 ) = trim(var_name(1))//'_w_s0 '
  ovarname(7 ) = trim(var_name(1))//'_e_a1 '
  ovarname(8 ) = trim(var_name(1))//'_w_a1 '
  ovarname(9 ) = trim(var_name(1))//'_e_s0x'
  ovarname(10) = trim(var_name(1))//'_w_s0x'
  ovarname(11) = trim(var_name(1))//'_e_a1x'
  ovarname(12) = trim(var_name(1))//'_w_a1x'

! GET AXES AND INITIALIZE ARRAYS

  call initialize

  ! get mean U
  call switch_para_in
  allocate( um(ny2,nz2) )
  um(:,:) = 0.

  deallocate( lon, lat, ht, ht_th, dim4 )
  iv_i = 1
  file_i(iv_i) = get_ifilename()
  call getdim(file_i(iv_i),var_i_name(iv_i))
  call get_iouter(lat    ,lat_rng, iy1)
  call get_iouter(ht/1.e3,z_rng  , iz1)
  if ( any(lat(iy1(1):iy1(2)) /= lat0) .or. &
       any(ht (iz1(1):iz1(2)) /= ht0 ) ) then
    print*, 'Check lat0, ht0'  ;  STOP
  end if

  do imon=1, nmon
    iv_i = 1
    file_i(iv_i) = get_ifilename()  ;  print*, trim(file_i(1))
    um(:,:) = um(:,:) + reshape(                                         &
       sum(get_ivara4d(1,nx,iy1(1),ny2,iz1(1),nz2,1,1),dim=1)/float(nx), &
       (/ny2,nz2/) )
    mon = mon + 1
    if (mon == 13) then
      year = year + 1  ;  mon = 1
    end if
  enddo
  um(:,:) = um(:,:)/float(nmon)

  ! rewind
  year = yyyy
  mon  = mm(1)

  call switch_para_in
  nx = nx2

  allocate( rcs(nx,ny2,nt,2), arcs(nx,nt,2,0:1) )
  allocate( o_hat(nk+1,nome*2,ny2) )

  L_LEV:  DO k=1, nz2

  do j=1, ny2
  do n=1, nome*2
    o_hat(:,n,j) = or(n) - um(j,k)*kr(:,j,k)
  enddo
  enddo

  call get_1var

  do j=1, ny2/2+1
    varo(:,j,k,:,1:2) = 0.5*(rcs(:,ny2/2+j,:,:)+rcs(:,ny2/2+2-j,:,:))
    varo(:,j,k,:,3:4) = 0.5*(rcs(:,ny2/2+j,:,:)-rcs(:,ny2/2+2-j,:,:))
  enddo

  do l=0, 1
  do idir=1, 2
  do n=1, nt
  do i=1, nx
    arcs(i,n,idir,l) = sum(rcs(i,:,n,idir)*pcf(:,l))*inv_pcfsum(l)
  enddo
  enddo
  enddo
  enddo
  do j=1, ny2/2+1
    varo(:,j,k,:,5:6) = arcs(:,:,:,0)*pcf(ny2/2+j,0)
    varo(:,j,k,:,7:8) = arcs(:,:,:,1)*pcf(ny2/2+j,1)
  enddo

  ENDDO  L_LEV

  varo(:,:,:,:,9:12) = varo(:,:,:,:,1:4) - varo(:,:,:,:,5:8)

  deallocate( rcs, arcs, o_hat )

  nd1a = NX
  nd2a = NY2/2+1
  nd3a = NZ2
  nd4a = NT

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = varo(:,:,:,:,iv)
  enddo

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'Reconstructed variable (large and '//  &
             'long scale)')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  real, parameter ::  a_earth = 6371229.
  real, parameter ::  pi = 3.14159265358979323846
!  real, parameter ::  ome_earth = 7.292116e-5
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
  call getdim(file_i(iv_i),'fcr_'//trim(var_name(iv_i)))

  call get_iouter(ht       ,lat_rng, iy2)
  call get_iouter(dim4/1.e3,z_rng  , iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  allocate( lat0(ny2), ht0(nz2) )
  lat0(:) = ht(iy2(1):iy2(2))  ;  ht0(:) = dim4(iz2(1):iz2(2))

  nta = ny
  if ( nta /= (nmon+2*nmon_patch)*30*nhour .or. opt_30d == 0 ) then
    print*, 'not yet programmed (for opt_30d = 0)'  ;  STOP
  end if

  nk = nx - 1
  if (k_largest /= -999)  nk = k_largest
  nx = nk*2

  tmpi = 30*nhour
  if (period_smallest /= -999)  tmpi = 30/period_smallest*2
  nt = nmon*tmpi
  nt0 = nmon_patch*tmpi
  dt = 1./(float(tmpi)/30.)
  nome = nta*tmpi/(30*nhour)/2

  allocate( kwn(nk+1), ome(nome*2) )
  kwn(:) = lon(1:nk+1)
  ome(1:nome+1) = lat(1:nome+1)
  ome(nome+2:nome*2) = ome(nome:2:-1)*(-1.)

  allocate( kr(nk+1,ny2,nz2), or(nome*2) )
  do k=1, nz2
  do j=1, ny2
    kr(:,j,k) = kwn(:)/((a_earth+ht0(k))*cos(lat0(j)*deg2rad))
  enddo 
  enddo
  or(:) = (-1.)*ome(:)*(2.*pi/(float(nmon+2*nmon_patch)*30.*86400.))

  allocate( lon0(nx), time0(nt) )

  do i=1, nx
    lon0(i) = float(i-1)*(360./float(nx))
  enddo

  time0(1) = (year-refdate(1))*360.+(mon-refdate(2))*30.+(date-refdate(3))+hour/24.
  do n=2, nt
    time0(n) = time0(1) + float(n-1)*dt
  enddo

  allocate( lat0s(ny2), pcf(ny2,0:1) )
  lat0s(:) = lat0(:)/6.
  allocate( tmp1d(ny2), tmp1d2(ny2) )
  tmp1d (:) = lat0s(:)*lat0s(:)
  tmp1d2(:) = exp(-0.25*tmp1d(:))
  pcf(:,0) = 1.
  pcf(:,1) = lat0s(:)
!  pcf(:,2) = tmp1d(:) - 1.
!  pcf(:,3) = lat0s(:)*(tmp1d(:) - 3.)
!  pcf(:,4) = tmp1d(:)*(tmp1d(:) - 6.) + 3.
!  pcf(:,5) = lat0s(:)*(tmp1d(:)*(tmp1d(:) - 10.) + 15.)
  pcf(:,:) = pcf(:,:)*spread(tmp1d2(:),2,1+1)
!  inv_pcfsum = (/1.,1.,2.,6.,24.,120./)
!  inv_pcfsum(:) = (lat0s(2)-lat0s(1))/(inv_pcfsum(:)*sqrt(2.*pi))
  inv_pcfsum(:) = 1./sum(pcf(:,:)*pcf(:,:), dim=1)
  deallocate( tmp1d, tmp1d2 )

  allocate( varo(nx,ny2/2+1,nz2,nt,nv) )

  nx2 = nx

  END subroutine initialize

  SUBROUTINE get_1var

  allocate( fc3d(nx,ny2,nome*2), fc4d(nk+1,nome*2,ny2,2) )

  iv_i = 1

  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  print*, trim(file_i(iv_i))

  allocate( var4d(nk+1,nome*2,ny2,2) )
  var_i_name(iv_i) = 'fcr_'//trim(var_name(iv_i))
  var4d(:,:nome,:,1:1) = get_ivara4d(1,nk+1,1,nome, &
                           iy2(1),ny2,iz2(1)-1+k,1)
  var4d(:,nome+2:,:,1:1) = get_ivara4d(1,nk+1,nta-nome+2,nome-1, &
                           iy2(1),ny2,iz2(1)-1+k,1)
  var_i_name(iv_i) = 'fci_'//trim(var_name(iv_i))
  var4d(:,:nome,:,2:2) = get_ivara4d(1,nk+1,1,nome, &
                           iy2(1),ny2,iz2(1)-1+k,1)
  var4d(:,nome+2:,:,2:2) = get_ivara4d(1,nk+1,nta-nome+2,nome-1, &
                           iy2(1),ny2,iz2(1)-1+k,1)
  if (nk == nres)  var4d(nk+1,:,:,:) = 0.
  var4d(:,nome+1,:,:) = 0.
  fc4d(:,:,:,1) = cmplx(var4d(:,:,:,1),var4d(:,:,:,2))
  fc4d(2:,:,:,1) = fc4d(2:,:,:,1)*2.
  fc4d(:,:,:,2) = fc4d(:,:,:,1)
  deallocate( var4d )

  print*, ' Calculating...', k, '/', nz2

  where (o_hat > 0.)
    fc4d(:,:,:,2) = 0.
  else where
    fc4d(:,:,:,1) = 0.
  end where

  fc3d(:,:,:) = 0.
  do idir=1, 2
  do j=1, ny2
    do i=2, nk+1
      call fft1d_b(nome*2,fc4d(i,:,j,idir),fc3d(i,j,:))
    enddo
    do n=1, nt
      call fft1d_b(nx,fc3d(:,j,n+nt0),rcs(:,j,n,idir))
    enddo
  enddo
  enddo
  print*, ' .'

  deallocate( fc3d, fc4d )

  END subroutine get_1var

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lon','lat','z','t'/)
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lon0
  set(iv)%axis2 = lat0(ny2/2+1:ny2)
  set(iv)%axis3 = ht0
  set(iv)%axis4 = time0

  END subroutine setdim

  SUBROUTINE finalize

  deallocate( varo, um )
  deallocate( lon, lat, ht, lat0, ht0, lon0, time0, lat0s, pcf )
  deallocate( kwn, ome, kr, or )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program RC_VAR_UM_z


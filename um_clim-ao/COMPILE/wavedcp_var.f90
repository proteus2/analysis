PROGRAM RC_VAR_UM_z

  use hadgem
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 8

  integer ::  k_largest, period_smallest, nmon_patch

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D, NRES
  namelist /PARAM/ VAR_NAME, K_LARGEST, PERIOD_SMALLEST, LAT_RNG, Z_RNG, &
                   NMON_PATCH
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FILE_I_HEAD, FILE_I_FORM,       &
                    FILE_I_XXXX, FILE_I_HEAD2, FILE_I_FORM2,             &
                    FILE_I_XXXX2, VAR_I_NAME2, FILE_O

  integer ::  imon, ii, l
  integer ::  iy1(2), iz1(2), iy2(2), iz2(2), ny2, nz2, nta, tmpi, nk, nome
  integer ::  i,j,k,n
  real    ::  inv_pcfsum(0:1)
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  varo
  real, dimension(:,:,:,:),   allocatable ::  var4d, arcs
  real, dimension(:,:,:),     allocatable ::  kr
  real, dimension(:,:),       allocatable ::  pcf
  real, dimension(:),         allocatable ::  lat0, ht0, kwn, ome, ome0, &
                                              or, lat0s, tmp1d, tmp1d2

  type(vset), dimension(nv) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  ovarname(1 ) = 'fcr_'//trim(var_name(1))//'_s  '
  ovarname(2 ) = 'fci_'//trim(var_name(1))//'_s  '
  ovarname(3 ) = 'fcr_'//trim(var_name(1))//'_a  '
  ovarname(4 ) = 'fci_'//trim(var_name(1))//'_a  '
  ovarname(5 ) = 'fcr_'//trim(var_name(1))//'_s0 '
  ovarname(6 ) = 'fci_'//trim(var_name(1))//'_s0 '
  ovarname(7 ) = 'fcr_'//trim(var_name(1))//'_a1 '
  ovarname(8 ) = 'fci_'//trim(var_name(1))//'_a1 '

! GET AXES AND INITIALIZE ARRAYS

  call initialize

  allocate( arcs(nk+1,nome*2,2,0:1) )

  L_LEV:  DO k=1, nz2

  call get_1var

  print*, ' Calculating...', k, '/', nz2

  do j=1, (ny2+1)/2
    varo(:,:,j,k,1:2) = 0.5*(var4d(:,:,ny2/2+j,:)+var4d(:,:,(ny2+3)/2-j,:))
    varo(:,:,j,k,3:4) = 0.5*(var4d(:,:,ny2/2+j,:)-var4d(:,:,(ny2+3)/2-j,:))
  enddo

  do l=0, 1
  do ii=1, 2
  do n=1, nome*2
  do i=1, nk+1
    arcs(i,n,ii,l) = sum(var4d(i,n,:,ii)*pcf(:,l))*inv_pcfsum(l)
  enddo
  enddo
  enddo
  enddo
  do j=1, (ny2+1)/2
    varo(:,:,j,k,5:6) = arcs(:,:,:,0)*pcf(ny2/2+j,0)
    varo(:,:,j,k,7:8) = arcs(:,:,:,1)*pcf(ny2/2+j,1)
  enddo

  ENDDO  L_LEV

  deallocate( arcs )

  nd1a = NK+1
  nd2a = NOME*2
  nd3a = (NY2+1)/2
  nd4a = NZ2

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

  nome = nta/2
  if (period_smallest /= -999)  nome = nta/(period_smallest*nhour)

  allocate( kwn(nk+1), ome(nome*2), ome0(nome*2) )
  kwn(:) = lon(1:nk+1)
  ome(1:nome+1) = lat(1:nome+1)
  ome(nome+2:nome*2) = ome(nome:2:-1)*(-1.)
  ome0(:) = lat(1:nome*2)

  allocate( kr(nk+1,ny2,nz2), or(nome*2) )
  do k=1, nz2
  do j=1, ny2
    kr(:,j,k) = kwn(:)/((a_earth+ht0(k))*cos(lat0(j)*deg2rad))
  enddo 
  enddo
  or(:) = (-1.)*ome(:)*(2.*pi/(float(nmon+2*nmon_patch)*30.*86400.))

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

  allocate( var4d(nk+1,nome*2,ny2,2) )
  allocate( varo(nk+1,nome*2,(ny2+1)/2,nz2,nv) )

  END subroutine initialize

  SUBROUTINE get_1var

  iv_i = 1

  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  print*, trim(file_i(iv_i))

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
  var4d(:,nome+1,:,:) = 0.
  if (nk == nres)  var4d(nk+1,:,:,:) = 0.

  END subroutine get_1var

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'k_wn    ','ome_fr','lat  ','z '/)
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = kwn
  set(iv)%axis2 = ome0
  set(iv)%axis3 = lat0(ny2/2+1:ny2)
  set(iv)%axis4 = ht0

  END subroutine setdim

  SUBROUTINE finalize

  deallocate( varo )
  deallocate( var4d )
  deallocate( lon, lat, ht, lat0, ht0, lat0s, pcf )
  deallocate( kwn, ome, ome0, kr, or )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program RC_VAR_UM_z


PROGRAM EPF_UM_z

  use epf
  use hadgem
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 8, nv1 = 4
  real,    parameter ::  zb_out = 8.e3

  integer ::  k_largest, period_smallest, nmon_patch
  character(len=128) ::  file_o2

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ K_LARGEST, PERIOD_SMALLEST, LAT_RNG, Z_RNG, NMON_PATCH
  namelist /FILEIO/ DAY1, NDAY_I, FID, FILE_I_HEAD, FILE_I_FORM,         &
                    FILE_I_XXXX, VAR_I_NAME, FILE_I_HEAD2, FILE_I_FORM2, &
                    FILE_I_XXXX2, VAR_I_NAME2, FILE_O, FILE_O2

  integer ::  iz, imon, ihour, i_time
  integer ::  nk, nome
  integer ::  iy2(2), iz2(2), ny2, nz2, nta
  integer ::  k
  real    ::  wgt
  character(len=32), dimension(8)  ::  ivarname
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:,:,:), allocatable ::  var6d
  real, dimension(:,:,:,:,:), allocatable ::  var5d
  real, dimension(:,:,:),   allocatable ::  var3d
  real, dimension(:,:,:),   allocatable ::  um, ptm, rhom, tem3d
  real, dimension(:),       allocatable ::  lat0, ht0, ht_th0, kwn, ome

  type(vset), dimension(nv) ::  set

  ovarname(:) = varname_epf(:)

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  call initialize

  i_time = 0

  ! get variable
  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------

  if (opt_30d == 0)  ndate = get_ndate()

  do date=1, ndate
    hour = hh(1)
    do ihour=1, nhour

      i_time = i_time + 1

      day_from_ref = get_dayfromref(year,mon,date,hour)

      call get_3var

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

  year = yyyy
  mon  = mm(1)
  day1 = -999

  call switch_para_in

  day_from_ref = get_dayfromref(year,mon,16,0)

  call get_8var

  allocate( var5d(nk,nome,ny2,nz2,nv) )
  var5d(:,:,:,:,:) = 0.

  print*, ' Calculating...'
  call epf_z_cp_fc2(                                                              &
     nk,nome,ny2,nz2,lat0,ht0,ht_th0,nt,um(:,:,1:nt),ptm(:,:,1:nt),rhom(:,:,1:nt),&
     1.e32, &
     var6d(:,:,:,:,1,1),var6d(:,:,:,:,1,2),var6d(:,:,:,:,1,3),var6d(:,:,:,:,1,4), &
     var6d(:,:,:,:,2,1),var6d(:,:,:,:,2,2),var6d(:,:,:,:,2,3),var6d(:,:,:,:,2,4), &
     var5d(:,:,:,:,1),var5d(:,:,:,:,2),var5d(:,:,:,:,3),var5d(:,:,:,:,4),         &
     var5d(:,:,:,:,5),var5d(:,:,:,:,6),var5d(:,:,:,:,7),var5d(:,:,:,:,8) )

  wgt = float(nmon+2*nmon_patch)/float(nmon+nmon_patch)
  where ( var5d /= 1.e32 )
    var5d = var5d*wgt
  end where

  deallocate( var6d )
  deallocate( um, ptm, rhom )


  nd1a = NK
  nd2a = NOME
  nd3a = NY2
  nd4a = NZ2

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = var5d(:,:,:,:,iv)
  enddo

  if (zb_out > ht0(1)) then
    iz = nd4a + 1
    do k=2, nd4a
      if (ht0(k) > zb_out) then  ;  iz = k - 1  ;  EXIT  ;  end if
    enddo
    do iv=1, nv
      set(iv)%var_out(:,:,:,:iz-1) = 1.e32
    enddo
  end if

  print*, 'Causion:  The dimension size of fc_w_zt (and fc_pt_zt) should' 
  print*, '          be the same as that of fc_u (and fc_v).'

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv1,set(1:nv1),'EP flux analysis in spectral domain')

  write(6,*)  ;  write(6,*) trim(file_o2)  ;  write(6,*)

  call outnc(trim(file_o2),nv-nv1,set(nv1+1:nv),'EP flux analysis in spectral domain')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  year = yyyy
  mon  = mm(1)

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 3  ! rho
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  call get_iouter(lat    ,lat_rng, iy2)
  call get_iouter(ht/1.e3,z_rng  , iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  allocate( lat0(ny2), ht0(nz2), ht_th0(nz2) )
  lat0(:) = lat(iy2(1):iy2(2))  ;  ht0(:) = ht(iz2(1):iz2(2))
  ht_th0(:) = ht_th(iz2(1):iz2(2))

  allocate( um(ny2,nz2,nmon*31*nhour), ptm(ny2,nz2,nmon*31*nhour),       &
            rhom(ny2,nz2,nmon*31*nhour) )

  print*, 'Causion:  The dimension size of fc_w_zt (and fc_pt_zt) should'
  print*, '          be the same as that of fc_u (and fc_v).'

  END subroutine initialize

  SUBROUTINE get_3var

  ex0 = .TRUE.
  do iv_i=1, 3
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  allocate( var3d(nx,ny2,nz2) )

  ! read 3 var.s
  print*, trim(file_i(1))
  iv_i = 1  ;  var3d(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  um(:,:,i_time) = sum(var3d(:,:,:), dim=1)/float(nx)
  print*, trim(file_i(2))
  iv_i = 2  ;  var3d(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  ptm(:,:,i_time) = sum(var3d(:,:,:), dim=1)/float(nx)
  print*, trim(file_i(3))
  iv_i = 3  ;  var3d(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  rhom(:,:,i_time) = sum(var3d(:,:,:), dim=1)/float(nx)
  print*, 'time index :', it_i(1)

  deallocate( var3d )

  END subroutine get_3var

  SUBROUTINE get_8var

  ex0 = .TRUE.
  do iv_i=1, 4
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  deallocate( lon, lat, ht, ht_th, dim4 )
  call getdim(file_i(1),'fcr_'//trim(var_i_name(1)))

  nk = k_largest+1
  if (k_largest == -999)  nk = nx
  allocate( kwn(nk) )
  kwn(:) = lon(:)

  nta = (nmon+2*nmon_patch)*30*nhour
  if ( opt_30d == 0 .and. nmon /= 1 ) then
    print*, 'not yet programmed (for opt_30d = 0)'  ;  STOP
  end if

  nome = nta/(period_smallest*nhour)*2+1
  if (period_smallest == -999)  nome = ny
  allocate( ome(nome) )
  ome(:) = lat(:)

  if ( any(lat0 /= ht) ) then
    print*, 'Check lat0'  ;  STOP
  end if

  allocate( var6d(nk,nome,ny2,nz2,2,4) )
  var6d(:,:,:,:,:,:) = 0.

  ! read 8 var.s
  do iv_i=1, 4
    print*, trim(file_i(iv_i))
    ivarname(iv_i) = trim(var_i_name(iv_i))
    var_i_name(iv_i) = 'fcr_'//trim(ivarname(iv_i))
    var6d(:,:,:,:,1,iv_i) = get_ivara4d(1,nk,1,nome,1,ny2,1,nz2)
    var_i_name(iv_i) = 'fci_'//trim(ivarname(iv_i))
    var6d(:,:,:,:,2,iv_i) = get_ivara4d(1,nk,1,nome,1,ny2,1,nz2)
  enddo

  END subroutine get_8var

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'k_wn    ','ome_fr','lat  ','z '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = kwn
  set(iv)%axis2 = ome
  set(iv)%axis3 = lat0
  set(iv)%axis4 = ht0

  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var5d )
  deallocate( lon, lat, ht, lat0, ht0, ht_th0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program EPF_UM_z


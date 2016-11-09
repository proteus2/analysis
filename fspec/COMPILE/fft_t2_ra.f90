PROGRAM FFT_t2nd

  use fft
  use reanal
  use netio

  implicit none

  integer ::  period_smallest, k_largest, nmon_patch

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE
  namelist /PARAM/ VAR_NAME, PERIOD_SMALLEST, K_LARGEST, LAT_RNG, P_RNG, &
                   NMON_PATCH
  namelist /FILEIO/ NT_F4, MISSV, FILE_I_HEAD, FILE_I_FORM, FILE_I_XXXX, &
                    FILE_O

  integer ::  imon, ihour, i_time, npad
  integer ::  nk, nome, iy2(2), iz2(2), ny2, nz2, nt0, it1, it2
  integer ::  i,j,k
  character(len=32), dimension(2) ::  ovarname

  real,    dimension(:,:,:,:,:), allocatable ::  var5d
  real,    dimension(:,:,:,:),   allocatable ::  var_k
  real,    dimension(:),         allocatable ::  kwn, ome, wgt
  complex, dimension(:,:,:,:),   allocatable ::  fc_var_k
  complex, dimension(:),         allocatable ::  fc_var

  type(vset), dimension(2) ::  set

  real, parameter ::  pi05 = 0.5*3.14159265358979323846

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

  L_MON:  DO imon=1, nmon+1
  !---------------------------------------------------------------------
  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------
  hour = 0

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
  if ( i_time >= nt ) then
    if ( imon == nmon+1 .and. hour >= hh(1) )  EXIT L_MON
  end if

  if ( i_time == 0 .and. hour < hh(1) ) then
    hour = hour + 24/nhour  ;  CYCLE L_HOUR
  end if

  i_time = i_time + 1
print*, year, mon, date, hour

  iv_i = 1
  file_i(1) = get_ifilename()
  file_i(2) = file_i(1)
  inquire(file=trim(file_i(1)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(1)),' not found.'  ;  STOP
  end if

  ! get variable
  call get_var

  t(i_time) = get_dayfromref(year,mon,date,hour)

  hour = hour + 24/nhour
  !---------------------------------------------------------------------
  ENDDO  L_HOUR

  !---------------------------------------------------------------------
  ENDDO  L_DATE

  mon = mon + 1
  if (mon == 13) then
    year = year + 1  ;  mon = 1
  end if
  !---------------------------------------------------------------------
  ENDDO  L_MON

  it1 = 1 + (i_time - nt)/2
  it2 = it1 - 1 + nt

  do k=1, nz2
  do j=1, ny2
  do i=2, nk
    fc_var_k(it1:it2,i,j,k) = wgt(1:nt)*fc_var_k(it1:it2,i,j,k)
  enddo
  enddo
  enddo

! CALCULATE FFT

  do k=1, nz2
  do j=1, ny2
  do i=1, nk
    call fft1d_f(nt,fc_var_k(it1:it2,i,j,k),fc_var)
    var5d(i,1:nome/2,j,k,1) = real (fc_var(1:nome/2))
    var5d(i,1:nome/2,j,k,2) = aimag(fc_var(1:nome/2))
    var5d(i,nome/2+1:nome,j,k,1) = real (fc_var(nt+1-nome/2:nt))
    var5d(i,nome/2+1:nome,j,k,2) = aimag(fc_var(nt+1-nome/2:nt))
  enddo
  enddo
  enddo
  if (nome /= nt)  var5d(:,nome/2+1,:,:,:) = 0.

  nd1a = NK
  nd2a = NOME
  nd3a = NY2
  nd4a = NZ2

  call setdim
  do iv=1, 2
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = var5d(:,:,:,:,iv)
  enddo

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),2,set,'FFT2D_k_omega',4)

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = get_ndate()  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename

  iv_i = 1
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  var_i_name(1) = 'fcr_'//trim(var_name(1))
  var_i_name(2) = 'fci_'//trim(var_name(1))
  call getdim(file_i(iv_i),var_i_name(iv_i))
  nk = k_largest+1
  if (k_largest == -999)  nk = nx
  allocate( kwn(nk) )
  kwn(:) = lon(:)

  call get_iouter(lat    ,lat_rng    , iy2)
  call get_iouter(p*(-1.),p_rng*(-1.), iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  nt0 = nmon*31*nhour
  nt = nmon*30*nhour

  nome = nt/(period_smallest*nhour)*2
  if (period_smallest == -999)  nome = nt
  allocate( ome(nome) )
  do i=1, nome
    ome(i) = float(i-1)
  enddo

  ovarname(1) = trim(var_i_name(1))
  ovarname(2) = trim(var_i_name(2))

  allocate( var_k(nk,ny2,nz2,2), fc_var_k(nt0,nk,ny2,nz2), fc_var(nt) )
  allocate( var5d(nk,nome,ny2,nz2,2), t(nt0), wgt(nt) )
  var5d(:,:,:,:,:) = 0.  ;  t(:) = 0.  ;  wgt(:) = 1.

  if (nmon_patch /= 0) then
    npad = nmon_patch*30*nhour
    do i=1, npad
      wgt(i) = sin(pi05*float(i)/float(30*nhour+1))
      wgt(nt+1-i) = wgt(i)
    enddo
  end if

  END subroutine initialize

  SUBROUTINE get_var

  iv_i = 1
  file_i(1) = get_ifilename()
  file_i(2) = file_i(1)
  inquire(file=trim(file_i(1)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(1)),' not found.'  ;  STOP
  end if

  print*, trim(file_i(1))
  iv_i = 1  ;  var_k(:,:,:,1) = get_ivara3d(1,nk,iy2(1),ny2,iz2(1),nz2)
  iv_i = 2  ;  var_k(:,:,:,2) = get_ivara3d(1,nk,iy2(1),ny2,iz2(1),nz2)
  fc_var_k(i_time,:,:,:) = cmplx(var_k(:,:,:,1),var_k(:,:,:,2))
  print*, 'time index :', it_i(1)

  END subroutine get_var

  SUBROUTINE setdim

  do iv=1, 2
    set(iv)%vname = trim(ovarname(iv))
    set(iv)%axis = (/'k_wn     ','ome_fr','latitude','p'/) 
    set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = kwn
    set(iv)%axis2 = ome
    set(iv)%axis3 = lat(iy2(1):iy2(2))
    set(iv)%axis4 = p  (iz2(1):iz2(2))
  enddo
    
  END subroutine setdim

  SUBROUTINE finalize

  if ( allocated(wgt) )  deallocate(wgt)
  deallocate( var5d, var_k, fc_var_k, fc_var )
  deallocate( kwn , ome )
  deallocate( lon, lat, p, t )
  do iv=1, 2
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program FFT_t2nd


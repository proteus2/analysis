PROGRAM FFT_UM_t2nd_int

  use fft
  use hadgem
  use netio

  implicit none

  include 'c_math.inc'
  include 'c_phys.inc'

  integer ::  period_smallest, k_largest, nmon_patch
  character(len=32) ::  u_intp

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ VAR_NAME, U_INTP, PERIOD_SMALLEST, K_LARGEST,         &
                   LAT_RNG, Z_RNG, NMON_PATCH
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FILE_I_HEAD, FILE_I_FORM,       &
                    FILE_I_XXXX, VAR_I_NAME, FILE_I_HEAD2, FILE_I_FORM2, &
                    FILE_I_XXXX2, VAR_I_NAME2, FILE_O

  integer ::  imon, ihour, i_time, npad
  integer ::  tag_uintp
  integer ::  nk, nome, iy2(2), iz2(2), ny2, nz2, iz3
  integer ::  i,j,k
  character(len=32), dimension(2) ::  ovarname

  real,    dimension(:,:,:,:,:), allocatable ::  var5d
  real,    dimension(:,:,:,:),   allocatable ::  var_k
  real,    dimension(:,:,:),     allocatable ::  ub
  real,    dimension(:,:),       allocatable ::  dl, s_rcosphi, uba
  real,    dimension(:),         allocatable ::  kwn, ome, wgt
  complex, dimension(:,:,:,:),   allocatable ::  fc_var_k
  complex, dimension(:,:,:),     allocatable ::  expdphi
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
  if (opt_30d == 0)  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------
  hour = 0

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
  if ( imon == nmon+1 .and. hour >= hh(1) )  EXIT L_MON

  if ( i_time == 0 .and. hour < hh(1) ) then
    hour = hour + 24/nhour  ;  CYCLE L_HOUR
  end if

  i_time = i_time + 1

  day_from_ref = get_dayfromref(year,mon,date,hour)
  if (hh(3) == 1)  day_from_ref = day_from_ref - 0.5*(24/nhour)/24.

  ! get variable
  call get_var

  t(i_time) = day_from_ref

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

! CALCULATE FFT

  do k=1, nz2
  do j=1, ny2
  do i=1, nk
    call fft1d_f(nt,fc_var_k(:,i,j,k),fc_var)
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

  call outnc(trim(file_o),2,set,'FFT2D_k_omega')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)
  if (hh(3) == 1)  day_from_ref = day_from_ref - 0.5*(24/nhour)/24.

  iv_i = 1
  file_i(1) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))
  nk = k_largest+1
  if (k_largest == -999)  nk = nx
  allocate( kwn(nk) )
  kwn(:) = lon(:)

  call get_iouter(lat    ,lat_rng, iy2)
  call get_iouter(ht/1.e3,z_rng  , iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  nt = nmon*ndate*nhour
  if ( opt_30d == 0 .and. nmon /= 1 ) then
    print*, 'not yet programmed (for opt_30d = 0)'  ;  STOP
  end if

  nome = nt/(period_smallest*nhour)*2
  if (period_smallest == -999)  nome = nt
  allocate( ome(nome) )
  do i=1, nome
    ome(i) = float(i-1)
  enddo

  ovarname(1) = trim(var_i_name(1))
  ovarname(2) = trim(var_i_name(2))

  allocate( var_k(nk,ny2,nz2,2), fc_var_k(nt,nk,ny2,nz2), fc_var(nt) )
  allocate( var5d(nk,nome,ny2,nz2,2), t(nt), wgt(nt) )
  var5d(:,:,:,:,:) = 0.  ;  t(:) = 0.  ;  wgt(:) = 1.
  allocate( ub(ny2,nz2,2), s_rcosphi(nk,ny2), dl(ny2,nz2) )
  allocate( expdphi(nk,ny2,nz2) )
  ub = 0.  ;  dl = 0.
  expdphi(:,:,:) = 0.

  if (nmon_patch /= 0) then
    npad = nmon_patch*30*nhour
    do i=1, npad
      wgt(i) = sin(pi05*float(i)/float(30*nhour+1))
      wgt(nt+1-i) = wgt(i)
    enddo
  end if

  do j=iy2(1), iy2(2)
  do i=1, nk
    s_rcosphi(i,j) = kwn(i)/(r_earth*cos(lat(j)*deg2rad))
  enddo
  enddo

  tag_uintp = 0
  if ( .not. ( trim(u_intp) == '' .or. trim(u_intp) == '-999' .or. &
       trim(u_intp) == '0' ) ) then
    tag_uintp = 1
    allocate( uba(ny2,nz2) )
    iz3 = iz2(2) + 1
    if (iz2(2) == nz)  iz3 = nz
    if (trim(var_name(iv_i)) == 'theta')  tag_uintp = tag_uintp + 10
    if ( trim(u_intp) == 'zr' .or. trim(u_intp) == 'zt' )                &
       tag_uintp = tag_uintp + 100
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
  fc_var_k(i_time,2:,:,:) = wgt(i_time)*fc_var_k(i_time,2:,:,:)

  call switch_para_in

  iv_i = 1
  file_i(1) = get_ifilename()
  inquire(file=trim(file_i(1)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(1)),' not found.'  ;  STOP
  end if

  print*, trim(file_i(1))
  ub(:,:,2) = reshape( get_ivara3d(1,1,iy2(1),ny2,iz2(1),nz2),           &
                       (/ny2,nz2/) )
  if (tag_uintp /= 0) then
    if (tag_uintp < 100) then
      print*, 'interpolate'
      uba(:,1:nz2-1) = ub(:,2:nz2,2)
      uba(:,nz2) = reshape( get_ivara3d(1,1,iy2(1),ny2,iz3,1), (/ny2/) )
      ub(:,:,2) = 0.5*(ub(:,:,2)+uba(:,:))
    else
      print*, 'Not yet coded.'
    end if
  end if

  call switch_para_in

  if (i_time > 1) then
    dl(:,:) = dl(:,:) + 0.5*(ub(:,:,1) + ub(:,:,2))*float(86400/nhour)
    expdphi(:,:,:) = exp( spread(s_rcosphi(:,:),3,nz2)* &
                          spread(dl(:,:),1,nk)*(0.,1.) )
    fc_var_k(i_time,2:,:,:) = fc_var_k(i_time,2:,:,:)*expdphi(2:,:,:)
  end if
  print*, 'time index, weight :', it_i(1), wgt(i_time)

  ub(:,:,1) = ub(:,:,2)

  END subroutine get_var

  SUBROUTINE setdim

  do iv=1, 2
    set(iv)%vname = trim(ovarname(iv))
    set(iv)%axis = (/'k_wn     ','ome_fr','latitude','ht'/) 
    set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = kwn
    set(iv)%axis2 = ome
    set(iv)%axis3 = lat(iy2(1):iy2(2))
    set(iv)%axis4 = ht (iz2(1):iz2(2))
  enddo
    
  END subroutine setdim

  SUBROUTINE finalize

  if ( allocated(wgt) )  deallocate(wgt)
  deallocate( var5d, var_k, fc_var_k, fc_var )
  deallocate( kwn , ome )
  deallocate( lon, lat, ht, t )
  do iv=1, 2
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program FFT_UM_t2nd_int


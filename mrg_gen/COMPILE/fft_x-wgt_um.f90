PROGRAM FFT_UM_x

  use fft
  use hadgem
  use netio

  implicit none

  include 'c_phys.inc'

  integer ::  k_largest
  character(len=32) ::  v_intp

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ VAR_NAME, V_INTP, K_LARGEST, LAT_RNG, Z_RNG
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, FILE_O

  integer ::  imon, ihour, i_time
  integer ::  tag_intp
  integer ::  nk, iy2(2), iz2(2), ny2, nz2, nv
  integer ::  i,j,k
  character(len=32), dimension(2) ::  ovarname

  real,    dimension(:,:,:,:,:), allocatable ::  var5d
  real,    dimension(:,:,:),     allocatable ::  var, wgt, var_th, wgt_r
  real,    dimension(:),         allocatable ::  kwn
  complex, dimension(:),         allocatable ::  fc_var

  type(vset), dimension(2) ::  set

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

  ! get variable
  call get_var

  ! calculate FFT
  do k=1, nz2
  do j=1, ny2
    call fft1d_f(nx,var(:,j,k),fc_var)
    var5d(:,j,k,i_time,1) = real (fc_var(1:nk))
    var5d(:,j,k,i_time,2) = aimag(fc_var(1:nk))
  enddo
  enddo

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

  nt = i_time

  nd1a = NK
  nd2a = NY2
  nd3a = NZ2
  nd4a = NT

  call setdim
  do iv=1, 2
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = var5d(:,:,:,1:nd4a,iv)
  enddo

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),2,set,'FFT1D_k')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 1
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  tag_intp = 0
  if ( trim(v_intp) == '' .or. trim(v_intp) == '-999' .or. &
       trim(v_intp) == '0' ) then
    call getdim(file_i(iv_i),var_name(iv_i))
  else
    if (trim(v_intp) == trim(var_name(1))) then
      tag_intp = 1
      if (trim(var_name(2)) == 'rho'  )  tag_intp = tag_intp + 10
    else
      tag_intp = 2
      if (trim(var_name(1)) == 'theta')  tag_intp = tag_intp + 10
    end if
    call getdim(file_i(mod(tag_intp,10)),trim(v_intp))
  end if
  nk = k_largest+1
  if (k_largest == -999)  nk = nx/2+1
  allocate( kwn(nk) )
  do i=1, nk
    kwn(i) = float(i-1)
  enddo

  call get_iouter(lat    ,lat_rng, iy2)
  call get_iouter(ht/1.e3,z_rng  , iz2)

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  var_i_name(1) = trim(var_name(1))
  ovarname(1) = trim(var_name(1))

  var_i_name(2) = trim(var_name(2))
  ovarname(1) = trim(var_name(2))//trim(var_name(1))
  allocate( wgt(nx,ny2,nz2) )

  ovarname(2) = 'fci_'//trim(ovarname(1))
  ovarname(1) = 'fcr_'//trim(ovarname(1))

  allocate( var(nx,ny2,nz2), fc_var(nx) )
  if ( mod(tag_intp,10) == 1 )  allocate( wgt_r (nx,ny2,1:nz2+1) )
  if ( mod(tag_intp,10) == 2 )  allocate( var_th(nx,ny2,0:nz2  ) )
  allocate( var5d(nk,ny2,nz2,nmon*31*nhour,2), t(nmon*31*nhour) )
  var5d(:,:,:,:,:) = 0.  ;  t(:) = 0.

  END subroutine initialize

  SUBROUTINE get_var

  nv = 2

  do iv_i=1, 2  ! should be "iv_i" for get_ifilename
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 ) then
      print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
    end if
  enddo

  print*, trim(file_i(1))
  iv_i = 1
  if ( mod(tag_intp,10) /= 2 ) then
    var(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  else
    print*, 'interpolate ', trim(var_i_name(iv_i))
    if (iz2(1) > 1) then
      var_th(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1)-1,nz2+1)
    else
      var_th(:,:,1:nz2) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
      var_th(:,:,0) = var_th(:,:,1)
    end if
    if (tag_intp > 10)  var_th(:,:,:) = log(var_th(:,:,:))
    do k=1, nz2
      var(:,:,k) = 0.5*(var_th(:,:,k-1)+var_th(:,:,k))
    enddo
    if (tag_intp > 10)  var(:,:,:) = exp(var(:,:,:))
  end if

  print*, trim(file_i(2))
  iv_i = 2
  if ( mod(tag_intp,10) /= 1 ) then
    wgt(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  else
    print*, 'interpolate ', trim(var_i_name(iv_i))
    if (z_rng(2) /= -999.) then
      wgt_r(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2+1)
    else
      wgt_r(:,:,1:nz2) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
      wgt_r(:,:,nz2+1) = wgt_r(:,:,nz2)
    end if
    if (tag_intp > 10)  wgt_r(:,:,:) = log(wgt_r(:,:,:))
    do k=1, nz2
      wgt(:,:,k) = 0.5*(wgt_r(:,:,k)+wgt_r(:,:,k+1))
    enddo
    if (tag_intp > 10)  wgt(:,:,:) = exp(wgt(:,:,:))
  end if
  var(:,:,:) = var(:,:,:)*wgt(:,:,:)

  print*, 'time index :', it_i(1:nv)

  END subroutine get_var

  SUBROUTINE setdim

  do iv=1, 2
    set(iv)%vname = trim(ovarname(iv))
    set(iv)%axis = (/'k_wn     ','latitude','ht','t'/) 
    set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = kwn
    set(iv)%axis2 = lat(iy2(1):iy2(2))
    set(iv)%axis3 = ht (iz2(1):iz2(2))
    set(iv)%axis4 = t(1:nd4a)
  enddo
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( wgt )
  deallocate( var5d, var, fc_var )
  if ( allocated(wgt_r ) )  deallocate( wgt_r  )
  if ( allocated(var_th) )  deallocate( var_th )
  deallocate( kwn )
  deallocate( lon, lat, ht, t )
  do iv=1, 2
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program FFT_UM_x


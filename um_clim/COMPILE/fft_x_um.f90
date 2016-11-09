PROGRAM FFT_UM_x

  use fft
  use hadgem
  use netio

  implicit none

  include 'c_phys.inc'

  integer ::  k_largest, h_intp
  character(len=32) ::  v_intp

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ VAR_NAME, V_INTP, H_INTP, K_LARGEST, LAT_RNG, Z_RNG
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, FILE_ALT, VAR_ALT, FILE_O

  integer ::  imon, ihour, i_time
  integer ::  tag_vintp
  integer ::  nk, iy2(2), iz2(2), ny2, nz2, iy3, ny3, nv
  integer ::  i,j,k
  character(len=32), dimension(2) ::  ovarname

  real,    dimension(:,:,:,:,:), allocatable ::  var5d
  real,    dimension(:,:,:),     allocatable ::  var, var_in
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
  if (hh(3) == 1)  day_from_ref = day_from_ref - 0.5*(24/nhour)/24.

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
  if (hh(3) == 1)  day_from_ref = day_from_ref - 0.5*(24/nhour)/24.

  iv_i = 1
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  tag_vintp = 0
  if ( trim(v_intp) == '' .or. trim(v_intp) == '-999' .or. &
       trim(v_intp) == '0' ) then
    call getdim(file_i(iv_i),var_name(iv_i))
  else
    tag_vintp = 1
    if (trim(var_name(iv_i)) == 'theta')  tag_vintp = tag_vintp + 10
    if ( trim(v_intp) == 'zr' .or. trim(v_intp) == 'zt' )                &
       tag_vintp = tag_vintp + 100
    if (tag_vintp < 100) then
      iv_i = 2
      file_i(iv_i) = get_ifilename()
      call getdim(file_i(iv_i),trim(v_intp))
      iv_i = 1
    else
      call getdim(file_i(iv_i),var_name(iv_i))
    end if
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

  if (tag_vintp > 100) then
    x1_i = 1   ;  y1_i = iy2(1)   ;  z1_i = iz2(1)   ! for getalt
    nx_i = nx  ;  ny_i = ny2      ;  nz_i = nz2
    if (h_intp /= 0)  then
      y1_i = y1_i - 1  ;  ny_i = ny_i + 2
    end if
    call getalt
  end if

  if ( tag_vintp > 100 .and. iz2(1) /= 1 ) then
    print*, 'not yet programmed.'  ;  STOP
  end if

  var_i_name(1) = trim(var_name(1))
  ovarname(1) = trim(var_name(1))
  ovarname(2) = 'fci_'//trim(ovarname(1))
  ovarname(1) = 'fcr_'//trim(ovarname(1))

  allocate( var(nx,ny2,nz2), fc_var(nx) )
  if (h_intp == 0) then
    iy3 = iy2(1)  ;  ny3 = ny2
    allocate( var_in(nx,ny2,0:nz2) )
  else
    iy3 = iy2(1) - 1  ;  ny3 = ny2 + 2
    allocate( var_in(0:nx+1,0:ny2+1,0:nz2) )
  end if
  allocate( var5d(nk,ny2,nz2,nmon*31*nhour,2), t(nmon*31*nhour) )
  var5d(:,:,:,:,:) = 0.  ;  t(:) = 0.

  END subroutine initialize

  SUBROUTINE get_var

  nv = 1
  iv_i = 1

  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if

  print*, trim(file_i(iv_i))
  if (tag_vintp == 0) then
    var_in(1:nx,:,1:) = get_ivara3d(1,nx,iy3,ny3,iz2(1),nz2)
  else if (tag_vintp < 100) then
    print*, 'interpolate'
    if (iz2(1) > 1) then
      var_in(1:nx,:,:) = get_ivara3d(1,nx,iy3,ny3,iz2(1)-1,nz2+1)
    else
      var_in(1:nx,:,1:nz2) = get_ivara3d(1,nx,iy3,ny3,iz2(1),nz2)
      var_in(1:nx,:,0) = var_in(1:nx,:,1)
    end if
    if (tag_vintp == 11)  var_in(1:nx,:,:) = log(var_in(1:nx,:,:))
    do k=nz2, 1, -1
      var_in(1:nx,:,k) = 0.5*(var_in(1:nx,:,k-1)+var_in(1:nx,:,k))
    enddo
    if (tag_vintp == 11)  var_in(1:nx,:,:) = exp(var_in(1:nx,:,:))
  else
    var_in(1:nx,:,1:) = get_ivara3d(1,nx,iy3,ny3,iz2(1),nz2)
    ! iz2(1) should be 1 for use of the subroutine below.
    nx_i = nx  ;  ny_i = ny3  ;  nz_i = min(k_const_rho+1,nz2)
    if (trim(v_intp) == 'zr') then
      print*, 'interpolate (rho2z)'
      call v_linintp_rho2z(var_in(1:nx,:,1:nz_i))
    else
      print*, 'interpolate (th2zth)'
      call v_linintp_th2zth(var_in(1:nx,:,1:nz_i))
    end if
  end if

  if (h_intp == 0) then
    var(:,:,:) = var_in(:,:,1:)
  else
    var_in(0   ,:,:) = var_in(nx,:,:)
    var_in(nx+1,:,:) = var_in(1 ,:,:)
    var(:,:,:) = 0.5*var_in(1:nx,1:ny2,1:) + 0.125* &
         ( var_in(0:nx-1,1:ny2  ,1:)+var_in(2:nx+1,1:ny2  ,1:)+ &
           var_in(1:nx  ,0:ny2-1,1:)+var_in(1:nx  ,2:ny2+1,1:) )
  end if

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

  deallocate( var5d, var, fc_var )
  deallocate( var_in )
  deallocate( kwn )
  deallocate( lon, lat, ht, t )
  do iv=1, 2
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program FFT_UM_x


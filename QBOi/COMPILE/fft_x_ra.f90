PROGRAM FFT_x

  use fft
  use reanal
  use netio

  implicit none

  integer ::  k_largest, h_intp

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE
  namelist /PARAM/ VAR_NAME, H_INTP, K_LARGEST, LAT_RNG, P_RNG
  namelist /FILEIO/ NT_F4, P_PREDEF, MISSV, FILE_I_HEAD, FILE_I_FORM,    &
                    FILE_I_XXXX, VAR_I, FILE_O

  integer ::  imon, ihour, i_time
  integer ::  nk, iy2(2), iz2(2), ny2, nz2, iy3, ny3, nv
  integer ::  iy2o(2), iz2o(2)
  integer ::  i,j,k, tmpi
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
  ndate = get_ndate()

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
print*, year, mon, date, hour

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

  ndate = get_ndate()  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename

  iv_i = 1
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_name(iv_i))
  nk = k_largest+1
  if (k_largest == -999)  nk = nx/2+1
  allocate( kwn(nk) )
  do i=1, nk
    kwn(i) = float(i-1)
  enddo

  call get_iouter(lat,lat_rng, iy2o)
  iy2(:) = iy2o(:)
  if ( l_rev(2) ) then
    iy2(:) = ny + 1 - iy2o(:)
    tmpi = iy2(1)
    iy2(1) = iy2(2)  ;  iy2(2) = tmpi
  end if
  call get_iouter(p*(-1.),p_rng*(-1.), iz2o)
  iz2(:) = iz2o(:)
  if ( l_rev(3) ) then
    iz2(:) = nz + 1 - iz2o(:)
    tmpi = iz2(1)
    iz2(1) = iz2(2)  ;  iz2(2) = tmpi
  end if

  ny2 = iy2(2) - iy2(1) + 1  ;  nz2 = iz2(2) - iz2(1) + 1

  var_i_name(1) = trim(var_name(1))
  ovarname(1) = trim(var_name(1))
  ovarname(2) = 'fci_'//trim(ovarname(1))
  ovarname(1) = 'fcr_'//trim(ovarname(1))

  allocate( var(nx,ny2,nz2), fc_var(nx) )
  if (h_intp == 0) then
    iy3 = iy2(1)  ;  ny3 = ny2
    allocate( var_in(nx,ny2,nz2) )
  else
    iy3 = iy2(1) - 1  ;  ny3 = ny2 + 2
    allocate( var_in(0:nx+1,0:ny2+1,nz2) )
  end if
  allocate( var5d(nk,ny2,nz2,nmon*31*nhour,2), t(nmon*31*nhour) )
  var5d(:,:,:,:,:) = 0.  ;  t(:) = 0.

  END subroutine initialize

  SUBROUTINE get_var

  nv = 1

  iv_i=1     ! should be "iv_i" for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if

  print*, trim(file_i(1))
  iv_i = 1
  var_in(1:nx,:,:) = get_ivara3d(1,nx,iy3,ny3,iz2(1),nz2)

  if (h_intp == 0) then
    var(:,:,:) = var_in(:,:,:)
  else
    var_in(0   ,:,:) = var_in(nx,:,:)
    var_in(nx+1,:,:) = var_in(1 ,:,:)
    var(:,:,:) = 0.5*var_in(1:nx,1:ny2,:) + 0.125* &
         ( var_in(0:nx-1,1:ny2  ,:)+var_in(2:nx+1,1:ny2  ,:)+ &
           var_in(1:nx  ,0:ny2-1,:)+var_in(1:nx  ,2:ny2+1,:) )
  end if

  print*, 'time index :', it_i(1:nv)

  END subroutine get_var

  SUBROUTINE setdim

  do iv=1, 2
    set(iv)%vname = trim(ovarname(iv))
    set(iv)%axis = (/'k_wn     ','latitude','p','t'/) 
    set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = kwn
    set(iv)%axis2 = lat(iy2o(1):iy2o(2))
    set(iv)%axis3 = p  (iz2o(1):iz2o(2))
    set(iv)%axis4 = t(1:nd4a)
  enddo
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var5d, var, fc_var )
  deallocate( var_in )
  deallocate( kwn )
  deallocate( lon, lat, p, t, dim4, t2pt )
  do iv=1, 2
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program FFT_x


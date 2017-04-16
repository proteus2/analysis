PROGRAM PRECIP_PDF

  use hadgem
  use netio

  implicit none

  integer, parameter ::  nv = 4

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
!  namelist /PARAM/ 
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, VAR_I_NAME, FILE_O

  integer ::  imon, ihour, i_time
  integer ::  j,k
  real    ::  wgt
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:)  , allocatable ::  prcp
  real, dimension(:,:,:), allocatable ::  var3d
  real, dimension(:,:)  , allocatable ::  p_zm, dp2_zm

  type(vset), dimension(nv) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
!  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  read(10, ANALCASE)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  year = yyyy
  mon  = mm(1)

  call initialize

  var3d(:,:,:) = 0.
  i_time = 0

  L_MON:  DO imon=1, nmon+1
  !---------------------------------------------------------------------
  if (opt_30d == 0)  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
  i_time = i_time + 1

  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 1  ! for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
 
  ! get variable
  call get_1var

  ! calculate zonal mean
!  p_zm(:,:) = sum(p3d, dim=1)/float(nx)
!  do k=1, nz
!  do j=1, ny
!    p3d(:,j,k) = p3d(:,j,k) - p_zm(j,k)
!  enddo
!  enddo
!  dp2_zm(:,:) = sum(p3d(:,:,:)*p3d(:,:,:), dim=1)/float(nx)

  wgt = 1.
  if ( i_time == 1 .or. imon == nmon+1 )  wgt = 0.5

  var3d(:,:,1) = var3d(:,:,1) + p_zm(:,:)*(1.e-2*wgt)            ! [hPa]
  var3d(:,:,2) = var3d(:,:,2) + dp2_zm(:,:)*(1.e-4*wgt)          ! [hPa^2]
  var3d(:,:,3) = var3d(:,:,3) + p_zm(:,:)*p_zm(:,:)*(1.e-4*wgt)  ! [hPa^2]

  hour = hour + 24/nhour

  if (imon == nmon+1)  EXIT L_MON
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

  nt = i_time - 1

  deallocate( prcp, p_zm, dp2_zm )

  var3d(:,:,3) = var3d(:,:,2) + var3d(:,:,3)

  ! calculate time mean
  var3d(:,:,1:3) = var3d(:,:,1:3)/float(nt)

  ! additional diagnostics
  var3d(:,:,4) = sqrt(var3d(:,:,3)/(var3d(:,:,1)*var3d(:,:,1)) - 1.)

  ! TOTAL_VAR = avg_yt(p2_m) - [avg_yt(p_m)]^2


  nd1a = NY
  nd2a = NZ
  nd3a = 1
  nd4a = 1

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,1,1) = var3d(:,:,iv)
  enddo


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'mean and variances of p at model levels')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 1  ! for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  allocate( var3d(ny,nz,nv) )
  var3d(:,:,:) = 0.
  allocate( prcp(nx,ny), p_zm(ny,nz), dp2_zm(ny,nz) )

  ovarname = (/'p_m   ','dp2_m','p2_m','sd_loc'/)

  END subroutine initialize

  SUBROUTINE get_1var

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz   ! for get_ivar3d

  ! read var.
  print*, trim(file_i(1))
  iv_i = 1  ;  prcp(:,:) = get_ivar3d()
  print*, 'time index :', it_i(1)

  END subroutine get_1var

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lat  ','z ',' ',' '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat
  set(iv)%axis2 = ht
  set(iv)%axis3 = -999.
  set(iv)%axis4 = -999.
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var3d )
  deallocate( lon, lat, ht, ht_th )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program PRECIP_PDF


PROGRAM DCHM_PDF_MERGE

  use hadgem
  use netio

  implicit none

  integer, parameter ::  nstat = 3
  ! 1: mean  /  2: mean of square  /  3: variance

  integer ::  yyyy2(2)

  namelist /ANALCASE/ EXPNAME, YYYY2, MM
!  namelist /PARAM/ PDF_RNG_LN10, LAT_RNG
  namelist /FILEIO/ FILE_I_HEAD, FILE_I_FORM, FILE_I_XXXX, VAR_I_NAME, FILE_O

  integer ::  nv, nyr, iyr, imon, i_time, npdfx, ncid
  real    ::  dpdfx
  integer ::  ib
  logical ::  tmpl
  character(len=128) ::  ftitle

  real, dimension(:,:,:), allocatable ::  var_m, var_m2
  real, dimension(:)    , allocatable ::  pdfx, npop, npop2

  type(vset), dimension(:), allocatable ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
!  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  read(10, ANALCASE)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  call initialize

  npop2(:) = 0.
  var_m2(:,:,:) = 0.
 
  i_time = 0

  L_YEAR:  DO iyr=1, nyr
  !---------------------------------------------------------------------
  year = yyyy2(1) + (iyr-1)
  mon = mm(1)

  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------
  i_time = i_time + 1

  ! get variable
  call get_vars

  npop2(:) = npop2(:) + npop(:)
  do iv=1, nv
    var_m2(:,iv,1) = var_m2(:,iv,1) + var_m(:,iv,1)*npop(:)
    var_m2(:,iv,2) = var_m2(:,iv,2) + var_m(:,iv,2)*npop(:)
  enddo

  mon = mon + 1
  if (mon == 13) then
    year = year + 1  ;  mon = 1
  end if
  !---------------------------------------------------------------------
  ENDDO  L_MON

  !---------------------------------------------------------------------
  ENDDO  L_YEAR

  nt = i_time

  deallocate( npop, var_m )

  do ib=1, npdfx
    if (npop2(ib) /= 0.)  var_m2(ib,:,:2) = var_m2(ib,:,:2)/npop2(ib)
  enddo
  var_m2(:,:,3) = var_m2(:,:,2) - var_m2(:,:,1)*var_m2(:,:,1)

  write(6,*) 'total columns  : ', int(sum(npop2(:)))
  write(6,*) 'inside the PDF : ', int(sum(npop2(2:npdfx-1)))
  write(6,*) 'outside        : ', int(npop2(1)), int(npop2(npdfx))


  nd1a = NPDFX
  nd2a = NSTAT
  nd3a = 1
  nd4a = 1

  allocate( set(nv+1) )

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,1,1) = var_m2(:,iv,:)
  enddo

  iv = nv+1
  call setdim_npop
  allocate( set(iv)%var_out(nd1a,1,1,1) )
  set(iv)%var_out(:,1,1,1) = npop2(:)


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv+1,set,trim(ftitle))

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nyr = yyyy2(2) - yyyy2(1) + 1
  nmon = mm(2)

  do iv=2, 99
    if ( trim(var_i_name(iv)) == '-999' )  nv = iv - 2  ! excluding N_pop
  enddo

  year = yyyy2(1)  ;  mon = mm(1)

  day1 = -999  ;  date = -999  ! for get_ifilename
  iv_i = 1  ! for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call opennc(file_i(1),ncid)
  call dilen(ncid,'x_pdf', npdfx,tmpl)
  allocate( pdfx(npdfx) )
  call get1d(ncid,'x_pdf',npdfx, pdfx)
  call gettitle(ncid,ftitle)
  call closenc(ncid)

  dpdfx = pdfx(5) - pdfx(4)

  allocate( npop(npdfx), var_m(npdfx,nv,nstat) )
  npop(:) = 0.  ;  var_m(:,:,:) = 0.

  allocate( npop2(npdfx), var_m2(npdfx,nv,nstat) )
  npop2(:) = 0.  ;  var_m2(:,:,:) = 0.

  END subroutine initialize

  SUBROUTINE get_vars

  day1 = -999  ;  date = -999  ! for get_ifilename
  iv_i = 1  ! for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
 
  ! read var.
  print*, trim(file_i(1))

  call opennc(file_i(1),ncid)
  iv = 1  ;  call get1d(ncid,var_i_name(iv),npdfx, npop)
  do iv=2, nv+1
    call get2d(ncid,var_i_name(iv),npdfx,nstat, var_m(:,iv-1,:))
  enddo
  call closenc(ncid)

  END subroutine get_vars

  SUBROUTINE setdim

  set(iv)%vname = trim(var_i_name(iv+1))
  set(iv)%axis = (/'x_pdf','stat',' ',' '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = pdfx(:)
  set(iv)%axis2 = (/1.,2.,3/)
  set(iv)%axis3 = -999.
  set(iv)%axis4 = -999.
    
  END subroutine setdim

  SUBROUTINE setdim_npop

  set(iv)%vname = trim(var_i_name(1))
  set(iv)%axis = (/'x_pdf',' ',' ',' '/) 
  set(iv)%nd(:) = (/npdfx,1,1,1/)
  allocate( set(iv)%axis1(npdfx) )
  allocate( set(iv)%axis2(1) )
  allocate( set(iv)%axis3(1) )
  allocate( set(iv)%axis4(1) )
  set(iv)%axis1 = pdfx(:)
  set(iv)%axis2 = -999.
  set(iv)%axis3 = -999.
  set(iv)%axis4 = -999.
    
  END subroutine setdim_npop

  SUBROUTINE finalize

  deallocate( var_m2, npop2, pdfx )
  do iv=1, nv+1
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo
  deallocate( set )

  END subroutine finalize


END program DCHM_PDF_MERGE


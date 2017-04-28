PROGRAM DCHD_PDF

  use hadgem
  use netio

  implicit none

  integer, parameter ::  nv = 15, nstat = 3, nzu = 33
  ! 1: mean  /  2: mean of square  /  3: variance

  real ::  pdf_rng(3)
  character(len=128) ::  file_o2

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ PDF_RNG, LAT_RNG
  namelist /FILEIO/ DAY1, NDAY_I, MISSV, FID, FILE_I_HEAD, FILE_I_FORM,  &
                    FILE_I_XXXX, VAR_I_NAME, FILE_ALT, VAR_ALT, FILE_O,  &
                    FILE_O2

  integer ::  imon, ihour, i_time, npdfx, iy2(2), ny2, ncol(2), ncolm
  real    ::  dpdfx
  integer ::  i,j,k,l,ib,ii, l2(2)
  character(len=32), dimension(nv) ::  ovarname

  real   , dimension(:,:,:,:), allocatable ::  var_m
  real   , dimension(:,:,:)  , allocatable ::  vari, var_u, var_v, vars, &
                                               ind_mc
  real   , dimension(:,:)    , allocatable ::  varx
  real   , dimension(:)      , allocatable ::  pdfx
  integer, dimension(:,:)    , allocatable ::  pdfi, npop

  type(vset), dimension(nv+1,2) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  year = yyyy
  mon  = mm(1)

  call initialize

  npop(:,:) = 0
  var_m(:,:,:,:) = 0.
 
  i_time = 0

  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------
  if (opt_30d == 0)  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
  i_time = i_time + 1

  day_from_ref = get_dayfromref(year,mon,date,hour)

  ! problem found in the data for 00 UTC 1 each month
  if ( date == 1 .and. hour == 0 ) then
    hour = hour + 24/nhour
    i_time = i_time - 1
    CYCLE
  end if

  ! get variable
  call get_vars

  if ( all( ncol(:) == 0 ) ) then
    hour = hour + 24/nhour
    CYCLE
  end if

  do ii=1, 2
  pdfi(:ncol(ii),ii) = (int( (varx(:ncol(ii),ii) - pdf_rng(1))/dpdfx + 1.5 )-1)+2
  do l=1, ncol(ii)
    pdfi(l,ii) = min(npdfx,max(1,pdfi(l,ii)))
  enddo
  enddo

  do ii=1, 2
  do l=1, ncol(ii)
    ib = pdfi(l,ii)
    npop(ib,ii) = npop(ib,ii) + 1
    var_m(ib,:,1,ii) = var_m(ib,:,1,ii) + vars(l,:,ii)
    var_m(ib,:,2,ii) = var_m(ib,:,2,ii) + vars(l,:,ii)*vars(l,:,ii)
  enddo
  enddo

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

  deallocate( varx, vars, vari, pdfi, ind_mc )

  do ii=1, 2
  do ib=1, npdfx
    if (npop(ib,ii) /= 0)  var_m(ib,:,:2,ii) =                           &
                           var_m(ib,:,:2,ii)/float(npop(ib,ii))
  enddo
  enddo
  var_m(:,:,3,:) = var_m(:,:,2,:) - var_m(:,:,1,:)*var_m(:,:,1,:)

  write(6,*) 'total columns  : ', int(sum(npop(:,:))), '/', nx*ny2*nt
  write(6,*) 'inside the PDF : ', int(sum(npop(2:npdfx-1,:)))
  write(6,*) 'outside        : ', int(sum(npop(1,:))), int(sum(npop(npdfx,:)))


  nd1a = NPDFX
  nd2a = NSTAT
  nd3a = 1
  nd4a = 1

  pdfx(1    ) = pdfx(2      ) - dpdfx*10.
  pdfx(npdfx) = pdfx(npdfx-1) + dpdfx*10.
!l  pdfx(:) = 10.**pdfx(:)

  do ii=1, 2
  do iv=1, nv
    call setdim
    allocate( set(iv,ii)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv,ii)%var_out(:,:,:,:) = 1.e32

    set(iv,ii)%var_out(:,:,1,1) = var_m(:,iv,:,ii)
  enddo
  enddo

  do ii=1, 2
  iv = nv+1
  call setdim_npop
  allocate( set(iv,ii)%var_out(nd1a,1,1,1) )
  set(iv,ii)%var_out(:,1,1,1) = float(npop(:,ii))
  enddo


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  ii = 1
  call outnc(trim(file_o),nv+1,set(:,ii),                                &
             'PDF for '//trim(var_i_name(1))//' for midlevel conv.')

  ii = 2
  call outnc(trim(file_o2),nv+1,set(:,ii),                               &
             'PDF for '//trim(var_i_name(1))//' for non-midlevel conv.')


! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nmon = mm(2)  ;  nhour = hh(2)

!  problem found in the data for 00 UTC 1 each month
!  ndate = 30  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename
  ndate = 30  ;  date = 2  ;  hour = hh(1)  ! for get_ifilename
  if (opt_30d == 0)  ndate = get_ndate()
  day_from_ref = get_dayfromref(year,mon,date,hour)

  iv_i = 3  ! for get_ifilename
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  call get_iouter(lat,lat_rng, iy2)

  ny2 = iy2(2) - iy2(1) + 1

  npdfx = int(pdf_rng(3))+2
  dpdfx = (pdf_rng(2) - pdf_rng(1))/(pdf_rng(3)-1)
  allocate( pdfx(npdfx) )
  pdfx(1) = -999.
  do ib=2, npdfx-1
    pdfx(ib) = pdf_rng(1) + (ib-2)*dpdfx
  enddo
  pdfx(npdfx) = 999.

  allocate( npop(npdfx,2), var_m(npdfx,nv,nstat,2) )
  npop(:,:) = 0  ;  var_m(:,:,:,:) = 0.

  x1_i = 1   ;  y1_i = iy2(1)   ;  z1_i = 1       ! for getalt
  nx_i = nx  ;  ny_i = ny2      ;  nz_i = nzu-1
  call getalt
  do k=1, nzu-1
    z_th(:,:,k) = z_th(:,:,k) - z_th(:,:,0)
  enddo
  z_th(:,:,0) = 0.

  ncolm = nx*ny2

  allocate( vari(nx,ny2,nv), varx(ncolm,2), vars(ncolm,nv,2), pdfi(ncolm,2) )
  allocate( ind_mc(nx,ny2,1) )
  allocate( var_u(nx,ny2,nzu), var_v(nx,ny2,nzu) )

  ovarname = (/'dchmax','zcba  ','zcta  ','rho_ct','n_q   ','n_ct  ',    &
               't_ct  ','cq_x  ','cq_y  ','u_ct  ','v_ct  ','u_sfc ',    &
               'v_sfc ','u_cb  ','v_cb  '/)

  END subroutine initialize

  SUBROUTINE get_vars

  ex0 = .TRUE.
  do iv_i=1, 5  ! for get_ifilename
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP
 
  ! read var.
  print*, trim(file_i(1))

  iv_i = 1  ;  vari(:,:,1:1) = get_ivara3d(1,nx,iy2(1),ny2,1,1)

  iv_i = 2  ;  vari(:,:,2:11) = get_ivara3d(1,nx,iy2(1),ny2,1,10)

  iv_i = 3  ;  var_u(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,1,nzu)
  iv_i = 4  ;  var_v(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,1,nzu)

  iv_i = 5  ;  ind_mc(:,:,1:1) = get_ivara3d(1,nx,iy2(1),ny2,1,1)

  print*, 'time index :', it_i(1)

  vari(:,:,1) = vari(:,:,1)*3600.  ! [K/hour]

  varx(:,:) = 0.  ;  vars(:,:,:) = 0.

  l2(:) = 0
  do j=1, ny2
  do i=1, nx
    if (vari(i,j,1) == missv)  CYCLE
!    if (ind_mc(i,j,1) /= 0.)  CYCLE
    if ( ind_mc(i,j,1) /= 0. .and. vari(i,j,2) .gt. 4.e3 ) then
      ii = 1   ! midlevel convection dominated
    else
      ii = 2
    end if
    l2(ii) = l2(ii) + 1
    varx(l2(ii),ii) = vari(i,j,3) - vari(i,j,2)  ! depth = zcta - zcba
    vars(l2(ii),:11,ii) = vari(i,j,:)
    vars(l2(ii),12,ii) = var_u(i,j,1)  ! u_sfc
    vars(l2(ii),13,ii) = var_v(i,j,1)  ! v_sfc
    k = minloc(abs(vari(i,j,2) - z_th(i,j,:)),1)-1  ! zcba: vari(:,:,2)
    vars(l2(ii),14,ii) = 0.5*(var_u(i,j,k) + var_u(i,j,k+1))  ! u_cb
    vars(l2(ii),15,ii) = 0.5*(var_v(i,j,k) + var_v(i,j,k+1))  ! v_cb
  enddo
  enddo
  ncol(:) = l2(:)
!l  do ii=1, 2
!l    if (ncol(ii) /= 0)  varx(:ncol(ii),ii) = log10(varx(:ncol(ii),ii))
!l  enddo

  END subroutine get_vars

  SUBROUTINE setdim

  set(iv,ii)%vname = trim(ovarname(iv))
  set(iv,ii)%axis = (/'x_pdf','stat',' ',' '/) 
  set(iv,ii)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv,ii)%axis1(set(iv,ii)%nd(1)) )
  allocate( set(iv,ii)%axis2(set(iv,ii)%nd(2)) )
  allocate( set(iv,ii)%axis3(set(iv,ii)%nd(3)) )
  allocate( set(iv,ii)%axis4(set(iv,ii)%nd(4)) )
  set(iv,ii)%axis1 = pdfx(:)
  set(iv,ii)%axis2 = (/1.,2.,3./)
  set(iv,ii)%axis3 = -999.
  set(iv,ii)%axis4 = -999.
    
  END subroutine setdim

  SUBROUTINE setdim_npop

  set(iv,ii)%vname = 'N_pop'
  set(iv,ii)%axis = (/'x_pdf',' ',' ',' '/) 
  set(iv,ii)%nd(:) = (/npdfx,1,1,1/)
  allocate( set(iv,ii)%axis1(npdfx) )
  allocate( set(iv,ii)%axis2(1) )
  allocate( set(iv,ii)%axis3(1) )
  allocate( set(iv,ii)%axis4(1) )
  set(iv,ii)%axis1 = pdfx(:)
  set(iv,ii)%axis2 = -999.
  set(iv,ii)%axis3 = -999.
  set(iv,ii)%axis4 = -999.
    
  END subroutine setdim_npop

  SUBROUTINE finalize

  deallocate( var_m, npop, pdfx )
  deallocate( lon, lat, ht, ht_th )
  do ii=1, 2
  do iv=1, nv+1
    deallocate( set(iv,ii)%axis1, set(iv,ii)%axis2, set(iv,ii)%axis3,    &
                set(iv,ii)%axis4 )
    deallocate( set(iv,ii)%var_out )
  enddo
  enddo

  END subroutine finalize


END program DCHD_PDF


PROGRAM EPFSA1_RA

  use epf
  use reanal
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 8, nv1 = 4
  real,    parameter ::  pb_out = 500.
  character(len=8), dimension(2) ::  ivarhead = (/'fcr_','fci_'/)

  integer ::  opt_30d, k_largest, period_smallest, nmon_patch
  character(len=128) ::  wavgrp(10), file_o2

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ K_LARGEST, PERIOD_SMALLEST, LAT_RNG, P_RNG, NMON_PATCH
  namelist /FILEIO/ NT_F4, P_PREDEF, WAVGRP, FILE_I_HEAD, FILE_I_FORM,   &
                    FILE_I_XXXX, VAR_I, VAR_I_NAME, FILE_I_HEAD2,        &
                    FILE_I_FORM2, FILE_I_XXXX2, VAR_I2, VAR_I_NAME2,     &
                    FILE_O, FILE_O2

  integer ::  iz, imon, ihour, i_time
  integer ::  nk, nome, nwavgrp
  integer ::  iy2(2), iz2(2), ny2, nz2, nta
  integer ::  iy2o(2), iz2o(2)
  integer ::  i,j,k, iw, iw2, iv0, iv2, iri, tmpi
  real    ::  wgt
  character(len=32), dimension(8)  ::  ivarname

  real, dimension(:,:,:,:,:), allocatable ::  var5d
  real, dimension(:,:,:,:),   allocatable ::  vpt, vu, wpt, wu
  real, dimension(:,:,:),   allocatable ::  var3d
  real, dimension(:,:,:),   allocatable ::  um, ptm
  real, dimension(:,:),     allocatable ::  t2pt0
  real, dimension(:),       allocatable ::  lat0, p0a, kwn, ome

  type(vset),  dimension(:), allocatable ::  set
  type(vset0), dimension(:), allocatable ::  set0, set0b

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  do i=1, 10
    if (wavgrp(i) == '-999')  nwavgrp = i - 1
  enddo

! GET AXES AND INITIALIZE ARRAYS

  call initialize

  i_time = 0

  ! get variable
  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------

  ndate = get_ndate()

  do date=1, ndate
    hour = hh(1)
    do ihour=1, nhour

      i_time = i_time + 1

      call get_2var

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

  call switch_para_in

  call get_8var_init


  wgt = float(nmon+2*nmon_patch)/float(nmon+nmon_patch)

  if (pb_out /= 1000.) then
    iz = 0
    do k=2, nd4a
      if (p0a(k) < pb_out) then  ;  iz = k - 1  ;  EXIT  ;  end if
    enddo
  end if

  print*, ' Creating output file...'
  call outnc_blank(trim(file_o),nv1*nwavgrp,set0,  &
             'EP flux analysis in spectral domain', 4)
  call outnc_blank(trim(file_o2),(nv-nv1)*nwavgrp,set0b,  &
             'EP flux analysis in spectral domain', 4)

  WAVE_GRP:  DO iw=1, nwavgrp
  !---------------------------------------------------------------------

  print*, iw, '/', nwavgrp

  print*, ' Reading...'
  call get_8var

  print*, ' Calculating...'
  call epf_hydro_p_fc2(                                                  &
     nk,nome,ny2,nz2,lat0,p0a*100.,nt,um(:,:,1:nt),ptm(:,:,1:nt),        &
     h_scale,vpt,vu,wpt,wu,1.e32,                                        &
     set(1)%var_out,set(2)%var_out,set(3)%var_out,set(4)%var_out,        &
     set(5)%var_out,set(6)%var_out,set(7)%var_out,set(8)%var_out )

  do iv=1, nv
    set(iv)%var_out(:,:,:,:) = set(iv)%var_out(:,:,:,:)*wgt
    if (jk_bnd(1,iv) > 0)  &
        set(iv)%var_out(:,:,1:jk_bnd(1,iv),:) = 1.e32
    if (jk_bnd(2,iv) > 0)  &
        set(iv)%var_out(:,:,ny2+1-jk_bnd(2,iv):ny2,:) = 1.e32
    if (jk_bnd(3,iv) > 0)  &
        set(iv)%var_out(:,:,:,1:jk_bnd(3,iv)) = 1.e32
    if (jk_bnd(4,iv) > 0)  &
        set(iv)%var_out(:,:,:,nz2+1-jk_bnd(4,iv):nz2) = 1.e32
  enddo

  if (iz > 1) then
    do iv=1, nv
      set(iv)%var_out(:,:,:,:iz-1) = 1.e32
    enddo
  end if

  print*, ' Dump...'
  do iv=1, nv
    set(iv)%vname = trim(varname_epf(iv))//'_'//trim(wavgrp(iw))
  enddo
  call outnc_ss(trim(file_o),nv1,set(1:nv1))
  call outnc_ss(trim(file_o2),nv-nv1,set(nv1+1:nv))

  !---------------------------------------------------------------------
  ENDDO  WAVE_GRP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*) trim(file_o2)

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  year = yyyy
  mon  = mm(1)

  nmon = mm(2)  ;  nhour = hh(2)

  ndate = get_ndate()  ;  date = 1  ;  hour = hh(1)  ! for get_ifilename

  iv_i = 1  ! u
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if

  call getdim(file_i(iv_i),var_i_name(iv_i))

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

  allocate( lat0(ny2), p0a(nz2), t2pt0(ny2,nz2) )
  lat0(:) = lat(iy2o(1):iy2o(2))  ;  p0a(:) = p(iz2o(1):iz2o(2))
  t2pt0(:,:) = spread(t2pt(iz2o(1):iz2o(2)),1,ny2)

  allocate( um(ny2,nz2,nmon*31*nhour), ptm(ny2,nz2,nmon*31*nhour) )

  END subroutine initialize

  SUBROUTINE get_2var

  ex0 = .TRUE.
  do iv_i=1, 2
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  allocate( var3d(nx,ny2,nz2) )

  ! read 2 var.s
  print*, trim(file_i(1))
  iv_i = 1  ;  var3d(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  um(:,:,i_time) = sum(var3d(:,:,:), dim=1)/float(nx)
  print*, trim(file_i(2))
  iv_i = 2  ;  var3d(:,:,:) = get_ivara3d(1,nx,iy2(1),ny2,iz2(1),nz2)
  ptm(:,:,i_time) = sum(var3d(:,:,:), dim=1)/float(nx)
  ptm(:,:,i_time) = ptm(:,:,i_time)*t2pt0(:,:)
  print*, 'time index :', it_i(1)

  deallocate( var3d )

  END subroutine get_2var

  SUBROUTINE get_8var_init

  ex0 = .TRUE.
  do iv_i=1, 4
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  if (.not. ex0)  STOP

  deallocate( lon, lat, p, dim4, t2pt )
  p_predef(1) = -999.
  call getdim(file_i(1),trim(ivarhead(1))//trim(var_i_name(1))//'_'//trim(wavgrp(1)))
  l_rev(1:3) = .False.

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

  if ( any(lat0(ny2/2+1:ny2) /= p((ny2+1)/2:1:-1)) ) then
    print*, 'Check lat0'  ;  STOP
  end if

  if ( any(p0a /= dim4) ) then
    print*, 'Mismatch between the two pressure coordinates.'
    print*, '  File_I_1:', p0a(1), p0a(nz2)
    print*, '  File_I_2:', dim4(1), dim4(size(dim4))
    STOP
  end if

  do iv_i=1, 4
    ivarname(iv_i) = trim(var_i_name(iv_i))
  enddo

  allocate( var5d(nk,nome,ny2,nz2,4) )

  allocate( vpt(nk,nome,ny2,nz2), vu(nk,nome,ny2,nz2) )
  allocate( wpt(nk,nome,ny2,nz2), wu(nk,nome,ny2,nz2) )

  allocate( set0(nv1*nwavgrp), set0b((nv-nv1)*nwavgrp) )
  allocate( set(nv) )

  nd1a = NK
  nd2a = NOME
  nd3a = NY2
  nd4a = NZ2

  do iw=1, nwavgrp
  do iv0=1, nv
    iv = (iw-1)*nv1 + iv0
    iv2 = (iw-1)*(nv-nv1) + (iv0-nv1)
    call setdim
  enddo
  enddo

  print*, ' Allocating...'
  do iv=1, nv
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32
  enddo

  END subroutine get_8var_init

  SUBROUTINE get_8var

  vpt(:,:,:,:) = 0.  ;  vu(:,:,:,:) = 0.
  wpt(:,:,:,:) = 0.  ;  wu(:,:,:,:) = 0.

  REAL_IMAG:  DO iri=1, 2

  do iv_i=1, 4
    iw2 = iw
    ! assume iw = 1 or 2, i.e., (/'s','a'/) or (/'a','s'/)
    if (iv_i == 4)  iw2 = 3-iw
    print*, iri, trim(file_i(iv_i))
    var_i_name(iv_i) = trim(ivarhead(iri))//trim(ivarname(iv_i))//'_'//  &
                       trim(wavgrp(iw2))
    var5d(:,:,ny2/2+1:ny2,:,iv_i) =  &
       get_ivara4d(1,nk,1,nome,1,(ny2+1)/2,1,nz2)
    if (trim(wavgrp(iw2)) == 's') then
      var5d(:,:,1:ny2/2,:,iv_i) = var5d(:,:,ny2:(ny2+3)/2:-1,:,iv_i)
    else if (trim(wavgrp(iw2)) == 'a') then
      var5d(:,:,1:ny2/2,:,iv_i) = var5d(:,:,ny2:(ny2+3)/2:-1,:,iv_i)*(-1.)
    end if
  enddo
  vpt(:,:,:,:) = vpt(:,:,:,:) + var5d(:,:,:,:,3)*var5d(:,:,:,:,2)
  vu (:,:,:,:) = vu (:,:,:,:) + var5d(:,:,:,:,3)*var5d(:,:,:,:,1)
  wpt(:,:,:,:) = wpt(:,:,:,:) + var5d(:,:,:,:,4)*var5d(:,:,:,:,2)
  wu (:,:,:,:) = wu (:,:,:,:) + var5d(:,:,:,:,4)*var5d(:,:,:,:,1)

  ENDDO  REAL_IMAG

  do k=1, nz2
  do j=1, ny2
    vpt(:,:,j,k) = vpt(:,:,j,k)*t2pt0(j,k)
    wpt(:,:,j,k) = wpt(:,:,j,k)*t2pt0(j,k)
  enddo
  enddo
  do k=1, nz2
    wpt(:,:,:,k) = -wpt(:,:,:,k)/(p0a(k)*100.)*h_scale
    wu (:,:,:,k) = -wu (:,:,:,k)/(p0a(k)*100.)*h_scale
  enddo

  END subroutine get_8var

  SUBROUTINE setdim

  if (iv0 <= nv1) then
    set0(iv)%vname = trim(varname_epf(iv0))//'_'//trim(wavgrp(iw))
    set0(iv)%axis = (/'k_wn    ','ome_fr','lat  ','p '/) 
    set0(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set0(iv)%axis1(set0(iv)%nd(1)) )
    allocate( set0(iv)%axis2(set0(iv)%nd(2)) )
    allocate( set0(iv)%axis3(set0(iv)%nd(3)) )
    allocate( set0(iv)%axis4(set0(iv)%nd(4)) )
    set0(iv)%axis1 = kwn
    set0(iv)%axis2 = ome
    set0(iv)%axis3 = lat0
    set0(iv)%axis4 = p0a
    if (iw == 1) then
      set(iv)%vname = ''
      set(iv)%axis = set0(iv)%axis
      set(iv)%nd(:) = set0(iv)%nd(:)
      allocate( set(iv)%axis1(set(iv)%nd(1)) )
      allocate( set(iv)%axis2(set(iv)%nd(2)) )
      allocate( set(iv)%axis3(set(iv)%nd(3)) )
      allocate( set(iv)%axis4(set(iv)%nd(4)) )
      set(iv)%axis1 = set0(iv)%axis1
      set(iv)%axis2 = set0(iv)%axis2
      set(iv)%axis3 = set0(iv)%axis3
      set(iv)%axis4 = set0(iv)%axis4
    end if
  else
    set0b(iv2)%vname = trim(varname_epf(iv0))//'_'//trim(wavgrp(iw))
    set0b(iv2)%axis = (/'k_wn    ','ome_fr','lat  ','p '/) 
    set0b(iv2)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set0b(iv2)%axis1(set0b(iv2)%nd(1)) )
    allocate( set0b(iv2)%axis2(set0b(iv2)%nd(2)) )
    allocate( set0b(iv2)%axis3(set0b(iv2)%nd(3)) )
    allocate( set0b(iv2)%axis4(set0b(iv2)%nd(4)) )
    set0b(iv2)%axis1 = kwn
    set0b(iv2)%axis2 = ome
    set0b(iv2)%axis3 = lat0
    set0b(iv2)%axis4 = p0a
    if (iw == 1) then
      set(iv)%vname = ''
      set(iv)%axis = set0b(iv2)%axis
      set(iv)%nd(:) = set0b(iv2)%nd(:)
      allocate( set(iv)%axis1(set(iv)%nd(1)) )
      allocate( set(iv)%axis2(set(iv)%nd(2)) )
      allocate( set(iv)%axis3(set(iv)%nd(3)) )
      allocate( set(iv)%axis4(set(iv)%nd(4)) )
      set(iv)%axis1 = set0b(iv2)%axis1
      set(iv)%axis2 = set0b(iv2)%axis2
      set(iv)%axis3 = set0b(iv2)%axis3
      set(iv)%axis4 = set0b(iv2)%axis4
    end if
  end if

  END subroutine setdim

  SUBROUTINE finalize

  deallocate( lon, lat, p, lat0, p0a, dim4, t2pt, t2pt0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo
  do iv=1, nv1*nwavgrp
    deallocate( set0(iv)%axis1, set0(iv)%axis2, set0(iv)%axis3,          &
                set0(iv)%axis4 )
  enddo
  do iv=1, (nv-nv1)*nwavgrp
    deallocate( set0b(iv)%axis1, set0b(iv)%axis2, set0b(iv)%axis3,       &
                set0b(iv)%axis4 )
  enddo
  deallocate( set, set0, set0b )

  deallocate( var5d )
  deallocate( vpt, vu, wpt, wu )
  deallocate( um, ptm )

  END subroutine finalize


END program EPFSA1_RA


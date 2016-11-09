PROGRAM EPF_RA

  use epf
  use reanal
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 8, nv1 = 4, max_read_ome = 60  ! 62
  real,    parameter ::  pb_out = 500.
  character(len=8), dimension(2) ::  ivarhead = (/'fcr_','fci_'/)

  integer ::  opt_30d, k_largest, period_smallest, nmon_patch, k_sum
  real    ::  lat_smpl0
  character(len=128) ::  file_o2

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, REFDATE, OPT_30D
  namelist /PARAM/ K_LARGEST, PERIOD_SMALLEST, LAT_RNG, P_RNG,           &
                   NMON_PATCH, K_SUM, LAT_SMPL0
  namelist /FILEIO/ NT_F4, P_PREDEF, FILE_I_HEAD, FILE_I_FORM,           &
                    FILE_I_XXXX, VAR_I, VAR_I_NAME, FILE_I_HEAD2,        &
                    FILE_I_FORM2, FILE_I_XXXX2, VAR_I2, VAR_I_NAME2,     &
                    FILE_O, FILE_O2

  integer ::  iz, imon, ihour, i_time, n0r, n9r
  integer ::  nk, nome, nk0, nksum, nkcal, nn, nread, nread_ome
  integer ::  iy2(2), iz2(2), ny2, nz2, nta
  integer ::  iy2o(2), iz2o(2)
  integer ::  j,k,n, iri, tmpi
  real    ::  wgt
  character(len=32), dimension(8)  ::  ivarname
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:,:), allocatable ::  var5d
  real, dimension(:,:,:,:), allocatable ::  vpt, vu, wpt, wu
  real, dimension(:,:,:),   allocatable ::  var3d
  real, dimension(:,:,:),   allocatable ::  um, ptm
  real, dimension(:,:),     allocatable ::  t2pt0
  real, dimension(:),       allocatable ::  lat0, p0a, kwn, ome

  type(vset),  dimension(nv) ::  set
  type(vset0), dimension(nv) ::  set0

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
  call outnc_blank(trim(file_o),nv1,set0(1:nv1),  &
       'EP flux analysis in spectral domain', 4)
  call outnc_blank(trim(file_o2),nv-nv1,set0(nv1+1:nv),  &
       'EP flux analysis in spectral domain', 4)

  SPLIT_OME:  DO n=1, nn
  !---------------------------------------------------------------------

  print*, n, '/', nn

  print*, ' Reading...'
  if (n == 1) then
    n0r = 1
  else
    n0r = n9r + 1
  end if
  call get_8var

  print*, ' Calculating...'
  call epf_hydro_p_fc2(                                                  &
     nkcal,nread,ny2,nz2,lat0,p0a*100.,nt,um(:,:,1:nt),ptm(:,:,1:nt),    &
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
  call outnc_ss(trim(file_o),nv1,set(1:nv1),(/2,n0r,nread/))
  call outnc_ss(trim(file_o2),nv-nv1,set(nv1+1:nv),(/2,n0r,nread/))

  !---------------------------------------------------------------------
  ENDDO  SPLIT_OME

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
  call getdim(file_i(1),trim(ivarhead(1))//trim(var_i_name(1)))
  l_rev(1:3) = .False.

  nk = k_largest+1
  if (k_largest == -999)  nk = nx
  allocate( kwn(nk) )
  kwn(:) = lon(:)

  nk0 = k_sum+1
  if (k_sum == -999)  nk0 = nk
  nksum = nk - nk0
  nkcal = min(nk,nk0+1)

  nta = (nmon+2*nmon_patch)*30*nhour
  if ( opt_30d == 0 .and. nmon /= 1 ) then
    print*, 'not yet programmed (for opt_30d = 0)'  ;  STOP
  end if

  nome = nta/(period_smallest*nhour)*2+1
  if (period_smallest == -999)  nome = ny
  allocate( ome(nome) )
  ome(:) = lat(:)

  if ( any(lat0 /= p(ny2:1:-1)) ) then
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

  nread_ome = min(max_read_ome, nome)

  allocate( var5d(nk,nread_ome,ny2,nz2,3) )

  allocate( vpt(nkcal,nread_ome,ny2,nz2), vu(nkcal,nread_ome,ny2,nz2) )
  allocate( wpt(nkcal,nread_ome,ny2,nz2), wu(nkcal,nread_ome,ny2,nz2) )

  nd1a = NKCAL
  nd2a = NOME
  nd3a = NY2
  nd4a = NZ2

  print*, ' Allocating...'
  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nread_ome,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32
  enddo

  nn = nome/nread_ome + 1
  if (mod(nome,nread_ome) == 0)  nn = nome/nread_ome

  END subroutine get_8var_init

  SUBROUTINE get_8var

  nread = nread_ome
  n9r = n0r + nread_ome - 1
  if (n9r > nome) then
    n9r = nome
    nread = n9r - n0r + 1

    deallocate( var5d, vpt, vu, wpt, wu )
    do iv=1, nv
      deallocate( set(iv)%var_out )
    end do

    allocate( var5d(nk,nread,ny2,nz2,3) )
    allocate( vpt(nkcal,nread,ny2,nz2), vu(nkcal,nread,ny2,nz2) )
    allocate( wpt(nkcal,nread,ny2,nz2), wu(nkcal,nread,ny2,nz2) )

    do iv=1, nv
      allocate( set(iv)%var_out(nd1a,nread,nd3a,nd4a) )
      set(iv)%var_out(:,:,:,:) = 1.e32
    enddo
  end if

  vpt(:,:,:,:) = 0.  ;  vu(:,:,:,:) = 0.
  wpt(:,:,:,:) = 0.  ;  wu(:,:,:,:) = 0.

  REAL_IMAG:  DO iri=1, 2

  do iv_i=1, 3
    print*, iri, trim(file_i(iv_i))
    var_i_name(iv_i) = trim(ivarhead(iri))//trim(ivarname(iv_i))
    var5d(:,:,:,:,iv_i) = get_ivara4d(1,nk,n0r,nread,1,ny2,1,nz2)
  enddo

  vpt(1:nk0,:,:,:) = vpt(1:nk0,:,:,:) +  &
                     var5d(1:nk0,:,:,:,3)*var5d(1:nk0,:,:,:,2)
  vu (1:nk0,:,:,:) = vu (1:nk0,:,:,:) +  &
                     var5d(1:nk0,:,:,:,3)*var5d(1:nk0,:,:,:,1)
  if (nksum > 0) then
    vpt(nk0+1,:,:,:) = vpt(nk0+1,:,:,:) +  &
       sum(var5d(nk0+1:,:,:,:,3)*var5d(nk0+1:,:,:,:,2), dim=1)
    vu (nk0+1,:,:,:) = vu (nk0+1,:,:,:) +  &
       sum(var5d(nk0+1:,:,:,:,3)*var5d(nk0+1:,:,:,:,1), dim=1)
  end if    

  iv_i = 4
  print*, iri, trim(file_i(iv_i))
  var_i_name(iv_i) = trim(ivarhead(iri))//trim(ivarname(iv_i))
  var5d(:,:,:,:,3) = get_ivara4d(1,nk,n0r,nread,1,ny2,1,nz2)

  wpt(1:nk0,:,:,:) = wpt(1:nk0,:,:,:) +  &
                     var5d(1:nk0,:,:,:,3)*var5d(1:nk0,:,:,:,2)
  wu (1:nk0,:,:,:) = wu (1:nk0,:,:,:) +  &
                     var5d(1:nk0,:,:,:,3)*var5d(1:nk0,:,:,:,1)
  if (nksum > 0) then
    wpt(nk0+1,:,:,:) = wpt(nk0+1,:,:,:) +  &
       sum(var5d(nk0+1:,:,:,:,3)*var5d(nk0+1:,:,:,:,2), dim=1)
    wu (nk0+1,:,:,:) = wu (nk0+1,:,:,:) +  &
       sum(var5d(nk0+1:,:,:,:,3)*var5d(nk0+1:,:,:,:,1), dim=1)
  end if    

  ENDDO  REAL_IMAG

  ! corrections for pt and w

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

  set0(iv)%vname = trim(ovarname(iv))
  set0(iv)%axis = (/'k_wn    ','ome_fr','lat  ','p '/) 
  set0(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set0(iv)%axis1(set0(iv)%nd(1)) )
  allocate( set0(iv)%axis2(set0(iv)%nd(2)) )
  allocate( set0(iv)%axis3(set0(iv)%nd(3)) )
  allocate( set0(iv)%axis4(set0(iv)%nd(4)) )
  set0(iv)%axis1(1:nkcal-1) = kwn(1:nkcal-1)
  set0(iv)%axis1(nkcal) = kwn(nk)
  set0(iv)%axis2 = ome
  set0(iv)%axis3 = lat0
  set0(iv)%axis4 = p0a

  set(iv)%vname = set0(iv)%vname
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

  END subroutine setdim

  SUBROUTINE finalize

  deallocate( lon, lat, p, lat0, p0a, dim4, t2pt, t2pt0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
    deallocate( set0(iv)%axis1, set0(iv)%axis2, set0(iv)%axis3,          &
                set0(iv)%axis4 )
  enddo

  deallocate( var5d )
  deallocate( vpt, vu, wpt, wu )
  deallocate( um, ptm )

  END subroutine finalize


END program EPF_RA


PROGRAM TEM_REANALYSIS
! attributes (scale facter, add_offset, "_FillValue" or "missing_value") must be considered.
! m_avg must be completed.

  use tem
  use reanal
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 12
  real,    parameter ::  pb_out = 500.

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH
  namelist /PARAM/ OPT_AVRG
  namelist /FILEIO/ NT_F4, MISSV, P_PREDEF, FILE_I_HEAD, FILE_I_FORM,    &
                    FILE_I_XXXX, VAR_I, VAR_I_NAME, FILE_O

  integer ::  iz, imon, ihour, i_time, i_time_last
  integer ::  i,j,k,n
  integer ::  tag_exit
  real    ::  ntavg
  character(len=32), dimension(nv) ::  ovarname

  real,    dimension(:,:,:,:), allocatable ::  var4d
  real,    dimension(:,:,:),   allocatable ::  var3d
  real,    dimension(:,:,:),   allocatable ::  u, v, w, pt
  logical, dimension(:,:),     allocatable ::  l_tag2d

  type(vset), dimension(nv) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  if ( nt_f4(4) /= 1 .or. ( nt_f4(3) /= 1 .and. nt_f4(3) /= 12 ) ) then
    print*, 'NT_I(3:4) should be (/1,1/) or (/12,1/).'  ;  STOP
  end if

  nmon = mm(2)  ;  nhour = hh(2)

  year = yyyy
  mon  = mm(1)

  tag_exit = 0

  i_time_last = 0
  i_time      = 0

  L_MON:  DO imon=1, nmon+1
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

  ndate = enddate(mon)
  if ( mon == 2 .and. mod(year,4) == 0 )  ndate = 29

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

  if (opt_avrg /= 0)  i_time = i_time + 1

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------
                                                if (tag_exit == 1)  EXIT

  if (opt_avrg == 0)  i_time = i_time + 1

  ex0 = .TRUE.
  do iv_i=1, 4
    file_i(iv_i) = get_ifilename()
    inquire(file=trim(file_i(iv_i)), exist=ex1)
    if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
    ex0 = ( ex0 .and. ex1 )
  enddo
  print*, trim(file_i(4))

  if (.not. ex0)  STOP
  if (.not. allocated(lon)) then
    call getdim(file_i(4),var_i_name(4))  ! temperature
    ovarname(1:11) = varname_tem_p(1:11)
    ovarname(12) = 'u_tend'
    if (opt_avrg == 0) then
      allocate( var4d(ny,nz,nmon*31*nhour,nv) )
    else
      allocate( var4d(ny,nz,nmon*31,nv) )
    end if
    var4d(:,:,:,:) = 0.
    allocate( var3d(ny,nz,nv) )
  end if

  ! get variable
  allocate( u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), pt(nx,ny,nz) )
  allocate( l_tag2d(ny,nz) )
  call get_4var
  if (missv /= 1.0) then
do j=1, ny
do i=1, nx
do k=nz, 1, -1
if (u(i,j,k) == missv) then
u(i,j,1:k) = u(i,j,k+1)
v(i,j,1:k) = v(i,j,k+1)
w(i,j,1:k) = w(i,j,k+1)
do n=k, 1, -1
pt(i,j,n) = pt(i,j,n+1)*exp( -log(pt(i,j,n+2)/pt(i,j,n+1))*log(p(n)/p(n+1))/log(p(n+1)/p(n+2)) )
enddo
exit
endif
enddo
enddo
enddo
  end if

  ! calculate zonal mean
  call tem_hydro_p(                                                      &
      nx,ny,nz,lat,p*100.,u,v,w,pt,h_scale,1.e32,                        &
      var3d(:,:,1 ),var3d(:,:,2 ),var3d(:,:,3 ),var3d(:,:,4 ),           &
      var3d(:,:,5 ),var3d(:,:,6 ),var3d(:,:,7 ),var3d(:,:,8 ),           &
      var3d(:,:,9 ),var3d(:,:,10),var3d(:,:,11) )

  ntavg = float(nhour)
  if (opt_avrg == 0)  ntavg = 1.

  if (i_time == i_time_last) then
    var4d(:,:,i_time,:11) = var4d(:,:,i_time,:11) +                      &
                            var3d(:,:,:11)/ntavg
  else
    var3d(:,:,12) = sum(u, dim=1)/float(nx)
    if (opt_avrg == 0)  var3d(:,:,12) = var3d(:,:,12)*float(nhour)

    if (i_time /= 1) then
      var4d(:,:,i_time_last,:11) = var4d(:,:,i_time_last,:11) +          &
                                   var3d(:,:,:11)/(ntavg*2.)
      var4d(:,:,i_time_last,12) = var4d(:,:,i_time_last,12) +            &
                                  var3d(:,:,12)
    end if
    if (imon /= nmon+1) then
      var4d(:,:,i_time,:11) = var4d(:,:,i_time,:11) +                    &
                              var3d(:,:,:11)/(ntavg*2.)
      var4d(:,:,i_time,12) = var4d(:,:,i_time,12) - var3d(:,:,12)
    end if
  end if

  deallocate( l_tag2d )
  deallocate( u, v, w, pt )

  hour = hour + 24/nhour

  i_time_last = i_time

  if (imon == nmon+1)  tag_exit = 1
  !---------------------------------------------------------------------
  ENDDO  L_HOUR

  !---------------------------------------------------------------------
  ENDDO  L_DATE

  mon = mon + 1
  if (mon == 13) then
    year = year + 1
    mon = 1
  end if

  !---------------------------------------------------------------------
  ENDDO  L_MON

  nt = i_time - 1

  ! pt2t
  ! missing

  nd1a = NY
  nd2a = NZ
  nd3a = NT
  nd4a = 1

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,1) = var4d(:,:,1:nd3a,iv)
  enddo

  if (pb_out /= 1000.) then
    iz = nd2a + 1
    do k=2, nd2a
      if (p(k) < pb_out) then  ;  iz = k - 1  ;  EXIT  ;  end if
    enddo
    do iv=1, nv
      set(iv)%var_out(:,:iz-1,:,1) = 1.e32
    enddo
  end if


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'transformed Eulerian mean eqn.')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE get_4var

  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz   ! for get_ivar3d

  ! read 4 var.s
  iv_i = 1  ;  u (:,:,:) = get_ivar3d()
  iv_i = 2  ;  v (:,:,:) = get_ivar3d()
  iv_i = 3  ;  w (:,:,:) = get_ivar3d()
  iv_i = 4  ;  pt(:,:,:) = get_ivar3d()

  ! omega to w(zp)  ;  t to pt
  do k=1, nz
  do j=1, ny
    l_tag2d(j,k) = .TRUE.
    do i=1, nx
      if (pt(i,j,k) == missv) then
        l_tag2d(j,k) = .FALSE.
      else
        w (i,j,k) = -w(i,j,k)/(p(k)*100.)*h_scale
        pt(i,j,k) = pt(i,j,k)*t2pt(k)
      end if
    enddo
  enddo
  enddo

  END subroutine get_4var

  SUBROUTINE setdim

  if ( .not. allocated(t) ) then
    allocate( t(nt) )
    do n=1, nt
      t(n) = float(n-1)+0.5
    enddo
    if (opt_avrg == 0)  t(:) = t(:)/nhour
  end if

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lat ','p   ','time','    '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat
  set(iv)%axis2 = p
  set(iv)%axis3 = t
  set(iv)%axis4 = -999.
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var4d, var3d )
  deallocate( lon, lat, p, t, t2pt )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program TEM_REANALYSIS


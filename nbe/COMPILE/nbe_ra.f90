PROGRAM NBE_REANALYSIS
! attributes (scale facter, add_offset, "_FillValue" or "missing_value") must be considered.

  use nbe
  use reanal
  use netio

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 5

  integer             ::  gp_discont_lon
  real, dimension(20) ::  plev

  namelist /ANALCASE/ EXPNAME, YYYY, MM, HH, PLEV
  namelist /PARAM/ GP_DISCONT_LON
  namelist /FILEIO/ NT_F4, MISSV, FILE_I_HEAD, FILE_I_FORM, FILE_I_XXXX, &
                    VAR_I, VAR_I_NAME, FILE_O

  integer ::  imon, ihour, i_time, np0, nx_tmp
  integer ::  i,j,k,n
  character(len=32), dimension(nv) ::  ovarname

  real,    dimension(:,:,:,:,:), allocatable ::  var5d
  real,    dimension(:,:,:),     allocatable ::  var3d
  real,    dimension(:,:),       allocatable ::  u, v, gp, omega, us, vs
  real,    dimension(:,:),       allocatable ::  gp6
  real,    dimension(:),         allocatable ::  tmp_gp, tmp_x, worka
  complex, dimension(:),         allocatable ::  fc
  integer, dimension(:),         allocatable ::  iz

  real ::  tmp

  type(vset), dimension(nv) ::  set

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  if ( nt_f4(4) /= 1 .or. ( nt_f4(3) /= 1 .and. nt_f4(3) /= 12 ) ) then
    print*, 'NT_I(3:4) should be (/1,1/) or (/12,1/).'  ;  STOP
  end if

  do i=1, 20
    if (plev(i) <= 0.) then
      np0 = i - 1  ;  EXIT
    end if
  enddo

  year = yyyy
  mon  = mm(1)

  call initialize

  i_time = 0

  L_MON:  DO imon=1, nmon+1
  !---------------------------------------------------------------------

  ndate = get_ndate()

  L_DATE:  DO date=1, ndate
  !---------------------------------------------------------------------

  hour = hh(1)

  L_HOUR:  DO ihour=1, nhour
  !---------------------------------------------------------------------

  i_time = i_time + 1

  do k=1, np0

    ! get variables
    call get_6var

    ! calculate NBE
    if (gp_discont_lon /= 0)  l_gp_discont_lon = .True.
    call nbe_p(nx,ny,lat,u,v,gp,omega,us,vs,1.e32,                       &
               var3d(:,:,1),var3d(:,:,2),var3d(:,:,3),var3d(:,:,4),      &
               var3d(:,:,5) )

    var3d(:,:,2) = var3d(:,:,2) + var3d(:,:,4)  ! NBE major + Curv.
    var3d(:,1:2,2) = 1.e32  ;  var3d(:,ny-1:ny,2) = 1.e32
    var5d(:,:,k,i_time,:) = var3d(:,:,:)

  enddo

  hour = hour + 24/nhour

  if ( imon == nmon .and. date == ndate .and. ihour == nhour )           &
      EXIT L_MON

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

  nt = i_time

  deallocate( fc )
  deallocate( gp6, tmp_gp, tmp_x, worka )
  deallocate( omega, us, vs )
  deallocate( u, v, gp )


  nd1a = NX
  nd2a = NY
  nd3a = NP0
  nd4a = NT

  do iv=1, nv
    call setdim
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = 1.e32

    set(iv)%var_out(:,:,:,:) = var5d(:,:,:,1:nd4a,iv)
  enddo


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'transformed Eulerian mean eqn.')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  nmon = mm(2)  ;  nhour = hh(2)

  date = 1  ;  hour = hh(1)  ! for get_ifilename
  ndate = get_ndate()

  iv_i = 3
  file_i(iv_i) = get_ifilename()
  inquire(file=trim(file_i(iv_i)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(iv_i)),' not found.'  ;  STOP
  end if
  call getdim(file_i(iv_i),var_i_name(iv_i))

  allocate( iz(np0) )
  do k=1, np0
    iz(k) = minloc( abs(log(p(:)/plev(k))), 1 )
  enddo
  print*, 'pressure levels: ', p(iz(:))
  if ( l_rev(3) )  iz(:) = nz + 1 - iz(:)

  ovarname(1:nv) = varname_nbe_p(1:nv)

  allocate( var5d(nx,ny,np0,nmon*31*nhour,nv) )
  var5d(:,:,:,:,:) = 0.
  allocate( var3d(nx,ny,nv) )

  allocate( u(nx,ny), v(nx,ny), gp(nx,ny) )
  allocate( omega(nx,ny), us(nx,ny), vs(nx,ny) )
  allocate( gp6(nx,6), tmp_gp(nx*3), tmp_x(nx*3), worka(0:nx*3) )
  allocate( fc(nx) )

  END subroutine initialize

  SUBROUTINE get_6var

  if (k == 1) then
    ex0 = .TRUE.
    do iv_i=1, 4
      file_i(iv_i) = get_ifilename()
      inquire(file=trim(file_i(iv_i)), exist=ex1)
      if ( .not. ex1 )  print*, '    ',trim(file_i(iv_i)),' not found.'
      ex0 = ( ex0 .and. ex1 )
    enddo
    if (.not. ex0)  STOP
    print*, trim(file_i(3))
  end if

  iv_i = 1  ;  u (:,:) = reshape(get_ivara3d(1,nx,1,ny,iz(k),1),(/nx,ny/))
  iv_i = 2  ;  v (:,:) = reshape(get_ivara3d(1,nx,1,ny,iz(k),1),(/nx,ny/))
  iv_i = 3  ;  gp(:,:) = reshape(get_ivara3d(1,nx,1,ny,iz(k),1),(/nx,ny/))

  iv_i = 4  ;  omega(:,:) = reshape(get_ivara3d(1,nx,1,ny,iz(k),1),      &
                                    (/nx,ny/))
  iv_i = 1  ;  us(:,:) = reshape( get_ivara3d(1,nx,1,ny,iz(k)+1,1) -     &
                                  get_ivara3d(1,nx,1,ny,iz(k)-1,1),      &
                                  (/nx,ny/) )
  iv_i = 2  ;  vs(:,:) = reshape( get_ivara3d(1,nx,1,ny,iz(k)+1,1) -     &
                                  get_ivara3d(1,nx,1,ny,iz(k)-1,1),      &
                                  (/nx,ny/) )
  tmp = (p(iz(k)+1) - p(iz(k)-1))*1.e2
  if ( l_rev(3) )  tmp = tmp*(-1.)
  us(:,:) = us(:,:)/tmp
  vs(:,:) = vs(:,:)/tmp

  ! gph to gp
  tmp = h_scale*log(p0/(plev(k)*1.e2))
  if ( abs(gp(nx/2,ny/2)/g-tmp) > abs(gp(nx/2,ny/2)-tmp) ) then
    print*, 'converting GPH to GP'
    gp(:,:) = gp(:,:)*g
  else
    print*, 'This RA may contain GP (cf. GPH). Check this.'
  end if

  END subroutine get_6var

  SUBROUTINE setdim

  if ( .not. allocated(t) ) then
    allocate( t(nt) )
    do n=1, nt
      t(n) = float(n-1)
    enddo
    t(:) = t(:)/nhour
  end if

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lon  ','lat  ','p ','time'/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lon
  set(iv)%axis2 = lat
  set(iv)%axis3 = plev(1:np0)
  set(iv)%axis4 = t
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( var5d, var3d, iz )
  deallocate( lon, lat, p, t )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program NBE_REANALYSIS


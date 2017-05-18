MODULE reanal

  implicit none

  real, private, parameter ::  kappa = 287./1004.

! variables in namelists

  integer ::  yyyy, mm(2), hh(2), refdate(3), opt_avrg, nt_f4(4), nres
  real    ::  missv, lat_rng(2), p_rng(2)

  character(len=32)  ::  expname, var_name(99)
  character(len=32)  ::  var_i(99), var_i_name(99), file_i_xxxx(10)
  character(len=32)  ::  var_i2(99), var_i_name2(99), file_i_xxxx2(10)
  character(len=128) ::  file_i_head, file_i_form(10), file_o
  character(len=128) ::  file_i_head2, file_i_form2(10)

! parameters for namelists

  character(len=8), parameter ::  xblank = 'XXXX'

! common variables

  integer ::  year, mon, date, hour
  integer ::  nmon, ndate, nhour
  logical ::  l_rev(3), ex0, ex1, l_gpheight
  character(len=256) ::  file_i(20), f_namelist

  integer ::  nt, nx, ny, nz, nl, nd1a, nd2a, nd3a, nd4a, nx_neg
  integer ::  nx_i, ny_i, nz_i, nt_i
  integer ::  iv_i, iv, it_i(20)
  integer, private ::  i,j,k,n
  data  it_i /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

  integer, dimension(12), parameter ::                                   &
     enddate = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  real, dimension(:), allocatable ::  lon, lat, p, t2pt, dim4
  real, dimension(:), allocatable ::  t

  character(len=32)  ::  var_i0(99), var_i_name0(99), file_i_xxxx0(10)
  character(len=128) ::  file_i_head0, file_i_form0(10)

! output variable type

!  type ::  vset
!    real, dimension(:,:,:,:), allocatable ::  var_out
!    real, dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
!    integer                               ::  nd(4)
!    character(len=32)                     ::  vname, axis(4)
!  end type vset


  contains


!::::  FUNCTIONS  ::::::::::::::::::::::::::::::::::::::::::::::::::::::

FUNCTION get_ifilename()

  character(len=256) ::  get_ifilename

  integer           ::  nform, lblank, iblank, ic
  character(len=32) ::  vblank

  lblank = len_trim(xblank)
  write(get_ifilename,'(a)')  trim(file_i_head)

  do i=1, 10
    if (file_i_form(i) == '-999')  nform = i - 1
  enddo

  iblank = 0
  do i=1, nform
    write(get_ifilename,'(a)')  trim(get_ifilename)//'/'
    ic = 1
    do while ( ic <= len_trim(file_i_form(i)) )
      if ( file_i_form(i)(ic:ic+lblank-1) == trim(xblank) ) then
        iblank = iblank + 1
        vblank = str2var(file_i_xxxx(iblank))
        write(get_ifilename,'(a)')  trim(get_ifilename)//trim(vblank)
        ic = ic + lblank
      else
        write(get_ifilename,'(a)')  trim(get_ifilename)//              &
                                    file_i_form(i)(ic:ic)
        ic = ic + 1
      end if
    enddo
  enddo

  RETURN

END function get_ifilename


FUNCTION str2var(str)

  character(len=32) ::  str2var

  character(len=*), intent(in)  ::  str

  SELECT case ( trim(str) )

  case ( 'EXPNAME', 'expname' )
    write(str2var,'(a)')  trim(expname)
  case ( 'YYYY', 'yyyy' )
    write(str2var,'(i4.4)')  year
  case ( 'MM', 'mm' )
    write(str2var,'(i2.2)')  mon
  case ( 'DD', 'dd' )
    write(str2var,'(i2.2)')  date
  case ( 'HH', 'hh' )
    write(str2var,'(i2.2)')  hour
  case ( 'VAR_I', 'var_i' )
    write(str2var,'(a)')  trim(var_i(iv_i))
  case ( 'VAR_I_N', 'var_i_n' )
    write(str2var,'(a)')  trim(var_i_name(iv_i))
  case default
    print*, 'Check FILE_I_XXXX in the shell script and function '//    &
            'str2var !!!'
    STOP

  END select

  RETURN

END function str2var


FUNCTION get_ndate()

  integer ::  get_ndate

  get_ndate = enddate(mon)
  if ( mon == 2 .and. mod(year,4) == 0 )  get_ndate = 29

END function get_ndate


FUNCTION get_dayfromref_30d(y,m,d,h)

  integer, intent(in) ::  y, m, d, h
  
  integer ::  get_dayfromref_30d
  
  get_dayfromref_30d = (y-refdate(1))*360.+(m-refdate(2))*30.+ &
                       (d-refdate(3))+h/24.
                   
END function get_dayfromref_30d


FUNCTION get_ivar3d()

  real, dimension(nx_i,ny_i,nz_i) ::  get_ivar3d

  get_ivar3d(:,:,:) = get_ivara3d(1,nx_i,1,ny_i,1,nz_i)

END function get_ivar3d


FUNCTION get_ivara3d(ix_i0,nx_i0,iy_i0,ny_i0,iz_i0,nz_i0)

  integer, intent(in) ::  ix_i0, nx_i0, iy_i0, ny_i0, iz_i0, nz_i0

  real, dimension(nx_i0,ny_i0,nz_i0) ::  get_ivara3d

  include 'netcdf.inc'

  integer ::  st, ncid, ncvarid, ncdimid(4), ncvtype, ncvstart(4),       &
              ncvcount(4), ncdimvarid
  integer ::  it_i0
  character(len=32) ::  dimname

  real,             dimension(:,:,:), allocatable ::  ncvar_r
  integer,          dimension(:,:,:), allocatable ::  ncvar_i
  integer*2,        dimension(:,:,:), allocatable ::  ncvar_s
  double precision, dimension(:,:,:), allocatable ::  ncvar_d

  it_i0 = 0
  if ( nt_f4(3) == 12 .and. mon /= 1 ) then
    it_i0 = sum(enddate(1:mon-1))
    if ( mon > 2 .and. mod(year,4) == 0 )  it_i0 = it_i0 + 1
  end if
  if (nt_f4(2) /= 1)  it_i0 = it_i0 + date - 1
  it_i0 = it_i0*nt_f4(1)
  if (nt_f4(1) /= 1)  it_i0 = it_i0 + hour/(24/nt_f4(1))
  it_i0 = it_i0 + 1

  ncvstart = (/ix_i0,iy_i0,iz_i0,it_i0/)
  ncvcount = (/nx_i0,ny_i0,nz_i0,1    /)

  st = nf_open(trim(file_i(iv_i)), NF_NOWRITE, ncid)

  st = nf_inq_varid(ncid,trim(var_i_name(iv_i)), ncvarid)

  allocate( ncvar_r(nx_i0,ny_i0,nz_i0) )

  st = nf_inq_vartype(ncid,ncvarid, ncvtype)

  SELECT case ( ncvtype )

  case ( NF_FLOAT )
    st = nf_get_vara_real(ncid,ncvarid,ncvstart,ncvcount, ncvar_r)
  case ( NF_DOUBLE )
    allocate( ncvar_d(nx_i0,ny_i0,nz_i0) )
    st = nf_get_vara_double(ncid,ncvarid,ncvstart,ncvcount, ncvar_d)
    ncvar_r(:,:,:) = real(ncvar_d(:,:,:))
    deallocate( ncvar_d )
  case ( NF_INT )
    allocate( ncvar_i(nx_i0,ny_i0,nz_i0) )
    st = nf_get_vara_int(ncid,ncvarid,ncvstart,ncvcount, ncvar_i)
    ncvar_r(:,:,:) = real(ncvar_i(:,:,:))
    deallocate( ncvar_i )
  case ( NF_SHORT )
    allocate( ncvar_s(nx_i0,ny_i0,nz_i0) )
    st = nf_get_vara_int2(ncid,ncvarid,ncvstart,ncvcount, ncvar_s)
    ncvar_r(:,:,:) = real(ncvar_s(:,:,:))
    deallocate( ncvar_s )
  case default
    print*, 'Check the type of variables in NETCDF file'
    STOP

  END select

  st = nf_close(ncid)

  get_ivara3d(:,:,:) = ncvar_r(:,:,:)

  if ( l_rev(1) ) then
    ncvar_r(:,:,:) = get_ivara3d(:,:,:)
    get_ivara3d(1:nx_i0-nx_neg,:,:) = ncvar_r(nx_neg+1:nx_i0,:,:)
    get_ivara3d(nx_i0-nx_neg+1:nx_i0,:,:) = ncvar_r(1:nx_neg,:,:)
  end if
  if ( l_rev(2) ) then
    ncvar_r(:,:,:) = get_ivara3d(:,:,:)
    do j=1, ny_i0
      get_ivara3d(:,j,:) = ncvar_r(:,ny_i0+1-j,:)
    enddo
  end if
  if ( l_rev(3) ) then
    ncvar_r(:,:,:) = get_ivara3d(:,:,:)
    do k=1, nz_i0
      get_ivara3d(:,:,k) = ncvar_r(:,:,nz_i0+1-k)
    enddo
  end if

  deallocate( ncvar_r )

  it_i(iv_i) = it_i0

  RETURN

END function get_ivara3d


FUNCTION get_ivara4d(ix_i0,nx_i0,iy_i0,ny_i0,iz_i0,nz_i0,it_i0,nt_i0)

  integer, intent(in) ::  ix_i0, nx_i0, iy_i0, ny_i0, iz_i0, nz_i0,      &
                          it_i0, nt_i0

  real, dimension(nx_i0,ny_i0,nz_i0,nt_i0) ::  get_ivara4d

  include 'netcdf.inc'

  integer ::  st, ncid, ncvarid, ncdimid(4), ncvtype, ncvstart(4),       &
              ncvcount(4), ncdimvarid
  character(len=32) ::  dimname

  real,             dimension(:,:,:,:), allocatable ::  ncvar_r
  integer,          dimension(:,:,:,:), allocatable ::  ncvar_i
  integer*2,        dimension(:,:,:,:), allocatable ::  ncvar_s
  double precision, dimension(:,:,:,:), allocatable ::  ncvar_d

  ncvstart = (/ix_i0,iy_i0,iz_i0,it_i0/)
  ncvcount = (/nx_i0,ny_i0,nz_i0,nt_i0/)

  st = nf_open(trim(file_i(iv_i)), NF_NOWRITE, ncid)

  st = nf_inq_varid(ncid,trim(var_i_name(iv_i)), ncvarid)

  allocate( ncvar_r(nx_i0,ny_i0,nz_i0,nt_i0) )

  st = nf_inq_vartype(ncid,ncvarid, ncvtype)

  SELECT case ( ncvtype )

  case ( NF_FLOAT )
    st = nf_get_vara_real(ncid,ncvarid,ncvstart,ncvcount, ncvar_r)
  case ( NF_DOUBLE )
    allocate( ncvar_d(nx_i0,ny_i0,nz_i0,nt_i0) )
    st = nf_get_vara_double(ncid,ncvarid,ncvstart,ncvcount, ncvar_d)
    ncvar_r(:,:,:,:) = real(ncvar_d(:,:,:,:))
    deallocate( ncvar_d )
  case ( NF_INT )
    allocate( ncvar_i(nx_i0,ny_i0,nz_i0,nt_i0) )
    st = nf_get_vara_int(ncid,ncvarid,ncvstart,ncvcount, ncvar_i)
    ncvar_r(:,:,:,:) = real(ncvar_i(:,:,:,:))
    deallocate( ncvar_i )
  case ( NF_SHORT )
    allocate( ncvar_s(nx_i0,ny_i0,nz_i0,nt_i0) )
    st = nf_get_vara_int2(ncid,ncvarid,ncvstart,ncvcount, ncvar_s)
    ncvar_r(:,:,:,:) = real(ncvar_s(:,:,:,:))
    deallocate( ncvar_s )
  case default
    print*, 'Check the type of variables in NETCDF file'
    STOP

  END select

  st = nf_close(ncid)

  get_ivara4d(:,:,:,:) = ncvar_r(:,:,:,:)

  if ( l_rev(1) ) then
    ncvar_r(:,:,:,:) = get_ivara4d(:,:,:,:)
    get_ivara4d(1:nx_i0-nx_neg,:,:,:) = ncvar_r(nx_neg+1:nx_i0,:,:,:)
    get_ivara4d(nx_i0-nx_neg+1:nx_i0,:,:,:) = ncvar_r(1:nx_neg,:,:,:)
  end if
  if ( l_rev(2) ) then
    ncvar_r(:,:,:,:) = get_ivara4d(:,:,:,:)
    do j=1, ny_i0
      get_ivara4d(:,j,:,:) = ncvar_r(:,ny_i0+1-j,:,:)
    enddo
  end if
  if ( l_rev(3) ) then
    ncvar_r(:,:,:,:) = get_ivara4d(:,:,:,:)
    do k=1, nz_i0
      get_ivara4d(:,:,k,:) = ncvar_r(:,:,nz_i0+1-k,:)
    enddo
  end if

  deallocate( ncvar_r )

  RETURN

END function get_ivara4d

!::::  SUBROUTINES  ::::::::::::::::::::::::::::::::::::::::::::::::::::

SUBROUTINE getdim(file_i0,var_i_name0)

  include 'netcdf.inc'

  character(len=*) ::  file_i0, var_i_name0

  integer ::  st, ncid, ncvarid, ncdimid(4), ncvtype(4), ncdimvarid
  integer ::  varndim, dimlen(4), nti0, idm, ii0
  character(len=32) ::  dimname

  real,             dimension(:,:), allocatable ::  axis_r
  integer,          dimension(:),   allocatable ::  axis_i
  integer*2,        dimension(:),   allocatable ::  axis_s
  double precision, dimension(:),   allocatable ::  axis_d

  st = nf_open(trim(file_i0), NF_NOWRITE, ncid)

  st = nf_inq_varid(ncid,trim(var_i_name0), ncvarid)

  st = nf_inq_varndims(ncid,ncvarid, varndim)
  if (varndim /= 4)  print*, 'Warning: No. of dimensions in input /= 4.'

  st = nf_inq_vardimid(ncid,ncvarid, ncdimid)

  do idm=1, 4
    st = nf_inq_dimlen(ncid,ncdimid(idm), dimlen(idm))
  enddo
  nx = dimlen(1)  ;  ny = dimlen(2)  ;  nz = dimlen(3)
  nl = dimlen(4)

  nti0 = nt_f4(4)*nt_f4(3)*nt_f4(1)
  if (nt_f4(2) /= 1)  nti0 = nti0*ndate
!! applicable only when the 4th dim. is time.
!  if (nti0 /= dimlen(4)) then
!    print*, 'Check NT_I in the shell script !!!'
!    STOP
!  end if

  allocate( lon(nx), lat(ny), p(nz), dim4(nl) )

  allocate( axis_r(maxval(dimlen(1:4)),4) )

  do idm=1, 4

    st = nf_inq_dimname(ncid,ncdimid(idm), dimname)

    st = nf_inq_varid(ncid,trim(dimname), ncdimvarid)

    st = nf_inq_vartype(ncid,ncdimvarid, ncvtype(idm))

    SELECT case ( ncvtype(idm) )

    case ( NF_FLOAT )
      st = nf_get_var_real(ncid,ncdimvarid, axis_r(1:dimlen(idm),idm))
    case ( NF_DOUBLE )
      allocate( axis_d(dimlen(idm)) )
      st = nf_get_var_double(ncid,ncdimvarid, axis_d)
      axis_r(1:dimlen(idm),idm) = real(axis_d(:))
      deallocate( axis_d )
    case ( NF_INT )
      allocate( axis_i(dimlen(idm)) )
      st = nf_get_var_int(ncid,ncdimvarid, axis_i)
      axis_r(1:dimlen(idm),idm) = real(axis_i(:))
      deallocate( axis_i )
    case ( NF_SHORT )
      allocate( axis_s(dimlen(idm)) )
      st = nf_get_var_int2(ncid,ncdimvarid, axis_s)
      axis_r(1:dimlen(idm),idm) = real(axis_s(:))
      deallocate( axis_s )
    case default
      print*, 'Check the type of variables in NETCDF file'
      STOP

    END select

  enddo

  st = nf_close(ncid)

  l_rev(:) = .False.
  if (axis_r(1,1) < 0.          )  l_rev(1) = .True.
  if (axis_r(1,2) > axis_r(ny,2))  l_rev(2) = .True.
  if (axis_r(1,3) < axis_r(nz,3))  l_rev(3) = .True.

  lon(:) = axis_r(1:nx,1)
  lat(:) = axis_r(1:ny,2)
  p  (:) = axis_r(1:nz,3)
  if ( l_rev(1) ) then
    do i=1, nx
      if (lon(i) >= 0.) then
        ii0 = i - 1
        EXIT
      end if
    enddo
    nx_neg = ii0
    lon(1:nx-nx_neg) = axis_r(nx_neg+1:nx,1)
    lon(nx-nx_neg+1:nx) = axis_r(1:nx_neg,1) + 360.
  end if
  if ( l_rev(2) ) then
    do j=1, ny
      lat(j) = axis_r(ny+1-j,2)
    enddo
  end if
  if ( l_rev(3) ) then
    do k=1, nz
      p(k) = axis_r(nz+1-k,3)
    enddo
  end if

  dim4(:) = axis_r(1:nl,4)

  deallocate( axis_r )

  allocate( t2pt(nz) )
  t2pt(:) = (1000./p(:))**kappa

  RETURN

END subroutine getdim


SUBROUTINE get_iouter(arr,arr_rng, ind2)

  real,    dimension(:), intent(in)  ::  arr
  real,    dimension(2), intent(in)  ::  arr_rng
  integer, dimension(2), intent(out) ::  ind2
  
  integer ::  narr

  narr = size(arr)

  ind2 = (/1,narr/)

  if (abs(arr_rng(1)) /= 999.) then
    do i=2, narr
      if (arr(i) > arr_rng(1)) then
        ind2(1) = i - 1  ;  EXIT
      end if
    enddo
  end if
  if (abs(arr_rng(2)) /= 999.) then
    do i=1, narr
      if (arr(i) >= arr_rng(2)) then
        ind2(2) = i  ;  EXIT
      end if
    enddo
  end if

END subroutine get_iouter


SUBROUTINE switch_para_in

  file_i_head0    = file_i_head
  file_i_form0(:) = file_i_form(:)
  file_i_xxxx0(:) = file_i_xxxx(:)
  var_i0      (:) = var_i      (:)
  var_i_name0 (:) = var_i_name (:)

  file_i_head    = file_i_head2
  file_i_form(:) = file_i_form2(:)
  file_i_xxxx(:) = file_i_xxxx2(:)
  var_i      (:) = var_i2      (:)
  var_i_name (:) = var_i_name2 (:)

  file_i_head2    = file_i_head0
  file_i_form2(:) = file_i_form0(:)
  file_i_xxxx2(:) = file_i_xxxx0(:)
  var_i2      (:) = var_i0      (:)
  var_i_name2 (:) = var_i_name0 (:)

END subroutine switch_para_in


END module reanal


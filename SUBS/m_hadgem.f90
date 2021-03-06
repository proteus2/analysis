MODULE hadgem

  implicit none

! variables in namelists

  integer ::  yyyy, mm(2), dd(2), hh(3), refdate(3), opt_30d,            &
              nres, day1, nday_i
  real    ::  missv, lat_rng(2), z_rng(2)

  character(len=2)   ::  fid(99)
  character(len=32)  ::  expname, var_name(99)
  character(len=32)  ::  var_i_name(99), file_i_xxxx(99), var_alt(2)
  character(len=32)  ::  var_i_name2(99), file_i_xxxx2(99)
  character(len=128) ::  file_i_head, file_i_form(99), file_alt, file_o
  character(len=128) ::  file_i_head2, file_i_form2(99)
  character(len=32)  ::  nl_aux(99)

! parameters for namelists

  character(len=8), parameter ::  xblank = 'XXXX'

! common variables

  integer ::  year, mon, date, hour
  integer ::  nmon, ndate, nhour
  real    ::  day_from_ref
  logical ::  ex0, ex1
  character(len=256) ::  file_i(99), f_namelist

  integer ::  nt, nx, ny, nz, nl, nd1a, nd2a, nd3a, nd4a
  integer ::  x1_i, y1_i, z1_i, t1_i
  integer ::  nx_i, ny_i, nz_i, nt_i, n_intp
  integer ::  iv_i, iv, it_i(99)
  integer ::  k_const_th, k_const_rho
  data  it_i /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

  integer, dimension(12), parameter ::                                   &
     enddate = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  real, dimension(:),     allocatable ::  lon, lat, ht, ht_th, dim4
  real, dimension(:),     allocatable ::  t, t_i
  real, dimension(:,:,:), allocatable ::  z_th, z_rho, wt_th, wt_rho,    &
                                          wt_th2zth

  integer, dimension(:,:),   allocatable ::  ij
  integer, dimension(:,:,:), allocatable ::  kk_th, kk_rho, kk_th2zth

  character(len=32)  ::  var_i_name0(99), file_i_xxxx0(99)
  character(len=128) ::  file_i_head0, file_i_form0(99)

  integer, private ::  i,j,k,n

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

  do i=1, 99
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
        write(get_ifilename,'(a)')  trim(get_ifilename)//                &
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

  integer ::  fyear, fmon, fdate, day1_0, datpnd

  fyear = year  ;  fmon = mon

  if (day1 /= -999) then
    day1_0 = mod(day1,nday_i)
    if ( mod(date,nday_i) == day1_0 .and. hour == 0 ) then
      fdate = date
    else
      datpnd = date+nday_i
      fdate = mod( datpnd-mod(datpnd-day1_0,nday_i), ndate )
      if (fdate < date) then
        fmon = fmon + 1
        if (fmon == 13) then
          fyear = fyear + 1  ;  fmon = 1
        end if
      end if
    end if
  else
    if ( date == 1 .and. hour < hh(1) ) then
      fmon = fmon - 1
      if (fmon == 0) then
        fyear = fyear - 1  ;  fmon = 12
      end if
    end if
  end if

  SELECT case ( trim(str) )

  case ( 'FID', 'fid' )
    write(str2var,'(a)'   )  fid(iv_i)
  case ( 'VAR_I_N', 'var_i_n' )
    write(str2var,'(a)'   )  var_i_name(iv_i)
  case ( 'YYYY', 'yyyy' )
    write(str2var,'(i4.4)')  fyear
  case ( 'MM', 'mm' )
    write(str2var,'(i2.2)')  fmon
  case ( 'DD', 'dd' )
    write(str2var,'(i2.2)')  fdate
  case ( 'AUX', 'aux' )
    write(str2var,'(a)')  nl_aux(iv_i)
  case default
    print*, 'Check FILE_I_XXXX in the shell script and function '//      &
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


FUNCTION get_dayfromref(y,m,d,h)

  integer, intent(in) ::  y, m, d, h

  real ::  get_dayfromref

  get_dayfromref = (y-refdate(1))*360.+(m-refdate(2))*30.+ &
                   (d-refdate(3))+h/24.

END function get_dayfromref


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
  character(len=32) ::  dimname

  integer,          dimension(:,:,:), allocatable ::  ncvar_i
  integer*2,        dimension(:,:,:), allocatable ::  ncvar_s
  double precision, dimension(:,:,:), allocatable ::  ncvar_d

  integer,          dimension(:), allocatable ::  axis_i
  integer*2,        dimension(:), allocatable ::  axis_s
  double precision, dimension(:), allocatable ::  axis_d

  st = nf_open(trim(file_i(iv_i)), NF_NOWRITE, ncid)

  st = nf_inq_varid(ncid,trim(var_i_name(iv_i)), ncvarid)

  st = nf_inq_vardimid(ncid,ncvarid, ncdimid)

  st = nf_inq_dimlen(ncid,ncdimid(4), nt_i)

  st = nf_inq_dimname(ncid,ncdimid(4), dimname)

  st = nf_inq_varid(ncid,trim(dimname), ncdimvarid)

  st = nf_inq_vartype(ncid,ncdimvarid, ncvtype)
  allocate( t_i(nt_i) )

  SELECT case ( ncvtype )

  case ( NF_FLOAT )
    st = nf_get_var_real(ncid,ncdimvarid, t_i)
  case ( NF_DOUBLE )
    allocate( axis_d(nt_i) )
    st = nf_get_var_double(ncid,ncdimvarid, axis_d)
    t_i = real(axis_d)
    deallocate( axis_d )
  case ( NF_INT )
    allocate( axis_i(nt_i) )
    st = nf_get_var_int(ncid,ncdimvarid, axis_i)
    t_i = real(axis_i)
    deallocate( axis_i )
  case ( NF_SHORT )
    allocate( axis_s(nt_i) )
    st = nf_get_var_int2(ncid,ncdimvarid, axis_s)
    t_i = real(axis_s)
    deallocate( axis_s )
  case default
    print*, 'Check the type of variables in NETCDF file'
    STOP

  END select

  it_i(iv_i) = minloc( abs(t_i(:) - day_from_ref), dim=1 )
  if ( abs(t_i(it_i(iv_i)) - day_from_ref) > 0.0104 ) then  ! ~15 min.
    print*, 'Data at the required time is not found.'  ;  STOP
  end if

  deallocate( t_i )

  ncvstart = (/ix_i0,iy_i0,iz_i0,it_i(iv_i)/)
  ncvcount = (/nx_i0,ny_i0,nz_i0,1         /)

  st = nf_inq_vartype(ncid,ncvarid, ncvtype)

  SELECT case ( ncvtype )

  case ( NF_FLOAT )
    st = nf_get_vara_real(ncid,ncvarid,ncvstart,ncvcount, get_ivara3d)
  case ( NF_DOUBLE )
    allocate( ncvar_d(nx_i0,ny_i0,nz_i0) )
    st = nf_get_vara_double(ncid,ncvarid,ncvstart,ncvcount, ncvar_d)
    get_ivara3d(:,:,:) = real(ncvar_d(:,:,:))
    deallocate( ncvar_d )
  case ( NF_INT )
    allocate( ncvar_i(nx_i0,ny_i0,nz_i0) )
    st = nf_get_vara_int(ncid,ncvarid,ncvstart,ncvcount, ncvar_i)
    get_ivara3d(:,:,:) = real(ncvar_i(:,:,:))
    deallocate( ncvar_i )
  case ( NF_SHORT )
    allocate( ncvar_s(nx_i0,ny_i0,nz_i0) )
    st = nf_get_vara_int2(ncid,ncvarid,ncvstart,ncvcount, ncvar_s)
    get_ivara3d(:,:,:) = real(ncvar_s(:,:,:))
    deallocate( ncvar_s )
  case default
    print*, 'Check the type of variables in NETCDF file'
    STOP

  END select

  st = nf_close(ncid)

  RETURN

END function get_ivara3d


FUNCTION get_ivara4d(ix_i0,nx_i0,iy_i0,ny_i0,iz_i0,nz_i0,it_i0,nt_i0)

  integer, intent(in) ::  ix_i0, nx_i0, iy_i0, ny_i0, iz_i0, nz_i0,      &
                          it_i0, nt_i0

  real, dimension(nx_i0,ny_i0,nz_i0,nt_i0) ::  get_ivara4d

  include 'netcdf.inc'

  integer ::  st, ncid, ncvarid, ncvtype, ncvstart(4), ncvcount(4),      &
              ncdimvarid
  character(len=32) ::  dimname

  integer,          dimension(:,:,:,:), allocatable ::  ncvar_i
  integer*2,        dimension(:,:,:,:), allocatable ::  ncvar_s
  double precision, dimension(:,:,:,:), allocatable ::  ncvar_d

  integer,          dimension(:), allocatable ::  axis_i
  integer*2,        dimension(:), allocatable ::  axis_s
  double precision, dimension(:), allocatable ::  axis_d

  st = nf_open(trim(file_i(iv_i)), NF_NOWRITE, ncid)

  st = nf_inq_varid(ncid,trim(var_i_name(iv_i)), ncvarid)

  ncvstart = (/ix_i0,iy_i0,iz_i0,it_i0/)
  ncvcount = (/nx_i0,ny_i0,nz_i0,nt_i0/)

  st = nf_inq_vartype(ncid,ncvarid, ncvtype)

  SELECT case ( ncvtype )

  case ( NF_FLOAT )
    st = nf_get_vara_real(ncid,ncvarid,ncvstart,ncvcount, get_ivara4d)
  case ( NF_DOUBLE )
    allocate( ncvar_d(nx_i0,ny_i0,nz_i0,nt_i0) )
    st = nf_get_vara_double(ncid,ncvarid,ncvstart,ncvcount, ncvar_d)
    get_ivara4d(:,:,:,:) = real(ncvar_d(:,:,:,:))
    deallocate( ncvar_d )
  case ( NF_INT )
    allocate( ncvar_i(nx_i0,ny_i0,nz_i0,nt_i0) )
    st = nf_get_vara_int(ncid,ncvarid,ncvstart,ncvcount, ncvar_i)
    get_ivara4d(:,:,:,:) = real(ncvar_i(:,:,:,:))
    deallocate( ncvar_i )
  case ( NF_SHORT )
    allocate( ncvar_s(nx_i0,ny_i0,nz_i0,nt_i0) )
    st = nf_get_vara_int2(ncid,ncvarid,ncvstart,ncvcount, ncvar_s)
    get_ivara4d(:,:,:,:) = real(ncvar_s(:,:,:,:))
    deallocate( ncvar_s )
  case default
    print*, 'Check the type of variables in NETCDF file'
    STOP

  END select

  st = nf_close(ncid)

  RETURN

END function get_ivara4d

!::::  SUBROUTINES  ::::::::::::::::::::::::::::::::::::::::::::::::::::

SUBROUTINE getdim(file_i0,var_i_name0)

  include 'netcdf.inc'

  character(len=*), intent(in) ::  file_i0, var_i_name0

  integer ::  st, ncid, ncvarid, ncdimid(4), ncvtype(4), ncdimvarid
  integer ::  varndim, dimlen(4), idm
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

  allocate( lon(nx), lat(ny), ht(nz), ht_th(0:nz), dim4(nl) )

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

  lon(:) = axis_r(1:nx,1)
  lat(:) = axis_r(1:ny,2)
  ht (:) = axis_r(1:nz,3)

  dim4(:) = axis_r(1:nl,4)

  deallocate( axis_r )

  ht_th(0) = 0.
  do k=1, nz
    ht_th(k) = 2.*ht(k) - ht_th(k-1)
  enddo

  RETURN

END subroutine getdim


SUBROUTINE getalt

  include 'netcdf.inc'

  integer ::  st, ncid, ncvarid(2), ncvtype
  integer ::  l, kk, kkk, is3(4), in3(4), jn3(4)

  integer,          dimension(:,:,:), allocatable ::  ncvar_i
  integer*2,        dimension(:,:,:), allocatable ::  ncvar_s
  double precision, dimension(:,:,:), allocatable ::  ncvar_d

  is3 = (/x1_i,y1_i,z1_i,1/)  ;  in3 = (/nx_i,ny_i,nz_i,1/)
  jn3 = (/nx_i,ny_i,nz_i+1,1/)

  st = nf_open(trim(file_alt), NF_NOWRITE, ncid)

  allocate( z_th(1:nx_i,1:ny_i,0:nz_i), z_rho(nx_i,ny_i,nz_i) )

  st = nf_inq_varid(ncid,trim(var_alt(1)), ncvarid(1))
  st = nf_inq_varid(ncid,trim(var_alt(2)), ncvarid(2))
  st = nf_inq_vartype(ncid,ncvarid(1), ncvtype)

  SELECT case ( ncvtype )

  case ( NF_FLOAT )
    st = nf_get_vara_real(ncid,ncvarid(1),is3,jn3, z_th )
    st = nf_get_vara_real(ncid,ncvarid(2),is3,in3, z_rho)
  case ( NF_DOUBLE )
    allocate( ncvar_d(1:nx_i,1:ny_i,0:nz_i) )
    st = nf_get_vara_double(ncid,ncvarid(1),is3,jn3, ncvar_d)
    z_th (:,:,:) = real(ncvar_d(:,:,:))
    st = nf_get_vara_double(ncid,ncvarid(2),is3,in3, ncvar_d(:,:,1:nz_i))
    z_rho(:,:,:) = real(ncvar_d(:,:,1:nz_i))
    deallocate( ncvar_d )
  case ( NF_INT )
    allocate( ncvar_i(1:nx_i,1:ny_i,0:nz_i) )
    st = nf_get_vara_int(ncid,ncvarid(1),is3,jn3, ncvar_i)
    z_th (:,:,:) = real(ncvar_i(:,:,:))
    st = nf_get_vara_int(ncid,ncvarid(2),is3,in3, ncvar_i(:,:,1:nz_i))
    z_rho(:,:,:) = real(ncvar_i(:,:,1:nz_i))
    deallocate( ncvar_i )
  case ( NF_SHORT )
    allocate( ncvar_s(1:nx_i,1:ny_i,0:nz_i) )
    st = nf_get_vara_int2(ncid,ncvarid(1),is3,jn3, ncvar_s)
    z_th (:,:,:) = real(ncvar_s(:,:,:))
    st = nf_get_vara_int2(ncid,ncvarid(2),is3,in3, ncvar_s(:,:,1:nz_i))
    z_rho(:,:,:) = real(ncvar_s(:,:,1:nz_i))
    deallocate( ncvar_s )
  case default
    print*, 'Check the type of variables in NETCDF file'
    STOP

  END select

  st = nf_close(ncid)

  k_const_th = nz_i + 1  ;  k_const_rho = nz_i + 1
  do k=0, nz_i
    if ( all(z_th(:,:,k) == z_th(1,1,k)) ) then
      k_const_th = max(1, k)  ;  EXIT
    end if
  enddo
  do k=1, nz_i
    if ( all(z_rho(:,:,k) == z_rho(1,1,k)) ) then
      k_const_rho = k  ;  EXIT
    end if
  enddo

  l = 0
  do j=1, ny_i
  do i=1, nx_i
    if (z_th(i,j,0) /= 0.)  l = l + 1
  enddo
  enddo
  n_intp = l

  allocate( ij(n_intp,2) )
  allocate( kk_th(n_intp,nz_i,2), kk_rho(n_intp,nz_i,2), kk_th2zth(n_intp,nz_i,2) )
  allocate( wt_th(n_intp,nz_i,2), wt_rho(n_intp,nz_i,2), wt_th2zth(n_intp,nz_i,2) )

  kk_th(:,:,:) = 1   ;  kk_rho(:,:,:) = 1   ;  kk_th2zth(:,:,:) = 1
  wt_th(:,:,:) = 0.  ;  wt_rho(:,:,:) = 0.  ;  wt_th2zth(:,:,:) = 0.

  l = 0
  do j=1, ny_i
  do i=1, nx_i
    if (z_th(i,j,0) /= 0.) then
      l = l + 1 
      ij(l,:) = (/i,j/)
      kkk = 0
      do k=1, nz_i
      do kk=kkk, nz_i-1
        if ( z_th(i,j,kk) < ht(k) .and. z_th(i,j,kk+1) >= ht(k) ) then
          kk_th(l,k,:) = (/kk,kk+1/)
          wt_th(l,k,:) = (/z_th(i,j,kk+1)-ht(k),ht(k)-z_th(i,j,kk)/)
          wt_th(l,k,:) = wt_th(l,k,:)/(z_th(i,j,kk+1)-z_th(i,j,kk))
          kkk = kk  ;  EXIT
        end if
      enddo
      enddo
      kkk = 1
      do k=1, nz_i
      do kk=kkk, nz_i-1
        if ( z_rho(i,j,kk) < ht(k) .and. z_rho(i,j,kk+1) >= ht(k) ) then
          kk_rho(l,k,:) = (/kk,kk+1/)
          wt_rho(l,k,:) = (/z_rho(i,j,kk+1)-ht(k),ht(k)-z_rho(i,j,kk)/)
          wt_rho(l,k,:) = wt_rho(l,k,:)/(z_rho(i,j,kk+1)-z_rho(i,j,kk))
          kkk = kk  ;  EXIT
        end if
      enddo
      enddo
      kkk = 0
      do k=1, nz_i
      do kk=kkk, nz_i-1
        if ( z_th(i,j,kk) < ht_th(k) .and. z_th(i,j,kk+1) >= ht_th(k) ) then
          kk_th2zth(l,k,:) = (/kk,kk+1/)
          wt_th2zth(l,k,:) = (/z_th(i,j,kk+1)-ht_th(k),ht_th(k)-z_th(i,j,kk)/)
          wt_th2zth(l,k,:) = wt_th2zth(l,k,:)/(z_th(i,j,kk+1)-z_th(i,j,kk))
          kkk = kk  ;  EXIT
        end if
      enddo
      enddo
    end if
  enddo
  enddo

  RETURN

END subroutine getalt


SUBROUTINE get_iouter(arr,arr_rng, ind2)

  real,    dimension(:), intent(in)  ::  arr
  real,    dimension(2), intent(in)  ::  arr_rng
  integer, dimension(2), intent(out) ::  ind2

  integer ::  narr

  narr = size(arr)

  ind2 = (/1,narr/)

  if (arr_rng(1) /= -999.) then
    do i=2, narr
      if (arr(i) > arr_rng(1)) then
        ind2(1) = i - 1  ;  EXIT
      end if
    enddo
  end if
  if (arr_rng(2) /= -999.) then
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
  var_i_name0 (:) = var_i_name (:)

  file_i_head    = file_i_head2
  file_i_form(:) = file_i_form2(:)
  file_i_xxxx(:) = file_i_xxxx2(:)
  var_i_name (:) = var_i_name2 (:)

  file_i_head2    = file_i_head0
  file_i_form2(:) = file_i_form0(:)
  file_i_xxxx2(:) = file_i_xxxx0(:)
  var_i_name2 (:) = var_i_name0 (:)

END subroutine switch_para_in


SUBROUTINE v_cubintp_rho2z(var_intp)

  use nr, only: spline, splint

  real, dimension(nx_i,ny_i,nz_i), intent(inout) ::  var_intp

  real, dimension(nz_i,nx_i,ny_i) ::  zz, vv, work

  do j=1, ny_i
  do i=1, nx_i
    zz(:,i,j) = z_rho   (i,j,1:nz_i)
    vv(:,i,j) = var_intp(i,j,:     )
  enddo
  enddo

  do j=1, ny_i
  do i=1, nx_i
    call spline(zz(:,i,j),vv(:,i,j),1.e32,1.e32,work(:,i,j))
  enddo
  enddo

  do k=1, nz_i
  do j=1, ny_i
  do i=1, nx_i
    var_intp(i,j,k) = splint(zz(:,i,j),vv(:,i,j),work(:,i,j),ht(k))
  enddo
  enddo
  enddo

  RETURN

END subroutine v_cubintp_rho2z


SUBROUTINE v_cubintp_0th2z(var_intp,var_zend)

  use nr, only: spline, splint

  real, dimension(nx_i,ny_i),            intent(in)    ::  var_zend
  real, dimension(1:nx_i,1:ny_i,0:nz_i), intent(inout) ::  var_intp

  real, dimension(0:nz_i+1,1:nx_i,1:ny_i) ::  zz, vv, work

  do j=1, ny_i
  do i=1, nx_i
    zz(0:nz_i,i,j) = z_th    (i,j,0:nz_i)
    vv(0:nz_i,i,j) = var_intp(i,j, :    )
    zz(nz_i+1,i,j) = z_rho   (i,j,nz_i+1)
    vv(nz_i+1,i,j) = var_zend(i,j       )
  enddo
  enddo

  var_intp(:,:,0) = 0.

  do j=1, ny_i
  do i=1, nx_i
    call spline(zz(:,i,j),vv(:,i,j),1.e32,1.e32,work(:,i,j))
  enddo
  enddo

  do k=1, nz_i
  do j=1, ny_i
  do i=1, nx_i
    var_intp(i,j,k) = splint(zz(:,i,j),vv(:,i,j),work(:,i,j),ht(k))
  enddo
  enddo
  enddo

  RETURN

END subroutine v_cubintp_0th2z


SUBROUTINE v_cubintp_th2z(var_intp,var_zend)

  use nr, only: spline, splint

  real, dimension(nx_i,ny_i),      intent(in)    ::  var_zend
  real, dimension(nx_i,ny_i,nz_i), intent(inout) ::  var_intp

  real, dimension(nz_i+1,nx_i,ny_i) ::  zz, vv, work

  do j=1, ny_i
  do i=1, nx_i
    zz(:nz_i ,i,j) = z_th    (i,j,1:nz_i)
    vv(:nz_i ,i,j) = var_intp(i,j, :    )
    zz(nz_i+1,i,j) = z_rho   (i,j,nz_i+1)
    vv(nz_i+1,i,j) = var_zend(i,j       )
  enddo
  enddo

  do j=1, ny_i
  do i=1, nx_i
    call spline(zz(:,i,j),vv(:,i,j),1.e32,1.e32,work(:,i,j))
  enddo
  enddo

  do k=1, nz_i
  do j=1, ny_i
  do i=1, nx_i
    var_intp(i,j,k) = splint(zz(:,i,j),vv(:,i,j),work(:,i,j),ht(k))
  enddo
  enddo
  enddo

  RETURN

END subroutine v_cubintp_th2z


SUBROUTINE v_cubintp_th2zth(var_intp)

  use nr, only: spline, splint

  real, dimension(nx_i,ny_i,nz_i), intent(inout) ::  var_intp

  real, dimension(nz_i,nx_i,ny_i) ::  zz, vv, work

  do j=1, ny_i
  do i=1, nx_i
    zz(:,i,j) = z_th    (i,j,1:nz_i)
    vv(:,i,j) = var_intp(i,j,:     )
  enddo
  enddo

  do j=1, ny_i
  do i=1, nx_i
    call spline(zz(:,i,j),vv(:,i,j),1.e32,1.e32,work(:,i,j))
  enddo
  enddo

  do k=1, nz_i
  do j=1, ny_i
  do i=1, nx_i
    var_intp(i,j,k) = splint(zz(:,i,j),vv(:,i,j),work(:,i,j),ht_th(k))
  enddo
  enddo
  enddo

  RETURN

END subroutine v_cubintp_th2zth


SUBROUTINE v_linintp_rho2z(var_intp)

  real, dimension(nx_i,ny_i,nz_i), intent(inout) ::  var_intp

  real, dimension(n_intp,k_const_rho) ::  vv

  integer ::  l, kl

  kl = min(nz_i, k_const_rho)

  do l=1, n_intp
    vv(l,1:kl) = var_intp(ij(l,1),ij(l,2),1:kl)
  enddo

  do k=1, kl
  do l=1, n_intp
    var_intp(ij(l,1),ij(l,2),k) = sum(vv(l,kk_rho(l,k,:))*wt_rho(l,k,:))
  enddo
  enddo

  RETURN

END subroutine v_linintp_rho2z


SUBROUTINE v_linintp_0th2z(var_intp)

  real, dimension(1:nx_i,1:ny_i,0:nz_i), intent(inout) ::  var_intp

  real, dimension(1:n_intp,0:k_const_th) ::  vv
  real, dimension(1:nx_i,1:ny_i,k_const_th:nz_i) ::  vv2

  integer ::  l, kl

  kl = min(nz_i, k_const_th)

!  do l=1, n_intp
!    vv(l,0:kl) = var_intp(ij(l,1),ij(l,2),0:kl)
!  enddo
  if (k_const_th < nz_i)  vv2(:,:,:) = var_intp(:,:,k_const_th:nz_i)
!
!  var_intp(:,:,0) = 0.
!  do k=1, kl
!  do l=1, n_intp
!    var_intp(ij(l,1),ij(l,2),k) = sum(vv(l,kk_th(l,k,:))*wt_th(l,k,:))
!  enddo
!  enddo

  do k=k_const_th+1, nz_i
    var_intp(:,:,k) = 0.5*(vv2(:,:,k-1)+vv2(:,:,k))
  enddo

  RETURN

END subroutine v_linintp_0th2z


SUBROUTINE v_linintp_th2z(var_intp)

  real, dimension(nx_i,ny_i,nz_i), intent(inout) ::  var_intp

  real, dimension(1:n_intp,0:k_const_th) ::  vv
  real, dimension(1:nx_i,1:ny_i,k_const_th:nz_i) ::  vv2

  integer ::  l, kl

  kl = min(nz_i, k_const_th)

!  do l=1, n_intp
!    vv(l,1:kl) = var_intp(ij(l,1),ij(l,2),1:kl)
!  enddo
!  vv(:,0) = vv(:,1)
  if (k_const_th < nz_i)  vv2(:,:,:) = var_intp(:,:,k_const_th:nz_i)
!
!  do k=1, kl
!  do l=1, n_intp
!    var_intp(ij(l,1),ij(l,2),k) = sum(vv(l,kk_th(l,k,:))*wt_th(l,k,:))
!  enddo
!  enddo

  do k=k_const_th+1, nz_i
    var_intp(:,:,k) = 0.5*(vv2(:,:,k-1)+vv2(:,:,k))
  enddo

  RETURN

END subroutine v_linintp_th2z


SUBROUTINE v_linintp_th2zth(var_intp)

  real, dimension(nx_i,ny_i,nz_i), intent(inout) ::  var_intp

  real, dimension(1:n_intp,0:k_const_th) ::  vv

  integer ::  l, kl

  kl = min(nz_i, k_const_th)

  do l=1, n_intp
    vv(l,1:kl) = var_intp(ij(l,1),ij(l,2),1:kl)
  enddo
  vv(:,0) = vv(:,1)

  do k=1, kl
  do l=1, n_intp
    var_intp(ij(l,1),ij(l,2),k) = sum(vv(l,kk_th2zth(l,k,:))*wt_th2zth(l,k,:))
  enddo
  enddo

  RETURN

END subroutine v_linintp_th2zth


END module hadgem


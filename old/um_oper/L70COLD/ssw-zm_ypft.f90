! UM_OPER / L70COLD

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  character(len=128) ::  f_namelist

  integer ::  yyyymm, date(2), utc(12)
  real    ::  fct(2), fct_itv
  character(len=10)  ::  expname, exception(99)
  character(len=128) ::  file_i, file_o, fnamet(2), vname_i(99), vname_o(99)

  namelist /ANALCASE/ EXPNAME, YYYYMM
  namelist /PARAM/ DATE, UTC, FCT, FCT_ITV, EXCEPTION
  namelist /FILEIO/ FILE_I, FILE_O, FNAMET, VNAME_I, VNAME_O

  ! input variables
  data  f_namelist  /'./namelist/nl.input'/

  integer ::  nutc, nt, nfct, nx, ny, nz, nxu, nyv, dlen(4)
  integer ::  nxr, nyr, nzr, ntr, nd1a, nd2a, nd3a, nd4a
  integer ::  i,j,k,n,nv, iv, idat, iutc, ifct, i_time, it, ncid
  integer ::  year, month
  logical ::  l_ftail, l_getdim, ex1
  character(len=10)  ::  timec
  character(len=32)  ::  c_axis(3,2)
  character(len=256) ::  fname1
  character(len=3)   ::  hhh1

  real, parameter ::  missv = 1.e32

  real, dimension(:,:,:,:), allocatable ::  var, var_za
  real, dimension(:),       allocatable ::  lon, lat, p, lonu, latv
  real, dimension(:),       allocatable ::  t, t_fct

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(:), allocatable ::  set

! READ NAMELISTS

  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  l_ftail = .TRUE.
  if ( len_trim(fnamet(2)) == len(fnamet(2)) .or. len_trim(fnamet(2)) == 0 )  &
     l_ftail = .FALSE.
  do i=1, 99
    if ( len_trim(vname_i(i)) == len(vname_i(i)) .or. len_trim(vname_i(i)) == 0 ) then
      nv = i - 1  ;  EXIT
    end if
  enddo

  allocate( set(nv) )

! DEFINE AXIS

  do n=1, 12
    if (utc(n) < 0) then  ;  nutc = n - 1  ;  EXIT  ;  end if
  enddo

  call gettaxis

  l_getdim = .TRUE.

! LOOP

  L_VAR:  DO iv=1, nv

  i_time = 0


  L_DAT:  DO idat=date(1), date(2)
  L_UTC:  DO iutc=1, nutc


  year  = int(yyyymm/100)
  month = yyyymm - year*100

  i_time = i_time + 1


  L_FCT:  DO ifct=1, nfct


  write(timec ,'(i6.6,2i2.2)') yyyymm, idat, utc(iutc)

  write(hhh1,'(i3.3)') int(t_fct(ifct)*24.)

  write(6,*)
  write(6,'(a)') ' '//timec//' + '//hhh1//' h'

  if (l_ftail) then
    fname1 = trim(file_i)//'/'//timec//'/'//trim(expname)//'.'//trim(fnamet(1))// &
             '.'//timec//'+'//hhh1//trim(fnamet(2))//'.nc'
  else
    fname1 = trim(file_i)//'/'//timec//'/'//trim(expname)//'.'//trim(fnamet(1))// &
             '.'//timec//'+'//hhh1//'.nc'
  end if

  call check_ex(fname1, ex1)
  if (.not. ex1)  CYCLE

  if (l_getdim)  call getdim

  ! allocate output
  nd1a = NY
  nd2a = NZ
  nd3a = NFCT
  nd4a = NT

  if (.not. allocated(set(iv)%var_out)) then
    allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )
    set(iv)%var_out(:,:,:,:) = missv
  end if

  ! get var.s
  call opennc(trim(fname1), ncid)

  if (.not. allocated(var)) then
    call varlen(ncid,trim(vname_i(iv)),4, dlen)
    nxr = dlen(1)
    nyr = dlen(2)
    nzr = dlen(3)
    ntr = 1
    allocate( var(nxr,nyr,nzr,ntr), var_za(1,nyr,nzr,ntr) )
  end if

  call get4d(ncid,trim(vname_i(iv)),nxr,nyr,nzr,ntr, var)

  call closenc(ncid)

  ! calculate zonal mean
  call zonal_avg(nxr,nyr,nzr,ntr,0.,var, var_za)

  if ( nyr == nd1a ) then
    set(iv)%var_out(:,:,ifct,i_time) = var_za(1,:,:,1)
  else if ( nyr == nd1a-1 ) then
    do k=1, nzr
    do j=2, nyr
      set(iv)%var_out(j,k,ifct,i_time) = 0.5*(var_za(1,j-1,k,1)+var_za(1,j,k,1))
    enddo
    enddo
    set(iv)%var_out(1   ,:,ifct,i_time) = 0.
    set(iv)%var_out(nd1a,:,ifct,i_time) = 0.
  else
    print*, 'Unexpected length of y-dimension !!!'
    STOP
  end if


  ENDDO  L_FCT


  ENDDO  L_UTC
  ENDDO  L_DAT


  deallocate( var, var_za )


  ENDDO  L_VAR


  do iv=1, nv
    set(iv)%vname = vname_o(iv)
    set(iv)%axis = (/'lat  ','p ','fcst','t'/)
    set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
    allocate( set(iv)%axis1(set(iv)%nd(1)) )
    allocate( set(iv)%axis2(set(iv)%nd(2)) )
    allocate( set(iv)%axis3(set(iv)%nd(3)) )
    allocate( set(iv)%axis4(set(iv)%nd(4)) )
    set(iv)%axis1 = lat
    set(iv)%axis2 = p
    set(iv)%axis3 = t_fct
    set(iv)%axis4 = t
  enddo

! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'Zonal Mean')

! END

  deallocate( t_fct, t )
  if ( allocated(lon) )  deallocate( lon )
  if ( allocated(lat) )  deallocate( lat )
  if ( allocated(p  ) )  deallocate( p   )
  if ( allocated(lonu) )  deallocate( lonu )
  if ( allocated(latv) )  deallocate( latv )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3, set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo
  deallocate( set )

  STOP


  CONTAINS

  SUBROUTINE gettaxis

    ! t
    nt = (date(2)-date(1)+1)*nutc
    allocate( t(nt) )
    n = 0
    do idat=date(1), date(2)
    do iutc=1, nutc
      n = n + 1
      t(n) = float(idat) + float(utc(iutc))/24.
    enddo
    enddo

    ! fct
    nfct = int((fct(2)-fct(1))/fct_itv)+1
    allocate( t_fct(nfct) )
    do ifct=1, nfct
      t_fct(ifct) = fct(1) + fct_itv*(ifct-1)
    enddo

  END subroutine gettaxis

  SUBROUTINE getdim

    call opennc(fname1, ncid)

    call diminfop(ncid,.TRUE., nx,ny,nz,nxu,nyv,c_axis)
    allocate( lat(ny) )
    if (trim(c_axis(1,1)) /= empty)  allocate( lon(nx) )
    if (trim(c_axis(3,1)) /= empty)  allocate( p  (nz) )
    if (trim(c_axis(1,2)) /= empty)  allocate( lonu(nxu) )
    if (trim(c_axis(2,2)) /= empty)  allocate( latv(nyv) )
    call axisinfo(ncid,nx,ny,nz,nxu,nyv,1,c_axis, lon,lat,p,lonu,latv,p)
    if (trim(c_axis(2,1)) == empty) then
      do j=1, ny
        lat(j) = float(j-1-(ny-1)/2)*180./(ny-1)
      enddo
    end if

    call closenc(ncid)

    l_getdim = .FALSE.

  END subroutine getdim


END program analysis_UM


SUBROUTINE check_ex(fname,existence)

  implicit none

  character(len=*), intent(in)  ::  fname
  logical,          intent(out) ::  existence

  inquire(file=trim(fname), exist=existence)
  if (.not. existence)  print*, '    ',trim(fname),' not found. - passed'

  RETURN

END subroutine check_ex


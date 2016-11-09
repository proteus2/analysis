PROGRAM STREAMFTN_RCIRC_UM

  use hadgem
  use netio
  use nr, only: svbksb

  implicit none

  include 'c_phys.inc'

  integer, parameter ::  nv = 1
!  real,    parameter ::  zb_out = 8.e3

  real ::  dz_intp

  namelist /ANALCASE/ EXPNAME, YYYY, MM
  namelist /PARAM/ DZ_INTP, Z_RNG
  namelist /FILEIO/ MISSV, FILE_I_HEAD, FILE_I_FORM, FILE_I_XXXX, &
                    VAR_I_NAME, FILE_O

  integer ::  iz, imon, i_time1, i_time9
  integer ::  i,k, ii, nz0, nt0, nt0m, nn, nyn, nzn
  character(len=32) ::  tname
  character(len=32), dimension(nv) ::  ovarname

  real, dimension(:,:,:,:), allocatable ::  var4d
  real, dimension(:,:,:),   allocatable ::  var3d
  real, dimension(:,:,:),   allocatable ::  v, w
  real, dimension(:,:),     allocatable ::  mat, mat2
  real, dimension(:),       allocatable ::  ht0
  real, dimension(:),       allocatable ::  sval, b, sol

  type(vset), dimension(nv) ::  set

  ovarname = (/'strf'/)

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! GET AXES AND INITIALIZE ARRAYS

  year = yyyy
  mon  = mm(1)

  call initialize

  call calc_matrix

  i_time1 = 1

  L_MON:  DO imon=1, nmon
  !---------------------------------------------------------------------
  call get_2var
!  call v_interpol2z_5var

  ! calculate zonal mean
b = 0.d0
  call svbksb(mat,sval,mat2,b,sol)

  i_time9 = i_time1 + nt0 - 1

  var4d(:,:,i_time1:i_time9,1) = w(:,:,:)

  mon = mon + 1
  if (mon == 13) then
    year = year + 1  ;  mon = 1
  end if

  i_time1 = i_time9 + 1
  !---------------------------------------------------------------------
  ENDDO  L_MON

  nt = i_time9

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

!  if (zb_out /= 0.) then
!    iz = nd2a + 1
!    do k=2, nd2a
!      if (ht(k) > zb_out) then  ;  iz = k - 1  ;  EXIT  ;  end if
!    enddo
!    do iv=1, nv
!      set(iv)%var_out(:,:iz-1,:,1) = 1.e32
!    enddo
!  end if


! DUMP

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  call outnc(trim(file_o),nv,set,'stream function in TEM eqn.')

! END

  call finalize

  STOP


  CONTAINS


  SUBROUTINE initialize

  include 'netcdf.inc'

  integer ::  st, ncid, ncvarid, ncdimid(3), ncdimvarid, idm
  character(len=32) ::  dimname(3)

  nmon = mm(2)
  day1 = -999  ;  date = 15  ! for get_ifilename

  iv_i = 1  ! for get_ifilename
  file_i(1) = get_ifilename()
  inquire(file=trim(file_i(1)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(1)),' not found.'  ;  STOP
  end if

  st = nf_open(trim(file_i(1)),NF_NOWRITE, ncid)
  st = nf_inq_varid(ncid,trim(var_i_name(1)), ncvarid)
  st = nf_inq_vardimid(ncid,ncvarid, ncdimid)
  st = nf_inq_dimlen(ncid,ncdimid(1), ny )
  st = nf_inq_dimlen(ncid,ncdimid(2), nz0)
  st = nf_inq_dimlen(ncid,ncdimid(3), nt0)
  allocate( lat(ny), ht0(nz0) )
  do idm=1, 3
    st = nf_inq_dimname(ncid,ncdimid(idm), dimname(idm))
  enddo
  st = nf_inq_varid(ncid,trim(dimname(1)), ncdimvarid)
  st = nf_get_var_real(ncid,ncdimvarid, lat)
  st = nf_inq_varid(ncid,trim(dimname(2)), ncdimvarid)
  st = nf_get_var_real(ncid,ncdimvarid, ht0)
  st = nf_close(ncid)
  tname = trim(dimname(3))

  nt0m = int(nt0*(31./28.))
  nz = nz0
  allocate( ht(nz) )
  ht(:) = ht0

  allocate( var4d(ny,nz,nmon*nt0m,nv), t(nmon*nt0m) )
  var4d(:,:,:,:) = 0.  ;  t(:) = 0.
  allocate( var3d(ny,nz,nv) )

  allocate( v(ny,nz0,nt0m), w(ny,nz0,nt0m) )

  nyn = ny - 4
  nzn = 100  ! nz
  nn = nyn*nzn
  allocate( mat(nn,nn), mat2(nn,nn), sval(nn), b(nn), sol(nn) )

  END subroutine initialize

  SUBROUTINE calc_matrix

  use nr, only: svdcmp

  double precision, dimension(nn,nn) ::  mat_8, mat2_8
  double precision, dimension(nn)    ::  sval_8
  real,             dimension(nn)    ::  tmp1d

  real    ::  tmpm, sval_m
  integer ::  ii

  mat(:,:) = 0.
  do i=nyn+2, nn-nyn-1
    mat(i,i) = -4.
  enddo
  do i=1, nn-1
    mat(i,i+1) = 1.  ;  mat(i+1,i) = 1.
  enddo
  do i=1, nn-nyn
    mat(i,i+nyn) = 1.  ;  mat(i+nyn,i) = 1.
  enddo
  ! for boundary conditions
  do ii=1, nzn-2
    i = ii*nyn+1
    mat(i,i) = -3.  ;  mat(i-1,i) = 0.
    i = (ii+1)*nyn
    mat(i,i) = -3.  ;  mat(i+1,i) = 0.
  enddo
  do i=1, nyn
    mat(i,i) = -3.
  enddo
  mat(1,1) = -2.  ;  mat(nyn,nyn) = -2.  ;  mat(nyn+1,nyn) = 0.
  do i=nn-nyn+1, nn
    mat(i,i) = -3.
  enddo
  mat(nn-nyn+1,nn-nyn+1) = -2.  ;  mat(nn,nn) = -2.  ;  mat(nn-1,nn) = 0.
 
  mat_8 = dble(mat)
  call svdcmp(mat_8,sval_8,mat2_8)
  sval = real(sval_8)
  mat2 = real(mat2_8)

  do i=1, nn
    tmp1d(i) = float(i)
  enddo

  call out2d1d('/data11/kyh/analy/um_clim-ao/mat.nc', &
               2,(/'mat','mat2'/),(/mat,mat2/),       &
               'dim1',nn,tmp1d,'dim2',nn,tmp1d,       &
               1,(/'sval'/),(/sval/),                 &
               'dim',nn,tmp1d,                        &
               'SVD matrices')
!  call opennc('/data11/kyh/analy/um_clim-ao/mat.nc',ncid)
!  call get2d(ncid,'mat' ,nn,nn,mat )
!  call get2d(ncid,'mat2',nn,nn,mat2)
!  call get1d(ncid,'sval',nn,sval)
!  call closenc(ncid)
print*,sval

  sval_m = maxval(sval)*1.e-6
  where ( sval < sval_m )
    sval = 0.
  end where 

  END subroutine calc_matrix

  SUBROUTINE get_2var

  include 'netcdf.inc'

  integer ::  st, ncid, ncvarid, ncdimid(3), ncdimvarid
  character(len=32) ::  dimname

  iv_i = 1  ! for get_ifilename
  file_i(1) = get_ifilename()
  inquire(file=trim(file_i(1)), exist=ex1)
  if ( .not. ex1 ) then
    print*, '    ',trim(file_i(1)),' not found.'  ;  STOP
  end if

  ! read 2 var.s and t
  st = nf_open(trim(file_i(1)),NF_NOWRITE, ncid)
  st = nf_inq_varid(ncid,trim(var_i_name(1)), ncvarid)
  st = nf_inq_vardimid(ncid,ncvarid, ncdimid)
  st = nf_inq_dimlen(ncid,ncdimid(3), nt0)
  st = nf_inq_varid(ncid,trim(tname), ncdimvarid)
  st = nf_get_var_real(ncid,ncdimvarid, t(i_time1:i_time1+nt0-1))
  st = nf_get_var_real(ncid,ncvarid, v(:,:,:nt0))
  st = nf_inq_varid(ncid,trim(var_i_name(2)), ncvarid)
  st = nf_get_var_real(ncid,ncvarid, w(:,:,:nt0))
  st = nf_close(ncid)

  END subroutine get_2var

!  SUBROUTINE v_interpol2z_5var
!
!  real, dimension(:,:,:), allocatable ::  lnvar
!
!  nx_i = nx  ;  ny_i = ny  ;  nz_i = k_const_rho+1   ! for v_cubintp_*
!  call v_cubintp_rho2z (u(:,:,:nz_i))
!  call v_cubintp_rho2z (v(:,:,:nz_i))
!  call v_cubintp_th2zth(w(:,:,:nz_i))
!
!  allocate( lnvar(nx_i,ny_i,nz_i) )
!
!  lnvar(:,:,:) = log(rho(:,:,:nz_i))
!  call v_cubintp_rho2z(lnvar)
!  rho(:,:,:nz_i) = exp(lnvar(:,:,:))
!
!  lnvar(:,:,:) = log(pt(:,:,:nz_i))
!  call v_cubintp_th2zth(lnvar)
!  pt(:,:,:nz_i) = exp(lnvar(:,:,:))
!
!  deallocate( lnvar )
!
!  nx_i = nx  ;  ny_i = ny  ;  nz_i = nz
!
!  END subroutine v_interpol2z_5var

  SUBROUTINE setdim

  set(iv)%vname = trim(ovarname(iv))
  set(iv)%axis = (/'lat  ','z ','time',' '/) 
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat
  set(iv)%axis2 = ht
  set(iv)%axis3 = t
  set(iv)%axis4 = -999.
    
  END subroutine setdim

  SUBROUTINE finalize

  deallocate( v, w )
  deallocate( var4d, var3d )
  deallocate( lat, ht, t, ht0 )
  do iv=1, nv
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3,             &
                set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  END subroutine finalize


END program STREAMFTN_RCIRC_UM


PROGRAM FFT_x

  use fft
  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2)
  character(len=16)  ::  var_name
  character(len=256) ::  f_namelist, file_i, file_o

  namelist /PARAM/ VAR_NAME, NNN, JJ, II
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  i1, j1, n1, nx, ny, nt, nk
  integer ::  i,j,n

  real,    dimension(:,:,:,:), allocatable ::  var4d
  real,    dimension(:,:,:),   allocatable ::  var_in
  real,    dimension(:),       allocatable ::  kwn, lat, t
  complex, dimension(:),       allocatable ::  fc_var

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  ! calculate FFT
  allocate( fc_var(nx) )
  do n=1, nt
  do j=1, ny
    call fft1d_f(nx,var_in(:,j,n),fc_var)
    var4d(:,j,n,1) = real (fc_var(1:nk))
    var4d(:,j,n,2) = aimag(fc_var(1:nk))
  enddo
  enddo
  deallocate( fc_var )

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  i1 = ii(1)  ;  j1 = jj(1)  ;  n1 = nnn(1)
  nx = ii(2) - i1 + 1  ;  ny = jj(2) - j1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1
  allocate( lat(ny), t(nt) )
  do j=1, ny
    lat(j) = float(j1 + j - 1)
  enddo
  do n=1, nt
    t(n) = float(n1 + (n - 1)*nnn(3))
  enddo

  nk = nx/2 + 1
  allocate( kwn(nk) )
  do i=1, nk
    kwn(i) = float(i-1)
  enddo

  allocate( var_in(nx,ny,nt) )

  call get_var(file_i,var_name,var_in,start=(/i1,j1,n1/), count=(/nx,ny,nt/), &
               stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))

  allocate( var4d(nk,ny,nt,2) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'k_wn',kwn,axis='k_wn')
  call put_var('append'   ,file_o,'lat' ,lat,axis='lat' )
  call put_var('append'   ,file_o,'t'   ,t  ,axis='t'   )

  call put_var('append'   ,file_o,'fc_'//trim(var_name)//'_r',var4d(:,:,:,1), &
               axes=(/'k_wn','lat','t'/))
  call put_var('append'   ,file_o,'fc_'//trim(var_name)//'_i',var4d(:,:,:,2), &
               axes=(/'k_wn','lat','t'/))

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( var4d, var_in )
  deallocate( kwn )
  deallocate( lat, t )

  END subroutine finalize

END program FFT_x


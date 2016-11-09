PROGRAM FFT_y

  use fft
  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2)
  character(len=16)  ::  var_name
  character(len=256) ::  f_namelist, file_i, file_o

  namelist /PARAM/ VAR_NAME, NNN, JJ, II
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  i1, j1, n1, nx, ny, nt, nl
  integer ::  i,j,n

  real,    dimension(:,:,:,:), allocatable ::  var4d
  real,    dimension(:,:,:),   allocatable ::  var_in
  real,    dimension(:),       allocatable ::  lwn, lon, t
  complex, dimension(:),       allocatable ::  fc_var

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  ! subtract zonal mean
  var_in(:,:,:) = var_in(:,:,:) - spread( sum(var_in, dim=1)/float(nx), 1, nx )

  ! calculate FFT
  allocate( fc_var(ny) )
  do n=1, nt
  do i=1, nx
    call fft1d_f(ny,var_in(i,:,n),fc_var)
    var4d(i,:,n,1) = real (fc_var(1:nl))
    var4d(i,:,n,2) = aimag(fc_var(1:nl))
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
  allocate( lon(nx), t(nt) )
  do i=1, nx
    lon(i) = float(i1 + i - 1)
  enddo
  do n=1, nt
    t(n) = float(n1 + (n - 1)*nnn(3))
  enddo

  nl = ny/2 + 1
  allocate( lwn(nl) )
  do j=1, nl
    lwn(j) = float(j-1)
  enddo

  allocate( var_in(nx,ny,nt) )

  call get_var(file_i,var_name,var_in,start=(/i1,j1,n1/), count=(/nx,ny,nt/), &
               stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))

  allocate( var4d(nx,nl,nt,2) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'lon' ,lon,axis='lon' )
  call put_var('append'   ,file_o,'l_wn',lwn,axis='l_wn')
  call put_var('append'   ,file_o,'t'   ,t  ,axis='t'   )

  call put_var('append'   ,file_o,'fc_'//trim(var_name)//'_r',var4d(:,:,:,1), &
               axes=(/'lon','l_wn','t'/))
  call put_var('append'   ,file_o,'fc_'//trim(var_name)//'_i',var4d(:,:,:,2), &
               axes=(/'lon','l_wn','t'/))

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( var4d, var_in )
  deallocate( lwn )
  deallocate( lon, t )

  END subroutine finalize

END program FFT_y


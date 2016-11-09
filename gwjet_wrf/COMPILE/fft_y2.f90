PROGRAM FFT_y2

  use fft
  use easync

  implicit none

  integer            ::  nnn(3), jj(2), k_rng(2)
  character(len=16)  ::  var_name(2)
  character(len=256) ::  f_namelist, file_i, file_o

  namelist /PARAM/ VAR_NAME, NNN, JJ, K_RNG
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  i1, j1, n1, ny, nt, nk, nl
  integer ::  i,j,n

  real,    dimension(:,:,:,:), allocatable ::  var4d
  real,    dimension(:,:,:,:), allocatable ::  var_in
  real,    dimension(:),       allocatable ::  kwn, lwn, t
  complex, dimension(:,:,:),   allocatable ::  var_in_c
  complex, dimension(:),       allocatable ::  fc_var

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  ! calculate FFT
  allocate( fc_var(ny) )
  do n=1, nt
  do i=1, nk
    call fft1d_f(ny,var_in_c(i,:,n),fc_var)
    var4d(i,:,n,1) = real (fc_var(:))
    var4d(i,:,n,2) = aimag(fc_var(:))
  enddo
  enddo
  deallocate( fc_var )

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  i1 = k_rng(1)+1  ;  j1 = jj(1)  ;  n1 = nnn(1)
  nk = k_rng(2)+1 - i1 + 1  ;  ny = jj(2) - j1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1
  allocate( t(nt) )
  do n=1, nt
    t(n) = float(n1 + (n - 1)*nnn(3))
  enddo

  allocate( kwn(nk) )
  call get_var(file_i,'k_wn',kwn,start=(/i1/), count=(/nk/))
  if ( kwn(1) /= float(k_rng(1)) ) then
    print*, 'k_wn miss match'  ;  STOP
  end if

  nl = ny
  allocate( lwn(nl) )
  do j=1, nl
    lwn(j) = float(j-1)
  enddo

  allocate( var_in_c(nk,ny,nt) )
  allocate( var_in(nk,ny,nt,2) )

  call get_var(file_i,var_name(1),var_in(:,:,:,1),start=(/i1,j1,n1/), &
               count=(/nk,ny,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*ny/))
  call get_var(file_i,var_name(2),var_in(:,:,:,2),start=(/i1,j1,n1/), &
               count=(/nk,ny,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*ny/))
  var_in_c = cmplx(var_in(:,:,:,1),var_in(:,:,:,2))

  deallocate( var_in )

  allocate( var4d(nk,nl,nt,2) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'k_wn',kwn,axis='k_wn')
  call put_var('append'   ,file_o,'l_wn',lwn,axis='l_wn')
  call put_var('append'   ,file_o,'t'   ,t  ,axis='t'   )

  call put_var('append'   ,file_o,var_name(1),var4d(:,:,:,1), &
               axes=(/'k_wn','l_wn','t'/))
  call put_var('append'   ,file_o,var_name(2),var4d(:,:,:,2), &
               axes=(/'k_wn','l_wn','t'/))

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( var4d, var_in_c )
  deallocate( kwn, lwn )
  deallocate( t )

  END subroutine finalize

END program FFT_y2

